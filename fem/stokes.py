#!/usr/bin/env python

from __future__ import division
import ap.mesh.meshes as msh
import fem.Functions as fn
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spl
import sys


def f1(x, y):
    fval = y * (y - 1)
    return fval


def f2(x, y):
    fval = 0
    return fval


def get_coordinates(element, nodes):
    """
    INPUT:
        element - numpy array of vertices that make up the elements
        nodes   - list of (x, y) points by element

    OUTPUT:
        coordinates - numpy array of the (x, y) coordinates that make up the
                    elements
    """
    coordinates = np.zeros((len(element), 2))
    for index, vertex in enumerate(element):
        coordinates[index] = nodes[vertex - 1]

    return coordinates


def calculate_B(coordinates):
    # unpack coordinates
    x1, y1 = coordinates[0]
    x2, y2 = coordinates[1]
    x3, y3 = coordinates[2]

    return np.linalg.inv(np.array([[x1 - x3, x2 - x3],
                                   [y1 - y3, y2 - y3]]))


if __name__ == '__main__':
    eps = sys.float_info.epsilon
    # USER SET PARAMETERS
    reynolds = 1
    eps = 0.01
    mesh_file = './files/box.mesh'

    print "Loading and parsing the mesh..."
    domain = msh.mesh_factory(mesh_file)
    elements = domain.elements
    nodes = np.unique(elements)
    pnodes = np.unique(elements.T[:3])
    interior_nodes = domain.interior_nodes
    coordinates = domain.nodes
    n, m = len(nodes), len(pnodes)
    k = len(interior_nodes)

    weight = 1 / 6.0
    gauss_pts = np.array([[0.5, 0.0],
                          [0.0, 0.5],
                          [0.5, 0.5]])

    A = sp.lil_matrix((n, n))
    Bup = sp.lil_matrix((n, m))
    Bvp = sp.lil_matrix((n, m))
    T = sp.lil_matrix((m, m))
    Fx = np.zeros(n)
    Fy = np.zeros(n)
    uvp = np.zeros(2 * k + m)
    F = np.zeros(2 * k + m)

    quad_basis = fn.Quadratic_Basis_Function()
    lin_basis = fn.Linear_Basis_Function()

    counter = 1
    total = len(elements)
    print "Constructing system..."
    for element in elements:
        if np.mod(counter, 100) == 1:
            print "Element %d of %d..." % (counter, total)

        counter += 1

        # Precalculate some stuff
        element_coords = get_coordinates(element, domain.nodes)
        B = calculate_B(element_coords[:3])
        detJ = 1 / np.abs(np.linalg.det(B))
        weight_scaled = weight * detJ

        # Allocate local matrices
        local_A = np.zeros((6, 6))
        local_Bup = np.zeros((6, 3))
        local_Bvp = np.zeros((6, 3))
        local_Fx = np.zeros(6)
        local_Fy = np.zeros(6)
        local_T = np.zeros((3, 3))

        for point in gauss_pts:
            x_g, y_g = point

            for i in xrange(6):
                # Assemble load vectors
                local_Fx[i] += weight_scaled * f1(x_g, y_g)
                local_Fy[i] += weight_scaled * f2(x_g, y_g)

                Fx[element[i] - 1] += local_Fx[i]
                Fy[element[i] - 1] += local_Fy[i]

                for j in xrange(6):

                    local_A[i, j] += (weight_scaled / reynolds) * np.dot(
                        np.dot(quad_basis.grad(j, x_g, y_g), B),
                        np.dot(B.T, quad_basis.grad(i, x_g, y_g)))

                    A[element[i] - 1, element[j] - 1] += local_A[i, j]

                    if j < 3:
                        pnode_j = np.where(pnodes == element[j])[0][0]

                        local_Bup[i, j] += weight_scaled * (
                            lin_basis(j, x_g, y_g) *
                            np.dot(quad_basis.grad(i, x_g, y_g), B[0]))

                        local_Bvp[i, j] += weight_scaled * (
                            lin_basis(j, x_g, y_g) *
                            np.dot(quad_basis.grad(i, x_g, y_g), B[1]))

                        Bup[element[i] - 1, pnode_j] += local_Bup[i, j]
                        Bvp[element[i] - 1, pnode_j] += local_Bvp[i, j]

                        if i < 3:
                            pnode_i = np.where(pnodes == element[i])[0][0]

                            local_T[i, j] += weight_scaled * \
                                lin_basis(j, x_g, y_g) * lin_basis(i, x_g, y_g)

                            T[pnode_i, pnode_j] += local_T[i, j]

    A = A.tocsc()
    A = A[interior_nodes - 1, :]
    A = A[:, interior_nodes - 1]

    Bup = Bup.tocsc()
    Bup = Bup[interior_nodes - 1, :]

    Bvp = Bvp.tocsr()
    Bvp = Bvp[interior_nodes - 1, :]

    Fx = Fx[interior_nodes - 1, :]
    Fy = Fy[interior_nodes - 1, :]

    T = T.tocsc()

    # Update dimensions since we nuked the Dirchlet BCs, all = 0

    M = sp.bmat([[A, sp.csc_matrix((k, k)), -Bup],
                 [sp.csc_matrix((k, k)), A, -Bvp],
                 [-Bup.T, -Bvp.T, eps * T]])
    M = M.tocsc()

    F[:k] = Fx[:]
    F[k:2 * k] = Fy[:]

    print "Solving the linear system..."
    M_lu = spl.splu(M)
    uvp = M_lu.solve(F)

    u_coeffs = np.zeros(n)
    v_coeffs = np.zeros(n)
    p_coeffs = np.zeros(n)
    u_coeffs[interior_nodes - 1] = uvp[k:]
    v_coeffs[interior_nodes - 1] = uvp[k:2 * k]
    p_coeffs[pnodes - 1] = uvp[2 * k:]

    np.save('./files/u0.npy', uvp)
    np.savetxt('./files/u0.txt', uvp)
    np.savetxt('./files/tri.txt', elements)
    np.savetxt('./files/u_s.txt', u_coeffs)
    np.savetxt('./files/v_s.txt', v_coeffs)
    np.savetxt('./files/p_s.txt', p_coeffs)
    np.savetxt('./files/x.txt', coordinates)

    print "All done! Thank you, come again!\n"
