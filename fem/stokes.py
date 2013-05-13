#!/usr/bin/env python

from __future__ import division
import ap.mesh.meshes as msh
import fem.Functions as fn
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spl
import sys


def f1(x, y):
    fval = np.power(x, 2) * np.power(y, 3) * np.cos(y * np.pi)
    return fval


def f2(x, y):
    fval = 1e-4 * np.abs((y - 0.5))
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
    reynolds = 1e0
    o = 2
    eps = 1e-4 / reynolds
    mesh_file = '2box-circle.0.025.mesh'

    root_dir = './files/'
    print "Loading and parsing the mesh..."
    domain = msh.mesh_factory(root_dir + mesh_file)
    elements = domain.elements
    nodes = np.unique(elements)
    pnodes = np.unique(elements.T[:3])
    interior_nodes = domain.interior_nodes
    coordinates = domain.nodes
    n, m = len(nodes), len(pnodes)

    weight = 1 / 6.0
    gauss_pts = np.array([[0.5, 0.0],
                          [0.0, 0.5],
                          [0.5, 0.5]])

    A = sp.lil_matrix((n, n))
    B_up = sp.lil_matrix((n, m))
    B_vp = sp.lil_matrix((n, m))
    F_x = np.zeros(n)
    F_y = np.zeros(n)

    quad_basis = fn.Quadratic_Basis_Function()
    lin_basis = fn.Linear_Basis_Function()

    counter = 1
    total = len(elements)
    print "Constructing system..."
    for element in elements:
        if np.mod(counter, 100) == 1:
            print "Element %d of %d..." % (counter, total)
        # Precalculate some stuff
        element_coords = get_coordinates(element, domain.nodes)
        B = calculate_B(element_coords[:3])
        detJ = 1 / np.abs(np.linalg.det(B))
        weight_scaled = weight * detJ

        # Allocate local matrices
        local_A = sp.lil_matrix((6, 6))
        local_B_up = sp.lil_matrix((6, 3))
        local_B_vp = sp.lil_matrix((6, 3))
        local_F_x = np.zeros(6)
        local_F_y = np.zeros(6)

        for i in xrange(6):
            # Assemble load vectors
            for point in gauss_pts:
                x_g, y_g = point
                local_F_x[i] += weight_scaled * f1(x_g, y_g)
                local_F_y[i] += weight_scaled * f2(x_g, y_g)

            F_x[element[i] - 1] += local_F_x[i]
            F_y[element[i] - 1] += local_F_y[i]

            for j in xrange(6):
                # Assemble S
                for point in gauss_pts:
                    x_g, y_g = point
                    local_A[i, j] += (weight_scaled / reynolds) * np.dot(
                        np.dot(quad_basis.grad(j, x_g, y_g), B),
                        np.dot(B.T, quad_basis.grad(i, x_g, y_g)))

                A[element[i] - 1, element[j] - 1] += local_A[i, j]

                if j < 3:
                    # Assemble Gx and Gy
                    for point in gauss_pts:
                        x_g, y_g = point
                        local_B_up[i, j] += weight_scaled * (
                            lin_basis(j, x_g, y_g) *
                            np.dot(quad_basis.grad(i, x_g, y_g), B[0]))
                        local_B_vp[i, j] += weight_scaled * (
                            lin_basis(j, x_g, y_g) *
                            np.dot(quad_basis.grad(i, x_g, y_g), B[1]))

                    pnode_j = np.where(pnodes == element[j])[0][0]
                    B_up[element[i] - 1, pnode_j] += local_B_up[i, j]
                    B_vp[element[i] - 1, pnode_j] += local_B_vp[i, j]

        counter += 1

    A = A.tocsr()
    A = A[interior_nodes - 1, :]
    A = A.tocsc()
    A = A[:, interior_nodes - 1]

    B_up = B_up.tocsr()
    B_up = B_up[interior_nodes - 1, :]

    B_vp = B_vp.tocsr()
    B_vp = B_vp[interior_nodes - 1, :]

    F_x = F_x[interior_nodes - 1, :]
    F_y = F_y[interior_nodes - 1, :]

    # Update dimensions since we nuked the Dirchlet BCs, all = 0
    k = len(interior_nodes)

    M = sp.bmat([[A, sp.csc_matrix((k, k))],
                 [sp.csc_matrix((k, k)), A]]) + \
        sp.bmat([[B_up.dot(B_up.T), B_up.dot(B_vp.T)],
                 [B_vp.dot(B_up.T), B_vp.dot(B_vp.T)]]) / eps
    M = M.tocsc()

    print "Solving the linear system..."
    UV = np.zeros(2 * k)
    LU = spl.splu(M)
    UV = LU.solve(np.append(F_x, F_y))

    P = np.zeros(m)
    P = -sp.bmat([[B_up.T, B_vp.T]]).dot(UV) / eps
    UVP = np.zeros(2 * k + m)
    UVP[:k] = UV[:k]
    UVP[k:2 * k] = UV[k:2 * k]
    UVP[2 * k:] = P[:]

    U = np.zeros(n)
    U[interior_nodes - 1] = UVP[:k]
    V = np.zeros(n)
    V[interior_nodes - 1] = UVP[k:2 * k]

    # POST PROCESSING for system one
    x, y = coordinates[pnodes - 1].T
    tri = np.zeros(np.shape(elements[:, :3]))
    for i, triangle in enumerate(elements[:, :3]):
        for j, vert in enumerate(triangle):
            tri[i, j] = np.where(pnodes == vert)[0][0]

    F = np.zeros(2 * k + m)
    F[:k] = F_x[:]
    F[k:2 * k] = F_y[:]
    print "Saving some data..."
    u = np.zeros(n)
    u[interior_nodes - 1] = UVP[k:]
    v = np.zeros(n)
    v[interior_nodes - 1] = UVP[k:2 * k]
    np.save('./files/UVP.npy', UVP)
    np.save('./files/F.npy', F)
    np.savetxt('./files/UVP.txt', UVP)
    np.savetxt('./files/tri.txt', tri)
    np.savetxt('./files/x.txt', x)
    np.savetxt('./files/y.txt', y)
    np.savetxt('./files/p.txt', P)
    np.savetxt('./files/u.txt', u[pnodes - 1])
    np.savetxt('./files/v.txt', v[pnodes - 1])
    print "All done! Thank you, come again!\n\n"
