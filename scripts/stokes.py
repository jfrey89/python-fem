#!/usr/bin/env python

from __future__ import division
import ap.mesh.meshes as m
import fem.Functions as fn
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spl
import sys
import matplotlib.pyplot as pyplot


def f1(x, y):
    #fval = np.sin(np.pi * x) * np.sin(np.pi * y) * np.exp(-(x * x + y * y))
    fval = -y
    return fval


def f2(x, y):
    #fval = np.sin(np.pi * x) * np.sin(np.pi * y) * np.exp(-(x * x + y * y))
    fval = x
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
    eps = 1e-4 / reynolds
    mesh_file = 'box-circle.mesh'

    root_dir = './files/'
    print "Loading and parsing the mesh...\t",
    domain = m.mesh_factory(root_dir + mesh_file)
    print "Done!\n"
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

    print "Allocating individual matrices...\t",
    S = sp.lil_matrix((n, n))
    T = sp.lil_matrix((m, m))
    Gx = sp.lil_matrix((n, m))
    Gy = sp.lil_matrix((n, m))
    fx = np.zeros(n)
    fy = fx.copy()
    print "Done!\n"

    quad_basis = fn.Quadratic_Basis_Function()
    lin_basis = fn.Linear_Basis_Function()

    counter = 1
    total = len(elements)
    print "Looping over elements to construct local matrices..."
    print "Element 1 of %d..." % total

    for element in elements:
        if np.mod(counter, 100) == 0:
            print "Element %d of %d..." % (counter, total)
        # Precalculate some stuff
        element_coords = get_coordinates(element, domain.nodes)
        B = calculate_B(element_coords[:3])
        detJ = 1 / np.abs(np.linalg.det(B))
        weight_scaled = weight * detJ

        # Allocate local matrices
        local_S = sp.lil_matrix((6, 6))
        local_T = sp.lil_matrix((3, 3))
        local_Gx = sp.lil_matrix((6, 3))
        local_Gy = sp.lil_matrix((6, 3))
        local_fx = np.zeros(n)
        local_fy = np.zeros(n)

        for i in xrange(6):
            # Assemble load vectors
            for point in gauss_pts:
                x_g, y_g = point
                local_fx[i] += weight_scaled * f1(x_g, y_g)
                local_fy[i] += weight_scaled * f2(x_g, y_g)
            for j in xrange(6):
                # Assemble S
                for point in gauss_pts:
                    x_g, y_g = point
                    local_S[i, j] += (weight_scaled / reynolds) * np.dot(
                        np.dot(quad_basis.grad(j, x_g, y_g), B),
                        np.dot(B.T, quad_basis.grad(i, x_g, y_g)))

                if j < 3:
                    # Assemble Gx and Gy
                    for point in gauss_pts:
                        x_g, y_g = point
                        local_Gx[i, j] += weight_scaled * (
                            lin_basis(j, x_g, y_g) *
                            np.dot(quad_basis.grad(i, x_g, y_g), B[0]))
                        local_Gy[i, j] += weight_scaled * (
                            lin_basis(j, x_g, y_g) *
                            np.dot(quad_basis.grad(i, x_g, y_g), B[1]))

                    if i < 3:
                        for point in gauss_pts:
                            x_g, y_g = point
                            local_T[i, j] += weight_scaled * (
                                lin_basis(j, x_g, y_g) *
                                lin_basis(i, x_g, y_g))

                    #for point in gauss_pts:
                        #x_g, y_g = point
                        #local_Hx[i, j] += weight_scaled * (
                            #np.dot(quad_basis.grad(j, x_g, y_g), B[0]) *
                            #lin_basis(i, x_g, y_g))
                        #local_Hy[i, j] += weight_scaled * (
                            #np.dot(quad_basis.grad(j, x_g, y_g), B[1]) *
                            #lin_basis(i, x_g, y_g))

        counter += 1
        # scatter the local matrices
        for i in xrange(6):
            fx[element[i] - 1] += local_fx[i]
            fy[element[i] - 1] += local_fy[i]
            for j in xrange(6):
                # Scatter S
                S[element[i] - 1, element[j] - 1] += local_S[i, j]
                # Scatter Gx, Gy
                if j < 3:
                    pnode_j = np.where(pnodes == element[j])[0][0]
                    Gx[element[i] - 1, pnode_j] += local_Gx[i, j]
                    Gy[element[i] - 1, pnode_j] += local_Gy[i, j]
                    if i < 3:
                        pnode_i = np.where(pnodes == element[i])[0][0]
                        T[pnode_i, pnode_j] += local_T[i, j]

    S = S.tocsr()
    T = T.tocsr()
    Gx = Gx.tocsr()
    Gy = Gy.tocsr()

    print "Applying homogenous boundary conditions...\t",
    S = S[interior_nodes - 1, :]
    S = S[:, interior_nodes - 1]
    Gx = Gx[interior_nodes - 1, :]
    Gy = Gy[interior_nodes - 1, :]
    fx = fx[interior_nodes - 1, :]
    fy = fy[interior_nodes - 1, :]
    print "Done!\n"

    # Make a BIG MATRIX YAY
    k = len(interior_nodes)
    # But preallocate it first...
    print "Constructing stiffness matrix with -eps...\t",
    Am = sp.bmat([[S, sp.csr_matrix((k, k)), -Gx],
                 [sp.csr_matrix((k, k)), S, -Gy],
                 [-Gx.T, -Gy.T, -eps * T]])
    Am = Am.tocsr()
    print "Done!\n"

    print "Constructing stiffness matrix with +eps...\t",
    Ap = sp.bmat([[S, sp.csr_matrix((k, k)), -Gx],
                 [sp.csr_matrix((k, k)), S, -Gy],
                 [-Gx.T, -Gy.T, eps * T]])
    Ap = Ap.tocsr()
    print "Done!\n"

    print "Constructing alternative system...\t",
    Ar = sp.bmat([[S, np.zeros((k, k))],
                  [np.zeros((k, k)), S]])
    Br = sp.bmat([[Gx.dot(Gx.T), Gx.dot(Gy.T)],
                  [Gy.dot(Gx.T), Gy.dot(Gy.T)]])

    Ar = Ar.tocsr()
    Br = Br.tocsr()
    print "Done!\n"

    # Change data type
    f = np.zeros(2 * k + m)
    print "Constructing f...\t",
    f[:k] = fx[:]
    f[k:2 * k] = fy[:]
    print "Done!\n"
    fr = f[:2 * k]

    print "Solving the linear systems...\t",
    cm = spl.spsolve(Am, f)
    print "1...",
    cp = spl.spsolve(Ap, f)
    print "2...",
    cr = spl.spsolve(Ar - Br / eps, fr)
    print "3...!\nDone!"

    # POST PROCESSING for system one
    x_p, y_p = coordinates[pnodes - 1].T
    tri_p = np.zeros(np.shape(elements[:, :3]))
    for i, tri in enumerate(elements[:, :3]):
        for j, vert in enumerate(tri):
            tri_p[i, j] = np.where(pnodes == vert)[0][0]

    np.savetxt('./files/tri.txt', tri_p)
    np.savetxt('./files/x_p.txt', x_p)
    np.savetxt('./files/y_p.txt', y_p)
    np.savetxt('./files/pm.txt', cm[2 * k:])
    np.savetxt('./files/pp.txt', cp[2 * k:])

    np.savetxt('./files/x.txt', coordinates[:, 0])
    np.savetxt('./files/y.txt', coordinates[:, 1])

    np.savetxt('./files/um.txt', cm[:k])
    np.savetxt('./files/vm.txt', cm[k:2 * k])
    np.savetxt('./files/up.txt', cp[:k])
    np.savetxt('./files/vp.txt', cp[k:2 * k])
    np.savetxt('./files/ur.txt', cr[:k])
    np.savetxt('./files/vr.txt', cr[k:2 * k])
