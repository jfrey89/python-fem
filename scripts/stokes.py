#!/usr/bin/env python

from __future__ import division
import ap.mesh.meshes as m
import fem.Functions as fn
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spl
import sys


def f1(x, y):
    #fval = np.sin(np.pi * x) * np.sin(np.pi * y) * np.exp(-(x * x + y * y))
    fval = 0
    return fval


def f2(x, y):
    #fval = np.sin(np.pi * x) * np.sin(np.pi * y) * np.exp(-(x * x + y * y))
    fval = -1
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
    # USER SET PARAMETERS
    reynolds = 1e0
    perturb = 1e-2
    mesh_file = 'strip_smaller.mesh'

    eps = sys.float_info.epsilon
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
    Gx = sp.lil_matrix((n, m))
    Gy = sp.lil_matrix((n, m))
    Hx = sp.lil_matrix((m, n))
    Hy = sp.lil_matrix((m, n))
    print "Done!\n"

    fx = np.zeros(n)
    fy = fx.copy()

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
        local_Hx = sp.lil_matrix((3, 6))
        local_Hy = sp.lil_matrix((3, 6))
        local_Gx = sp.lil_matrix((6, 3))
        local_Gy = sp.lil_matrix((6, 3))

        local_fx = np.zeros(n)
        local_fy = local_fx.copy()

        for i in xrange(6):

            # load vectors
            for point in gauss_pts:
                x_g, y_g = point
                local_fx[i] += weight_scaled * f1(x_g, y_g)
                local_fy[i] += weight_scaled * f2(x_g, y_g)

            for j in xrange(6):

                # S(i, j) = (grad psi_j, grad psi_i)
                for point in gauss_pts:
                    x_g, y_g = point
                    local_S[i, j] += (weight_scaled / reynolds) * np.dot(
                        np.dot(quad_basis.grad(j, x_g, y_g), B),
                        np.dot(B.T, quad_basis.grad(i, x_g, y_g)))

                # H_x(i, j) = (d/dx psi_j, phi_i)
                # H_y(i, j) = (d/dy psi_j, phi_i)
                if i < 3:
                    for point in gauss_pts:
                        x_g, y_g = point
                        local_Hx[i, j] += weight_scaled * (
                            np.dot(quad_basis.grad(j, x_g, y_g), B[0]) *
                            lin_basis(i, x_g, y_g))

                        local_Hy[i, j] += weight_scaled * (
                            np.dot(quad_basis.grad(j, x_g, y_g), B[1]) *
                            lin_basis(i, x_g, y_g))

                # G_x(i, j) = (phi_j, d/dx psi_i)
                # G_y(i, j) = (phi_j, d/dy psi_i)
                if j < 3:
                    for point in gauss_pts:
                        x_g, y_g = point
                        local_Gx[i, j] += weight_scaled * (
                            lin_basis(j, x_g, y_g) *
                            np.dot(quad_basis.grad(i, x_g, y_g), B[0]))

                        local_Gy[i, j] += weight_scaled * (
                            lin_basis(j, x_g, y_g) *
                            np.dot(quad_basis.grad(i, x_g, y_g), B[1]))
        counter += 1
        # scatter the local matrices
        for i in xrange(6):
            fx[element[i] - 1] += local_fx[i]
            fy[element[i] - 1] += local_fy[i]

            for j in xrange(6):
                # Scatter S
                S[element[i] - 1, element[j] - 1] += local_S[i, j]

                # Scatter Hx, Hy
                # Hx, Hy is m x n
                if i < 3:
                    pnode_i = np.where(pnodes == element[i])[0][0]
                    Hx[pnode_i, element[j] - 1] += local_Hx[i, j]
                    Hy[pnode_i, element[j] - 1] += local_Hy[i, j]

                # Scatter Gx, Gy
                # Gx, Gy is n x m
                if j < 3:
                    pnode_j = np.where(pnodes == element[j])[0][0]
                    Gx[element[i] - 1, pnode_j] += local_Gx[i, j]
                    Gy[element[i] - 1, pnode_j] += local_Gy[i, j]

    print "Applying homogenous boundary conditions...\t",
    # Apply the boundary conditions
    # Assuming that velocity (u, v) = (0, 0) on the boundary
    # That implies that the mean pressure is zero over the domain
    # The mean pressure is taken into account by allowing the space for
    # (u, v) to not necessarily be divergence-free. It integrates to zeros
    # against the pressure space.
    S = S.tocsr()
    Hx = Hx.tocsr()
    Hy = Hy.tocsr()
    Gx = Gx.tocsr()
    Gy = Gy.tocsr()

    S = S[interior_nodes - 1, :]
    S = S[:, interior_nodes - 1]
    Hx = Hx[:, interior_nodes - 1]
    Hy = Hy[:, interior_nodes - 1]
    Gx = Gx[interior_nodes - 1, :]
    Gy = Gy[interior_nodes - 1, :]
    fx = fx[interior_nodes - 1, :]
    fy = fy[interior_nodes - 1, :]
    print "Done!\n"

    # Make a BIG MATRIX YAY
    k = len(interior_nodes)
    print "Constructing stiffness matrix and load vector.\n"
    A = sp.bmat([[S, sp.csr_matrix((k, k)), -Gx],
                 [sp.csr_matrix((k, k)), S, -Gy],
                 [Hx, Hy, perturb * sp.eye(m, m, format='csr')]],
                format='csr')

    f = np.concatenate((fx, fy, np.zeros(m)))

    print "Solving the linear system...\t",
    c = spl.spsolve(A, f)
    print "Done!\n"

    #print "Saving files...\t",
    #np.savetxt('./files/A.txt', A))
    #np.savetxt('./files/f.txt', f)
    np.savetxt('./files/c.txt', c)
    print 'Done!'

    u_coef = np.zeros(n)
    v_coef = u_coef.copy()

    u_coef[interior_nodes - 1] = c[:k]
    v_coef[interior_nodes - 1] = c[k: 2 * k]
    p_coef = c[2 * k:]

    # POST PROCESSING
    pcoords = coordinates[pnodes - 1]
    x_p, y_p = pcoords.T
    z_p = p_coef
    tri_p = np.zeros(np.shape(elements[:, :3]))

    for i, tri in enumerate(elements[:, :3]):
        for j, vert in enumerate(tri):
            tri_p[i, j] = np.where(pnodes == vert)[0][0]

    np.savetxt('./files/x_p.txt', x_p)
    np.savetxt('./files/y_p.txt', y_p)
    np.savetxt('./files/z_p.txt', z_p)
    np.savetxt('./files/tri_p.txt', tri_p)
