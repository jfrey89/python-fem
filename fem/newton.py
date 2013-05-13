#!/usr/bin/env python

from __future__ import division
import ap.mesh.meshes as msh
import fem.Functions as fn
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spl


def f1(x, y):
    fval = -1e-4 * y
    return fval


def f2(x, y):
    fval = 1e-4 * x
    return fval


def calculate_B(coordinates):
    # unpack coordinates
    x1, y1 = coordinates[0]
    x2, y2 = coordinates[1]
    x3, y3 = coordinates[2]

    return np.linalg.inv(np.array([[x1 - x3, x2 - x3],
                                   [y1 - y3, y2 - y3]]))


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

if __name__ == "__main__":
    reynolds = 1e0
    mesh_file = 'box0.1.mesh'
    root_dir = './files/'
    eps = 1e-4 / reynolds

    uvp = np.load(root_dir + 'UVP.npy')
    z = uvp.copy()

    domain = msh.mesh_factory(root_dir + mesh_file)
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

    quad_basis = fn.Quadratic_Basis_Function()
    lin_basis = fn.Linear_Basis_Function()

    Fx = np.zeros(n)
    Fy = np.zeros(n)

    total = len(elements)

    print "Starting loop 5evr.\n"
    while True:
        P = np.zeros(m)
        U = np.zeros(n)
        V = np.zeros(n)
        U[interior_nodes - 1] = uvp[:k]
        V[interior_nodes - 1] = uvp[k:2 * k]
        P[:] = uvp[2 * k:]

        Auu = sp.lil_matrix((n, n))
        Auv = sp.lil_matrix((n, n))
        Avu = sp.lil_matrix((n, n))
        Avv = sp.lil_matrix((n, n))
        Bup = sp.lil_matrix((m, n))
        Bvp = sp.lil_matrix((m, n))
        Cup = sp.lil_matrix((n, m))
        Cvp = sp.lil_matrix((n, m))

        counter = 1
        for element in elements:
            if np.mod(counter, 100) == 1:
                print "Element %d of %d..." % (counter, total)

            ue = U[element - 1]
            ve = V[element - 1]

            pnode_0 = np.where(pnodes == element[0])[0][0]
            pnode_1 = np.where(pnodes == element[1])[0][0]
            pnode_2 = np.where(pnodes == element[2])[0][0]

            pe = np.zeros(3)
            pe[0] = P[pnode_0]
            pe[1] = P[pnode_1]
            pe[2] = P[pnode_2]

            counter += 1

            local_Bup = np.zeros((3, 6))
            local_Bvp = np.zeros((3, 6))
            local_Cup = np.zeros((6, 3))
            local_Cvp = np.zeros((6, 3))
            local_Auu = np.zeros((6, 6))
            local_Auv = np.zeros((6, 6))
            local_Avu = np.zeros((6, 6))
            local_Avv = np.zeros((6, 6))

            element_coords = get_coordinates(element, domain.nodes)
            B = calculate_B(element_coords[:3])
            detJ = 1 / np.abs(np.linalg.det(B))
            weight_scaled = weight * detJ

            u_f = lambda x, y: ue[0] * quad_basis(0, x, y) + ue[1] * \
                quad_basis(1, x, y) + ue[2] * quad_basis(2, x, y) + \
                ue[3] * quad_basis(3, x, y) + \
                ue[4] * quad_basis(4, x, y) + ue[5] * quad_basis(5, x, y)

            v_f = lambda x, y: ve[0] * quad_basis(0, x, y) + ve[1] * \
                quad_basis(1, x, y) + ve[2] * quad_basis(2, x, y) + \
                ve[3] * quad_basis(3, x, y) + \
                ve[4] * quad_basis(4, x, y) + ve[5] * quad_basis(5, x, y)

            du_dx = lambda x, y: ue[0] * np.dot(
                quad_basis.grad(0, x, y), B[0]) + ue[1] * np.dot(
                    quad_basis.grad(1, x, y), B[0]) + \
                ue[2] * np.dot(quad_basis.grad(2, x, y), B[0]) + \
                ue[3] * np.dot(quad_basis.grad(3, x, y), B[0]) + \
                ue[4] * np.dot(quad_basis.grad(4, x, y), B[0]) + \
                ue[5] * np.dot(quad_basis.grad(5, x, y), B[0])

            du_dy = lambda x, y: ue[0] * np.dot(
                quad_basis.grad(0, x, y), B[1]) + ue[1] * np.dot(
                    quad_basis.grad(1, x, y), B[1]) + \
                ue[2] * np.dot(quad_basis.grad(2, x, y), B[1]) + \
                ue[3] * np.dot(quad_basis.grad(3, x, y), B[1]) + \
                ue[4] * np.dot(quad_basis.grad(4, x, y), B[1]) + \
                ue[5] * np.dot(quad_basis.grad(5, x, y), B[1])

            dv_dx = lambda x, y: ve[0] * np.dot(
                quad_basis.grad(0, x, y), B[0]) + ve[1] * np.dot(
                    quad_basis.grad(1, x, y), B[0]) + \
                ve[2] * np.dot(quad_basis.grad(2, x, y), B[0]) + \
                ve[3] * np.dot(quad_basis.grad(3, x, y), B[0]) + \
                ve[4] * np.dot(quad_basis.grad(4, x, y), B[0]) + \
                ve[5] * np.dot(quad_basis.grad(5, x, y), B[0])

            dv_dy = lambda x, y: ve[0] * np.dot(
                quad_basis.grad(0, x, y), B[1]) + ve[1] * np.dot(
                    quad_basis.grad(1, x, y), B[1]) + \
                ve[2] * np.dot(quad_basis.grad(2, x, y), B[1]) + \
                ve[3] * np.dot(quad_basis.grad(3, x, y), B[1]) + \
                ve[4] * np.dot(quad_basis.grad(4, x, y), B[1]) + \
                ve[5] * np.dot(quad_basis.grad(5, x, y), B[1])

            p_f = lambda x, y: pe[0] * lin_basis(0, x, y) + \
                pe[1] * lin_basis(1, x, y) + pe[2] * lin_basis(2, x, y)

            dp_dx = lambda x, y: pe[0] * np.dot(
                lin_basis.grad(0), B[0]) + pe[1] * np.dot(
                    lin_basis.grad(1), B[0]) + pe[2] * np.dot(
                        lin_basis.grad(2), B[0])

            dp_dy = lambda x, y: pe[0] * np.dot(
                lin_basis.grad(0), B[1]) + pe[1] * np.dot(
                    lin_basis.grad(1), B[1]) + pe[2] * np.dot(
                        lin_basis.grad(2), B[1])

            for i in xrange(6):



                for j in xrange(6):

                    u = U[element[j] - 1]
                    v = V[element[j] - 1]

                    for point in gauss_pts:
                        x_g, y_g = point

                        local_Auu[i, j] += (weight_scaled / reynolds) * \
                            np.dot(
                                np.dot(quad_basis.grad(j, x_g, y_g), B),
                                np.dot(B.T, quad_basis.grad(i, x_g, y_g))) + \
                            weight_scaled * (u_f(x_g, y_g) +
                                             np.dot(
                                             quad_basis.grad(j, x_g, y_g),
                                             B[0]) + v_f(x_g, y_g) *
                                             np.dot(
                                             quad_basis.grad(j, x_g, y_g),
                                             B[1]) +
                                             quad_basis(j, x_g, y_g) *
                                             du_dx(x_g, y_g)) * \
                            quad_basis(i, x_g, y_g)

                        local_Avv[i, j] += (weight_scaled / reynolds) * \
                            np.dot(
                                np.dot(quad_basis.grad(j, x_g, y_g), B),
                                np.dot(B.T, quad_basis.grad(i, x_g, y_g))) + \
                            weight_scaled * (u_f(x_g, y_g) *
                                             np.dot(
                                             quad_basis.grad(j, x_g, y_g),
                                             B[0]) + v_f(x_g, y_g) *
                                             np.dot(
                                             quad_basis.grad(j, x_g, y_g),
                                             B[1]) +
                                             quad_basis(j, x_g, y_g) *
                                             dv_dy(x_g, y_g)) * \
                            quad_basis(i, x_g, y_g)

                        local_Auv[i, j] += weight_scaled * (
                            quad_basis(j, x_g, y_g) * du_dy(x_g, y_g)) * \
                            quad_basis(i, x_g, y_g)

                        local_Avu[i, j] += weight_scaled * (
                            quad_basis(j, x_g, y_g) * dv_dx(x_g, y_g)) * \
                            quad_basis(i, x_g, y_g)

                    Auu[element[i] - 1, element[j] - 1] += local_Auu[i, j]
                    Auv[element[i] - 1, element[j] - 1] += local_Auv[i, j]
                    Avu[element[i] - 1, element[j] - 1] += local_Avu[i, j]
                    Avv[element[i] - 1, element[j] - 1] += local_Avv[i, j]

            for i in xrange(6):
                for j in xrange(3):
                    pnode_j = np.where(pnodes == element[j])[0][0]
                    for point in gauss_pts:
                        x_g, y_g = point

                        local_Cup[i, j] += weight_scaled * np.dot(
                            lin_basis.grad(j), B[0]) * \
                            quad_basis(i, x_g, y_g)

                        local_Cvp[i, j] += weight_scaled * np.dot(
                            lin_basis.grad(j), B[1]) * \
                            quad_basis(i, x_g, y_g)

                    Cup[element[i] - 1, pnode_j] += local_Cup[i, j]
                    Cvp[element[i] - 1, pnode_j] += local_Cvp[i, j]

            for i in xrange(3):
                for j in xrange(6):
                    pnode_i = np.where(pnodes == element[i])[0][0]
                    for point in gauss_pts:
                        x_g, y_g = point

                        local_Bup[i, j] += weight_scaled * (np.dot(
                            quad_basis.grad(j, x_g, y_g), B[0]) *
                            lin_basis(i, x_g, y_g))

                        local_Bvp[i, j] += weight_scaled * (np.dot(
                            quad_basis.grad(j, x_g, y_g), B[1]) *
                            lin_basis(i, x_g, y_g))

                    Bup[pnode_i, element[j] - 1] += local_Bup[i, j]
                    Bvp[pnode_i, element[j] - 1] += local_Bvp[i, j]

        Auu = Auu.tocsr()
        Auu = Auu[interior_nodes - 1, :]
        Auu = Auu.tocsc()
        Auu = Auu[:, interior_nodes - 1]

        Auv = Auv.tocsr()
        Auv = Auv[interior_nodes - 1, :]
        Auv = Auv.tocsc()
        Auv = Auv[:, interior_nodes - 1]

        Avv = Avv.tocsr()
        Avv = Avv[interior_nodes - 1, :]
        Avv = Avv.tocsc()
        Avv = Avv[:, interior_nodes - 1]

        Avu = Avu.tocsr()
        Avu = Avu[interior_nodes - 1, :]
        Avu = Avu.tocsc()
        Avu = Avu[:, interior_nodes - 1]

        Bup = Bup.tocsc()
        Bup = Bup[:, interior_nodes - 1]
        Bvp = Bvp.tocsc()
        Bvp = Bvp[:, interior_nodes - 1]

        Cup = Cup.tocsc()
        Cup = Cup[interior_nodes - 1, :]
        Cvp = Cvp.tocsc()
        Cvp = Cvp[interior_nodes - 1, :]

        M = sp.bmat([[Auu, Auv, -Cup],
                    [Avu, Avv, -Cvp],
                    [Bup, Bvp, -eps * sp.eye(m, m)]])

        M = M.tocsc()

        d_k = np.zeros(len(uvp))
        d_k = spl.spsolve(M, -uvp)

        uvp = uvp - d_k

        print "Solving the system."

        print "Norm:\t\t%0.2f" % (np.linalg.norm(d_k))

        if np.linalg.norm(d_k) < 1e-1:

            break
        print 'AGAIN!'
