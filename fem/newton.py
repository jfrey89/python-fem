#!/usr/bin/env python

from __future__ import division
import ap.mesh.meshes as msh
import fem.Functions as fn
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spl


def f1(x, y):
    fval = y * (y - 1)
    return fval


def f2(x, y):
    fval = 0
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


def u_f(x, y, ue, quad_basis):
    u = ue[0] * quad_basis(0, x, y) + ue[1] * quad_basis(1, x, y) + \
        ue[2] * quad_basis(2, x, y) + ue[3] * quad_basis(3, x, y) + \
        ue[4] * quad_basis(4, x, y) + ue[5] * quad_basis(5, x, y)
    return u


def v_f(x, y, ve, quad_basis):
    v = ve[0] * quad_basis(0, x, y) + ve[1] * quad_basis(1, x, y) + \
        ve[2] * quad_basis(2, x, y) + ve[3] * quad_basis(3, x, y) + \
        ve[4] * quad_basis(4, x, y) + ve[5] * quad_basis(5, x, y)
    return v


def du_dx(x, y, B, ue, quad_basis):
    du_dx = ue[0] * np.dot(quad_basis.grad(0, x, y), B[0]) + \
        ue[1] * np.dot(quad_basis.grad(1, x, y), B[0]) + \
        ue[2] * np.dot(quad_basis.grad(2, x, y), B[0]) + \
        ue[3] * np.dot(quad_basis.grad(3, x, y), B[0]) + \
        ue[4] * np.dot(quad_basis.grad(4, x, y), B[0]) + \
        ue[5] * np.dot(quad_basis.grad(5, x, y), B[0])
    return du_dx


def du_dy(x, y, B, ue, quad_basis):
    du_dy = ue[0] * np.dot(quad_basis.grad(0, x, y), B[1]) + \
        ue[1] * np.dot(quad_basis.grad(1, x, y), B[1]) + \
        ue[2] * np.dot(quad_basis.grad(2, x, y), B[1]) + \
        ue[3] * np.dot(quad_basis.grad(3, x, y), B[1]) + \
        ue[4] * np.dot(quad_basis.grad(4, x, y), B[1]) + \
        ue[5] * np.dot(quad_basis.grad(5, x, y), B[1])
    return du_dy


def dv_dx(x, y, B, ve, quad_basis):
    dv_dx = ve[0] * np.dot(quad_basis.grad(0, x, y), B[0]) + \
        ve[1] * np.dot(quad_basis.grad(1, x, y), B[0]) + \
        ve[2] * np.dot(quad_basis.grad(2, x, y), B[0]) + \
        ve[3] * np.dot(quad_basis.grad(3, x, y), B[0]) + \
        ve[4] * np.dot(quad_basis.grad(4, x, y), B[0]) + \
        ve[5] * np.dot(quad_basis.grad(5, x, y), B[0])
    return dv_dx


def dv_dy(x, y, B, ve, quad_basis):
    dv_dy = ve[0] * np.dot(quad_basis.grad(0, x, y), B[1]) + \
        ve[1] * np.dot(quad_basis.grad(1, x, y), B[1]) + \
        ve[2] * np.dot(quad_basis.grad(2, x, y), B[1]) + \
        ve[3] * np.dot(quad_basis.grad(3, x, y), B[1]) + \
        ve[4] * np.dot(quad_basis.grad(4, x, y), B[1]) + \
        ve[5] * np.dot(quad_basis.grad(5, x, y), B[1])
    return dv_dy

if __name__ == "__main__":
    eps = 0.001
    reynolds = 1
    mesh_file = './files/box.mesh'

    domain = msh.mesh_factory(mesh_file)
    elements = domain.elements
    nodes = np.unique(elements)
    pnodes = np.unique(elements.T[:3])
    interior_nodes = domain.interior_nodes
    coordinates = domain.nodes
    n, m = len(nodes), len(pnodes)
    k = len(interior_nodes)

    u_0 = np.zeros(2 * k + m)
    u_0 = np.load('./files/u0.npy')

    weight = 1 / 6.0
    gauss_pts = np.array([[0.5, 0.0],
                          [0.0, 0.5],
                          [0.5, 0.5]])

    quad_basis = fn.Quadratic_Basis_Function()
    lin_basis = fn.Linear_Basis_Function()

    total = len(elements)
    newt_count = 1

    coeffs = u_0.copy()
    d_n = np.zeros(2 * k + m)

    print "Starting loop 5evr.\n"
    while True:
        print '*' * 40
        print '*' * 40
        print '\nNEWTON METHOD ITERATION NUMBER: %d\n' % newt_count
        print '*' * 40
        print '*' * 40
        print '\n'
        P = np.zeros(m)
        U = np.zeros(n)
        V = np.zeros(n)
        U[interior_nodes - 1] = coeffs[:k]
        V[interior_nodes - 1] = coeffs[k:2 * k]
        P[:] = coeffs[2 * k:]

        Auu = sp.lil_matrix((n, n))
        Auv = sp.lil_matrix((n, n))
        Avu = sp.lil_matrix((n, n))
        Avv = sp.lil_matrix((n, n))
        Bup = sp.lil_matrix((n, m))
        Bvp = sp.lil_matrix((n, m))

        rhsx = np.zeros(n)
        rhsy = np.zeros(n)
        rhs = np.zeros(2 * k + m)

        counter = 1
        for element in elements:
            if np.mod(counter, 100) == 1:
                print "Element %d of %d..." % (counter, total)

            counter += 1

            ue = np.zeros(6)
            ve = np.zeros(6)

            ue = U[element - 1]
            ve = V[element - 1]

            element_coords = get_coordinates(element, domain.nodes)

            B = calculate_B(element_coords[:3])
            detJ = 1 / np.abs(np.linalg.det(B))
            weight_scaled = weight * detJ

            for point in gauss_pts:
                x_g, y_g = point

                for i in xrange(6):
                    rhsy[element[i] - 1] += weight_scaled * (
                        u_f(x_g, y_g, ue, quad_basis) *
                        du_dx(x_g, y_g, B, ue, quad_basis) +
                        v_f(x_g, y_g, ve, quad_basis) *
                        du_dy(x_g, y_g, B, ue, quad_basis)) * \
                        quad_basis(i, x_g, y_g)

                    rhsy[element[i] - 1] += weight_scaled * (
                        u_f(x_g, y_g, ue, quad_basis) *
                        dv_dx(x_g, y_g, B, ve, quad_basis) +
                        v_f(x_g, y_g, ve, quad_basis) *
                        dv_dy(x_g, y_g, B, ve, quad_basis)) * \
                        quad_basis(i, x_g, y_g)

                    for j in xrange(6):

                        db_dx = np.dot(quad_basis.grad(j, x_g, y_g), B[0])
                        db_dy = np.dot(quad_basis.grad(j, x_g, y_g), B[1])

                        Auu[element[i] - 1, element[j] - 1] += \
                            (weight_scaled / reynolds) * np.dot(
                                np.dot(quad_basis.grad(j, x_g, y_g), B),
                                np.dot(B.T, quad_basis.grad(i, x_g, y_g))) + \
                            weight_scaled * (
                                u_f(x_g, y_g, ue, quad_basis) * db_dx +
                                v_f(x_g, y_g, ve, quad_basis) * db_dy) * \
                            quad_basis(i, x_g, y_g)

                        Avv[element[i] - 1, element[j] - 1] += \
                            Auu[element[i] - 1, element[j] - 1]

                        Auu[element[i] - 1, element[j] - 1] += \
                            weight_scaled * (
                                quad_basis(j, x_g, y_g) *
                                du_dx(x_g, y_g, B, ue, quad_basis) *
                                quad_basis(i, x_g, y_g))

                        Avv[element[i] - 1, element[j] - 1] += \
                            weight_scaled * (
                                quad_basis(j, x_g, y_g) *
                                dv_dy(x_g, y_g, B, ve, quad_basis) *
                                quad_basis(i, x_g, y_g))

                        Auv[element[i] - 1, element[j] - 1] += \
                            weight_scaled * (
                                quad_basis(j, x_g, y_g) *
                                du_dy(x_g, y_g, B, ue, quad_basis))

                        Avu[element[i] - 1, element[j] - 1] += \
                            weight_scaled * (
                                quad_basis(j, x_g, y_g) *
                                dv_dx(x_g, y_g, B, ve, quad_basis))

                        if j < 3:
                            pnode_j = np.where(pnodes == element[j])[0][0]

                            db_dx = np.dot(quad_basis.grad(i, x_g, y_g), B[0])
                            db_dy = np.dot(quad_basis.grad(i, x_g, y_g), B[1])

                            Bup[element[i] - 1, pnode_j] +=  \
                                weight_scaled * db_dx * \
                                lin_basis(j, x_g, y_g)

                            Bvp[element[i] - 1, pnode_j] +=  \
                                weight_scaled * db_dy * \
                                lin_basis(j, x_g, y_g)

        Auu = Auu.tocsc()
        Auu = Auu[interior_nodes - 1, :]
        Auu = Auu[:, interior_nodes - 1]

        Auv = Auv.tocsc()
        Auv = Auv[interior_nodes - 1, :]
        Auv = Auv[:, interior_nodes - 1]

        Avv = Avv.tocsc()
        Avv = Avv[interior_nodes - 1, :]
        Avv = Avv[:, interior_nodes - 1]

        Avu = Avu.tocsc()
        Avu = Avu[interior_nodes - 1, :]
        Avu = Avu[:, interior_nodes - 1]

        Bup = Bup.tocsc()
        Bup = Bup[interior_nodes - 1, :]

        Bvp = Bvp.tocsc()
        Bvp = Bvp[interior_nodes - 1, :]

        rhsx = rhsx[interior_nodes - 1]
        rhs[:k] = rhsx[:]

        rhsy = rhsy[interior_nodes - 1]
        rhs[k:2 * k] = rhsy[:]

        if np.allclose(np.linalg.norm(rhs), 0):
            print "||residual|| sufficiently close to 0!"
            print '*' * 60
            print "NEWTON'S METHOD WAS A SUCCESS."
            print "HOORAY" * 10
            print "IT TOOK:\t%d ITERATIONS" % newt_count
            print "*" * 60
            np.savetxt('./files/ns_soln.txt', coeffs)
            np.save('./files/ns_soln.npy', coeffs)
            print "Have a nice day!\n\n"
            break

        M = sp.bmat([[Auu, Auv, -Bup],
                     [Avu, Avv, -Bvp],
                     [-Bup.T, -Bvp.T, sp.csc_matrix((m, m))]])
        M = M.tocsc()

        M_lu = spl.splu(M)
        print "Solving the system."
        d_n = M_lu.solve(rhs)
        # Update
        u_n = np.zeros(2 * k + m)
        u_n = coeffs + d_n

        print "||d_n|| = %0.10f" % (np.linalg.norm(d_n))

        if np.allclose(np.linalg.norm(d_n), 0):
            print "||update vector|| sufficiently close to 0!"
            print '*' * 60
            print "NEWTON'S METHOD WAS A SUCCESS."
            print "HOORAY" * 10
            print "IT TOOK:\t%d ITERATIONS" % newt_count
            print "*" * 60
            np.savetxt('./files/ns_soln.txt', u_n)
            np.save('./files/ns_soln.npy', u_n)
            print "Have a nice day!\n\n"
            break

        coeffs = u_n.copy()
        newt_count += 1

        print 'AGAIN!\n\n'
