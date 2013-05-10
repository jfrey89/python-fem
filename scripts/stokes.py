#!/usr/bin/env python

from __future__ import division
import ap.mesh.meshes as m
import fem.Functions as fn
import numpy as np
import sys


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
    root_dir = './files/'
    mesh_file = 'unit-square.mesh'
    # load mesh
    domain = m.mesh_factory(root_dir + mesh_file)
    elements = domain.elements
    nodes = np.unique(elements)
    pnodes = np.unique(elements.T[:3])
    interior_nodes = domain.interior_nodes

    nodes_2_coords = domain.nodes

    n, m = len(nodes), len(pnodes)

    reynolds = 1.0
    weight = 1 / 6.0
    gauss_pts = np.array([[0.5, 0.0],
                          [0.0, 0.5],
                          [0.5, 0.5]])

    S = np.zeros((n, n))
    Gx = np.zeros((n, m))
    Hx = np.zeros((m, n))
    Gy = Gx.copy()
    Hy = Hx.copy()

    f = np.zeros(3 * len(S))
    f1 = np.zeros(len(S))
    f2 = f1.copy()

    quad_basis = fn.Quadratic_Basis_Function()
    lin_basis = fn.Linear_Basis_Function()

    for count, element in enumerate(elements):
        print '*' * 20
        print "Element %d of %d" % (count + 1, len(elements))
        print '#' * 20

        coordinates = get_coordinates(element, nodes_2_coords)
        B = calculate_B(coordinates[:3])
        detJ = 1 / np.abs(np.linalg.det(B))
        weight_scaled = weight * detJ

        print "Triangular element has vertices:"
        print coordinates[:3]

        local_S = np.zeros((6, 6))
        local_Hx = np.zeros((3, 6))
        local_Hy = local_Hx.copy()
        local_Gx = np.zeros((6, 3))
        local_Gy = local_Gx.copy()

        for i in xrange(6):
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

        # scatter the local matrices
        for i in xrange(6):
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

    # Apply the boundary conditions
    # Assuming that velocity (u, v) = (0, 0) on the boundary
    # That implies that the mean pressure is zero over the domain
    # The mean pressure is taken into account by allowing the space for
    # (u, v) to not necessarily be divergence-free. It integrates to zeros
    # against the pressure space.

    S = S[interior_nodes - 1]
    S = S[:, interior_nodes - 1]
    Hx = Hx[:, interior_nodes - 1]
    Hy = Hy[:, interior_nodes - 1]
    Gx = Gx[interior_nodes - 1]
    Gy = Gy[interior_nodes - 1]

    k = len(interior_nodes)
    A = np.bmat([[S, np.zeros((k, k)), -Gx],
                 [np.zeros((k, k)), S, -Gy],
                 [Hx, Hy, np.zeros((m, m))]])

    A += np.eye(len(A)) * eps

    np.savetxt('./files/A.txt', A)
    np.savetxt('./files/S.txt', S)
    np.savetxt('./files/Gx.txt', Gx)
    np.savetxt('./files/Gy.txt', Gy)
    np.savetxt('./files/Hx.txt', Hx)
    np.savetxt('./files/Hy.txt', Hy)
