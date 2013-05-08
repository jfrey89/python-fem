#!/usr/bin/env python

from __future__ import division
import ap.mesh.meshes as m
import fem.Functions as fn
import numpy as np


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
    root_dir = './files/'
    mesh_file = 'unit-square_h-0.5.mesh'
    # load mesh
    domain = m.mesh_factory(root_dir + mesh_file)
    elements = domain.elements
    nodes_2_coords = domain.nodes
    nodes = np.unique(elements)
    pnodes = np.unique(elements.T[:][:3])

    weight = 1 / 6.0
    gauss_pts = np.array([[0.5, 0.0],
                          [0.0, 0.5],
                          [0.5, 0.5]])

    A = np.zeros((nodes.max(), nodes.max()))
    Bx = np.zeros((pnodes.max(), nodes.max()))
    By = Bx.copy()
    Cx = np.zeros((nodes.max(), pnodes.max()))
    Cy = Cx.copy()

    quad_basis = fn.Quadratic_Basis_Function()
    lin_basis = fn.Linear_Basis_Function()

    for element in elements:
        coordinates = get_coordinates(element, nodes_2_coords)
        B = calculate_B(coordinates[:3])
        detJ = 1 / np.abs(np.linalg.det(B))
        weight_scaled = weight * detJ

        for i in xrange(6):
            for j in xrange(6):
                # global index for scattering
                node_i, node_j = element[i] - 1, element[j] - 1

                for point in gauss_pts:
                    x_g, y_g = point
                    A[node_i, node_j] += weight_scaled * np.dot(
                        np.dot(quad_basis.grad(i, x_g, y_g), B),
                        np.dot(B.T, quad_basis.grad(j, x_g, y_g)))

        for i in xrange(3):
            for j in xrange(6):
                node_i, node_j = element[i] - 1, element[j] - 1

                for point in gauss_pts:
                    x_g, y_g = point
                    Bx[node_i, node_j] += weight_scaled * (
                        np.dot(lin_basis.grad(i), B[0]) *
                        quad_basis(j, x_g, y_g))

                    By[node_i, node_j] += weight_scaled * (
                        np.dot(lin_basis.grad(i), B[1]) *
                        quad_basis(j, x_g, y_g))

        for i in xrange(6):
            for j in xrange(3):
                node_i, node_j = element[i] - 1, element[j] - 1

                for point in gauss_pts:
                    x_g, y_g = point
                    Cx[node_i, node_j] += weight_scaled * (
                        np.dot(quad_basis.grad(i, x_g, y_g), B[0]) *
                        lin_basis(j, x_g, y_g))

                    Cy[node_i, node_j] += weight_scaled * (
                        np.dot(quad_basis.grad(i, x_g, y_g), B[1]) *
                        lin_basis(j, x_g, y_g))

    np.savetxt('./files/A.txt', A)
    np.savetxt('./files/Bx.txt', Bx)
    np.savetxt('./files/By.txt', By)
    np.savetxt('./files/Cx.txt', Cx)
    np.savetxt('./files/Cy.txt', Cy)
