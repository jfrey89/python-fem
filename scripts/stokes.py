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


def calculate_jacobian(coordinates):
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
    nodes = domain.nodes

    gauss_points = np.array([[0.5, 0.5],
                             [0.0, 0.5],
                             [0.5, 0.0]])

    stiff_vx = np.zeros((len(nodes), len(nodes)))
    stiff_vy = np.zeros((len(nodes), len(nodes)))
    stiff_p = np.zeros((len(nodes), len(nodes)))

    quad_basis = fn.Quadratic_Basis_Function()
    lin_basis = fn.Linear_Basis_Function()

    for element in elements:

        coordinates = get_coordinates(element, nodes)
        corners = coordinates[:3]
        J = calculate_jacobian(corners)
        detJ = np.abs(np.linalg.det(J))

        for i in xrange(6):
            for j in xrange(i, 6):

                gauss_quad = 0.0
                node_i = element[i] - 1
                node_j = element[j] - 1

                for point in gauss_points:
                    gauss_quad = gauss_quad + np.dot(
                        np.dot(J.T, quad_basis.grad(i, point[0], point[1])),
                        np.dot(J.T, quad_basis.grad(j, point[0], point[1])))

                gauss_quad = gauss_quad / (6.0 * detJ)

                stiff_vx[node_i, node_j] = stiff_vx[node_i, node_j] + \
                    gauss_quad
                stiff_vy[node_i, node_j] = stiff_vy[node_i, node_j] + \
                    gauss_quad

                #if j < 3:
                    #pres_cont_x = lambda x, y: \
                        #grad_basis(i, x, y)[0] * lin_basis(j, x, y) / detJ

                    #stiff_vx[node_i, node_j] += \
                        #(pres_cont_x(0.5, 0.5) + pres_cont_x(0.0, 0.5) +
                         #pres_cont_x(0.5, 0.0)) / (6 * detJ)

                    #pres_cont_y = lambda x, y: \
                        #grad_basis(i, x, y)[1] * lin_basis(j, x, y) / detJ

                    #stiff_vy[node_i, node_j] += \
                        #(pres_cont_y(0.5, 0.5) + pres_cont_y(0.0, 0.5) +
                         #pres_cont_y(0.5, 0.0)) / (6 * detJ)

                stiff_vx[node_j, node_i] = stiff_vx[node_i, node_j]
                stiff_vy[node_j, node_i] = stiff_vy[node_i, node_j]

                pass
            pass
