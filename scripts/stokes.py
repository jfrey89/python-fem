#!/usr/bin/env python

from __future__ import division
import ap.mesh.meshes as m
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


def quad_basis(index, x, y):
    coefficients = np.array([[1, -3, -3, 4, 2, 2],
                             [0, -1, 0, 0, 2, 0],
                             [0, 0, -1, 0, 0, 2],
                             [0, 4, 0, -4, -4, 0],
                             [0, 0, 4, -4, 0, -4],
                             [0, 0, 0, 4, 0, 0]])
    return np.dot(coefficients[index],
                  np.array([1, x, y, x * y, x * x, y * y]))


def quad_basis_gradient(index, x, y):
    coefficients = np.array([[-3, 4, 4], [-3, 4, 4],
                             [-1, 4, 0], [0, 0, 0],
                             [0, 0, 0], [-1, 0, 4],
                             [4, -8, -4], [0, -4, 0],
                             [0, 0, -4], [4, -4, -8],
                             [0, 0, 4], [0, 4, 0]])

    return np.dot(coefficients[index:index + 2], np.array([1, x, y]))


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

    stiffness_velocity_x = np.zeros((len(nodes), len(nodes)))
    stiffness_velocity_y = np.zeros((len(nodes), len(nodes)))
    stiffness_pressure = np.zeros((len(nodes), len(nodes)))

    for element in elements:

        coordinates = get_coordinates(element, nodes)
        corners = coordinates[:3]
        J = calculate_jacobian(corners)

        for i in xrange(6):
            for j in xrange(i, 6):

                grad_basis = lambda k, x, y: \
                    np.dot(J.T, quad_basis_gradient(k, x, y))

                grad_dot_grad = lambda x, y: \
                    np.dot(grad_basis(i, x, y), grad_basis(j, x, y)) / \
                    np.abs(np.linalg.det(J))

                I_x = (grad_dot_grad(0.5, 0.5) + grad_dot_grad(0.0, 0.5) +
                       grad_dot_grad(0.5, 0.0)) / 6

                I_y = (grad_dot_grad(0.5, 0.5) + grad_dot_grad(0.0, 0.5) +
                       grad_dot_grad(0.5, 0.0)) / 6

                # scatter the local operator back
                stiffness_velocity_x[element[i] - 1, element[j] - 1] = \
                    stiffness_velocity_x[element[i] - 1, element[j] - 1] + I_x

                stiffness_velocity_y[element[i] - 1, element[j] - 1] = \
                    stiffness_velocity_y[element[i] - 1, element[j] - 1] + I_y

                stiffness_velocity_x[element[j] - 1, element[i] - 1] = \
                    stiffness_velocity_x[element[i] - 1, element[j] - 1]

                stiffness_velocity_y[element[j] - 1, element[i] - 1] = \
                    stiffness_velocity_y[element[i] - 1, element[j] - 1]

                pass
            pass
