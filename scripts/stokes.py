#!/usr/bin/env python

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


def calculate_jacobian(coordinates):
    # unpack coordinates
    x1, y1 = coordinates[0]
    x2, y2 = coordinates[1]
    x3, y3 = coordinates[2]

    return np.linalg.inv(np.array([[x1 - x3, x2 - x3],
                                   [y1 - y3, y2 - y3]]))

if __name__ == '__main__':
    root_dir = './files/'
    mesh_file = 'unit-square_h-0.3.mesh'
    # load mesh
    domain = m.mesh_factory(root_dir + mesh_file)
    elements = domain.elements
    nodes = domain.nodes

    A = np.zeros((len(nodes), len(nodes)))
    B = np.zeros((len(nodes), len(nodes)))
    C = np.zeros((len(nodes), len(nodes)))

    for element in elements:
        coordinates = get_coordinates(element, nodes)
        corners = coordinates[:4]
        J = calculate_jacobian(corners)

        for i in xrange(6):
            for j in xrange(i, 6):
                A(i, j) = A(i, j) +

