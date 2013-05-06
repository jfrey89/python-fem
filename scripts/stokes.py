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
    

def quadratic_basis(index, x, y):
    if index == 0:
        coefficients = np.array([1, -3, -3, 4, 2, 2])
    elif index == 1:
        coefficients = np.array([0, -1, 0, 0, 2, 0])
    elif index == 2:
        coefficients = np.array([0, 0, -1, 0, 0, 2])
    elif index == 3:
        coefficients = np.array([0, 4, 0, -4, -4, 0])
    elif index == 4:
        coefficients = np.array([0, 0, 4, -4, 0, -4])
    elif index == 5:
        coefficients = np.array([0, 0, 0, 4, 0, 0])
        
    return np.dot(np.array([1, x, y, x * y, x * x, y * y]), coefficients)


def quadratic_basis_gradient(index, x, y):
    if index == 0:
        gradient_coefficients = np.array([[-3, -3], [4, 4], [4, 4]])
    elif index == 1:
        gradient_coefficients = np.array([[-1, 0], [4, 0], [4, 4]])
    elif index == 2:
        gradient_coefficients = np.array([[0, -1], [0, 0], [0, 4]])
    elif index == 3:
        gradient_coefficients = np.array([[4, 0], [-8, 4], [4, 0]])
    elif index == 4:
        gradient_coefficients = np.array([[0, 4], [0, -4], [-4, -8]])
    elif index == 5:
        gradient_coefficients = np.array([[0, 0], [0, 4], [4, 0]])
        
    return np.dot(np.array([1, x, y]), gradient_coefficients)

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

    stiffness_velocity_x = np.zeros((len(nodes), len(nodes)))
    stiffness_velocity_y = np.zeros((len(nodes), len(nodes)))
    stiffness_pressure = np.zeros((len(nodes), len(nodes)))
    
    
    for element in elements:
        
        coordinates = get_coordinates(element, nodes)
        corners = coordinates[:3]
        J = calculate_jacobian(corners)
        
        for i in xrange(6):
            for j in xrange(i, 6):
                
                quad_grad_basis = lambda k, x, y: \
                np.dot(J.T, quadratic_basis_gradient(k, x, y))
                
                dot_product = lambda x, y: \
                np.dot(quad_grad_basis(i, x, y), quad_grad_basis(j, x, y)) * \
                np.linalg.det(J)
                
                I_x = (dot_product(0.5, 0.5) + dot_product(0.0, 0.5) + \
                dot_product(0.5, 0.0)) / 6
                
                I_y = (dot_product(0.5, 0.5) + dot_product(0.0, 0.5) + \
                dot_product(0.5, 0.0)) / 6
                
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