#!/usr/bin/env python

from __future__ import division
import ap.mesh.meshes as m
import fem.Functions as fn
import numpy as np
from scipy.sparse.linalg import spsolve
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
    mesh_file = 'unit-square_h-0.5.mesh'
    # load mesh
    domain = m.mesh_factory(root_dir + mesh_file)
    elements = domain.elements
    nodes_2_coords = domain.nodes
    nodes = np.unique(elements)
    pnodes = np.unique(elements.T[:][:3])

    reynolds = 1.0
    weight = 1 / 6.0
    gauss_pts = np.array([[0.5, 0.0],
                          [0.0, 0.5],
                          [0.5, 0.5]])

    S = np.zeros((nodes.max(), nodes.max()))
    Gx = np.zeros((nodes.max(), nodes.max()))
    Gy = Gx.copy()
    Hx = np.zeros((nodes.max(), nodes.max()))
    Hy = Hx.copy()

    f = np.zeros(3 * len(S))
    f1 = np.zeros(len(S))
    f2 = f1.copy()

    quad_basis = fn.Quadratic_Basis_Function()
    lin_basis = fn.Linear_Basis_Function()

    for element in elements:
        coordinates = get_coordinates(element, nodes_2_coords)
        B = calculate_B(coordinates[:3])
        detJ = 1 / np.abs(np.linalg.det(B))
        weight_scaled = weight * detJ

        for i in xrange(6):
            for j in xrange(6):
                node_i, node_j = element[i] - 1, element[j] - 1

                # (grad u_h, grad v)
                for point in gauss_pts:
                    x_g, y_g = point

                    S[node_i, node_j] += (weight_scaled / reynolds) * np.dot(
                        np.dot(quad_basis.grad(i, x_g, y_g), B),
                        np.dot(B.T, quad_basis.grad(j, x_g, y_g)))

                    f2[node_j] += -weight_scaled * quad_basis(j, x_g, y_g)

        # -(p_h, div v)
        for i in xrange(3):
            for j in xrange(6):
                node_i, node_j = element[i] - 1, element[j] - 1

                for point in gauss_pts:
                    x_g, y_g = point
                    Gx[node_i, node_j] += weight_scaled * (
                        lin_basis(i, x_g, y_g) *
                        np.dot(quad_basis.grad(j, x_g, y_g), B[0]))

                    Gy[node_i, node_j] += weight_scaled * (
                        lin_basis(i, x_g, y_g) *
                        np.dot(quad_basis.grad(j, x_g, y_g), B[1]))

        for i in xrange(6):
            for j in xrange(3):
                node_i, node_j = element[i] - 1, element[j] - 1

                for point in gauss_pts:
                    x_g, y_g = point
                    Hx[node_i, node_j] += weight_scaled * (
                        np.dot(quad_basis.grad(i, x_g, y_g), B[0]) *
                        lin_basis(j, x_g, y_g))

                    Hy[node_i, node_j] += weight_scaled * (
                        np.dot(quad_basis.grad(i, x_g, y_g), B[1]) *
                        lin_basis(j, x_g, y_g))

    # nuke the boundary nodes
    boundary = domain.boundary_nodes['land']
    for node in boundary:
        node = node - 1
        S[node] = np.zeros(np.shape(S[node]))
        S.T[node] = np.zeros(np.shape(S.T[node]))
        #Gx[node] = np.zeros(np.shape(Gx[node]))
        Gx.T[node] = np.zeros(np.shape(Gx.T[node]))
        #Gy[node] = np.zeros(np.shape(Gy[node]))
        Gy.T[node] = np.zeros(np.shape(Gy.T[node]))
        Hx[node] = np.zeros(np.shape(Hx[node]))
        #Hx.T[node] = np.zeros(np.shape(Hx.T[node]))
        Hy[node] = np.zeros(np.shape(Hy[node]))
        #Hy.T[node] = np.zeros(np.shape(Hy.T[node]))

    Z = np.zeros(np.shape(S))
    A = np.bmat([[S, Z, Gx], [Z, S, Gy], [Hx, Hy, Z]])
    # make it nonsingular, nearby matrix
    A += np.eye(np.shape(A)[0], np.shape(A)[1]) * eps

    #f[:len(f1)] += f1[:]
    #f[len(f1):len(f1) + len(f2)] += f2[:]

    #c = spsolve(A, f)

    np.savetxt('./files/A.txt', A)
    np.savetxt('./files/S.txt', S)
    np.savetxt('./files/Gx.txt', Gx)
    np.savetxt('./files/Gy.txt', Gy)
    np.savetxt('./files/Hx.txt', Hx)
    np.savetxt('./files/Hy.txt', Hy)
    #np.savetxt('./files/f.txt', f)
    #np.savetxt('./files/c.txt', c)
