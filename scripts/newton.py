#!/usr/bin/env python

from __future__ import division
import ap.mesh.meshes as m
import fem.Functions as fn
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spl
from fem.stokes import f1 as f1
from fem.stokes import f2 as f2
from fem.stokes import get_coordinates as get_coordinates
from fem.stokes import calculate_B


if __name__ == "__main__":
    reynolds = 1e0
    mesh_file = 'convergence/01.mesh'
    root_dir = '.files/'

    UVP = np.loadtxt(root_dir + 'UVP.txt')

    print "Loading and parsing the mesh..."
    domain = m.mesh_factory(root_dir + mesh_file)
    elements = domain.elements
    nodes = np.unique(elements)
    pnodes = np.unique(elements.T[:3])
    interior_nodes = domain.interior_nodes
    coordinates = domain.nodes
    m, n = len(nodes), len(pnodes)

    weight = 1 / 6.0
    gauss_pts = np.array([[0.5, 0.0],
                          [0.0, 0.5],
                          [0.5, 0.5]])

    quad_basis = fn.Quadratic_Basis_Function()
    lin_basis = fn.Linear_Basis_Function()

    # Allocate the Jacobian matrix
    J = sp.lil_matrix((2 * n + m, 2 * n + m))

    for element in elements:
        element_coords = get_coordinates(element, domain.nodes)
        B = calculate_B(element_coords[:3])
        detJ = 1 / np.abs(np.linalg.det(B))
        weight_scaled = weight * detJ

        # get the coefficients from the last solution
        u_c =



        pass
