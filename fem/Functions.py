#!/usr/bin/env/python

from __future__ import division
import numpy as np


class Quadratic_Basis_Function(object):

    def __init__(self):
        self._coefficients = np.array([[1, -3, -3, 4, 2, 2],
                                       [0, -1, 0, 0, 2, 0],
                                       [0, 0, -1, 0, 0, 2],
                                       [0, 4, 0, -4, -4, 0],
                                       [0, 0, 4, -4, 0, -4],
                                       [0, 0, 0, 4, 0, 0]])
        self._gradient_coefficients = np.array([[[-3, 4, 4], [-3, 4, 4]],
                                                [[-1, 4, 0], [0, 0, 0]],
                                                [[0, 0, 0], [-1, 0, 4]],
                                                [[4, -8, -4], [0, -4, 0]],
                                                [[0, 0, -4], [4, -4, -8]],
                                                [[0, 0, 4], [0, 4, 0]]])

    def grad(self, i, x, y):
        return np.dot(self._gradient_coefficients[i], np.array([1, x, y]))

    def diff(self, i, x, y, derivative_variable):
        if derivative_variable == 'x':
            return self.grad(i, x, y)[0]
        else:
            return self.grad(i, x, y)[1]

    def __call__(self, i, x, y):
        return np.dot(self._coefficients[i],
                      np.array([1, x, y, x * y, x * y, y * y]))


class Linear_Basis_Function(object):
    def __init__(self):
        self._coefficients = np.array([[1, -1, -1], [0, 1, 0], [0, 0, 1]])

        self._derivative_coefficients = np.array([[-1, -1], [1, 0], [0, 1]])

    def diff(self, i, x, y, derivative_variable):
        if derivative_variable == 'x':
            return self._derivative_coefficients[i, 0]
        else:
            return self._derivative_coefficients[i, 1]

    def __call__(self, i, x, y):
        return np.dot(self._coefficients[i], np.array([1, x, y]))
