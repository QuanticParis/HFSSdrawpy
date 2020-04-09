# -*- coding: utf-8 -*-
"""
Created on Thu Nov 14 13:46:42 2019

@author: Zaki
"""
import numpy
import variable_string
# from .variable_string import VariableString, parse_entry

class Vector(list):
    def __init__(self, vec, vec_y=None):
        if vec_y is not None:
            vec = [vec, vec_y]
        super().__init__(variable_string.parse_entry(vec))

    def check(self, elt):
        return isinstance(elt, (list, tuple, numpy.ndarray))

    def check_nb(self, nb):
        return isinstance(nb, float) or isinstance(nb, int) or isinstance(nb, variable_string.VariableString)

    def __add__(self, other):
        if self.check(other):
            return Vector([self[0]+other[0], self[1]+other[1]])
        else:
            raise TypeError('Could not perform add operation')

    def __radd__(self, other):
        return self + other

    def __sub__(self, other):
        if self.check(other):
            return Vector([self[0]-other[0], self[1]-other[1]])
        else:
            raise TypeError('Could not perform sub operation')

    def __neg__(self):
        return Vector([-self[0], -self[1]])

    def __rsub__(self, other):
        return -self + other

    def __mul__(self, other):
        if self.check(other):
            return Vector([self[0]*other[0], self[1]*other[1]])
        elif self.check_nb(other):
            return Vector([other*self[0], other*self[1]])
        else:
            raise TypeError('Could not perform mul operation')

    def __rmul__(self, other):
        return self * other

    def __truediv__(self, other):
        if self.check(other):
            return Vector([self[0]/other[0], self[1]/other[1]])
        elif self.check_nb(other):
            return Vector([self[0]/other, self[1]/other])
        else:
            raise TypeError('Could not perform div operation')

    def __rtruediv__(self, other):
        if self.check(other):
            return Vector([other[0]/self[0], other[1]/self[1]])
        elif self.check_nb(other):
            return Vector([other/self[0], other/self[1]])
        else:
            raise TypeError('Could not perform rdiv operation')

    def dot(self, other):
        if self.check(other):
            return self[0]*other[0]+self[1]*other[1]
        else:
            raise TypeError('Could not perform dot operation')

    def cross(self, other):
        if self.check(other):
            return self[0]*other[1]-self[1]*other[0]
        else:
            raise TypeError('Could not perform dot operation')

    def norm(self):
        return (self[0]**2+self[1]**2)**0.5

    def abs(self):
        return Vector([abs(self[0]), abs(self[1])])

    def unit(self):
        norm = self.norm()
        return Vector([self[0]/norm, self[1]/norm])

    def orth(self):
        return Vector([-self[1], self[0]])

    def rot(self, other):
        '''
        Inputs:
        -------
        other: vector

        Returns
        -------
        vector: rotation around z of self by an angle given by other w.r.t to x
        '''
        if self.check(other):
            unitOther = Vector(other).unit()
            return Vector([self.dot(unitOther.refx()), self.dot(unitOther.orth().refy())])
        else:
            raise TypeError('Could not perform rdiv operation')

    def px(self):
        return Vector([self[0], 0])

    def py(self):
        return Vector([0, self[1]])

    def refx(self, offset=0):
        return Vector([self[0], -self[1]+2*offset])

    def refy(self, offset=0):
        return Vector([-self[0]+2*offset, self[1]])
