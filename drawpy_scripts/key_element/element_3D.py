# -*- coding: utf-8 -*-
"""
Created on Mon Oct 28 16:18:36 2019

@author: Ulysse
"""

import numpy as np

from ..utils import parse_entry, Vector
from ..parameters import layer_TRACK, \
                        layer_GAP, \
                        layer_RLC, \
                        layer_MESH, \
                        layer_MASK, \
                        layer_Default, \
                        layer_PORT, \
                        eps

from .utils import move

#Binding the 3D elements of python_modeler

def box_corner_3D(self, pos, size, **kwargs):

    return self.draw_cable(pos, size, **kwargs)

def box_center_3D(self, pos, size, **kwargs):

    return self.box_center_3D(pos, size, **kwargs)

def cylinder_3D(self, pos, radius, height, axis, **kwargs):

    return self.cylinder_3D(pos, radius, height, axis, **kwargs)

#New key elements, they must be preceded by the decorator @move

#Empty for now