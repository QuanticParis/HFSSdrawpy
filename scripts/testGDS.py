# -*- coding: utf-8 -*-
"""
Created on Mon Nov  4 17:54:30 2019

@author: antho
"""

import gdspy

import numpy as np

poly = gdspy.Polygon([[0,0],[0,1],[1,1],[1,0]])
points = poly.polygons
polySet = gdspy.PolygonSet(points)
print(polySet.polygons)
poly.fillet([np.array([0,0.1,0,0])])
print(poly.polygons)
polySet.fillet([0.1])
#print(polySet.polygons[0])

