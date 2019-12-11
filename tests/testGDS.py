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

#%%

import PythonModeler
import traceback

PM = PythonModeler.PythonMdlr('gds')
PM.set_variable('track', '42um')
PM.set_variable('bond', '100um')
PM.set_variable('pad_spacing', '0.12mm')

#2 draw a cylinder of vacuum with perfect_E boundaries in the center of the global coordinate system
chip2, net2 = PM.body('chip2', 'Global')
chip2.set_current_coor(pos = ['0mm', '0mm','0mm'], ori=[0,1])
#cavity1 = chip2.cavity_3D_simple('cavity1', '3mm', '10mm', '0.5mm', '3mm')
#chip2.set_current_coor(pos = ['1mm', '0mm','0mm'], ori=[0,1])
#cavity2 = chip2.cavity_3D_simple('cavity2', '3mm', '10mm', '0.5mm', '3mm')
#
#a = chip2.rect_corner_2D([0,0],[0.5,0.5], name='rectangle1', layer ='layer1')
#b = chip2.rect_corner_2D([0,0.25],[0.7,0.7], name='rectangle2', layer ='layer1')
#chip2.unite([a,b])
#print("objets finaux", hfss.ModelEntity.dict_instances)
###
###3 Setup another body for the transmon and draw the transmon
#chip1, net1 = PM.body('chip1', "chip_1", [['0mm','3mm','2.5mm'], [0,0,-1], [0,-1,0]])
chip2.set_current_coor(pos = ['0mm', '0mm','0mm'], ori=[0,1])
right_quarter_up1 = chip2.draw_quarter_circle('right_quarter_up1', 'TRACK', 1)
#print("objets finaux", hfss.ModelEntity.instances_to_move)
