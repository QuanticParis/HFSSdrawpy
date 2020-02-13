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

import scripts
import gdspy

PM = PythonMdlr('gds')
PM.set_variable('track', '42um')
PM.set_variable('bond', '100um')
PM.set_variable('pad_spacing', '0.12mm')

#2 draw a cylinder of vacuum with perfect_E boundaries in the center of the global coordinate system
chip2= PM.body('chip2', 'Global')
chip2.set_current_coor(pos = ['0mm', '0mm','0mm'], ori=[0,1])

#chip2.rect_corner_2D([0,0,1],[1,1,1], name="ee", layer='Z')
right_quarter_up2 = chip2.polyline_2D([[1,1],[0.5,0.5],[1,0.5]], name="rect", layer=1)
#right_quarter_up1 = chip2.rect_center_2D([2,2,2],[0.5,0.5], name="rect2", layer=1)
right_quarter_up1 = chip2.draw_quarter_circle('right_quarter_up1', 4, 1)

tes = chip2.polyline_2D([[-0.01,0],[0.01,0]], name="toflex", layer=1)
tes2 = chip2._sweep_along_path(tes, right_quarter_up2)

chip2.unite([right_quarter_up1, right_quarter_up2], name='new_name')
#print("objets finaux", hfss.ModelEntity.instances_to_move)
print("final", chip2.interface.cell.polygons)
chip2.generate_gds("quarter_circle.gds")
gdspy.LayoutViewer(library=gdspy.current_library, pattern={'default': 8},background='#FFFFFF')


#%%

#import scripts
import gdspy

PM = PythonMdlr('gds')
PM.set_variable('track', '42um')
PM.set_variable('bond', '100um')
PM.set_variable('pad_spacing', '0.12mm')

#2 draw a cylinder of vacuum with perfect_E boundaries in the center of the global coordinate system
chip2= PM.body('chip2', 'Global')
chip2.set_current_coor(pos = ['0mm', '0mm','0mm'], ori=[0,1])

#chip2.rect_corner_2D([0,0,1],[1,1,1], name="ee", layer='Z')
#rect = chip2.rect_corner_2D([0,0,0], [1,1,0], name='nom', layer=2)
#chip2._fillet(0.1, [0,1], rect)
#print("objets finaux", hfss.ModelEntity.instances_to_move)
PM.set_variable('track', '42um')
PM.set_variable('gap', '25um')
track_right = '0mm'
track_left = '0mm'
PM.set_variable('Lj', '12nH')
PM.set_variable('trm_junction_width', '10um')
PM.set_variable('gap_mem', '50um')


chip2.set_current_coor([0, 0], [-1,0])
chip2.draw_ZR_transmon('trm', ['1.47mm', '0.75mm'], '0.12mm', ['0.5mm', '0.5mm'],
                        track_right, PM.gap, '0.30mm', '30um', '0mm', PM.trm_junction_width, '50um', PM.Lj,
                        pad_size_left=['0.5mm', '0.5mm'], track_left=track_left,
                        gap_left=PM.gap_mem, length_left='0.2mm', spacing_left='-30um', 
                        short_left='0um', fillet=None)


print("eee", chip2.interface.gds_cells)
#chip2.generate_gds("quarter_circle.gds")
gdspy.LayoutViewer(library=gdspy.current_library, pattern={'default': 8},background='#FFFFFF')
