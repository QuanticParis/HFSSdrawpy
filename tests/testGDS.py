# -*- coding: utf-8 -*-
"""
Created on Mon Nov  4 17:54:30 2019

@author: antho
"""

import gdspy
import os
import numpy as np

gdspy.current_library = gdspy.GdsLibrary()

# The GDSII file is called a library, which contains multiple cells.
lib = gdspy.GdsLibrary()

# Geometry must be placed in cells.
cell = lib.new_cell("FIRST")

### rect
if 0:
    rect = gdspy.Rectangle((0, 0), (2, 1))
    cell.add(rect)

###
if 1:
    # Path created with automatic bends of radius 5
    points = [[0, 0], [0, 10], [10, 10], [10, 20]]
    cable = gdspy.FlexPath(
        points,
        [1, 2, 1],
        offset=[-1, 3, 2],
        corners="circular bend",
        bend_radius=[5, 5],
        gdsii_path=False,
        layer=[0, 0, 0],
    )

    # Same path, generated with natural corners, for comparison
    # sp5 = gdspy.FlexPath(points, 1, layer=1, gdsii_path=True)

    # cell.add(sp4)
    # cell.add(sp4)

    polygons = cable.get_polygons()
    print("enumerate dict keys")
    for key in cable._polygon_dict.keys():
        print(key)
        print(len(cable._polygon_dict[key]))
        for elt in cable._polygon_dict[key]:
            print(elt)
    for ii, poly in enumerate(polygons):
        poly = gdspy.Polygon(poly)
        poly.layers = [ii]
        cell.add(poly)
    print("max_points")
    print(cable.max_points)


# Display all cells using the internal viewer.
display = 0
if display:
    gdspy.LayoutViewer()
else:
    cwd = os.getcwd()
    gdspy.write_gds(os.path.join(cwd, "test_path.gds"), unit=1.0e-6, precision=1e-8)

# print(polySet.polygons[0])

#%%

# import scripts
# import gdspy

# PM = PythonMdlr('gds')
# PM.set_variable('track', '42um')
# PM.set_variable('bond', '100um')
# PM.set_variable('pad_spacing', '0.12mm')

# #2 draw a cylinder of vacuum with perfect_E boundaries in the center of the global coordinate system
# chip2= PM.body('chip2', 'Global')
# chip2.set_current_coor(pos = ['0mm', '0mm','0mm'], ori=[0,1])

# #chip2.rect_corner_2D([0,0,1],[1,1,1], name="ee", layer='Z')
# right_quarter_up2 = chip2.polyline_2D([[1,1],[0.5,0.5],[1,0.5]], name="rect", layer=1)
# #right_quarter_up1 = chip2.rect_center_2D([2,2,2],[0.5,0.5], name="rect2", layer=1)
# right_quarter_up1 = chip2.draw_quarter_circle('right_quarter_up1', 4, 1)

# tes = chip2.polyline_2D([[-0.01,0],[0.01,0]], name="toflex", layer=1)
# tes2 = chip2._sweep_along_path(tes, right_quarter_up2)

# chip2.unite([right_quarter_up1, right_quarter_up2], name='new_name')
# #print("objets finaux", hfss.ModelEntity.instances_to_move)
# print("final", chip2.interface.cell.polygons)
# chip2.generate_gds("quarter_circle.gds")
# gdspy.LayoutViewer(library=gdspy.current_library, pattern={'default': 8},background='#FFFFFF')


# #%%

# #import scripts
# import gdspy

# PM = PythonMdlr('gds')
# PM.set_variable('track', '42um')
# PM.set_variable('bond', '100um')
# PM.set_variable('pad_spacing', '0.12mm')

# #2 draw a cylinder of vacuum with perfect_E boundaries in the center of the global coordinate system
# chip2= PM.body('chip2', 'Global')
# chip2.set_current_coor(pos = ['0mm', '0mm','0mm'], ori=[0,1])

# #chip2.rect_corner_2D([0,0,1],[1,1,1], name="ee", layer='Z')
# #rect = chip2.rect_corner_2D([0,0,0], [1,1,0], name='nom', layer=2)
# #chip2._fillet(0.1, [0,1], rect)
# #print("objets finaux", hfss.ModelEntity.instances_to_move)
# PM.set_variable('track', '42um')
# PM.set_variable('gap', '25um')
# track_right = '0mm'
# track_left = '0mm'
# PM.set_variable('Lj', '12nH')
# PM.set_variable('trm_junction_width', '10um')
# PM.set_variable('gap_mem', '50um')


# chip2.set_current_coor([0, 0], [-1,0])
# chip2.draw_ZR_transmon('trm', ['1.47mm', '0.75mm'], '0.12mm', ['0.5mm', '0.5mm'],
#                         track_right, PM.gap, '0.30mm', '30um', '0mm', PM.trm_junction_width, '50um', PM.Lj,
#                         pad_size_left=['0.5mm', '0.5mm'], track_left=track_left,
#                         gap_left=PM.gap_mem, length_left='0.2mm', spacing_left='-30um',
#                         short_left='0um', fillet=None)


# print("eee", chip2.interface.gds_cells)
# #chip2.generate_gds("quarter_circle.gds")
# gdspy.LayoutViewer(library=gdspy.current_library, pattern={'default': 8},background='#FFFFFF')
