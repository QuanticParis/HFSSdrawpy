# -*- coding: utf-8 -*-
"""
Created on Mon Nov  4 17:54:30 2019

@author: antho
"""

from PythonModeler import PythonModeler
import gdspy

global current_library
current_library = gdspy.GdsLibrary()
PM = PythonModeler('gds')
chip = PM.body('coord_chip1', "chip_4", [['1mm','0mm','0mm'], [1,0,0], [0,1,0]])
chip.polyline([(1,1),(2,2),(1,2)], 2, name = "polyline")
chip.polyline([(1,3),(1,2),(1,2)], 2, name = "polyline2")
chip.interface.generate_gds('test_gds.gds', chip.interface.cell)
chip.reset_cell()