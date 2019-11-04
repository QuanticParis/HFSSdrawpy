# -*- coding: utf-8 -*-
"""
Created on Mon Nov  4 17:54:30 2019

@author: antho
"""

from PythonModeler import PythonModeler

PM = PythonModeler('gds')
chip = PM.body()
chip.polyline([(1,1),(2,2),(1,2)], 2, name = "polyline")
chip.polyline([(1,3),(1,2),(1,2)], 2, name = "polyline2")
chip.generate_gds('test_gds.gds')