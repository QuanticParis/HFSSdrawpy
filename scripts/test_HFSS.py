# -*- coding: utf-8 -*-
"""
Created on Wed Oct 30 15:39:44 2019

@author: antho
"""

from PythonModeler import PythonModeler

PM = PythonModeler('hfss')

# A body is the shared coordinate system of several elements.
# If one wants to insert a chip into a 3D cavity, one should create
# One body for the cavity and one body for the chip. 
chip = PM.body('chip1', [['1mm','0mm','0mm'], [1,0,0], [0,1,0]])

# When drawing 2D object, the code assume we draw them in the plane z=0
info = 'con1', ['1mm','1mm'], 90
chip.draw_connector(*info, '42um', '25um', '100um')



