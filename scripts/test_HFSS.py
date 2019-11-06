# -*- coding: utf-8 -*-
"""
Created on Wed Oct 30 15:39:44 2019

@author: antho
"""

from PythonModeler import PythonModeler

PM = PythonModeler('hfss')

PM.set_variable('track', '42um')
PM.set_variable('bond', '100um')

# A body is the shared coordinate system of several elements.
# If one wants to insert a chip into a 3D cavity, one should create
# One body for the cavity and one body for the chip. 
chip = PM.body('coord_chip1', "chip_1", [['1mm','0mm','0mm'], [1,0,0], [0,1,0]])
        
# When drawing 2D object, the code assume we draw them in the plane z=0
info = 'con1', ['1mm','1mm'], 90

#def draw_connector(self, name, iTrack, iGap, iBondLength, iSlope=1, pcbTrack=None, pcbGap=None, tr_line=True):
#chip.draw_connector(*info, PM.track, '25um', PM.bond*2+'100um')

#draw_quarter_circle(self, name, coor, ori, size ?):
info1 = 'quarter1', ['0mm','0mm'], 90
#chip.draw_quarter_circle(*info1, 0.01)

##rect_corner_2D(self, pos, size, z=0, **kwargs):
#chip.rect_corner_2D([0,0],[5,5], name='rectangle1')
##rect_center_2D (self, pos, size, z=0, **kwargs):
#chip.rect_center_2D([0,0],[5,5], name='rectangle2')

#chip.mesh_zone(['0mm','0mm'], 90,[10,20],0.1)

#chip.cutout('cutout', [0,0],90, [1,2])
chip.draw_T('T', 0.3, 0.2)
