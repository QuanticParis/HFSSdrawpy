# -*- coding: utf-8 -*-
"""
Created on Wed Oct 30 15:39:44 2019

@author: antho
"""

from PythonModeler import PythonModeler
from CustomElement import Port
import traceback
from ConnectElement2 import ConnectElt2, Vector

PM = PythonModeler('hfss')

PM.set_variable('track', '42um')
PM.set_variable('bond', '100um')

# A body is the shared coordinate system of several elements.
# If one wants to insert a chip into a 3D cavity, one should create
# One body for the cavity and one body for the chip. 
chip = PM.body('coord_chip1', "chip_1", [['1mm','0mm','0mm'], [1,0,0], [0,1,0]])
connector = ConnectElt2(chip)      
# When drawing 2D object, the code assume we draw them in the plane z=0
info = 'con1', ['1mm','1mm'], 90

#def draw_connector(self, name, iTrack, iGap, iBondLength, iSlope=1, pcbTrack=None, pcbGap=None, tr_line=True):
#chip.draw_connector(*info, PM.track, '25um', PM.bond*2+'100um')

#draw_quarter_circle(self, name, coor, ori, size ?):
info1 = 'quarter1', ['0mm','0mm'], 90


# TEST BATCH 1
#try:
#    print("Quarter Cirle")
#    chip.draw_quarter_circle(*info1, 'track', 0.01)
#except Exception:
#    print("Quarter Cirle error")
#    traceback.print_exc()
#
#    
#try:
#    print("Rect Corner")
#    chip.rect_corner_2D([0,0],[5,5], name='rectangle1', layer ='layer1')
#except Exception:
#    print("Rectangle error")
#    traceback.print_exc()
#
#    
#try:
#    print("mesh_zone")
#
#    a = chip.rect_center_2D([0,0],[5,5], name='rectangle2', layer ='layer1')
#    chip.mesh_zone(a,0.1)
#
#except Exception:
#    print("Mesh Zone error")
#
#try:
#    print("cutout")
#    chip.cutout('cutout', [0,0],90, [1,2])
#except Exception:
#    print("Cutout error")
#
#try:
#    print("draw_T")
#    chip.draw_T('T', [0,0,0], 90, 0.3, 0.2)
#except Exception:
#    print("DrawT error")
#    traceback.print_exc()
#
#    
#try:
#    print("end_cable")
#    chip.draw_end_cable('end_cable', [1,1,1], 90, 5,2)
#except Exception:
#    print("Draw_end_cable error")
#    traceback.print_exc()
#    
#try:
#    print("draw_JJ")
#    chip.draw_JJ('JJ', 5, 2, 4, 2)
#except Exception:
#    print("Draw_JJ error")   
#    traceback.print_exc()

#TEST BATCH 2
print("connector")
chip.draw_connector(*info, PM.track, '25um', PM.bond*2+'100um')
print("rectangles")
chip.rect_corner_2D([0,0],[0.5,0.5], name='rectangle1', layer ='layer1')
chip.rect_corner_2D([0,1],[0.5,0.5], name='rectangle2', layer ='layer1')
P1 = Port('port1', [0.25,0.5], [0,1], 2, 3)
P2 = Port('port2', [0.25,1], [0,-1], 1, 2)
print("capa")
connector.draw_capa('capacite', 'port1', 'port2', 10, 5, 1)

#%% TEST BATCH 3
from PythonModeler import PythonModeler
from CustomElement import Port
from ConnectElement2 import ConnectElt2

connector = ConnectElt2(chip)      
PM = PythonModeler('hfss', connector)
PM.set_variable('track', '42um')
PM.set_variable('bond', '100um')
chip = PM.body('coord_chip1', "chip_1", [['1mm','0mm','0mm'], [1,0,0], [0,1,0]])
info = 'con1', ['1mm','1mm'], 90
chip.rect_corner_2D([0,0],[0.5,0.5], name='rectangle1', layer ='layer1')
chip.rect_corner_2D([0,1],[0.5,0.5], name='rectangle2', layer ='layer1')
P1 = Port('port1', Vector([0,0]), Vector([0,1]), '2mm', '3mm')
P2 = Port('port2', Vector([1,1]), Vector([0,1]), '1mm','2mm')
#SL_PTH = connector.find_slanted_path('slanded_path', 'port1', 'port2')
chip.draw_connector(*info, PM.track, '25um', PM.bond*2+'100um')

final_choice1 = PM.connector.find_path('path', 'port1', 'port2', 0.05, True, [0,2,0,0,0,0,0], 0.4, 0.1)
final_choice2 = PM.connector.find_path('path', 'port1', 'port2', 0.05, False, [0,2,0,0,0,0,0], 0.4, 0.1)
print("final_choice1",final_choice1)
print("final_choice1",final_choice2)

longueur2 = connector.length(final_choice2, 0, 3, 0.05)
longueur1 = connector.length(final_choice1, 0, 3, 0.05)

cable_starter(self, width = 'track', index=None, border=parse_entry('15um'))

#%% TEST BATCH 4
from PythonModeler import PythonModeler
from CustomElement import Port
import traceback
from ConnectElement2 import ConnectElt2

connector = ConnectElt2(chip)      
PM = PythonModeler('hfss', connector)
PM.set_variable('track', '42um')
PM.set_variable('bond', '100um')
chip = PM.body('coord_chip1', "chip_1", [['1mm','0mm','0mm'], [1,0,0], [0,1,0]])
info = 'con1', ['1mm','1mm'], 90

P1 = Port('port1', Vector([0,0]), Vector([0,1]), '30mm', '10mm')
P2 = Port('port2', Vector([1,1]), Vector([0,1]), '20mm','20mm')
#SL_PTH = connector._connect_JJ('jojo', 'port1', 'port2', 2)

#chip.draw_IBM_tansmon(['1.47', '0.75'],'0.12',['0.5', '0.5'],'30mm', '25mm','42mm' ,'25nH')


chip.draw_ZR_transmon(['1.47mm', '0.75mm'],'0.12mm',['0.5mm', '0.5mm'],'42um','25um', '0.30mm', '30um', '0mm', '0.01mm' ,'25nH', pad_size_left=['0.5mm','0.5mm'], track_left = '10um', gap_left='50um', length_left='0.2mm', spacing_left='50um', short_left='0um' ,fillet=True)

