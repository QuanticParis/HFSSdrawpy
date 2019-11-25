# -*- coding: utf-8 -*-
"""
Created on Fri Nov 22 17:18:25 2019

@author: Zaki
"""

import PythonModeler
import traceback
import hfss
from Vector import Vector

PythonModeler.Port.reset()
hfss.ModelEntity.reset()
PM = PythonModeler.PythonMdlr('hfss')
PM.set_variable('track', '42um')
PM.set_variable('bond', '100um')
PM.set_variable('pad_spacing', '0.12mm')


chip1 = PM.body('chip2', 'Global')
chip1.set_current_coor(pos = ['0mm', '10mm','0mm'], ori=[0,1])
P1 = chip1.port('port1', Vector(['10mm',0]), Vector([0,1]), '0.2mm', '0.1mm')
chip1.set_current_coor(pos = ['0mm', '10mm','0mm'], ori=[0,1])
P2 = chip1.port('port2', Vector(['10mm','10mm']), Vector([0,1]), '0.1mm','0.05mm')
SL_PTH = chip1._connect_JJ('jojo', 'port1', 'port2', "0.02mm")

#%%
import PythonModeler
import traceback
import hfss
from Vector import Vector

PythonModeler.Port.reset()
hfss.ModelEntity.reset()
PM = PythonModeler.PythonMdlr('hfss')
PM.set_variable('track', '42um')
PM.set_variable('bond', '100um')
PM.set_variable('pad_spacing', '0.12mm')

chip1 = PM.body('chip2', 'Global')
chip1.set_current_coor(pos = ['0mm', '0mm','0mm'], ori=[0,1])
P1 = chip1.port('port1', Vector(['10mm','0mm']), Vector([1,0]), '0.2mm', '0.1mm')
P2 = chip1.port('port2', Vector(['10mm','10mm']), Vector([0,1]), '0.1mm','0.05mm')
chip1.draw_cable("cable", "port1", "port2")
P3 = chip1.port('port3', Vector(['10mm','0mm']), Vector([1,0]), '0.2mm', '0.1mm')

chip1.cable_starter('CBSTRT', 'port3')
#%%
import PythonModeler
import traceback
import hfss
from Vector import Vector

PythonModeler.Port.reset()
hfss.ModelEntity.reset()
PM = PythonModeler.PythonMdlr('hfss')
PM.set_variable('track', '42um')
PM.set_variable('bond', '100um')
PM.set_variable('pad_spacing', '0.12mm')

chip1 = PM.body('chip2', 'Global')
P1 = chip1.port('port1', Vector(['10mm','0mm']), Vector([1,0]), '0.2mm', '0.1mm')
chip1.set_current_coor(pos = ['0mm', '0mm','0mm'], ori=[0,1])
#chip1.draw_capa_inline('CAPA_IN_LINE', '20mm', '10mm', '5mm', '1mm', n_pad=5, iTrack_capa='5mm', iGap_capa=None, premesh=True, tight=False)
chip1.draw_capa_interdigitated('CAPA_INTERDIGITATED', '20mm', '10mm', ['2mm','2mm'], ['20mm','20mm'], 10, '0.1mm')

#%%
import PythonModeler
import traceback
import hfss

PythonModeler.Port.reset()
hfss.ModelEntity.reset()

PM = PythonModeler.PythonMdlr('hfss')
PM.set_variable('track', '42um')
PM.set_variable('bond', '100um')
PM.set_variable('cylinder_height', '30mm')
PM.set_variable('antenna_height', '10mm')
PM.set_variable('cutout_x', '1.47mm')
PM.set_variable('cutout_y', '0.75mm')

PM.set_variable('cylinder_radius', '4mm')


chip1 = PM.body('chip1', "chip_1", [['0mm',PM.cylinder_radius,PM.antenna_height], [1,0,0], [0,1,0]])
chip1.set_current_coor(pos = ['0mm', '0mm','0mm'], ori=[0,1])
chip1.insert_transmon("insert", [PM.cutout_x, PM.cutout_y, '0.5mm'],'0.3mm',['0.2mm', '0.2mm'],'0.1mm', '0.35mm','0.42mm' ,'0.25um')

chip2 = PM.body('chip2', 'Global')
chip2.set_current_coor(pos = ['0mm', '0mm','0mm'], ori=[1,0])
cavity = chip2.cavity_3D_simple('cavity', PM.cylinder_radius,PM.cylinder_height, '3mm', PM.antenna_height, '1mm', '1mm')
#cavity = chip1.cavity_3D_simple('cavity', '3mm', '5mm', '0.5mm', '2.5mm', '2mm', '1mm')

#%%
from importlib import reload
import PythonModeler
import hfss
reload(PythonModeler)

PythonModeler.Port.reset()
hfss.ModelEntity.reset()

PM = PythonModeler.PythonMdlr('hfss')
PM.set_variable('track', '42um')
PM.set_variable('bond', '100um')
PM.set_variable('cylinder_height', '30mm')
PM.set_variable('antenna_height', '10mm')
PM.set_variable('cutout_x', '1.47mm')
PM.set_variable('cutout_y', '0.75mm')

PM.set_variable('cylinder_radius', '4mm')


chip1 = PM.body('chip1', "chip_1", [['0mm',PM.cylinder_radius,PM.antenna_height], [1,0,0], [0,1,0]])
chip1.set_current_coor(pos = ['0mm', '0mm','0mm'], ori=[1,0])
sapphire = chip1.box_center([0,0,'-0.25mm'], [PM.cutout_y, PM.cutout_x, '0.5mm'], "3D", name="sapphire")
chip1.make_material(sapphire, "\"sapphire\"")
chip1.rect_center_2D([PM.cutout_y/10,0,0], [PM.cutout_y/10, PM.cutout_x, 0], layer="TRACK", name="tracked1")
chip1.rect_center_2D([-PM.cutout_y/10,0,0], [PM.cutout_y/10, PM.cutout_x, 0], layer="TRACK", name="tracked2")
chip2 = PM.body('chip2', 'Global')
chip2.set_current_coor(pos = ['0mm', '0mm','0mm'], ori=[1,0])
cavity = chip2.cavity_3D_simple('cavity', PM.cylinder_radius,PM.cylinder_height, '3mm', PM.antenna_height, '1mm', '1mm')
#cavity = chip1.cavity_3D_simple('cavity', '3mm', '5mm', '0.5mm', '2.5mm', '2mm', '1mm')


