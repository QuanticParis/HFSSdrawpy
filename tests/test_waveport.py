import os
import numpy as np


from HFSSdrawpy import Modeler, Body, Entity
from HFSSdrawpy.utils import parse_entry, Vector
from HFSSdrawpy.libraries.base_elements import *

'''
testing script for waveport assignment
'''
modeler = 'hfss'
pm = Modeler(modeler)

pm.is_litho = False
pm.is_hfss = True

chip_body = Body(pm, 'chip')

# Drawing

track = pm.set_variable('0.32mm')
sub_h = pm.set_variable('0.43mm')
cover_H = pm.set_variable('0.57mm')
MSL_length = pm.set_variable('1mm')
width = pm.set_variable('3mm')

# define substrate + MSL + GND
chip_subs = box(chip_body,
                [-width/2, 0, -sub_h],
                [width, MSL_length, sub_h],
                name="chip_subs")
chip_subs.assign_material("sapphire")

MSL = rect(chip_body,
           [-track/2, 0, 0],
           [track, MSL_length, 0], 
           name="MSL")
MSL.assign_perfect_E('_perfE')

GND = rect(chip_body,
           [-width/2, 0, 0],
           [width, MSL_length, 0], 
           name="GND")
GND.assign_perfect_E('_perfE')

# define vacuum
cover = box(chip_body,
                [-width/2, 0, 0],
                [width, MSL_length, cover_H],
                name="air_top")
cover.assign_material("vacuum")

# define ports
port1 = rect(chip_body,
              [-width/2, 0, -sub_h], 
              [width, 0, cover_H+sub_h],
              name="1")
port1.assign_waveport(Nmodes=1)
port1.assign_terminal_auto(GND)

port2 = rect(chip_body,
              [-width/2, MSL_length, -sub_h], 
              [width, 0, cover_H+sub_h],
              name="2")
port2.assign_waveport(Nmodes=1)
port2.assign_terminal_auto(GND)