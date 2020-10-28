import os
import numpy as np


from HFSSdrawpy import Modeler, Body, Entity
from HFSSdrawpy.utils import parse_entry, Vector
from HFSSdrawpy.libraries.base_elements import *

'''
testing script for waveport assignment

generates a shielded microstrip line on sapphire including ground plane 
for terminal assignment

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
<<<<<<< HEAD
<<<<<<< HEAD
MSL_length = pm.set_variable('1mm')
width = pm.set_variable('3mm')

# define substrate + MSL + GND
chip_subs = chip_body.box([-width/2, 0, -sub_h],
                [width, MSL_length, sub_h],
                name="chip_subs")
chip_subs.assign_material("sapphire")

MSL = chip_body.rect([-track/2, 0, 0],
           [track, MSL_length, 0], 
           name="MSL")
MSL.assign_perfect_E('_perfE')

GND = chip_body.rect([-width/2, 0, -sub_h],
           [width, MSL_length, 0], 
           name="GND")
GND.assign_perfect_E('_perfE')

# define vacuum
cover = chip_body.box([-width/2, 0, 0],
                [width, MSL_length, cover_H],
                name="air_top")
cover.assign_material("vacuum")

# define ports
port1 = chip_body.rect([-width/2, 0, -sub_h], 
              [width, 0, cover_H+sub_h],
              name="1")
port1.assign_waveport(Nmodes=1)
port1.assign_terminal_auto(GND)

<<<<<<< HEAD
port2 = chip_body.rect([-width/2, MSL_length, -sub_h], 
              [width, 0, cover_H+sub_h],
              name="2")
port2.assign_waveport(Nmodes=1)
port2.assign_terminal_auto(GND)
=======
rect = chip.rect([0, -(gap+track/2)*10, -lower], 
                [0, (gap+track/2)*20, upper+lower],
                name='_waveport')
rect.assign_waveport(Nmodes=2, DoRenorm=True, RenormValue="50ohm", DoDeembed=True, DeembedDist="2mm")
>>>>>>> 633aff3 (update assign_waveport)
=======
length = pm.set_variable('1mm')
=======
MSL_length = pm.set_variable('1mm')
>>>>>>> e40649f (Update test_waveport.py)
width = pm.set_variable('3mm')

# define substrate + MSL
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
port1.assign_waveport(Nmodes=1, DoRenorm=True, RenormValue="50ohm")
<<<<<<< HEAD
<<<<<<< HEAD
port1 = chip.rect([-width/2, length, -sub_h], 
                  [width, length, cover_H],
                  name='2')
port1.assign_waveport(Nmodes=1, DoRenorm=True, RenormValue="50ohm")
>>>>>>> 62d226b (Update test_waveport.py)
=======
port1 = chip_body.rect([-width/2, length, -sub_h], 
                       [width, length, cover_H],
                       name='2')
port1.assign_waveport(Nmodes=1, DoRenorm=True, RenormValue="50ohm")
>>>>>>> 5655412 (Update test_waveport.py)
=======

port2 = rect(chip_body,
              [-width/2, MSL_length, -sub_h], 
              [width, 0, cover_H+sub_h],
<<<<<<< HEAD
              name="_2")
port2.assign_waveport(Nmodes=1, DoRenorm=True, RenormValue="50ohm")
>>>>>>> e40649f (Update test_waveport.py)
=======
              name="2")
port2.assign_waveport(Nmodes=1, DoRenorm=True, RenormValue="50ohm")
>>>>>>> 22869bd (update waveport)
