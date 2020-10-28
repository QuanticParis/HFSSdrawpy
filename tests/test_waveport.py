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
length = pm.set_variable('1mm')
width = pm.set_variable('3mm')

# define substrate + MSL
chip_subs = chip_body.box(chip_body,
                          [-width/2, 0, -sub_h],
                          [width, length, sub_h],
                          name="chip_subs")
chip_subs.assign_material("sapphire")

MSL = chip_body.rect(chip_body,
                     [-track/2, 0, 0],
                     [track, length, 0], 
                     name="MSL")
MSL.assign_perfect_E('_perfE')

# define vacuum
chip_subs = chip_body.box(chip_body,
                          [-width/2, 0, 0],
                          [width, length, cover_H],
                          name="chip_subs")
chip_subs.assign_material("sapphire")

# define ports
port1 = chip_body.rect([-width/2, 0, -sub_h], 
                       [width, 0, cover_H],
                       name='1')
port1.assign_waveport(Nmodes=1, DoRenorm=True, RenormValue="50ohm")
port1 = chip_body.rect([-width/2, length, -sub_h], 
                       [width, length, cover_H],
                       name='2')
port1.assign_waveport(Nmodes=1, DoRenorm=True, RenormValue="50ohm")