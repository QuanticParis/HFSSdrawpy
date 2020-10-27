import os
import numpy as np


from HFSSdrawpy import Modeler, Body, Entity
import drawpylib.cpw_elements as elt
from HFSSdrawpy.utils import parse_entry, Vector
from drawpylib.parameters import TRACK, GAP, RLC, MESH, MASK, DEFAULT, ELEC, \
    eps

'''
testing script for waveport assignment
'''
modeler = 'hfss'
pm = Modeler(modeler)

pm.is_litho = False
pm.is_hfss = True

chip = Body(pm, 'chip')

# Drawing

track = '42um'
gap = '25um'

with chip([0, 0], [-1, 0]):
    port, = elt.draw_end_cable(chip, track, gap, typeEnd='waveport', 
                               lower='0.3mm', upper='1mm', name='test')