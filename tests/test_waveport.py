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

track = pm.set_variable('19.5mm')
gap = pm.set_variable('19.5mm')
upper = pm.set_variable('19.5mm')
lower = pm.set_variable('19.5mm')

rect = chip.rect([0, -(gap+track/2)*10, -lower], 
                [0, (gap+track/2)*20, upper+lower],
                name='_waveport')
rect.assign_waveport(Nmodes=2, DoRenorm=True, RenormValue="50ohm", DoDeembed=True, DeembedDist="2mm")