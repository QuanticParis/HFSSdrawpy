# -*- coding: utf-8 -*-
"""
Created on Wed Mar 11 09:07:00 2020

@author: wcs
"""

from scripts.hfss import get_active_project, release
from scripts.designer import Circuit, KeyElt

'''
testing script for consolidation of junction functions in designer

changes:
    -
'''

project = get_active_project()
design = project.get_active_design()
modeler = design.modeler
modeler.set_units('mm')
modeler.delete_all_objects()

c = Circuit(design, modeler)

KeyElt.is_litho = False
KeyElt.is_hfss = True

# Drawing

c.set_variable('xx', '100um')  # vertical spacing
c.set_variable('yy', '50um')  # horizontal spacing


def pos(i, j):
    return [(i + 0.5)*c.xx, (j + 0.5)*c.yy]

# Finger junctions


gap = '20um'
track = '15um'
Jwidth = '1um'
Jlength = '0.8um'
Jinduc = '1nH'
left = '100um'
right = left
down = '100um'
up_left = '50um'
dn_left = '50um'
up_right = '50um'
dn_right = '50um'

c.key_elt('temp', pos(0, 0), [0, 1])

c.temp.draw_pHmon(track, gap, left, right, down, up_left, dn_left,
                  up_right, dn_right, Jwidth, Jlength, Jinduc)
release()
