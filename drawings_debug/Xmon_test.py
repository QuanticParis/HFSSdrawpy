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
KeyElt.is_hfss = False

# Drawing

c.set_variable('xx', '100um')  # vertical spacing
c.set_variable('yy', '50um')  # horizontal spacing


def pos(i, j):
    return [(i + 0.5)*c.xx, (j + 0.5)*c.yy]

# Finger junctions


gap_trm = '6um'
track_trm = '5um'
length = '20um'
gap_cpl = ['0.5um', '1um', '2um']
gnd_cpl = ['2um', '3um', '4um']
gap = ['1um', '2um', '3um']
track = ['1um', '2um', '3um']
Jwidth = '1um'
Jlength = '0.8um'
Jinduc = '1nH'

c.key_elt('temp', pos(0, 0), [0, 1])

c.temp.draw_Xmon(gap_trm, track_trm, length, gap_cpl, gnd_cpl,
                 gap, track, Jwidth, Jlength, Jinduc)
release()
