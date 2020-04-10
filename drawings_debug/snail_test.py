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


pad_spacing = '150um'
loop_width = '20um'
loop_length = '30um'
N = 3
length_island = '4um'
width_bridge_left = '0.8um'
width_bridge_right = '0.4um'
width_jct_left = '4um'
width_jct_right = '1um'
c.key_elt('temp', pos(0, 0), [0, 1])
c.temp.draw_test_snails(['10um', '10um'], pad_spacing, loop_width, loop_length,
                        N, length_island,
                        width_bridge_left, width_bridge_right,
                        width_jct_left, width_jct_right,
                        n_left=4, n_right=1,
                        spacing_bridge_left='3um', spacing_bridge_right='0',
                        yoffset=1, litho='elec')

release()
