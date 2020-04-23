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


iTrack = '10um'
iGap = '5um'
array_room = '100um'
array_offset = '0um'
iTrackPump = '5um'
iGapPump = '10um'

snail_dict = {'loop_width': '10um', 'loop_length': '20um', 'N': 3,
              'length_island': '4um', 'width_bridge_left': '0.8um',
              'width_bridge_right': '0.4um', 'width_jct_left': '4um',
              'width_jct_right': '1um', 'n_left': 1, 'n_right': 1,
              'spacing_bridge_left': '3um', 'spacing_bridge_right': '0um',
              'yoffset': 0, 'litho': 'elec', 'iInduct': '0nH'}

c.key_elt('temp', pos(0, 0), [0, 1])

c.temp.draw_snails_flux_line(iTrack, iGap, array_room, array_offset,
                             iTrackPump, iGapPump, snail_dict,
                             iTrackSnail=None, fillet=None, typePump='down',
                             doublePump=False)

release()
