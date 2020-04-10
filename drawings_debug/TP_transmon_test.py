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


gap_tr = '20um'
big_capa_width = '100um'
big_capa_length = '100um'
small_capa_width = '100um'
small_capa_length = '100um'
track_big = '3um'
gap_big = '5um'
length_big = '1um'
short_big = '0um'

short_small = '0um'
track_small = '3um'
gap_small = '5um'
capa_add_small = '0um'
                        
Jwidth = '1um'
Jlength = '0.8um'
Jinduc = '1nH'


c.key_elt('temp', pos(0, 0), [0, 1])

c.temp.draw_TP_transmon(gap_tr, big_capa_width, big_capa_length,
                        small_capa_width, small_capa_length, track_big, 
                        gap_big, length_big, short_big, short_small, 
                        track_small, gap_small, capa_add_small, Jwidth, Jlength,
                        Jinduc, small_type='Left')

release()
