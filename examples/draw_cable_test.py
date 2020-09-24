# -*- coding: utf-8 -*-
"""
Created on Tue Oct 29 14:05:00 2019

@author: wcs
"""
import os

from HFSSdrawpy import Modeler, Body
from HFSSdrawpy.parameters import TRACK, GAP
import HFSSdrawpy.libraries.example_elements as elt
# import HFSSdrawpy.libraries.base_elements as base

pm = Modeler('gds')

relative = pm.set_variable('1mm')

chip1 = Body(pm, 'chip1', rel_coor=[[0, 0, relative], [1, 0, 0], [0, 1, 0]])
chip2 = Body(pm, 'chip2', rel_coor=[[0,0,relative], [0,1,0], [-1,0,0]])
chip3 = Body(pm, 'chip3', rel_coor=[[0,relative,0], [1,0,0], [0,1,0]],
                ref_name='chip2')

track = pm.set_variable('20um')
gap = pm.set_variable('10um', name='gap')

track_big = pm.set_variable('25um')
gap_big = pm.set_variable('15um')

track_middle = pm.set_variable('22.5um')
gap_middle = pm.set_variable('12.5um')

offset = pm.set_variable('-50um')

# chip2
with chip2(['0.5mm', '0.5mm'], [1, 0]):
    port0, = elt.create_port(chip2, [track, track+2*gap], name='port0') # default is the widths of track and gap

    with chip2(['1.0mm', '0.1mm'], [-1, 0]):
        port1, = elt.create_port(chip2, name='port1') # default is the widths of track and gap
    with chip2(['2.0mm', '0.1mm'], [-1, 0]):
        port2, = elt.create_port(chip2, [track, track+2*gap], name='port2') # default is the widths of track and gap

chip2.draw_cable(port0, port1, port2, is_bond=False, fillet='200um',
                 reverse_adaptor=False, to_meander=[0, 0, 0],
                 meander_length=0)

# 3D
chip2.box([0, 0, 0], ['3mm', '3mm', '-1mm'], material='silicon')

ground_plane = chip2.rect([0, 0], ['3mm', '3mm'], layer=TRACK)

ground_plane.subtract(chip2.entities[GAP])
ground_plane.unite(chip2.entities[TRACK])
ground_plane.assign_perfect_E()

# generate gds file
pm.generate_gds(os.path.join(os.getcwd(), 'gds_files'), 'cable_test')

