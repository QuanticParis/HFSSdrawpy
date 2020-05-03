# -*- coding: utf-8 -*-
"""
Created on Tue Oct 29 14:05:00 2019

@author: wcs
"""
import os

from HFSSdrawpy import Modeler, Body
from HFSSdrawpy.parameters import layer_TRACK, layer_GAP, layer_RLC
import HFSSdrawpy.libraries.example_elements as elt
# import HFSSdrawpy.libraries.base_elements as base

pm = Modeler('gds')

relative = pm.set_variable('1mm')

chip1 = Body(pm, 'chip1')
chip2 = Body(pm, 'chip2', rel_coor=[[0,0,relative], [0,1,0], [-1,0,0]])
chip3 = Body(pm, 'chip3', rel_coor=[[0,relative,0], [1,0,0], [0,1,0]],
                ref_name='chip2')

track = pm.set_variable('20um')
gap = pm.set_variable('10um')

track_big = pm.set_variable('25um')
gap_big = pm.set_variable('15um')

track_middle = pm.set_variable('22.5um')
gap_middle = pm.set_variable('12.5um')

offset = pm.set_variable('-50um')

# chip2
in_port0 = elt.create_port(chip2, 'in0', [track, track+2*gap]) # default is the widths of track and gap

with chip2(['2.0mm', '0.0mm'], [1, 0]):
    in_port2 = elt.create_port(chip2, 'in2', [track, track+2*gap]) # default is the widths of track and gap

bond_length, bond_slope, pcb_track, pcb_gap = '200um', 0.5, '300um', '200um'

with chip2(['0.5mm', '0.5mm'], [0, 1]):
    elt.draw_connector(chip2, 'in_flux_top', track, gap, bond_length, pcb_track,
                   pcb_gap, 0.5)

    with chip2(['1.5mm', '-1.0mm'], [0, 1]):
        in_port2 = elt.create_port(chip2, 'test', [track, track+2*gap])

chip2.draw_cable('cable2', 'in_flux_top', 'test', is_bond=True, fillet='100um',
              reverse_adaptor=False, to_meander=[0, 0, 0], meander_length=0)#, is_mesh=True)

# 3D
chip2.box_corner_3D([0, 0, 0], ['3mm', '3mm', '-1mm'], name='box', material='silicon')
ground_plane = chip2.rect_corner_2D([0, 0], ['3mm', '3mm'], name='ground_plane', layer=layer_TRACK)

ground_plane.subtract(chip2.entities[layer_GAP])
ground_plane.unite(chip2.entities[layer_TRACK])
ground_plane.assign_perfect_E()

# generate gds file
pm.generate_gds(os.path.join(os.getcwd(), 'gds_files'), 'cable_test')

