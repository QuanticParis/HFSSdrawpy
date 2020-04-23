# -*- coding: utf-8 -*-
"""
Created on Tue Oct 29 14:05:00 2019

@author: wcs
"""
import os

from drawpy_scripts.python_modeler import PythonModeler, ModelEntity
from drawpy_scripts.variable_string import parse_entry, Vector
from drawpy_scripts.python_modeler import layer_TRACK, layer_GAP, layer_RLC

pm = PythonModeler('hfss')

relative = pm.set_variable('1mm')

chip1 = pm.body('chip1')
chip2 = pm.body('chip2', rel_coor=[[0,0,relative], [0,1,0], [-1,0,0]])
chip3 = pm.body('chip3', rel_coor=[[0,relative,0], [1,0,0], [0,1,0]],
                ref_name='chip2')

track = pm.set_variable('20um')
gap = pm.set_variable('10um')

track_big = pm.set_variable('25um')
gap_big = pm.set_variable('15um')

track_middle = pm.set_variable('22.5um')
gap_middle = pm.set_variable('12.5um')

offset = pm.set_variable('-50um')

# chip1
# port_design = {'widths':[track_big, track_big+2*gap_big, track_big], 'subnames':['track', 'gap', 'other'], 'offsets':[0, 0, offset], 'layers':[1, 3, 4]}
# chip1.set_current_coor(['-0.5mm', '0.5mm'], [1, 0])
# in_port = chip1.create_port('in', **port_design)

# chip1.set_current_coor(['0.5mm', '-0.5mm'], [-1, 0])
# middle_port = chip1.create_port('middle')#, [track_middle, track_middle+2*gap_middle, track_middle], subnames=['track', 'gap', 'other'], offsets=[0, 0, -offset], layers=[1, 3, 4])

# chip1.set_current_coor(['2mm', '2mm'], [-1, 0])
# out_port = chip1.create_port('out', [track, track+2*gap, track], subnames=['track', 'gap', 'other'], offsets=[0, 0, -offset], layers=[1, 3, 4])

# chip1.draw_cable('cable', 'in', 'middle', 'out', is_bond=True, fillet='100um', reverse_adaptor=False, to_meander=[0, 0, 0], meander_length=0)#, is_mesh=True)

# chip2
chip2.set_current_coor(['0.0mm', '0.0mm'], [1, 0])
in_port2 = chip2.create_port('in2', [track, track+2*gap]) # default is the widths of track and gap

bond_length, bond_slope, pcb_track, pcb_gap = '200um', 0.5, '300um', '200um'
chip2.set_current_coor(['0.5mm', '0.5mm'], [0, 1])
chip2.draw_connector('in_flux_top', track, gap, bond_length, pcb_track,
                     pcb_gap, 0.5)

chip2.draw_cable('cable2', 'in_flux_top', 'in2', is_bond=True, fillet='100um',
                 reverse_adaptor=False, to_meander=[0, 0, 0], meander_length=0)#, is_mesh=True)


# generate gds file
pm.generate_gds(os.getcwd(), 'test')

# ModelEntity.print_instances()
