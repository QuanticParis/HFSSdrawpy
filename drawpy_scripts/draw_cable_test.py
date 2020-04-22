# -*- coding: utf-8 -*-
"""
Created on Tue Oct 29 14:05:00 2019

@author: wcs
"""

from drawpy_scripts.python_modeler import PythonModeler, ModelEntity
from drawpy_scripts.variable_string import parse_entry, Vector
import os
import gdspy

pm = PythonModeler('gds')

# body(self, body_name, coor_name='Global', coor_sys=None):
chip1 = pm.body('chip1')  # create a chip

relative = pm.set_variable('20um')
chip2 = pm.body('chip2', rel_coor=[[0,0,relative], [1,0,0], [0,1,0]])

# use litho for optical
# use litho + mask for ground plane mesh
# use litho + ebeam for ebeam
# use nothing for hfss

litho = False
mask = False
ebeam = False

if litho:
    is_bond = False
    pm.is_overdev = False
    pm.is_litho = True
    # KeyElt.overdev = parse_entry('0.5um')
    if mask:
        pm.is_mask = True
        pm.gap_mask = parse_entry('20um')
else:
    mask = False
    is_bond = True

track = pm.set_variable('20um')
gap = pm.set_variable('10um')

track_big = pm.set_variable('25um')
gap_big = pm.set_variable('15um')

track_middle = pm.set_variable('22.5um')
gap_middle = pm.set_variable('12.5um')

offset = pm.set_variable('-50um')

CPW=0
if CPW:
    chip1.set_current_coor(['0.1mm', '0.1mm'], [1, 0])
    # chip1.create_port(name, widths, subnames=None, layers=layer_Default, offsets=0)
    in_port = chip1.create_port('in', [track_big, track_big+2*gap_big], subnames=['track', 'gap'])

    chip1.set_current_coor(['0.5mm', '0.5mm'], [0, 1])
    middle_port = chip1.create_port('middle', [track, track+2*gap], subnames=['track', 'gap'])

    chip1.set_current_coor(['0.7mm', '0.3mm'], [-1, 0])
    out_port = chip1.create_port('out', [track, track+2*gap], subnames=['track', 'gap'])
else:
    port_design = {'widths':[track_big, track_big+2*gap_big, track_big], 'subnames':['track', 'gap', 'other'], 'offsets':[0, 0, offset], 'layers':[1, 3, 4]}
    chip1.set_current_coor(['0.0mm', '0.0mm'], [1, 0])
    in_port = chip1.create_port('in', **port_design)

    chip1.set_current_coor(['0.5mm', '-0.5mm'], [-1, 0])
    middle_port = chip1.create_port('middle')#, [track_middle, track_middle+2*gap_middle, track_middle], subnames=['track', 'gap', 'other'], offsets=[0, 0, -offset], layers=[1, 3, 4])

    chip1.set_current_coor(['2mm', '2mm'], [-1, 0])
    out_port = chip1.create_port('out', [track, track+2*gap, track], subnames=['track', 'gap', 'other'], offsets=[0, 0, -offset], layers=[1, 3, 4])


if 1:
    chip1.draw_cable('cable', 'in', 'middle', 'out', is_bond=True, fillet='100um', reverse_adaptor=False, to_meander=[0, 0, 0], meander_length=0)#, is_mesh=True)

chip2.set_current_coor(['0.0mm', '0.0mm'], [1, 0])
in_port2 = chip2.create_port('in2', **port_design)

ModelEntity.print_instances()
pm.generate_gds(os.getcwd(), 'test')

#%%

#### Test junctions
#
#if litho and 1:
#    pm.set_variable('x_D', '4.6mm')
#    pm.set_variable('y_D', '3.6mm')
#    pm.set_variable('dose_pad', '100um')
#    pm.set_variable('dose_sep', '100um')
#    pm.set_variable('dose_cor', '10um')
#    pm.key_elt('dose_test', [pm.x_D, pm.y_D], [1,0])
#    dose_mat = [2,6]
#
#    pm.key_elt('align_dose_test', pm.dose_test.pos + 0.5*Vector([2*pm.dose_pad*dose_mat[0] + 3*pm.dose_sep*dose_mat[0] - pm.dose_cor*(2*dose_mat[0]-1), pm.dose_pad*dose_mat[1] + pm.dose_sep*(dose_mat[1]+1)]), pm.dose_test.ori)
#    pm.align_dose_test.draw_alignement_marks('20um', ['0.3mm', '0.7mm'])
#    pm.key_elt('align_chip', [0.5*pm.chip_width, 0.5*pm.chip_length], pm.dose_test.ori)
#    pm.align_chip.draw_alignement_marks('80um', ['3.5mm', '3.2mm'])
#    pm.dose_test.draw_dose_test_Nb([pm.dose_pad, pm.dose_pad], pm.dose_sep, dose_mat, correction=pm.dose_cor)
#
#    if ebeam:
#        xpos = pm.x_D - pm.dose_pad - 1.5*pm.dose_sep
#        ypos = pm.y_D + 0.5*pm.dose_pad + pm.dose_sep
#        for ii in range(dose_mat[0]):
#            xpos = xpos + 3*pm.dose_sep + 2*pm.dose_pad
#            pm.key_elt('temp', [xpos, ypos], [1,0])
#            pm.temp.draw_dose_test_junction(['5um', '5um'],
#                                            pm.dose_sep-pm.dose_cor, '0.172um', '400nm',
#                                            alternate_width=False, version=1)
#            pm.key_elt('temp', [xpos, ypos + pm.dose_sep + pm.dose_pad], [1,0])
#            pm.temp.draw_dose_test_junction(['5um','5um'],
#                                            pm.dose_sep-pm.dose_cor, '0.147um', '400nm',
#                                            alternate_width=False)
#            pm.key_elt('temp', [xpos, ypos + 2*pm.dose_sep + 2*pm.dose_pad], [1,0])
#            pm.temp.draw_dose_test_junction(['5um', '5um'],
#                                            pm.dose_sep-pm.dose_cor, '0.852um', '400nm',
#                                            alternate_width=False)
#            pm.key_elt('temp', [xpos, ypos + 3*pm.dose_sep + 3*pm.dose_pad], [1,0])
#            pm.temp.draw_dose_test_junction(['5um', '5um'],
#                                            pm.dose_sep-pm.dose_cor, '0.852um', '400nm',
#                                            n_bridge=26, spacing_bridge='1.047um', alternate_width=False)
#            pm.key_elt('temp', [xpos, ypos + 4*pm.dose_sep + 4*pm.dose_pad], [1,0])
#            pm.temp.draw_dose_test_junction(['5um', '5um'],
#                                            pm.dose_sep-pm.dose_cor, '0.852um', '400nm',
#                                            n_bridge=52, spacing_bridge='1.047um', alternate_width=False)
#            pm.key_elt('temp', [xpos, ypos + 5*pm.dose_sep + 5*pm.dose_pad], [0,1])
#            pm.temp.draw_meander_array(['5um', '5um'], '90um', '0.852um', '400nm',
#                                            ['60um', '80um'], [52, 26],
#                                            spacing_bridge='1.047um')
#
#### Other stuff
#%%
#
#bottom, top, left, right  = pm.get_extent(margin='500um')
#width = right-left
#height = top-bottom
#if not litho:
#    pm.draw_box('pcb', [left, bottom, -pm.chip_thickness], [width, height, -pm.pcb_thickness], "Rogers TMM 10i (tm)")
#    pm.draw_box('chip', [left, bottom, 0], [width, height, -pm.chip_thickness], 'silicon')
#    pm.draw_box('box', [left, bottom, 0], [width, height, pm.vaccuum_thickness], 'vacuum')
#    pm.draw_rect('ground_plane', [left, bottom], [width, height])
#
#if litho:
#    pm.draw_rect('ground_plane', [0, 0], [pm.chip_width, pm.chip_length])
#    pm.draw_rect('negatif', [0, 0], [pm.chip_width, pm.chip_length])
#
#pm.assign_perfE(pm.ground_plane)
#
#if mask:
#    maskObject = pm.unite(pm.maskObjects, 'mask')
#
#for obj in pm.trackObjects:
#    pm.assign_perfE(obj)
#
#if 0:
#    gapObject = pm.unite(pm.gapObjects)
#    pm.subtract(pm.ground_plane, [gapObject])
#
#if litho and 0:
#    pm.subtract(pm.negatif, [pm.ground_plane]+pm.trackObjects)
#
#release()