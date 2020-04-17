# -*- coding: utf-8 -*-
"""
Created on Tue Oct 29 14:05:00 2019

@author: wcs
"""

from drawpy_scripts.python_modeler import PythonModeler
from drawpy_scripts.variable_string import parse_entry, Vector
import os
import gdspy

PM = PythonModeler('hfss')

# body(self, body_name, coor_name='Global', coor_sys=None):
chip1 = PM.body('chip1', 'Global')  # create a chip

# use litho for optical
# use litho + mask for ground plane mesh
# use litho + ebeam for ebeam
# use nothing for hfss

litho = False
mask = False
ebeam = False

if litho:
    is_bond = False
    PM.is_overdev = False
    PM.is_litho = True
    # KeyElt.overdev = parse_entry('0.5um')
    if mask:
        PM.is_mask = True
        PM.gap_mask = parse_entry('20um')
else:
    mask = False
    is_bond = True

track = PM.set_variable('20um')
gap = PM.set_variable('10um')

track_big = PM.set_variable('25um')
gap_big = PM.set_variable('15um')

track_middle = PM.set_variable('22.5um')
gap_middle = PM.set_variable('12.5um')

offset = PM.set_variable('-50um')

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

cwd  = os.getcwd()
gdspy.write_gds(os.path.join(cwd, 'test.gds'), unit=1.0, precision=1e-9)
#PM.interface.generate_gds('test_parity.gds')

#%%

#### Test junctions
#
#if litho and 1:
#    PM.set_variable('x_D', '4.6mm')
#    PM.set_variable('y_D', '3.6mm')
#    PM.set_variable('dose_pad', '100um')
#    PM.set_variable('dose_sep', '100um')
#    PM.set_variable('dose_cor', '10um')
#    PM.key_elt('dose_test', [PM.x_D, PM.y_D], [1,0])
#    dose_mat = [2,6]
#
#    PM.key_elt('align_dose_test', PM.dose_test.pos + 0.5*Vector([2*PM.dose_pad*dose_mat[0] + 3*PM.dose_sep*dose_mat[0] - PM.dose_cor*(2*dose_mat[0]-1), PM.dose_pad*dose_mat[1] + PM.dose_sep*(dose_mat[1]+1)]), PM.dose_test.ori)
#    PM.align_dose_test.draw_alignement_marks('20um', ['0.3mm', '0.7mm'])
#    PM.key_elt('align_chip', [0.5*PM.chip_width, 0.5*PM.chip_length], PM.dose_test.ori)
#    PM.align_chip.draw_alignement_marks('80um', ['3.5mm', '3.2mm'])
#    PM.dose_test.draw_dose_test_Nb([PM.dose_pad, PM.dose_pad], PM.dose_sep, dose_mat, correction=PM.dose_cor)
#
#    if ebeam:
#        xpos = PM.x_D - PM.dose_pad - 1.5*PM.dose_sep
#        ypos = PM.y_D + 0.5*PM.dose_pad + PM.dose_sep
#        for ii in range(dose_mat[0]):
#            xpos = xpos + 3*PM.dose_sep + 2*PM.dose_pad
#            PM.key_elt('temp', [xpos, ypos], [1,0])
#            PM.temp.draw_dose_test_junction(['5um', '5um'],
#                                            PM.dose_sep-PM.dose_cor, '0.172um', '400nm',
#                                            alternate_width=False, version=1)
#            PM.key_elt('temp', [xpos, ypos + PM.dose_sep + PM.dose_pad], [1,0])
#            PM.temp.draw_dose_test_junction(['5um','5um'],
#                                            PM.dose_sep-PM.dose_cor, '0.147um', '400nm',
#                                            alternate_width=False)
#            PM.key_elt('temp', [xpos, ypos + 2*PM.dose_sep + 2*PM.dose_pad], [1,0])
#            PM.temp.draw_dose_test_junction(['5um', '5um'],
#                                            PM.dose_sep-PM.dose_cor, '0.852um', '400nm',
#                                            alternate_width=False)
#            PM.key_elt('temp', [xpos, ypos + 3*PM.dose_sep + 3*PM.dose_pad], [1,0])
#            PM.temp.draw_dose_test_junction(['5um', '5um'],
#                                            PM.dose_sep-PM.dose_cor, '0.852um', '400nm',
#                                            n_bridge=26, spacing_bridge='1.047um', alternate_width=False)
#            PM.key_elt('temp', [xpos, ypos + 4*PM.dose_sep + 4*PM.dose_pad], [1,0])
#            PM.temp.draw_dose_test_junction(['5um', '5um'],
#                                            PM.dose_sep-PM.dose_cor, '0.852um', '400nm',
#                                            n_bridge=52, spacing_bridge='1.047um', alternate_width=False)
#            PM.key_elt('temp', [xpos, ypos + 5*PM.dose_sep + 5*PM.dose_pad], [0,1])
#            PM.temp.draw_meander_array(['5um', '5um'], '90um', '0.852um', '400nm',
#                                            ['60um', '80um'], [52, 26],
#                                            spacing_bridge='1.047um')
#
#### Other stuff
#%%
#
#bottom, top, left, right  = PM.get_extent(margin='500um')
#width = right-left
#height = top-bottom
#if not litho:
#    PM.draw_box('pcb', [left, bottom, -PM.chip_thickness], [width, height, -PM.pcb_thickness], "Rogers TMM 10i (tm)")
#    PM.draw_box('chip', [left, bottom, 0], [width, height, -PM.chip_thickness], 'silicon')
#    PM.draw_box('box', [left, bottom, 0], [width, height, PM.vaccuum_thickness], 'vacuum')
#    PM.draw_rect('ground_plane', [left, bottom], [width, height])
#
#if litho:
#    PM.draw_rect('ground_plane', [0, 0], [PM.chip_width, PM.chip_length])
#    PM.draw_rect('negatif', [0, 0], [PM.chip_width, PM.chip_length])
#
#PM.assign_perfE(PM.ground_plane)
#
#if mask:
#    maskObject = PM.unite(PM.maskObjects, 'mask')
#
#for obj in PM.trackObjects:
#    PM.assign_perfE(obj)
#
#if 0:
#    gapObject = PM.unite(PM.gapObjects)
#    PM.subtract(PM.ground_plane, [gapObject])
#
#if litho and 0:
#    PM.subtract(PM.negatif, [PM.ground_plane]+PM.trackObjects)
#
#release()