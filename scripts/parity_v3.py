# -*- coding: utf-8 -*-
"""
Created on Tue Oct 29 14:05:00 2019

@author: wcs
"""

'''
Changes from v2 to v3:
 - trsmn junction width reduced from 0.237 to 0.147 um to reduce EJ from 21.9 to 13.6 GHz
 - small junction width reduced from 0.185 to 0.172 um to reduce EJ from 12.9 to 12.0 GHz
 - large junction width increased from 0.779 to 0.852 um to increase EJ from 45.7 to 50.0 GHz
 - half as many junctions in shunt inductance
 - both capacitances to parity circuit increased
 - interdigitated capacitances removed
'''  
  

import scripts
PM = scripts.PythonMdlr('gds')
from scripts.hfss import parse_entry as parse_entry
from scripts.Vector import Vector as Vector
chip1 = PM.body('chip1', 'Global')

# use litho for optical
# use litho + mask for ground plane mesh
# use litho + ebeam for ebeam
# use nothing for hfss

litho = True
mask = False
ebeam = True

if litho:
    is_bond = False
    PM.is_overdev = False
    PM.is_litho = True
    #KeyElt.overdev = parse_entry('0.5um')
    if mask:
        PM.is_mask = True
        PM.gap_mask = parse_entry('20um')
else:
    mask = False
    is_bond = True

qubit = True
memory = True
readout = True
drive_mem = True
parity = True
flux_lines = True

#### Standard dimensions

#TODO
chip_thickness = "12mm" # 6*PM.chip_thickness

PM.set_variable("chip_width", "8.67mm")
PM.set_variable("chip_length", "8.12mm")
PM.set_variable("chip_thickness", "280um")
PM.set_variable("pcb_thickness", "320um")
PM.set_variable('vaccuum_thickness', chip_thickness)

#### Connectors

con1, con2, con3, con4, con5 = '4.06mm', '3.01mm', '7.54mm', '3.01mm', '6.36mm'
PM.pcb_track, PM.pcb_gap = parse_entry('300um'), parse_entry('200um')
bond_length, bond_slope = '200um', '0.5'
PM.set_variable('track', '42um')
PM.set_variable('gap', '25um')
PM.set_variable('track_mem', '10um')
PM.set_variable('gap_mem', '50um')
PM.set_variable('track_flux', '50um')
PM.set_variable('min_track_flux', '4um')
PM.set_variable('gap_flux', '10um')
PM.set_variable('fillet', '200um')
#print([(k.name, k.ori) for k in PythonModeler.Port.dict_instances.values()])
 
chip1.set_current_coor([con2, PM.chip_length], [0, -1])
chip1.draw_connector(chip1.name+'in_flux_top', PM.track, PM.gap, bond_length, bond_slope)
chip1.set_current_coor([0, con1], [1, 0])
chip1.draw_connector(chip1.name+'in_readout', PM.track, PM.gap, bond_length, bond_slope)
chip1.set_current_coor([con5, 0], [0, 1])
chip1.draw_connector(chip1.name+'in_mem', PM.track, PM.gap, bond_length, bond_slope)
chip1.set_current_coor([con4, 0], [0, 1])
chip1.draw_connector(chip1.name+'in_unused', PM.track, PM.gap, bond_length, bond_slope)
chip1.set_current_coor([con3, PM.chip_length], [0, -1])
chip1.draw_connector(chip1.name+'in_other_flux', PM.track, PM.gap, bond_length, bond_slope)
#print([(k.name, k.ori) for k in PythonModeler.Port.dict_instances.values()])
#PM.set_variable('gap_capa_readout', '8um')
#chip1.set_current_coor(['0.9mm', con1], [1,0])
#chip1.draw_capa_inline('capa_readout', PM.track, PM.gap, '100um', PM.gap_capa_readout, n_pad=2)


#chip1.draw_cable('mem', chip1.name+'in_flux_topiOut', chip1.name+'in_unusediOut', is_bond=is_bond, fillet=PM.fillet, is_mesh=True)
#### Drawing
#
PM.set_variable('x_T', '4.685mm') # (con4 + con5)/2
PM.set_variable('y_T', '2.5mm')
#
if qubit:
    PM.set_variable('Lj', '12nH')
    PM.set_variable('trm_junction_width', '10um')
    track_right = '0mm'
    track_left = '0mm'
    if readout:
        track_right = PM.track
    if memory:
        track_left = PM.track_mem
        
    chip1.set_current_coor([PM.x_T, PM.y_T], [-1,0])
    chip1.draw_ZR_transmon('trm', ['1.47mm', '0.75mm'], '0.12mm', ['0.5mm', '0.5mm'],
                            track_right, PM.gap, '0.30mm', '30um', '0mm', PM.trm_junction_width, '50um', PM.Lj,
                            pad_size_left=['0.5mm', '0.5mm'], track_left=track_left,
                            gap_left=PM.gap_mem, length_left='0.2mm', spacing_left='-30um', 
                            short_left='0um', fillet=True)

    if litho:
        chip1.set_current_coor([PM.x_T, PM.y_T], [-1,0])
        chip1.draw_alignement_marks('align_trm', '20um', ['0.1mm', '0.43mm'])
        
        if ebeam:
            chip1.set_current_coor([PM.x_T, PM.y_T], [-1,0])
            chip1.draw_dose_test_junction('junction_trm', ['10um', PM.trm_junction_width], 
                                                    '50um', '0.147um', '400nm',
                                                    alternate_width=False)

if parity:
    PM.set_variable('cutout_depth', '0.26mm')
    PM.set_variable('cutout_width', '0.295mm')
    PM.set_variable('cap_gap', '4um')
    PM.set_variable('cap_depth', '0.0mm')
    PM.set_variable('ind_gap', '90um')
    PM.set_variable('jn_gap', '90um')
    PM.set_variable('buffer', '50um')
    PM.set_variable('par_track', '5um')
    PM.set_variable('x_I', con5-0.5*(-PM.track_mem+PM.par_track+PM.ind_gap)-PM.buffer)
    PM.set_variable('y_I', '2.9mm')
    chip1.set_current_coor([PM.x_I, PM.y_I], [0,1])
    #PM.key_elt('cap1', [PM.x_I - 0.5*PM.ind_gap - PM.buffer, PM.y_I - 0.5*PM.cutout_depth + 2*PM.buffer + 0.5*PM.cap_depth], [0,1])
    #PM.key_elt('cap2', [PM.x_I + 0.5*PM.ind_gap + PM.buffer, PM.y_I - 0.5*PM.cutout_depth + 2*PM.buffer + 0.5*PM.cap_depth], [0,1])
        
    chip1.draw_selfparity('parity', [PM.cutout_depth, PM.cutout_width], PM.cap_depth,
                             PM.ind_gap, PM.jn_gap, PM.buffer,
                             PM.track_mem, PM.gap, PM.track_mem, PM.gap_mem, PM.par_track,
                             fillet=PM.fillet)
    #PM.cap1.draw_capa_interdigitated(PM.par_track, PM.gap_mem, [0.5*PM.cap_depth-PM.par_track, PM.cap_gap], 0.5*PM.cap_gap, 0, '0um')
    #PM.cap2.draw_capa_interdigitated(PM.par_track, PM.gap_mem, [0.5*PM.cap_depth-PM.par_track, PM.cap_gap], 0.5*PM.cap_gap, 0, '0um')
    if litho:
#        chip1.set_current_coor([PM.x_I, PM.y_I - 0.5*PM.cutout_depth + 2*PM.buffer + 0.5*PM.cap_depth], [0,1])
        chip1.set_current_coor([PM.x_I, PM.y_I], [0,1])
        chip1.draw_alignement_marks('align_parity', '20um', ['0.15mm', '0.18mm'])
        
        if ebeam:
            if 1:
                chip1.set_current_coor([PM.x_I, PM.y_I - 0.5*PM.cutout_depth + 3*PM.buffer + PM.cap_depth], [0,1])
                chip1.draw_cos2phi('junction_parity', ['5um', '5um'], '90um', ['0.852um', '0.172um'], '400nm',
                                                  ['60um', '80um'], [46, 26],
                                                  spacing_bridge='1.047um')
            if 1:
                chip1.set_current_coor([PM.x_I, PM.y_I - 0.5*PM.cutout_depth + PM.buffer], Vector([0,1]).orth())
                chip1.draw_dose_test_junction('ind_parity', ['5um', '5um'], 
                                               '90um', '0.852um', '400nm',
                                               n_bridge=10, spacing_bridge='1.047um', alternate_width=False) #n_bridge = 52
                #PM.ind_parity.draw_meander_array(['5um', '5um'], '90um', '0.852um', '400nm',
                #                                  ['60um', '80um'], [52, 26],
                #                                  spacing_bridge='1.047um')

if memory and parity:
    chip1.draw_cable('mem', 'parity_portOut1', 'trm_portOut1', is_bond=is_bond, fillet=PM.fillet, is_mesh=True)
    
if flux_lines:
    PM.set_variable('approach', '40um')
    PM.set_variable('offset', '-40um')
    
    chip1.set_current_coor([PM.x_I + PM.offset, PM.y_I + 0.5*PM.cutout_depth - 1.75*PM.approach], [0,1])
    chip1.draw_fluxline('flux_top', PM.track, PM.gap, 10*PM.min_track_flux, PM.min_track_flux, sym='center', slope=0.2)

    chip1.set_current_coor([PM.x_I + 0.5*PM.cutout_width - PM.approach, PM.y_I - 0.55*PM.offset], [1,0])
    chip1.draw_fluxline('flux_bottom', PM.track, PM.gap, 10*PM.min_track_flux, 2*PM.min_track_flux, sym='center', slope=0.2)
#    print(PythonModeler.Port.dict_instances)
    chip1.draw_cable('flux_line_top', chip1.name+'in_flux_topiOut', 'flux_top_outPort', is_bond=is_bond, fillet = PM.fillet)
    chip1.draw_cable('flux_line_bottom', chip1.name+'in_other_fluxiOut', 'flux_bottom_outPort', is_bond=is_bond, fillet = PM.fillet)
    
    if 0:
        chip1.trackObjects[-2] = chip1.trackObjects[-2] + '_1'
        chip1.gapObjects[-2] = chip1.gapObjects[-2] + '_1'
        chip1.trackObjects[-1] = chip1.trackObjects[-1] + '_1'
        chip1.gapObjects[-1] = chip1.gapObjects[-1] + '_1'

if drive_mem:
    PM.set_variable('len_capa_drive', '100um')
    chip1.set_current_coor([con5, PM.y_I - 0.5*PM.cutout_depth - 0.5*PM.len_capa_drive], [0,1])
    chip1.draw_capa_inline('capa_drive', PM.track, PM.gap, PM.len_capa_drive, 0.5*(PM.track - 3*PM.track_mem), n_pad=3)
    
    chip1.draw_cable('chip1in_flux_topiOut', 'chip1in_memiOut', 'capa_drive_outPort2', is_bond=is_bond)

if not(readout):
    PM.set_variable('gap_capa_readout', '8um')
    chip1.set_current_coor(['0.9mm', con1], [1,0])
    chip1.draw_capa_inline('capa_readout', PM.track, PM.gap, '100um', PM.gap_capa_readout, n_pad=2)
    
    PM.set_variable('tune_ro', '3mm')
    chip1.double_port('constrain_readout1', [PM.x_T-PM.tune_ro, PM.y_T], [1,0], PM.track, PM.gap)
    chip1.double_port('constrain_readout2', [PM.x_T-0.3574*PM.tune_ro, 0.5*(PM.y_T+con1)], [0,-1] ,PM.track, PM.gap)

#    chip1.draw_cable('readout', 'trm_portOut2','constrain_readout1_front', 'constrain_readout1_back','constrain_readout2_front', 'constrain_readout2_back', 'capa_readout_outPort1', is_bond=is_bond, fillet=PM.fillet)
    chip1.draw_cable('bef_capa', 'capa_readout_outPort2', chip1.name+'in_readoutiOut', is_bond=is_bond, fillet=PM.fillet)
#

#%%
PM.interface.generate_gds("test_parity.gds")
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