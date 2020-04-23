# -*- coding: utf-8 -*-
"""
Created on Tue Oct 29 14:05:00 2019

@author: wcs
"""
from drawpy_scripts.python_modeler import PythonModeler, ModelEntity
from drawpy_scripts.variable_string import parse_entry, Vector
import os

pm = PythonModeler('gds')

chip1 = pm.body('chip1')

qubit = True
memory = True
readout = True
drive_mem = True
parity = True
flux_lines = True

### Standard dimensions

chip_width = pm.set_variable('8.67mm')
chip_length = pm.set_variable('8.12mm')
chip_thickness = pm.set_variable('280um')
pcb_thickness = pm.set_variable('320um')
vaccuum_thickness = pm.set_variable(6*chip_thickness)

### Connectors

con1, con2, con3, con4, con5 = '4.06mm', '3.01mm', '7.54mm', '3.01mm', '6.36mm'

bond_length, bond_slope, pcb_track, pcb_gap = '200um', '0.5', '300um', '200um'
track = pm.set_variable('42um')
gap = pm.set_variable('25um')
track_mem = pm.set_variable('10um')
gap_mem = pm.set_variable('50um')
track_flux = pm.set_variable('50um')
min_track_flux = pm.set_variable('4um')
gap_flux = pm.set_variable('10um')
fillet = pm.set_variable('200um')

chip1.set_current_coor([con2, chip_length], [0, -1])
chip1.draw_connector(chip1.name+'in_flux_top', track, gap, bond_length,
                     pcb_track, pcb_gap, bond_slope)

chip1.set_current_coor([0, con1], [1, 0])
chip1.draw_connector(chip1.name+'in_readout', track, gap, bond_length,
                     pcb_track, pcb_gap, bond_slope)
chip1.set_current_coor([con5, 0], [0, 1])
chip1.draw_connector(chip1.name+'in_mem', track, gap, bond_length,
                     pcb_track, pcb_gap, bond_slope)
chip1.set_current_coor([con4, 0], [0, 1])
chip1.draw_connector(chip1.name+'in_unused', track, gap, bond_length,
                     pcb_track, pcb_gap, bond_slope)
chip1.set_current_coor([con3, chip_length], [0, -1])
chip1.draw_connector(chip1.name+'in_other_flux', track, gap, bond_length,
                     pcb_track, pcb_gap, bond_slope)
# print([(k.name, k.ori) for k in PythonModeler.Port.dict_instances.values()])
# pm.set_variable('gap_capa_readout', '8um')
# chip1.set_current_coor(['0.9mm', con1], [1,0])
# chip1.draw_capa_inline('capa_readout', track, gap, '100um',
#                        gap_capa_readout, n_pad=2)


# chip1.draw_cable('mem', chip1.name+'in_flux_topiOut',
#                  chip1.name+'in_unusediOut', is_bond=is_bond,
#                  fillet=fillet, is_mesh=True)

### Drawing

x_T = pm.set_variable('4.685mm')  # (con4 + con5)/2
y_T = pm.set_variable('2.5mm')

#
if qubit:
    Lj = pm.set_variable('12nH')
    trm_junction_width = pm.set_variable('10um')
    track_right = '0mm'
    track_left = '0mm'
    if readout:
        track_right = track
    if memory:
        track_left = track_mem

    chip1.set_current_coor([x_T, y_T], [-1, 0])
    chip1.draw_ZR_transmon('trm', ['1.47mm', '0.75mm'], '0.12mm',
                           ['0.5mm', '0.5mm'], track_right, gap, '0.30mm',
                           '30um', '0mm', trm_junction_width, '50um', Lj,
                           pad_size_left=['0.5mm', '0.5mm'],
                           track_left=track_left, gap_left=gap_mem,
                           length_left='0.2mm', spacing_left='-30um',
                           short_left='0um', fillet=True)

    if litho:
        chip1.set_current_coor([x_T, y_T], [-1, 0])
        chip1.draw_alignement_marks('align_trm', '20um', ['0.1mm', '0.43mm'])

        if ebeam:
            chip1.set_current_coor([x_T, y_T], [-1, 0])
            chip1.draw_dose_test_junction('junction_trm',
                                          ['10um', trm_junction_width],
                                          '50um', '0.147um', '400nm',
                                          alternate_width=False)

if parity:
    cutout_depth = pm.set_variable('0.26mm')
    cutout_width = pm.set_variable('0.295mm')
    cap_gap = pm.set_variable('4um')
    cap_depth = pm.set_variable('0.0mm')
    ind_gap = pm.set_variable('90um')
    jn_gap = pm.set_variable('90um')
    buffer = pm.set_variable('50um')
    par_track = pm.set_variable('5um')
    x_I = pm.set_variable(con5-0.5*(-track_mem+par_track+ind_gap)-buffer)
    y_I = pm.set_variable('2.9mm')
    chip1.set_current_coor([x_I, y_I], [0, 1])
    chip1.draw_selfparity('parity', [cutout_depth, cutout_width],
                          cap_depth, ind_gap, jn_gap, buffer,
                          track_mem, gap, track_mem, gap_mem,
                          par_track, fillet=fillet)

    if litho:
        chip1.set_current_coor([x_I, y_I], [0, 1])
        chip1.draw_alignement_marks('align_parity', '20um',
                                    ['0.15mm', '0.18mm'])

        if ebeam:

            chip1.set_current_coor([x_I, y_I - 0.5*cutout_depth + 3*buffer + cap_depth], [0, 1])
            chip1.draw_cos2phi('junction_parity', ['5um', '5um'], '90um', ['0.852um', '0.172um'], '400nm',
                                              ['60um', '80um'], [46, 26],
                                              spacing_bridge='1.047um')

            chip1.set_current_coor([x_I, y_I - 0.5*cutout_depth + buffer], Vector([0,1]).orth())
            chip1.draw_dose_test_junction('ind_parity', ['5um', '5um'],
                                           '90um', '0.852um', '400nm',
                                           n_bridge=10, spacing_bridge='1.047um', alternate_width=False) #n_bridge = 52
                #pm.ind_parity.draw_meander_array(['5um', '5um'], '90um', '0.852um', '400nm',
                #                                  ['60um', '80um'], [52, 26],
                #                                  spacing_bridge='1.047um')

if memory and parity:
    chip1.draw_cable('mem', 'parity_portOut1', 'trm_portOut1', is_bond=is_bond, fillet='200um', is_mesh=True)

if flux_lines:
    approach = pm.set_variable('40um')
    offset = pm.set_variable('-40um')

    chip1.set_current_coor([x_I + offset, y_I + 0.5*cutout_depth - 1.75*approach], [0,1])
    chip1.draw_fluxline('flux_top', track, gap, 10*min_track_flux, min_track_flux, sym='center', slope=0.2)

    chip1.set_current_coor([x_I + 0.5*cutout_width - approach, y_I - 0.55*offset], [1,0])
    chip1.draw_fluxline('flux_bottom', track, gap, 10*min_track_flux, 2*min_track_flux, sym='center', slope=0.2)
#    print(PythonModeler.Port.dict_instances)
    chip1.draw_cable('flux_line_top', chip1.name+'in_flux_topiOut', 'flux_top_outPort', is_bond=is_bond, fillet = fillet)
    chip1.draw_cable('flux_line_bottom', chip1.name+'in_other_fluxiOut', 'flux_bottom_outPort', is_bond=is_bond, fillet = fillet)

    if 0:
        chip1.trackObjects[-2] = chip1.trackObjects[-2] + '_1'
        chip1.gapObjects[-2] = chip1.gapObjects[-2] + '_1'
        chip1.trackObjects[-1] = chip1.trackObjects[-1] + '_1'
        chip1.gapObjects[-1] = chip1.gapObjects[-1] + '_1'

if drive_mem:
    len_capa_drive = pm.set_variable('100um')
    chip1.set_current_coor([con5, y_I - 0.5*cutout_depth - 0.5*len_capa_drive], [0,1])
    chip1.draw_capa_inline('capa_drive', track, gap, len_capa_drive, 0.5*(track - 3*track_mem), n_pad=3)

#    chip1.draw_cable('chip1in_flux_topiOut', 'chip1in_memiOut', 'capa_drive_outPort2', is_bond=is_bond)

if not(readout):
    gap_capa_readout = pm.set_variable('8um')
    chip1.set_current_coor(['0.9mm', con1], [1,0])
    chip1.draw_capa_inline('capa_readout', track, gap, '100um', gap_capa_readout, n_pad=2)

    tune_ro = pm.set_variable('3mm')
    chip1.double_port('constrain_readout1', [x_T-tune_ro, y_T], [1,0], track, gap)
    chip1.double_port('constrain_readout2', [x_T-0.3574*tune_ro, 0.5*(y_T+con1)], [0,-1] ,track, gap)

#    chip1.draw_cable('readout', 'trm_portOut2','constrain_readout1_front', 'constrain_readout1_back','constrain_readout2_front', 'constrain_readout2_back', 'capa_readout_outPort1', is_bond=is_bond, fillet=fillet)
    chip1.draw_cable('bef_capa', 'capa_readout_outPort2', chip1.name+'in_readoutiOut', is_bond=is_bond, fillet=fillet)
#
#gdspy.write_gds('test.gds', unit=1.0, precision=1e-9)

cwd  = os.getcwd()
gdspy.write_gds(os.path.join(cwd, 'test.gds'), unit=1.0, precision=1e-9)
#pm.interface.generate_gds('test_parity.gds')

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