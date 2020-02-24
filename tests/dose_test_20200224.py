# -*- coding: utf-8 -*-
"""
Created on Tue Jan 21 09:54:09 2020

@author: wcs
"""

'''
Third dose test for more oxidation and cross junctions
    First dose test had collapsed bridges on the small junctions
    and the resist was thicker (675 nm) than expected (600 nm),
    so only the larger bridge spacing worked
    Second dose test worked, and cross junctions had small Ej (~2-5 GHz)
    but we need to know about reproducibility
'''  
  
import scripts
import gdspy
from scripts.hfss import parse_entry as parse_entry
from scripts.Vector import Vector as Vector

PM = scripts.PythonMdlr('gds')

chip1 = PM.body('chip1', 'Global')

is_bond = False
is_overdev = False
is_litho = False

### Test junctions

PM.set_variable('x_D', '0mm')
PM.set_variable('y_D', '0mm')
PM.set_variable('dose_pad', '100um') #100um
PM.set_variable('dose_sep', '100um') #100um
PM.set_variable('dose_cor', 0.5*PM.dose_sep)
PM.set_variable('overlap', '20um')
PM.set_variable('bridge_width','400nm') #or 350nm ?
chip1.set_current_coor([PM.x_D, PM.y_D], [1,0])#'dose_test',
dose_mat = [8,16] #8,16

junction_width_1= ['2.0um','2.25um','2.5um','2.75um','3.0um','3.25um','3.5um','3.75um','4.0um','4.25um','4.5um','4.75um','5um']
junction_width_2= ['0.2um','0.25um','0.3um','0.35um','0.4um']
junction_array_width=['3.75um','4um','4.25um'] #or '6um' ?
bs_large = ['1012.5nm'] #bridge spacing
bridge_spacing_tab=['0um','1.3um','1.4um','1.5um','1.6um','1.7um']

chip1.draw_dose_test_Nb('dose_test', [PM.dose_pad, PM.dose_pad], PM.dose_sep, dose_mat, correction=PM.dose_cor, alum=False)

width = 2*dose_mat[0]*PM.dose_pad + 3*dose_mat[0]*PM.dose_sep - PM.dose_cor
height = dose_mat[1]*PM.dose_pad + (1 + dose_mat[1])*PM.dose_sep
xpos = PM.x_D - 0.5*width - 1.5*PM.dose_sep - PM.dose_pad - 0.5*PM.dose_cor
ypos = PM.y_D - 0.5*height + PM.dose_sep

chip1.set_current_coor([PM.x_D, PM.y_D], PM.dose_test.ori)
chip1.draw_alignement_marks('align_dose_test', '80um',['1.6mm', '1.8mm'] )#[0.6*width, 0.4*height]) #1.6 1.8mm

if not is_litho:
    for ii in range(2): #8
#        xpos=  xpos +*PM.dose_pad+1.25*PM.dose_sep
        xpos = xpos + 3*PM.dose_sep + 2*PM.dose_pad
        for jj in range(0,12): #single junctions
            chip1.set_current_coor([xpos, ypos + jj*PM.dose_sep + jj*PM.dose_pad+0.5*PM.dose_pad], [1,0])
            chip1.draw_dose_test_junction_corrected('junction_temp', ['10um', '10um'], 
                                        PM.dose_sep-PM.dose_cor+PM.overlap, junction_width_1[jj], PM.bridge_width,n_bridge=1, spacing_bridge=0,
                                        alternate_width=False, version=1, dose=True)
    for ii in range(2,5):
# def draw_dose_test_junction(self, pad_size, pad_spacing, width, width_bridge, iInduct='0nH', n_bridge=1, spacing_bridge=0, alternate_width=True, version=0, override=False, dose=False, rot=False, rotspace=None):
        for jj in range(0,7): #array junctions bridge spacing 1
            indjj=jj
            indii=ii-2
            chip1.set_current_coor([xpos, ypos + jj*PM.dose_sep + jj*PM.dose_pad+0.5*PM.dose_pad], [1,0])
            chip1.draw_dose_test_junction_corrected('junction_temp', ['10um', '10um'], 
                                        PM.dose_sep-PM.dose_cor+PM.overlap, junction_array_width[indii], PM.bridge_width,n_bridge=3, spacing_bridge=bridge_spacing_tab[indjj],
                                        alternate_width=False, version=1, dose=True)
            
    for ii in range(5,7):
        for jj in range(0,5): #array junctions bridge spacing 2
            indjj=jj
            indii=ii-5
            chip1.set_current_coor('junction_temp',[xpos, ypos +jj*PM.dose_sep + jj*PM.dose_pad+0.5*PM.dose_pad], [1,0])
            chip1.draw_dose_test_junction_corrected('junction_temp', ['10um', '10um'], 
                                        PM.dose_sep-PM.dose_cor+PM.overlap, junction_width_2[indjj], PM.bridge_width,n_bridge=3, spacing_bridge=bridge_spacing_tab[1],
                                        alternate_width=False, version=1, dose=True)
#        for jj in range(10,13) : #array junctions bridge spacing 3
#            indjj=jj-10
#            PM.key_elt('junction_temp',[xpos, ypos + jj*PM.dose_sep + jj*PM.dose_pad+0.5*PM.dose_pad], [1,0])
#            PM.junction_temp.draw_dose_test_junction_corrected(['10um', '10um'], 
#                                        PM.dose_sep-PM.dose_cor+PM.overlap, junction_array_width[indjj], PM.bridge_width,n_bridge=3, spacing_bridge=bridge_spacing_tab[2],
#                                        alternate_width=False, version=1, dose=True)
#        for jj in range(13,14):
#            indjj=jj-13
#            PM.key_elt('junction_temp',[xpos, ypos + jj*PM.dose_sep + jj*PM.dose_pad+0.5*PM.dose_pad], [1,0])
#            PM.junction_temp.draw_dose_test_junction_corrected(['10um', '10um'], 
#                                        PM.dose_sep-PM.dose_cor+PM.overlap, '0.3um', PM.bridge_width,n_bridge=1, spacing_bridge=bridge_spacing_tab[0],
#                                        alternate_width=False, version=1, dose=True)
#            
#if not KeyElt.is_litho:
#    for ii in range(1):
#        xpos = xpos + 3*PM.dose_sep + 2*PM.dose_pad
#        for jj in range(len(width_cross)):
#            PM.key_elt('temp', [xpos, ypos + jj*PM.dose_sep + jj*PM.dose_pad], [1,0])
#            PM.temp.draw_dose_test_junction(['10um', '10um'], 
#                                        PM.dose_sep-PM.dose_cor+PM.overlap, width_cross[jj], '350nm',
#                                        alternate_width=False, version=1, dose=True)
#        for jj in range(len(width_large)):
#            PM.key_elt('temp', [xpos, ypos + (len(width_cross) + jj)*PM.dose_sep + (len(width_cross) + jj)*PM.dose_pad], [1,0])
#            PM.temp.draw_dose_test_junction(['10um', '10um'], 
#                                            PM.dose_sep-PM.dose_cor+PM.overlap, width_large[jj], '337.5nm',
#                                            n_bridge=10, spacing_bridge=bs_large[0], alternate_width=False, dose=True)
#        for jj in range(1):
#            PM.key_elt('temp', [xpos, ypos + (len(width_cross) + len(width_large) + jj)*PM.dose_sep + (len(width_cross) + len(width_large) + jj)*PM.dose_pad], [1,0])
#            PM.temp.draw_dose_test_gral(PM.dose_sep-PM.dose_cor+PM.overlap, '10um', '2um', PM.dose_sep-PM.dose_cor)
#            
### Other stuff

bottom, top, left, right  = PM.get_extent(margin='500um')
width = right-left
height = top-bottom

#if KeyElt.is_litho:
#    PM.draw_rect('ground_plane', [left, bottom], [width, height])
#    PM.draw_rect('negative', [left, bottom], [width, height])
#    gapObject = PM.unite(PM.gapObjects)
#    PM.subtract(PM.ground_plane, [gapObject])
#    PM.subtract(PM.negative, [PM.ground_plane]+PM.trackObjects)

gdspy.write_gds('test.gds', unit=1.0, precision=1e-9)
