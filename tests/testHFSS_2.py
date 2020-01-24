# -*- coding: utf-8 -*-
"""
Created on Fri Nov 22 17:18:25 2019

@author: Zaki
"""
#%% Connect JJ
import scripts # We import the package with gdsPy scripts.
from scripts import Vector as Vector
PM = scripts.PythonMdlr('hfss') # We set the interface to gds.

#Beginning of the instructions


#End of the instructions



PM.set_variable('track', '42um')
PM.set_variable('bond', '100um')
PM.set_variable('pad_spacing', '0.12mm')


chip1 = PM.body('chip2', 'Global')
chip1.set_current_coor(pos = ['0mm', '10mm','0mm'], ori=[0,1])
P1 = chip1.port('port1', Vector(['10mm',0]), Vector([0,1]), '0.2mm', '0.1mm')
chip1.set_current_coor(pos = ['0mm', '10mm','0mm'], ori=[0,1])
P2 = chip1.port('port2', Vector(['10mm','10mm']), Vector([0,1]), '0.1mm','0.05mm')
#SL_PTH = chip1._connect_JJ('jojo', 'port1', 'port2', "0.02mm")

#%% draw_cable
import PythonModeler
import traceback
import hfss
from Vector import Vector

PythonModeler.Port.reset()
hfss.ModelEntity.reset()
PM = PythonModeler.PythonMdlr('hfss')
PM.set_variable('track', '42um')
PM.set_variable('bond', '100um')
PM.set_variable('pad_spacing', '0.12mm')

chip1 = PM.body('chip2', 'Global')
chip1.set_current_coor(pos = ['0mm', '0mm','0mm'], ori=[0,1])
P1 = chip1.port('port1', Vector(['10mm','0mm']), Vector([1,0]), '0.2mm', '0.1mm')
P2 = chip1.port('port2', Vector(['10mm','10mm']), Vector([0,1]), '0.1mm','0.05mm')
chip1.draw_cable("cable", "port1", "port2")
P3 = chip1.port('port3', Vector(['10mm','0mm']), Vector([1,0]), '0.2mm', '0.1mm')

chip1.cable_starter('CBSTRT', 'port3')
#%% capa_interdigitated

import PythonModeler
import traceback
import hfss
from Vector import Vector

PythonModeler.Port.reset()
hfss.ModelEntity.reset()
PM = PythonModeler.PythonMdlr('hfss')
PM.set_variable('track', '42um')
PM.set_variable('bond', '100um')
PM.set_variable('pad_spacing', '0.12mm')

chip1 = PM.body('chip2', 'Global')
P1 = chip1.port('port1', Vector(['10mm','0mm']), Vector([1,0]), '0.2mm', '0.1mm')
chip1.set_current_coor(pos = ['0mm', '0mm','0mm'], ori=[0,1])
#chip1.draw_capa_inline('CAPA_IN_LINE', '20mm', '10mm', '5mm', '1mm', n_pad=5, iTrack_capa='5mm', iGap_capa=None, premesh=True, tight=False)
chip1.draw_capa_interdigitated('CAPA_INTERDIGITATED', '20mm', '10mm', ['2mm','2mm'], ['20mm','20mm'], 10, '0.1mm')

#%% 3D cavity + transmon

import PythonModeler
import traceback
import hfss

PythonModeler.Port.reset()
hfss.ModelEntity.reset()

PM = PythonModeler.PythonMdlr('hfss')
PM.set_variable('track', '42um')
PM.set_variable('bond', '100um')
PM.set_variable('cylinder_height', '30mm')
PM.set_variable('antenna_height', '10mm')
PM.set_variable('cutout_x', '1.47mm')
PM.set_variable('cutout_y', '0.75mm')

PM.set_variable('cylinder_radius', '4mm')


chip1 = PM.body('chip1', "chip_1", [['0mm',PM.cylinder_radius,PM.antenna_height], [1,0,0], [0,1,0]])
chip1.set_current_coor(pos = ['0mm', '0mm','0mm'], ori=[0,1])
chip1.insert_transmon("insert", [PM.cutout_x, PM.cutout_y, '0.5mm'],'0.05mm',['0.2mm', '0.2mm'],'0.1mm', '0.35mm','0.42mm' ,'0.25um')

chip2 = PM.body('chip2', 'Global')
chip2.set_current_coor(pos = ['0mm', '0mm','0mm'], ori=[1,0])
cavity = chip2.cavity_3D_simple('cavity', PM.cylinder_radius,PM.cylinder_height, '3mm', PM.antenna_height, '1mm', '1mm')
#cavity = chip1.cavity_3D_simple('cavity', '3mm', '5mm', '0.5mm', '2.5mm', '2mm', '1mm')

#%% Cavity + transmon + Box

from importlib import reload
import PythonModeler
import hfss
reload(PythonModeler)

PythonModeler.Port.reset()
hfss.ModelEntity.reset()

PM = PythonModeler.PythonMdlr('hfss')
PM.set_variable('track', '42um')
PM.set_variable('bond', '100um')
PM.set_variable('cylinder_height', '30mm')
PM.set_variable('antenna_height', '10mm')
PM.set_variable('cutout_x', '1.47mm')
PM.set_variable('cutout_y', '0.75mm')

PM.set_variable('cylinder_radius', '4mm')


chip1 = PM.body('chip1', "chip_1", [['0mm',PM.cylinder_radius,PM.antenna_height], [1,0,0], [0,1,0]])
chip1.set_current_coor(pos = ['0mm', '0mm','0mm'], ori=[1,0])
sapphire = chip1.box_center([0,0,'-0.25mm'], [PM.cutout_y, PM.cutout_x, '0.5mm'], "3D", name="sapphire")
chip1.make_material(sapphire, "\"sapphire\"")
chip1.insert_transmon("insert", [PM.cutout_x, PM.cutout_y, '0.5mm'],'0.05mm',['0.2mm', '0.2mm'],'0.1mm', '0.35mm','0.42mm' ,'0.25um')


chip2 = PM.body('chip2', 'Global')
chip2.set_current_coor(pos = ['0mm', '0mm','0mm'], ori=[1,0])
cavity = chip2.cavity_3D_simple('cavity', PM.cylinder_radius,PM.cylinder_height, PM.cylinder_radius/2, PM.antenna_height, '1mm', '1mm')
#cavity = chip1.cavity_3D_simple('cavity', '3mm', '5mm', '0.5mm', '2.5mm', '2mm', '1mm')

#%% 
from importlib import reload
import PythonModeler
import hfss
reload(PythonModeler)

PythonModeler.Port.reset()
hfss.ModelEntity.reset()

PM = PythonModeler.PythonMdlr('hfss')
PM.set_variable('track', '42um')
PM.set_variable('bond', '100um')
PM.set_variable('cylinder_height', '30mm')
PM.set_variable('antenna_height', '10mm')
PM.set_variable('cutout_x', '1.47mm')
PM.set_variable('cutout_y', '0.75mm')
PM.set_variable('cylinder_radius', '4mm')
PM.set_variable('cable_length', '1mm')

body1 = PM.body('chip1', "chip_1", [['0mm',PM.cylinder_radius,PM.antenna_height], [1,0,0], [0,1,0]])
body1.set_current_coor(pos = ['0mm', '0mm','-0mm'], ori=[1,0])
outer_cable = body1.cylinder([0,'-0.5mm',0], '0.5mm', PM.cable_length, 'Y', layer='3D', name='cable1')
inner_cable = body1.cylinder([0,'-0.5mm',0], '0.2mm', PM.cable_length, 'Y', layer='3D', name='cable2')
body1.assign_perfect_E(inner_cable)

chip2 = PM.body('chip2', 'Global')
chip2.set_current_coor(pos = ['0mm', '0mm','0mm'], ori=[1,0])
cavity = chip2.cavity_3D('cavity', PM.cylinder_radius,PM.cylinder_height, '3mm', PM.antenna_height)
#cavity = chip1.cavity_3D_simple('cavity', '3mm', '5mm', '0.5mm', '2.5mm', '2mm', '1mm')


#%% Only transmon
from importlib import reload
import PythonModeler
import hfss
#reload(PythonModeler)

PythonModeler.Port.reset()
hfss.ModelEntity.reset()

PM = PythonModeler.PythonMdlr('hfss')

PM.set_variable('box_height', '10mm')
PM.set_variable('box_size', '10mm')
PM.set_variable('cutout_x', '3mm')
PM.set_variable('cutout_y', '1.5mm')


chip1 = PM.body('chip1', "Global")
chip1.set_current_coor(pos = ['0mm', '0mm','0mm'], ori=[1,0])
chip1.insert_transmon("insert", [PM.cutout_x, PM.cutout_y, '0.2mm'],'0.1mm',['0.5mm', '0.5mm'],'0.1mm', '0.35mm','0.42mm' ,'0.25um')
vacuum = chip1.box_center([0,0,0], [PM.box_size, PM.box_size, PM.box_height], name='vacuum', layer='3D')
chip1.assign_perfect_E(vacuum, 'environment')

#%% Only Cavity

from importlib import reload
import PythonModeler
import hfss
reload(PythonModeler)

PythonModeler.Port.reset()
hfss.ModelEntity.reset()

PM = PythonModeler.PythonMdlr('hfss')
PM.set_variable('cylinder_height', '15mm')
PM.set_variable('antenna_height', '3mm')
PM.set_variable('cylinder_radius', '2mm')

chip2 = PM.body('chip2', 'Global')
chip2.set_current_coor(pos = ['0mm', '0mm','0mm'], ori=[1,0])
cylinder = chip2.cavity_3D('cavity', PM.cylinder_radius,PM.cylinder_height, PM.cylinder_radius/4, PM.antenna_height)

#%% Cavity + insert
from importlib import reload
import PythonModeler
import hfss
reload(PythonModeler)

PythonModeler.Port.reset()
hfss.ModelEntity.reset()

PM = PythonModeler.PythonMdlr('hfss')
PM.set_variable('cylinder_height', '15mm')
PM.set_variable('antenna_height', '3mm')
PM.set_variable('cylinder_radius', '2mm')

chip2 = PM.body('chip2', 'Global')
chip2.set_current_coor(pos = ['0mm', '0mm','0mm'], ori=[1,0])
cylinder = chip2.cavity_3D_simple('cavity', PM.cylinder_radius,PM.cylinder_height, PM.cylinder_radius/4, PM.antenna_height,'0.5mm', '0.2mm' )

#%% Transmon plus cylinder
import scripts


PM = scripts.PythonMdlr('hfss')
PM.set_variable('cylinder_height', '10mm')
PM.set_variable('antenna_height', '3mm')
PM.set_variable('cutout_x', '1.2mm')
PM.set_variable('cutout_y', '0.6mm')
PM.set_variable('cylinder_radius', '1mm')

chip2 = PM.body('chip2', 'Global')
chip2.set_current_coor(pos = ['0mm', '0mm','0mm'], ori=[1,0])
chip2.cavity_3D_simple('cavity', PM.cylinder_radius,PM.cylinder_height, PM.cylinder_radius/4, PM.antenna_height, '0.5mm', '0.2mm')
chip1 = PM.body('chip1', "chip_1", [['0mm',PM.cylinder_radius,PM.antenna_height], [0,1,0], [0,0,-1]])
chip1.set_current_coor(pos = ['0mm', '0mm','0mm'], ori=[1,0])
chip1.insert_transmon("insert", [PM.cutout_x, PM.cutout_y, '0.2mm'],'0.1mm',['0.5mm', '0.5mm'],'0.1mm', '0.35mm','0.42mm' ,'0.25um')

#poly = chip2.polyline_2D([[0,0], [0,1], [1,1]],closed = False, name='nom', layer='POLY')
#chip2._fillets(0.1, poly)

#%%
import scripts
PM = scripts.PythonMdlr('hfss')
PM.set_variable('rext', '6.5mm')
PM.set_variable('hint', '6mm')
PM.set_variable('rint', '2mm')
PM.set_variable('hext', '30mm')
chip = PM.body('chip', 'Global')
PM.set_variable('W_chip', '5mm')

cavity_param = [PM.rext, PM.hext, PM.rint, PM.hint]
transmons_param = []
ports_param = []
PM.set_variable('L_insert1', '10mm')
PM.set_variable('W_insert1', '6mm')
PM.set_variable('pene1', '1.5mm')
PM.set_variable('ep1', '0.2mm')
PM.set_variable('L_tr1', '1mm')
PM.set_variable('W_tr1', '0.5mm')
PM.set_variable('Lj_tr1', '0.1mm')
PM.set_variable('Wj_tr1', '0.05mm')
PM.set_variable('ecart_tr1', '0.2mm')
PM.set_variable('h_tr5', '1mm')
PM.set_variable('Induc_tr1', '5mm')

#transmons_param.append([0, -PM.rext+PM.pene1, PM.hint+PM.h_tr1, -90, ["tr1", [PM.L_insert1,PM.W_chip], PM.Lj_tr1, [PM.L_tr1, PM.W_tr1], PM.Wj_tr1, PM.ecart_tr1, PM.ep1, PM.Induc_tr1 ]])

PM.set_variable('L_insert2', '25mm')
PM.set_variable('W_insert2', '6mm')
PM.set_variable('pene2', '2mm')
PM.set_variable('ep2', '0.2mm')
PM.set_variable('L_tr2', '1mm')
PM.set_variable('W_tr2', '0.5mm')
PM.set_variable('Lj_tr2', '0.1mm')
PM.set_variable('Wj_tr2', '0.05mm')
PM.set_variable('ecart_tr2', '0.2mm')
PM.set_variable('Induc_tr2', '5mm')
PM.set_variable('h_tr2', '1mm')

PM.set_variable('dQRO', '2.7mm')
PM.set_variable('W_RO', '0.5mm')
PM.set_variable('L_RO', '14mm')
PM.set_variable('angle_rot', '0deg')

transmons_param.append([0, PM.rext-PM.pene2, PM.hint+PM.h_tr2, 90, ["tr2", [PM.L_insert2, PM.W_chip], PM.Lj_tr2, [PM.L_tr2, PM.W_tr2], PM.Wj_tr2, PM.ecart_tr2, PM.dQRO, PM.L_RO, PM.W_RO, PM.ep2, PM.Induc_tr2 ]])


PM.set_variable('pyPP', '0mm')
PM.set_variable('pzPP', '0mm')
PM.set_variable('rPP', '0mm')
PM.set_variable('hPP', '0mm')
PM.set_variable('hLTP', '0mm')
PM.set_variable('rdiElecP', '0mm')
PM.set_variable('rgaineP', '1mm')
PM.set_variable('rameP', '0mm')


PM.set_variable('x_port1', '0mm')
PM.set_variable('y_port1', '22mm')
PM.set_variable('z_port1', '-2.7mm')
PM.set_variable('h_port1', '5mm')
PM.set_variable('r_cable_port1', '2.5mm')
PM.set_variable('r_ame_port1', '0.64mm')
PM.set_variable('r_gaine_port1', '2.1mm')
PM.set_variable('h_connector_port1', '4mm')

ports_param.append([PM.x_port1, PM.y_port1, PM.z_port1, 90, ["port_tomo", PM.r_cable_port1, PM.r_gaine_port1, PM.r_ame_port1, PM.h_port1, PM.h_connector_port1, "Z"]])


PM.set_variable('rdiElecT', '0mm')
PM.set_variable('rgaineT', '0mm')
PM.set_variable('rameT', '0mm')
PM.set_variable('pxPM', '-13mm')
PM.set_variable('pyPM', '0mm')
PM.set_variable('pzPM', '25mm')

#ports_param.append([PM.pxPM, PM.pyPM, PM.pzPM, 90, ["port_purcell", PM.r_cable_port1, PM.r_gaine_port1, PM.r_ame_port1, PM.h_port1, PM.h_connector_port1, "Z"]])


chip.set_current_coor(pos = ['0mm', '0mm','0mm'], ori=[1,0])
union = chip.cavity_3D_with_ports("ports_cavity", cavity_param, transmons_param, ports_param)

chip1 = PM.body('chip1', "chip_1", [[PM.pxPM,PM.pyPM,PM.pzPM], [0,1,0], [0,0,1]])
chip1.set_current_coor(pos = ['0mm', '0mm','0mm'], ori=[1,0])
gaine_probe, ame_probe = chip1.cable_3D(*["port_purcell", PM.r_cable_port1, PM.r_gaine_port1, PM.r_ame_port1, PM.h_port1, PM.h_connector_port1, "Z"])
chip1.unite([union, gaine_probe], "union2")
chip1.subtract(union, [ame_probe])
chip1.assign_perfect_E(union)
#chip.cable_3D("port1", PM.r_cable_port1, PM.r_gaine_port1, PM.r_ame_port1, PM.h_port1, PM.h_connector_port1, "Z")