'''
TODO

- debug comments
- automate Bonds
- delete "mm", instead, select the units in HFSS
- define in classes
class CircuitElement(object):
    def __init__(self, gates, params):

    def get_gates(self):
        returns list of gates (inputs and outpts)

class Capacitor(CircuitElement):
    def draw(self):
        draws the element

- LitExp: check if variable exists, updates value if exists, create a new one if not.
- Create drawing script which can start from blank file
- Premesh

'''

'''
Assumption

- only assign lumpRLC to rectangles
- assume to have list and not tuples for polylines
- TODO, do not do Lj+'1nH' for now, but can do bx+'1mm'
'''



from scripts.hfss import get_active_project, release, parse_entry
from scripts.designer import Circuit, KeyElt, ConnectElt, Vector
import numpy as np

project = get_active_project()
design = project.get_active_design()
modeler = design.modeler

modeler.set_units('mm')

c = Circuit(design, modeler)

KeyElt.is_mask = False
KeyElt.gap_mask = parse_entry('20um')
#KeyElt.overdev = parse_entry('0.9um')

#######################
### DRAWING STARTS HERE 
#######################

### PCB, Chip and Ground Plane

#c.set_variable("chip_width", "8.67mm")
#c.set_variable("chip_length", "8.12mm")
#c.set_variable("chip_thickness", "280um")
#c.set_variable("pcb_thickness", "320um")
#c.set_variable('vaccuum_thickness', 6*c.chip_thickness)
#
#c.draw_box('pcb', [0, 0, -c.chip_thickness], [c.chip_width, c.chip_length, -c.pcb_thickness], "Rogers TMM 10i (tm) no loss")
#c.draw_box('chip', [0, 0, 0], [c.chip_width, c.chip_length, -c.chip_thickness], 'silicon')
#c.draw_box('box', [0, 0, 0], [c.chip_width, c.chip_length, c.vaccuum_thickness], 'vacuum')
#
#c.draw_rect('ground_plane', [0, 0], [c.chip_width, c.chip_length])
#c.assign_perfE(c.ground_plane)


### Connectors

#con1, con2, con3, con4, con5 = '4.06mm', '3.01mm', '7.54mm', '3.01mm', '6.36mm'
#KeyElt.pcb_track, KeyElt.pcb_gap = parse_entry('300um'), parse_entry('200um')
#bond_length, bond_slope = '200um', '0.5'
#c.set_variable('track', '84um')
#c.set_variable('gap', '50um')
#c.set_variable('track_mem', '3um')
#c.set_variable('gap_mem', '50um')
#c.set_variable('track_flux', '50um')
#c.set_variable('gap_flux', '30um')
#
#
#c.key_elt('in_mem', [con3, c.chip_length], [0, -1])
#c.key_elt('in_flux_top', [con2, c.chip_length], [0, -1])
#c.key_elt('in_buffer', [0, con1], [1, 0])
#c.key_elt('in_readout', [con5, 0], [0, 1])
#c.key_elt('in_flux_bot', [con4, 0], [0, 1])

#c.in_mem.draw_connector(c.track, c.gap, bond_length, bond_slope)
#c.in_flux_top.draw_connector(c.track_flux, c.gap_flux, bond_length, bond_slope)
#c.in_buffer.draw_connector(c.track, c.gap, bond_length, bond_slope)
#c.in_readout.draw_connector(c.track, c.gap, bond_length, bond_slope)
#c.in_flux_bot.draw_connector(c.track_flux, c.gap_flux, bond_length, bond_slope)

### Circuit
#

# for ,ulti junction dose test
#for ii, val in enumerate(np.linspace(0.3,1,8)):
#    width = str(val)+'um'
#    pos = str(ii*400)+'um'
#    c.key_elt('pads', ['000um',pos], [1,0])
#    c.pads.draw_dose_test2(['100um', '100um'], '50um', width,width, '400e-9')
#    c.key_elt('padz', ['400um',pos], [1,0])
#    c.padz.draw_dose_test2(['100um', '100um'], '50um', width,width, '0e-9')
#    
#
#for ii, val in enumerate(np.linspace(1,3.5,6)):
#    width = str(val)+'um'
#    pos = str(ii*400)+'um'
#    c.key_elt('pads', ['800um',pos], [1,0])
#    c.pads.draw_dose_test2(['100um', '100um'], '50um', width,width, '400e-9')
#    c.key_elt('padz', ['1200um',pos], [1,0])
#    c.padz.draw_dose_test2(['100um', '100um'], '50um', width,width, '0e-9')
#    

#for ii, val in enumerate(np.linspace(1.9,2.3,5)):
#    width = str(val)+'um'
#    pos = str(ii*400)+'um'
#    c.key_elt('pads', ['000um',pos], [1,0])
#    c.pads.draw_dose_test2(['100um', '100um'], '50um', '9um',width, '400e-9')
#    c.key_elt('padz', ['200um',pos], [1,0])
#    c.padz.draw_dose_test2(['100um', '100um'], '50um', '9um',width, '0e-9')


#c.key_elt('pads', ['400um','00um'], [1,0])
#c.pads.draw_dose_test2(['100um', '100um'], '50um', '9um', 0, '400e-9')

#for ii, val in enumerate(np.linspace(0.3,1,8)):
#    width = str(val)+'um'
#    pos = str(ii*400)+'um'
#    c.key_elt('pads', ['000um',pos], [1,0])
#    c.pads.draw_dose_test2(['100um', '100um'], '50um', width,width, '400e-9')
#    c.key_elt('padz', ['400um',pos], [1,0])
#    c.padz.draw_dose_test2(['100um', '100um'], '50um', width,width, '0e-9')

for ii, val in enumerate(np.arange(1, 10.1, 1)):
    pos_x = str(ii//10*400)+'um'
    pos_y = str(ii%10*200)+'um'
    c.key_elt('pads', [pos_x,pos_y], [1,0])
    width = str(val)+'um'
    c.pads.draw_dose_test_junction(['100um', '100um'], '50um', width, '500nm', n_bridge=1, spacing_bridge='2um')#(pad_size, pad_spacing, width, width_bridge, n_bridge=1, spacing_bridge=0)

for ii, val in enumerate(np.arange(1, 10.1, 1)):
    pos_x = str(ii//10*400+400)+'um'
    pos_y = str(ii%10*200)+'um'
    c.key_elt('pads', [pos_x,pos_y], [1,0])
    width = str(val)+'um'
    c.pads.draw_dose_test_junction(['100um', '100um'], '50um', width, '1000nm', n_bridge=1, spacing_bridge='2um')#(pad_size, pad_spacing, width, width_bridge, n_bridge=1, spacing_bridge=0)

for ii, val in enumerate(np.arange(1, 10.1, 1)):
    pos_x = str(ii//10*400+800)+'um'
    pos_y = str(ii%10*200)+'um'
    c.key_elt('pads', [pos_x,pos_y], [1,0])
    width = str(val)+'um'
    c.pads.draw_dose_test_junction(['100um', '100um'], '50um', width, '1500nm', n_bridge=1, spacing_bridge='2um')#(pad_size, pad_spacing, width, width_bridge, n_bridge=1, spacing_bridge=0)

for ii, val in enumerate(np.arange(1, 10.1, 1)):
    pos_x = str(ii//10*400+1200)+'um'
    pos_y = str(ii%10*200)+'um'
    c.key_elt('pads', [pos_x,pos_y], [1,0])
    width = str(val)+'um'
    c.pads.draw_dose_test_junction(['100um', '100um'], '50um', width, '2000nm', n_bridge=1, spacing_bridge='2um')#(pad_size, pad_spacing, width, width_bridge, n_bridge=1, spacing_bridge=0)

for ii, val in enumerate(np.arange(1, 10.1, 1)):
    pos_x = str(ii//10*400+1600)+'um'
    pos_y = str(ii%10*200)+'um'
    c.key_elt('pads', [pos_x,pos_y], [1,0])
    width = str(val)+'um'
    c.pads.draw_dose_test_junction(['100um', '100um'], '50um', width, '2500nm', n_bridge=1, spacing_bridge='2um')#(pad_size, pad_spacing, width, width_bridge, n_bridge=1, spacing_bridge=0)

#c.set_variable('x_T', '3.8mm')
#c.set_variable('y_T', '5.5mm')
#c.key_elt('T', [c.x_T,c.y_T], [1,0])
#c.T.draw_T(c.track_mem, c.gap_mem, is_overdev=True)

#c.set_variable('squid_width', '50um')
#c.set_variable('squid_length', '50um')
#c.key_elt('squid', [con2+c.squid_width/2+c.track_mem, '5.5mm'], [-1,0])
#c.squid.draw_squid_protect(c.track_mem, c.gap_mem, [c.squid_width, c.squid_length], c.track_flux, c.gap_flux, iTrackSquid=c.track_mem, iTrackJ=None, Lj_down='20nH', Lj_up=None,  typePump='down', doublePump=True, iSlope=1, iSlopePump=0.5, fillet=None, is_overdev=True)

#c.key_elt('end_mem', ['4.8mm','2.5mm'], [0,-1])
#c.end_mem.draw_end_cable(c.track_mem, c.gap_mem, typeEnd='short', is_overdev=True)

#c.key_elt('end_buffer', ['1.71mm','3.1mm'], [0,1])
#c.end_buffer.draw_end_cable(c.track_mem, c.gap_mem, typeEnd='short', is_overdev=True)


#c.key_elt('end_drive_buffer', c.squid.pos+Vector(['-0.15mm','0.08mm']), [-1,0])
#c.end_drive_buffer.draw_end_cable(c.track, c.gap, typeEnd='open')

#c.key_elt('end_drive_mem', [c.x_T-'0.4mm',c.y_T+'0.25mm'], [0,1])
#c.end_drive_mem.draw_end_cable(c.track, c.gap, typeEnd='open')

#c.set_variable('gap_capa_readout', '40um')
#c.key_elt('capa_readout', [con5,'0.9mm'], [0,1])
#c.capa_readout.draw_capa_inline(c.track, c.gap, '80um', c.gap_capa_readout , n_pad=1)




#c.connect_elt('squid_T', 'squid_2', 'T_2')
#c.squid_T.draw_cable(is_bond=False)
#
#c.connect_elt('T_qubit', 'T_1', 'trm_2')
#c.T_qubit.draw_cable(is_bond=False)
#
#c.connect_elt('flux_top', 'squid_pump1', 'in_flux_top')
#c.flux_top.draw_cable(is_bond=False)
#
#c.connect_elt('flux_bot', 'squid_pump2', 'in_flux_bot')
#c.flux_bot.draw_cable(is_bond=False)

#c.set_variable('fillet', '500um')
#c.connect_elt('buffer', 'squid_1', 'end_buffer')
#c.buffer.draw_cable(is_bond=False, fillet=c.fillet)
#c.connect_elt('mem', 'T_3', 'end_mem')
#c.mem.draw_cable(is_bond=False, fillet=c.fillet)

#c.connect_elt('drive_mem', 'in_mem', 'end_drive_mem')
#c.drive_mem.draw_cable(is_bond=False, fillet=c.fillet)
#
#c.key_elt('constr_readout', c.in_readout.pos+Vector(['1.9mm', '4mm']), [0,1])
#c.constr_readout.create_port(c.track, c.gap)
#c.connect_elt('readout', 'capa_readout_1', 'trm_1')
#c.readout.draw_cable(is_bond = False, fillet=c.fillet, constrains=['constr_readout'])
#c.connect_elt('bef_capa', 'capa_readout_2', 'in_readout')
#c.bef_capa.draw_cable(is_bond = False, fillet=c.fillet)
#
#c.key_elt('constr_drive_buffer', c.in_buffer.pos+Vector(['1.26mm', '0.75mm']), [0,1])
#c.constr_drive_buffer.create_port(c.track, c.gap)
#c.connect_elt('drive_buffer', 'in_buffer', 'end_drive_buffer')
#c.drive_buffer.draw_cable(is_bond = False, fillet=c.fillet, constrains=['constr_drive_buffer'])
##
#c.key_elt('array', ['0','0'], [1,0])
#c.array.draw_snail_array('600nm', '9um', 3, '1um', 1, 20, '400nm', ['10um', '10um'])
#    def draw_ZR_tansmon(self,
#                        cutout_size,
#                        pad_spacing,
#                        pad_size_left,
#                        track_left,
#                        gap_left,
#                        length_left,
#                        spacing_left,
#                        short_left,
#                        Jwidth,
#                        Jinduc,
#                        pad_size_right=None,
#                        track_right=None,
#                        gap_right=None,
#                        length_right=None,
#                        spacing_right=None,
#                        short_right=None,
#                        fillet=None,
#                        maskbox=False):
#
#c.set_variable('offsetX_capa_right', '2mm')
#c.set_variable('offsetY_capa_right', '3mm')
#c.set_variable("capa_right_spacing", '0.01mm')
#c.set_variable("capa_right_width", c.track/2)
#c.set_variable("capa_right_length", '0.45mm')
#c.key_elt('capa_right', [c.chip_width-c.offsetX_capa_right, c.chip_length-c.offsetY_capa_right], [-1,0])
#c.capa_right.draw_capa_inline(c.track, c.gap, c.capa_right_length, c.capa_right_spacing, n_pad=5)
#
#c.set_variable('offsetX_capa_left', '1mm')
#c.set_variable('offsetY_capa_left', '3mm')
#c.set_variable("capa_left_spacing", c.capa_right_spacing)
#c.set_variable("capa_left_width", c.track/2)
#c.set_variable("capa_left_length", '0.45mm')
#c.key_elt('capa_left', [c.chip_width-c.offsetX_capa_left, c.offsetY_capa_left], [-1,0])
#c.capa_left.draw_capa_inline(c.track, c.gap, c.capa_left_length, c.capa_left_spacing, n_pad=5)                                                 #nport
#
#c.set_variable('length_right', '5.5mm')
#c.key_elt('end_right', c.capa_right.pos-[c.length_right, 0], [1,0])
#c.end_right.draw_end_cable(c.track, c.gap, typeEnd = 'open')
#
#c.set_variable('length_left', '6.5mm')
#c.key_elt('end_left', c.capa_left.pos-[c.length_left, 0], [1,0])
#c.end_left.draw_end_cable(c.track, c.gap, typeEnd = 'open')
#
#
#c.set_variable('fillet', '0.2mm')
#
#c.connect_elt('in_cable_right', 'capa_right_2', 'in_right')
#c.in_cable_right.draw_cable(fillet=c.fillet)
#c.connect_elt('in_cable_left', 'capa_left_2', 'in_left')
#c.in_cable_left.draw_cable(fillet=c.fillet)
#
#c.set_variable('Lj_down', '0.16nH')
#c.set_variable('Lj_up', '0.138nH')
#
#c.key_elt('squid_right', (c.end_right.pos+c.capa_right.pos)/2, [-1,0])
#c.squid_right.draw_squid(c.track, c.gap, [2*c.track, c.track], c.track, c.gap, iTrackSquid=c.track/4, iTrackJ=c.track/10,  Lj_down=c.Lj_down, Lj_up=c.Lj_up, typePump='down', iSlopePump=0.25)
#
#c.key_elt('squid_left', (c.end_left.pos+c.capa_left.pos)/2, [-1,0])
#c.squid_left.draw_squid(c.track, c.gap, [2*c.track, c.track], c.track, c.gap, iTrackSquid=c.track/4, iTrackJ=c.track/10,  Lj_down=c.Lj_down, Lj_up=c.Lj_up, typePump='up', iSlopePump=0.25)
#
#c.connect_elt('flux_cable_right', 'flux_right', 'squid_right_pump')
#c.flux_cable_right.draw_cable(fillet=c.fillet)
#c.connect_elt('flux_cable_left', 'flux_left', 'squid_left_pump')
#c.flux_cable_left.draw_cable(fillet=c.fillet)
#
#c.set_variable('meander_length_right', '0.56mm')
#c.connect_elt('res1_right', 'capa_right_1', 'squid_right_2')
#c.res1_right.draw_cable(fillet=c.fillet, is_meander=True, to_meander=[1, 1], meander_length=c.meander_length_right)
#c.connect_elt('res2_right', 'end_right', 'squid_right_1')
#c.res2_right.draw_cable(fillet=c.fillet, is_meander=True, to_meander=[1, 1], meander_length=c.meander_length_right)
#
#c.set_variable('meander_length_left', '0.653mm')
#c.connect_elt('res1_left', 'capa_left_1', 'squid_left_2')
#c.res1_left.draw_cable(fillet=c.fillet, is_meander=True, to_meander=[1, 1], meander_length=c.meander_length_left)
#c.connect_elt('res2_left', 'end_left', 'squid_left_1')
#c.res2_left.draw_cable(fillet=c.fillet, is_meander=True, to_meander=[1, 1], meander_length=c.meander_length_left)
#

### Unite - Subtract - PerfE
#Finalisation globale du dessin

#if KeyElt.is_mask:
#    maskObject = c.unite(c.maskObjects, 'mask')
#
#for obj in c.trackObjects:
#    c.assign_perfE(obj)
#gapObject = c.unite(c.gapObjects)
#c.subtract(c.ground_plane, [gapObject])
#
#release()