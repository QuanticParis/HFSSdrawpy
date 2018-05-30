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

project = get_active_project()
design = project.get_active_design()
modeler = design.modeler

modeler.set_units('mm')

c = Circuit(design, modeler)

KeyElt.is_mask = False
KeyElt.gap_mask = parse_entry('20um')
KeyElt.is_overdev = False
KeyElt.overdev = parse_entry('0um')

#######################
### DRAWING STARTS HERE 
#######################

### PCB, Chip and Ground Plane

c.set_variable("chip_width", "8.67mm")
c.set_variable("chip_length", "8.12mm")
c.set_variable("chip_thickness", "280um")
c.set_variable("pcb_thickness", "320um")
c.set_variable('vaccuum_thickness', 6*c.chip_thickness)

c.draw_box('pcb', [0, 0, -c.chip_thickness], [c.chip_width, c.chip_length, -c.pcb_thickness], "Rogers TMM 10i (tm) no loss")
c.draw_box('chip', [0, 0, 0], [c.chip_width, c.chip_length, -c.chip_thickness], 'silicon')
c.draw_box('box', [0, 0, 0], [c.chip_width, c.chip_length, c.vaccuum_thickness], 'vacuum')

c.draw_rect('ground_plane', [0, 0], [c.chip_width, c.chip_length])
c.assign_perfE(c.ground_plane)


### Connectors

con1, con2, con3, con4, con5 = '4.06mm', '3.01mm', '7.54mm', '3.01mm', '6.36mm'
KeyElt.pcb_track, KeyElt.pcb_gap = parse_entry('300um'), parse_entry('200um')
bond_length, bond_slope = '200um', '0.5'
c.set_variable('track', '84um')
c.set_variable('gap', '50um')
c.set_variable('track_mem', '3um')
c.set_variable('gap_mem', '50um')
c.set_variable('track_flux', '50um')
c.set_variable('gap_flux', '30um')

c.key_elt('in_mem', [con3, c.chip_length], [0, -1])
c.key_elt('in_flux_top', [con2, c.chip_length], [0, -1])
c.key_elt('in_buffer', [0, con1], [1, 0])
c.key_elt('in_readout', [con5, 0], [0, 1])
c.key_elt('in_flux_bot', [con4, 0], [0, 1])

c.in_mem.draw_connector(c.track, c.gap, bond_length, bond_slope)
c.in_flux_top.draw_connector(c.track_flux, c.gap_flux, bond_length, bond_slope)
c.in_buffer.draw_connector(c.track, c.gap, bond_length, bond_slope)
c.in_readout.draw_connector(c.track, c.gap, bond_length, bond_slope)
c.in_flux_bot.draw_connector(c.track_flux, c.gap_flux, bond_length, bond_slope)

### Key ELements
c.set_variable('Lj', '12nH')
c.key_elt('trm', ['5mm','5.5mm'], [1,0])
c.trm.draw_ZR_transmon(['1.52mm', '0.8mm'], '0.12mm', ['0.5mm', '0.5mm'], '84um', '50um', '0.2mm', '30um', '0mm', '0.01mm', c.Lj, pad_size_left=['0.6mm', '0.5mm'], track_left='3um', gap_left='50um', length_left='0.2mm', spacing_left='50um', short_left='10um', fillet=True)#, pad_size_right=['0.5mm', '0.5mm'], track_right='84um', gap_right='50um', length_right'0.2mm', spacing_right='30um', short_right='0mm', 

c.set_variable('x_T', '3.8mm')
c.set_variable('y_T', '5.5mm')
c.key_elt('T', [c.x_T,c.y_T], [1,0])
c.T.draw_T(c.track_mem, c.gap_mem)
#

c.set_variable('squid_width', '50um')
c.set_variable('squid_length', '50um')
c.key_elt('squid', [con2+c.squid_width/2+c.track_mem, '5.5mm'], [-1,0])
c.squid.draw_squid_protect(c.track_mem, c.gap_mem, [c.squid_width, c.squid_length], c.track_flux, c.gap_flux, iTrackSquid=c.track_mem, iTrackJ=None, Lj_down='20nH', Lj_up=None,  typePump='down', doublePump=True, iSlope=1, iSlopePump=0.5, fillet=None)


c.key_elt('end_mem', ['4.8mm','2.5mm'], [0,-1])
c.end_mem.draw_end_cable(c.track_mem, c.gap_mem, typeEnd='short', fillet=True)
#
c.key_elt('end_buffer', ['1.71mm','3.1mm'], [0,1])
c.end_buffer.draw_end_cable(c.track_mem, c.gap_mem, typeEnd='short', fillet=True)
#

c.key_elt('end_drive_buffer', c.squid.pos+Vector(['-0.15mm','0.08mm']), [-1,0])
c.end_drive_buffer.draw_end_cable(c.track, c.gap, typeEnd='open', fillet=True)

c.key_elt('end_drive_mem', [c.x_T-'0.4mm',c.y_T+'0.25mm'], [0,1])
c.end_drive_mem.draw_end_cable(c.track, c.gap, typeEnd='open', fillet=True)

c.set_variable('gap_capa_readout', '10um')  
c.key_elt('capa_readout', [con5,'0.9mm'], [0,1])
c.capa_readout.draw_capa_inline(c.track, c.gap, '80um', c.gap_capa_readout , n_pad=1)

### Connect ELements

c.connect_elt('squid_T', 'squid_2', 'T_2')
c.squid_T.draw_cable(is_bond=True)

c.connect_elt('T_qubit', 'T_1', 'trm_2')
c.T_qubit.draw_cable(is_bond=True)

c.connect_elt('flux_top', 'squid_pump1', 'in_flux_top')
c.flux_top.draw_cable(is_bond=True)

c.connect_elt('flux_bot', 'squid_pump2', 'in_flux_bot')
c.flux_bot.draw_cable(is_bond=True)

c.set_variable('fillet', '500um')
c.connect_elt('buffer', 'squid_1', 'end_buffer')
c.buffer.draw_cable(is_bond=True, fillet=c.fillet)
c.connect_elt('mem', 'T_3', 'end_mem')
c.mem.draw_cable(is_bond=True, fillet=c.fillet)

c.connect_elt('drive_mem', 'in_mem', 'end_drive_mem')
c.drive_mem.draw_cable(is_bond=True, fillet=c.fillet)

c.set_variable('constrain_ro', '1.932mm')
c.key_elt('constr_readout', c.in_readout.pos+Vector([c.constrain_ro, '4mm']), [0,1])
c.constr_readout.create_port(c.track, c.gap)
c.connect_elt('readout', 'capa_readout_1', 'trm_1')
c.readout.draw_cable(is_bond = True, fillet=c.fillet, constrains=['constr_readout'])
c.connect_elt('bef_capa', 'capa_readout_2', 'in_readout')
c.bef_capa.draw_cable(is_bond = True, fillet=c.fillet)

c.key_elt('constr_drive_buffer', c.in_buffer.pos+Vector(['1.26mm', '0.75mm']), [0,1])
c.constr_drive_buffer.create_port(c.track, c.gap)
c.connect_elt('drive_buffer', 'in_buffer', 'end_drive_buffer')
c.drive_buffer.draw_cable(is_bond = True, fillet=c.fillet, constrains=['constr_drive_buffer'])

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


### Unite - Subtract - PerfE
#Finalisation globale du dessin

if KeyElt.is_mask:
    maskObject = c.unite(c.maskObjects, 'mask')

for obj in c.trackObjects:
    c.assign_perfE(obj)
gapObject = c.unite(c.gapObjects)
c.subtract(c.ground_plane, [gapObject])

release()