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
from scripts.designer import Circuit, KeyElt, ConnectElt

project = get_active_project()
design = project.get_active_design()
modeler = design.modeler

modeler.set_units('mm')

c = Circuit(design, modeler)

###########################################################
################# DRAWING STARTS HERE #####################
###########################################################

# Assumptions :
# The 2D circuit is in plane z=0 (can be made more general later
#
#


########## PCB, Chip and Ground Plane ###########
chip_width = c.set_variable("chip_width", "8.67mm")
chip_length = c.set_variable("chip_length", "8.12mm")
chip_thickness = c.set_variable("chip_thickness", "280um")
pcb_thickness = c.set_variable("pcb_thickness", "320um")
vaccuum_thickness = c.set_variable('vaccuum_thickness', 6*chip_thickness)

pcb = c.draw_box('pcb', [0, 0, -chip_thickness], [chip_width, chip_length, -pcb_thickness], "Rogers TMM 10i (tm) no loss")
chip = c.draw_box('chip', [0, 0, 0], [chip_width, chip_length, -chip_thickness], 'silicon')
box = c.draw_box('box', [0, 0, 0], [chip_width, chip_length, vaccuum_thickness], 'vacuum')

ground_plane = c.draw_rect('ground_plane', [0, 0], [chip_width, chip_length])
c.assign_perfE(ground_plane)


######### Connectors #######

con1, con2, con3, con4, con5 = '4.06mm', '3.01mm', '7.54mm', '3.01mm', '6.36mm'
# posistions of connectors: distance (mm) from left or bottom edge (custom parameters)
KeyElt.pcb_track, KeyElt.pcb_gap = parse_entry('300um'), parse_entry('200um')
# Made an attribute of class KeyElt, since it should be the same for all connectors
bond_length, bond_slope = '200um', '0.5'
# bond_length is the size of the region on which we solder the chip parpart of the bond
track, gap = c.set_variable('track', '42um'), c.set_variable('gap', '25um')
# track_width and gap_width of the CPW lines


buffer = KeyElt('buffer_con', [con4, 0], [0, 1])
# the port created along with the connector will be named 'buffer con'
# the connector will be placed at position [con4, 0] and with orientation [1,0] a.k.a. +X
# you can acces all ports names through command : c.available_ports()
# More generally ports are stored in dict c.ports

waste = KeyElt('waste_con', [con2, chip_length], [0, -1])
qubit = KeyElt('qubit_con', [con3, chip_length], [0, -1])
flux = KeyElt('flux_con', [0, con1], [1, 0])
memory = KeyElt('memory_con', [con5, 0], [0, 1])

buffer.draw_connector(2*track, 2*gap, bond_length, bond_slope)
# the connector is given necessary infos here and is drawn

waste.draw_connector(2*track, 2*gap, bond_length, bond_slope)
qubit.draw_connector(2*track, 2*gap, bond_length, bond_slope)
flux.draw_connector(track, gap, bond_length, bond_slope)
memory.draw_connector(track, gap, bond_length, bond_slope)

Lj = c.set_variable('Lj', '12nH')

T_cutout_Length = c.set_variable("T_cutout_Length", '1mm')
T_cutout_Width = c.set_variable("T_cutout_Width", '0.65mm')
T_pad_Width = c.set_variable("T_pad_Width", '0.1mm')
T_pad_Spacing = c.set_variable("T_pad_Spacing", '0.06mm')
T_pad_Length = c.set_variable("T_pad_Length", '0.75mm')
T_offset_width = c.set_variable("T_offset_width", '1.0mm')
T_offset_length = c.set_variable("T_offset_length", '0.0mm')


transmon = KeyElt('transmon', [chip_width/2+T_offset_width, chip_length/2+T_offset_length], [0,1])
transmon.draw_IBM_tansmon([T_cutout_Width, T_cutout_Length],        #cutout_size_vec
                           T_pad_Spacing,                           #pad_Spacing
                           [T_pad_Width, T_pad_Length],             #pad_size_vec
                           '5um',                                   #Jwidth
                           track,                                   #track
                           gap,                                     #gap
                           Lj,                                      #Jinduc                                                         #ori
                           nport=5)                                 #Should be corrected for now doesn't make much sense

capa_pad_spacing = c.set_variable("capa_pad_spacing", '0.05mm')
capa_pad_width = c.set_variable("capa_pad_width", '0.05mm')
capa_pad_length = c.set_variable("capa_pad_length", 2*track)

capa = KeyElt('capa', [con2, chip_length-'1mm'], [0,1])
capa.draw_capa(2*track, 2*gap, capa_pad_spacing, [capa_pad_width, capa_pad_length])                                                       #nport
# there is also a capa defined in ConnectElt and take an input port rather than
# another a position and orientation

adaptor = ConnectElt('adaptor', 'qubit_con', [track, gap])
adaptor.draw_adaptor()
# take and input ports and the track and gap of a new port which will be created

TW_capa_width = c.set_variable("TW_capa_width", '50um')
TB_capa_width = c.set_variable("TB_capa_width", '50um')
TC_capa_width = c.set_variable("TC_capa_width", '50um')
TW_capa_length = c.set_variable("TW_capa_length", 2*track)
TB_capa_length = c.set_variable("TB_capa_length", 2*track)
TC_capa_length = c.set_variable("TC_capa_length", track+'10um')

hcap_waste = ConnectElt('TransmonWasteCapa', 'transmon_1')
hcap_waste.draw_half_capa(TW_capa_length, TW_capa_width, '50um', add_gap=False)
# draw half a capa with or without and added gap with bool add_gap

hcap_buffer = ConnectElt('TransmonBufferCapa', 'transmon_2')
hcap_control = ConnectElt('TransmonControlCapa', 'transmon_3a')
hcap_buffer.draw_half_capa(TB_capa_length, TW_capa_width, '50um')
hcap_control.draw_half_capa(TC_capa_length, TW_capa_width, '5um')

fillet = c.set_variable('fillet', '0.15mm')

cable1 = ConnectElt('Control_cable', 'qubit_con_bis', 'transmon_3a')
cable1.draw_cable(fillet=fillet, is_meander=True, to_meander=[0, 1, 0, 1, 0, 1], meander_length='0.5mm') #for now does not support 2 meander in a row
# semi-smart connect elt which connects two ports.
# fillet is the curvature radius
# is_meander is a bool stating whether there should be meanders rather than a direct connection
# to_meander indicates what are the segments of the direct connection cable that need to be meandered
#       for now, only supports not consecutive meanders

cable2 = ConnectElt('capa_cable', 'waste_con', 'capa_1')
cable2.draw_cable(fillet=fillet, is_meander=False)
cable4 = ConnectElt('Waste_cable', 'capa_2', 'transmon_1')
cable4.draw_cable(fillet=fillet)





squid = KeyElt('squid', [(con4+chip_width/2+T_offset_width)/2,'2mm'], [1,0])
squid.draw_squid(track, gap, [2*track, track], 2*track, 2*gap, iTrackSquid=track/4, iTrackJ=track/10,  Lj_down='1nH', Lj_up='3nH', typePump='down', doublePump=False)

# first track and gap correspond to the track and gap of the line in which the squid is embedded
# the list tells the inside size of the rectangular loop forming the squid
# then the track and gap of the input port for current
# if Lj_up not mentionned draws a symetrical squid
# typePump is the position of the driving line
# doublePump is a boolean stating whether there should be two driving ports

cable4 = ConnectElt('pump_cable', 'memory_con', 'squid_pump1')
cable4.draw_cable(fillet=fillet)
cable5 = ConnectElt('bef_squid_cable', 'buffer_con', 'squid_2')
cable5.draw_cable(fillet=fillet, is_meander=True)
cable6 = ConnectElt('aft_squid_cable', 'transmon_2', 'squid_1')
cable6.draw_cable(fillet=fillet, is_meander=True)


end = KeyElt('end', ['2mm',con1], [-1,0])
end.draw_end_cable(track, gap, typeEnd = 'short')

cable7 = ConnectElt('end_cable', 'flux_con', 'end')
cable7.draw_cable(fillet=fillet)
# track and gap are the characteristic of the created port
# typeEnd is 'short' or 'open'




######## UNITE - SUBSTRACT - PerfE #####
#Finalisation globale du dessin

for obj in c.trackObjects:
    c.assign_perfE(obj)
gapObject = c.unite(c.gapObjects)
c.subtract(ground_plane, [gapObject])

release()