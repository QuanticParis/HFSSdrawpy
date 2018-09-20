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



from scripts.hfss  import get_active_project, release, parse_entry
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

############
is_mask=True

########## PCB, Chip and Ground Plane ###########
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






######### Connectors #######
con1, con2, con3, con4, con5 = '4.06mm', '3.01mm', '7.54mm', '3.01mm', '6.36mm'
#con1, con2, con3, con4, con5 = '3.06mm', '3.01mm', '7.54mm', '3.01mm', '6.36mm'

# posistions of connectors: distance (mm) from left or bottom edge (custom parameters)
KeyElt.pcb_track, KeyElt.pcb_gap = parse_entry('300um'), parse_entry('200um')
# Made an attribute of class KeyElt, since it should be the same for all connectors
bond_length, bond_slope = '200um', '0.5'
# bond_length is the size of the region on which we solder the chip parpart of the bond
c.set_variable('track', '42um')
c.set_variable('gap', '25um')
# track_width and gap_width of the CPW lines


######### DEFINE PORTS
# the port created along with the connector will be named 'buffer con'
# the connector will be placed at position [con4, 0] and with orientation [1,0] a.k.a. +X
# you can acces all ports names through command : c.available_ports()
# More generally ports are stored in dict c.ports
c.key_elt('buffer_con', [con4, 0], [0, 1])
c.key_elt('waste_con', [con2, c.chip_length], [0, -1])
#c.key_elt('qubit_con', [con3, c.chip_length], [0, -1])
c.key_elt('flux_con', [con5, 0], [0, 1])
c.key_elt('qubit_con', [0, con1], [1, 0])


######### DRAW CONNECTORS
# the connector is given necessary infos here and is drawn
c.buffer_con.draw_connector(c.track, c.gap, bond_length, bond_slope)
c.waste_con.draw_connector(c.track, c.gap, bond_length, bond_slope)
c.qubit_con.draw_connector(c.track/2, c.gap/2, bond_length, bond_slope)
c.flux_con.draw_connector(c.track/2, c.gap/2, bond_length, bond_slope)
#c.memory_con.draw_connector(c.track, c.gap, bond_length, bond_slope)


##########DRAW TRANSMON QUBIT IN THE IBM STYLE
c.set_variable('Lj', '12nH')
c.set_variable("T_cutout_Length", '1mm')
c.set_variable("T_cutout_Width", '0.65mm')
c.set_variable("T_pad_Width", '0.1mm')
c.set_variable("T_pad_Spacing", '0.06mm')
c.set_variable("T_pad_Length", '0.8mm')
c.set_variable("T_offset_width", '1.0mm')
c.set_variable("T_offset_length", '0.0mm')
c.set_variable("T_fillet", '0.015mm')


c.key_elt('transmon', [c.chip_width/2+c.T_offset_width, c.chip_length/2+c.T_offset_length], [0,-1])
c.transmon.draw_IBM_tansmon([c.T_cutout_Width, c.T_cutout_Length],     #cutout_size_vec
                           c.T_pad_Spacing,                            #pad_Spacing
                           [c.T_pad_Width, c.T_pad_Length],            #pad_size_vec
                           '5um',                                      #Jwidth
                           c.track,                                    #track
                           c.gap,                                      #gap
                           c.Lj,                                       #Jinduc                                                         #ori
                           nport=5,
                           fillet=c.T_fillet)
#if is_mask:
#    c.set_variable("markT_size", '0.02mm')
#    c.set_variable("markT_offset_x", '0.38mm')
#    c.set_variable("markT_offset_y", '0.15mm')
#    c.transmon.draw_alignement_mark(c.mark_size, c.mark_offset_x, c.mark_offset_y)
#    c.transmon.draw_alignement_mark(c.mark_size, -c.mark_offset_x, c.mark_offset_y)
#    c.transmon.draw_alignement_mark(c.mark_size, c.mark_offset_x, -c.mark_offset_y)
#    c.transmon.draw_alignement_mark(c.mark_size, -c.mark_offset_x, -c.mark_offset_y)
#
#if is_mask:
#    c.set_variable("T_dose_cutout_Length", '0.15mm')
#    c.set_variable("T_dose_cutout_Width", '0.25mm')
#    c.set_variable("T_dose_pad_Width", '0.1mm')
#    c.set_variable("T_dose_pad_Length", '0.1mm')
#    c.set_variable("T_dose_spacing_Width", '0.2mm')
#    c.set_variable("T_dose_spacing_Length", '-0.35mm')
#    c.set_variable("T_dose_offset_width", '3.0mm')
#    c.set_variable("T_dose_offset_length", '0.0mm')
#    for i in range(3):
#        for j in range(3):
#            transmon_dose_temp=c.key_elt('dose_transmon'+str(i)+str(j), [c.chip_width/2+c.T_dose_offset_width+i*c.T_dose_spacing_Width, c.chip_length/2+c.T_dose_offset_length+j*c.T_dose_spacing_Length], [0,-1])
#            transmon_dose_temp.draw_IBM_tansmon([c.T_dose_cutout_Width, c.T_dose_cutout_Length],        #cutout_size_vec
#                               c.T_pad_Spacing,                           #pad_Spacing
#                               [c.T_dose_pad_Width, c.T_dose_pad_Length],             #pad_size_vec
#                               '5um',                                   #Jwidth
#                               c.track,                                   #track
#                               c.gap,                                     #gap
#                               c.Lj,                                      #Jinduc                                                         #ori
#                               nport=1,
#                               fillet=c.T_fillet,
#                               is_mask=is_mask)  

                               #Should be corrected for now doesn't make much sense
##########DRAW CAPA IN THE QUBIT CUTOUT
c.set_variable("TW_capa_width", '180um')
c.set_variable("TB_capa_width", '84um')
c.set_variable("TC_capa_width", '22um')
c.set_variable("TW_capa_length", '130um')
c.set_variable("TB_capa_length", '90um')
c.set_variable("TC_capa_length", '5um')
c.set_variable("TW_capa_thickness", '50um')
c.set_variable("TB_capa_thickness", '50um')
c.set_variable("TC_capa_thickness", '5um')

# draw half a capa with or without and added gap with bool add_gap
c.connect_elt('TW_Capa', 'transmon_2')
c.connect_elt('TB_Capa', 'transmon_1')
c.connect_elt('TC_Capa', 'transmon_3a')
c.set_variable("HalfCapa_fillet", '0.01mm')

c.TW_Capa.draw_half_capa( c.TW_capa_width, c.TW_capa_thickness, c.TW_capa_length, add_gap=False,fillet=c.HalfCapa_fillet)
c.TB_Capa.draw_half_capa( c.TB_capa_width, c.TB_capa_thickness, c.TB_capa_length,fillet=c.HalfCapa_fillet)
c.TC_Capa.draw_half_capa( c.TC_capa_width, c.TC_capa_thickness, c.TC_capa_length)



#######CREATE CONTROL LINE
c.set_variable('fillet', '0.15mm')
c.set_variable('fillet_control', '0.05mm')

c.connect_elt('Control_cable', 'qubit_con', 'transmon_3a')
c.Control_cable.draw_cable(fillet=c.fillet_control), #is_meander=True, to_meander=[0, 1, 0, 1, 0, 1], meander_length='0.5mm') #for now does not support 2 meander in a row
# semi-smart connect elt which connects two ports.
# fillet is the curvature radius
# is_meander is a bool stating whether there should be meanders rather than a direct connection
# to_meander indicates what are the segments of the direct connection cable that need to be meandered
#       for now, only supports not consecutive meanders

######DEFINE WASTE PURCELL CAPA
#c.set_variable("capaPurW_spacing", '0.005mm')
#c.set_variable("capaPurW_width",  0.5*c.track)
#c.set_variable("capaPurW_length", '0.3mm')
#c.key_elt('capaPurW', [con2, c.chip_length-'0.8mm'], [0,1])
#c.capaPurW.draw_capa(c.track, c.gap, c.capaPurW_spacing, [c.capaPurW_width, c.capaPurW_length])                                                       #nport
## there is also a capa defined in ConnectElt and take an input port rather than
## another a position and orientation

c.set_variable("capaPurW_teethwidth", '0.02mm')
c.set_variable("capaPurW_teethlength",  '0.045mm')
c.set_variable("cappaPurW_fillet", '0.01mm')
c.set_variable("cappaPurW_gap", '0.003mm')


c.key_elt('capaPurW', [con2, c.chip_length-'0.8mm'], [0,1])

c.capaPurW.draw_capa_interdigitated(c.track, c.gap, [c.capaPurW_teethlength, c.capaPurW_teethwidth],c.cappaPurW_gap,3,c.cappaPurW_fillet)                                                       #nport



######CONNECT WASTE CAPA to CONNECTOR
c.connect_elt('cable_Wcon_capa', 'waste_con', 'capaPurW_1')
c.cable_Wcon_capa.draw_cable(fillet=c.fillet, is_meander=False)

#####CREATE WASTE CAPA
c.set_variable("capaW_spacing", '0.035mm')
c.set_variable("capaW_width",  0.5*c.track)
c.set_variable("capaW_length", c.track)
c.key_elt('capaW', (c.capaPurW.pos+c.TW_Capa.pos)/2, [1,0])
c.capaW.draw_capa(c.track, c.gap, c.capaW_spacing, [c.capaW_width, c.capaW_length])                                                       #nport



######CREATE WASTE RESONATOR
c.set_variable('meander_length_W','1.355mm')
c.connect_elt('res_Waste', 'capaW_1', 'transmon_2')
c.res_Waste.draw_cable(fillet=c.fillet, is_meander=True, to_meander=[0, -1, 0], meander_length=c.meander_length_W, is_mesh=True)

######CREATE PURCELL RESONATOR
c.set_variable('meander_length_PurW','1.29mm')
c.connect_elt('resPurcell_Waste', 'capaW_2', 'capaPurW_2')
c.resPurcell_Waste.draw_cable(fillet=c.fillet, is_meander=True, to_meander=[0, -1, 0], meander_length=c.meander_length_PurW, is_mesh=True)



######DEFINE BUFFER CAPA
c.set_variable("capaB_spacing", '0.015mm')
c.set_variable("capaB_width",  0.5*c.track)
c.set_variable("capaB_length", '0.2mm')
c.key_elt('capaB', [con4, '0.8mm'], [0,1])
c.capaB.draw_capa(c.track, c.gap, c.capaB_spacing, [c.capaB_width, c.capaB_length])                                                       #nport

######CONNECT BUFFER CAPA to CONNECTOR
c.connect_elt('cable_Bcon_capa', 'buffer_con', 'capaB_2')
c.cable_Bcon_capa.draw_cable(fillet=c.fillet, is_meander=False)

c.set_variable("L_SQUID", '0.1nH')


######DEFINE SQUID
c.key_elt('squid', (c.capaB.pos+c.TB_Capa.pos)/2, [1,0])
c.squid.draw_squid(c.track, c.gap, [2*c.track, c.track], c.track/2, c.gap/2, iTrackSquid=c.track/4, iTrackJ=c.track/10,  Lj_down=c.L_SQUID, Lj_up=c.L_SQUID, typePump='down', doublePump=False, fillet=None)#c.track/2)
# first track and gap correspond to the track and gap of the line in which the squid is embedded
# the list tells the inside size of the rectangular loop forming the squid
# then the track and gap of the input port for current
# if Lj_up not mentionned draws a symetrical squid
# typePump is the position of the driving line
# doublePump is a boolean stating whether there should be two driving ports

#if is_mask:
#    c.set_variable("markSQUID_size", '0.02mm')
#    c.set_variable("markSQUID_offset_x", '0.38mm')
#    c.set_variable("markSQUID_offset_y", '0.15mm')
#    c.squid.draw_alignement_mark(c.markSQUID_size, c.markSQUID_offset_x, c.markSQUID_offset_y)
#    c.squid.draw_alignement_mark(c.markSQUID_size, -c.markSQUID_offset_x, c.markSQUID_offset_y)
#    c.squid.draw_alignement_mark(c.markSQUID_size, c.markSQUID_offset_x, -c.markSQUID_offset_y)
#    c.squid.draw_alignement_mark(c.markSQUID_size, -c.markSQUID_offset_x, -c.markSQUID_offset_y)


######CREATE BUFFER RESONATOR
c.set_variable('meander_length_B','0.535mm')
c.connect_elt('res_Buffer_bef', 'capaB_1', 'squid_2')
c.res_Buffer_bef.draw_cable(fillet=c.fillet, is_meander=True, to_meander=[-1, 0], meander_length=c.meander_length_B, is_mesh=True)
c.connect_elt('res_Buffer_aft', 'transmon_1', 'squid_1')
c.res_Buffer_aft.draw_cable(fillet=c.fillet, is_meander=True, to_meander=[-1, 0], meander_length=c.meander_length_B, is_mesh=True)

######CREATE PUMP LINE
c.connect_elt('pump_line', 'flux_con', 'squid_pump')
c.pump_line.draw_cable(fillet=c.fillet)






######## UNITE - SUBSTRACT - PerfE #####
#Finalisation globale du dessin

for obj in c.trackObjects:
    c.assign_perfE(obj)
gapObject = c.unite(c.gapObjects)
c.subtract(c.ground_plane, [gapObject])

#if is_mask:
#    #draw_fluxtrappingholes(self, chip_length,chip_width,margin,hole_spacing,hole_size):
#    #c.key_elt('holes',[0,0],[0,1])
#    #holes=c.holes.draw_fluxtrappingholes(c.chip_length,c.chip_width,"0.5mm","0.5mm","10um")
#
#    maskObjects = c.unite(c.maskObjects)
#    #c.subtract(holes,[maskObjects])

release()