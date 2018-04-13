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
c.assign_perfE(box)

ground_plane = c.draw_rect('ground_plane', [0, 0], [chip_width, chip_length])


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

#waste.draw_connector(2*track, 2*gap, bond_length, bond_slope)
#qubit.draw_connector(2*track, 2*gap, bond_length, bond_slope)
##flux.draw_connector(track, gap, bond_length, bond_slope)
#memory.draw_connector(track, gap, bond_length, bond_slope)

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
capa_pad_length = c.set_variable("capa_pad_length", '0.3mm')

#capa = KeyElt('capa', [con2, chip_length-'1mm'], [0,1])
#capa.draw_capa(2*track, 2*gap, capa_pad_spacing, [capa_pad_width, capa_pad_length])                                                       #nport
# there is also a capa defined in ConnectElt and take an input port rather than
# another a position and orientation

#adaptor = ConnectElt('adaptor', 'qubit_con', [track, gap])
#adaptor.draw_adaptor()
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

#cable1 = ConnectElt('Control_cable', 'qubit_con', 'transmon_3a')
#cable1.draw_cable(fillet=fillet, is_meander = True, to_meander = [0, 1, 1, 0, 1, 0], meander_length='0.5mm') #for now does not support 2 meander in a row
# semi-smart connect elt which connects two ports.
# fillet is the curvature radius
# is_meander is a bool stating whether there should be meanders rather than a direct connection
# to_meander indicates what are the segments of the direct connection cable that need to be meandered
#       for now, only supports not consecutive meanders

#cable2 = ConnectElt('Qubit_cable', 'qubit_con', 'transmon_1')
#cable2.draw_cable(iMaxfillet=fillet, is_meander = True)
#cable3 = ConnectElt('Buffer_cable', 'buffer_con', 'transmon_2')
#cable3.draw_cable(iMaxfillet=fillet)

#cable5 = ConnectElt('Waste_cable', 'capa_2', 'transmon_2')
#cable5.draw_cable(iMaxfillet=fillet, is_meander=True, to_meander = [1, 0, 1, 0, 1, 1, 0], meander_length='0.5mm')



squid = KeyElt('squid', [0,0], [1,0])
squid.draw_squid(track, gap, [2*track, track], 2*track, 2*gap, iTrackSquid=track/4, iTrackJ=track/10,  Lj_down='1nH', Lj_up='3nH', typePump='down', doublePump=True)
# first track and gap correspond to the track and gap of the line in which the squid is embedded
# the list tells the inside size of the rectangular loop forming the squid
# then the track and gap of the input port for current
# if Lj_up not mentionned draws a symetrical squid
# typePump is the position of the driving line
# doublePump is a boolean stating whether there should be two driving ports

end = KeyElt('end', ['-1mm','0.5mm'], [1,0])
end.draw_end_cable(track, gap, typeEnd = 'short')
# track and gap are the characteristic of the created port
# typeEnd is 'short' or 'open'


cable4 = ConnectElt('squid_cable', 'buffer_con', 'squid_pump2')
cable4.draw_cable(fillet=fillet, is_meander=False)



cable6 = ConnectElt('end_cable', 'end', 'squid_2')
cable6.draw_cable(fillet=fillet, is_meander=False)

####################################################################################################################
####################################################################################################################
####################################################################################################################
####################################################################################################################
####################################################################################################################

#
#
#TW_capa_width=LitExp("$TW_capa_width", 2*track,VarDef=True)
#TB_capa_width=LitExp("$TB_capa_width", 2*track,VarDef=True)
#TC_capa_width=LitExp("$TC_capa_width", track+0.01,VarDef=True)
#
#TW_capa_length=LitExp("$TW_capa_length", 0.05,VarDef=True)
#TB_capa_length=LitExp("$TB_capa_length", 0.05,VarDef=True)
#TC_capa_length=LitExp("$TC_capa_length", 0.005,VarDef=True)
#
#drawHalfCapa("TransmonWasteCapa",TransmonIn_W,TW_capa_width,0.05,TW_capa_length)
#drawHalfCapa("TransmonBufferCapa",TransmonIn_B,TB_capa_width,0.05,TB_capa_length)
#drawHalfCapa("TransmonControlCapa",TransmonIn_C,TC_capa_width,0.05,TC_capa_length)
#
#
#WasteCapaLength, WasteCapaSize = LitExp("$WasteCapa_length", 0.1,VarDef=True), LitExp("$WasteCapa_size", 0.02,VarDef=True)
#BufferCapaLength, BufferCapaSize = LitExp("$BufferCapa_length", 0.1,VarDef=True), LitExp("$BufferCapa_size", 0.02,VarDef=True)
#
#WasteCapaIn, WasteCapaOut = WasteOut, [None, None, track, gap]
#BufferCapaIn, BufferCapaOut = BufferOut, [None, None, track, gap]
#
#WasteCapaIn, WasteCapaOut = drawCapa("WasteCapa", WasteCapaIn, WasteCapaOut, WasteCapaLength, track/2, WasteCapaSize)
#BufferCapaIn, BufferCapaOut = drawCapa("BufferCapa", BufferCapaIn, BufferCapaOut, BufferCapaLength, track/2, BufferCapaSize)
#
#W_resLength=LitExp("$W_resLength", 9,VarDef=True)
#W_resShift=LitExp("$W_resShift", 0.,VarDef=True)
#
#B_resLength=LitExp("$B_resLength", 9.5,VarDef=True)
#B_resShift=LitExp("$B_resShift", 0.,VarDef=True)
#
#
#
#
##drawWave2("BufferMeander", BufferCapaOut, TransmonIn_B, 3, 0.35,B_resShift, Lres=B_resLength, iStraight=0.3)
##Number of meanders, Shift, total length (shift is ignored if not None)
#
#
#pos1, ori = Vector(BufferCapaOut[POS]), Vector(BufferCapaOut[ORI])
#pos2 = Vector(TransmonIn_B[POS])
#_, _, inTrack, inGap = BufferCapaOut
#BufferSQUID_In = [0.5*(pos1+pos2), ori, inTrack, inGap]
#
#LSQUID=LitExp("$LSQUID", 0.1,VarDef=True)
#BufferSQUID_In,BufferSQUID_Out,FluxLine_In,_=drawSQUID("BufferSQUID", BufferSQUID_In,BufferSQUID_In,0.01,0.01,0.05,0.01,iInduct=LSQUID)
##(iName, iIn, iOut, iSize, iWidth, iLength, iSpace, iInduct=0.1)
#
#pos, ori = Vector(FluxLine_In[POS]), Vector(FluxLine_In[ORI])
#FluxLine_In = [pos, ori, track/4., gap/4.]
#FluxLine_In=drawFlux("fluxline", FluxLine_In, 0.02, 0.02)
## drawFlux(iName, iIn, iSize, iMargin, iGndGap=None):
#
#pos, ori = Vector(FluxLine_In[POS]), Vector(FluxLine_In[ORI])
#_, _, inTrack, inGap = FluxLine_In
#FluxLine_Out = [pos, ori, track, gap]
#
#FluxLine_In, FluxLine_Out = drawAdaptor("fluxline_adapt",FluxLine_In,FluxLine_Out)
##drawAdaptor(iName, iIn, iOut, iSlope=1):
#
#drawCable("Flux_cable", FluxLine_Out, memoryOut, extrabonds=3)
#
#drawWave2("BufferMeander1", BufferCapaOut,BufferSQUID_In, 1, 0.35,B_resShift, Lres=B_resLength/2., iStraight=0.3)
#drawWave2("BufferMeander2",  TransmonIn_B, BufferSQUID_Out, 1, 0.35,B_resShift, Lres=B_resLength/2., iStraight=0.3)
#
#drawWave2("WasteMeander", WasteCapaOut, TransmonIn_W, 3, 0.35, W_resShift,Lres=W_resLength, iStraight=0.3)









######## Connectors #######
# posistions of connectors: distance (mm) from left or bottom edge
# con1, con2, con3, con4, con5 = 4.06, 3.01, 7.537, 3.01, 6.36

# chipIn, chipOut = [[0, con1], [1, 0], pcbTrack, pcbGap], [None, None, track, gap]
# tomoIn, tomoOut = [[con2, chip_length], [0, -1], pcbTrack, pcbGap], [None, None, track, gap]
# qubitIn, qubitOut = [[con3, chip_length], [0, -1], pcbTrack, pcbGap], [None, None, track, gap]
# fluxIn, fluxOut = [[con4, 0], [0, 1], pcbTrack, pcbGap], [None, None, track, gap]
# memoryIn, memoryOut = [[con5, 0], [0, 1], pcbTrack, pcbGap], [None, None, track, gap]

# chipIn, chipOut = drawConnector("chip_con", chipIn, chipOut, BondLength, BondSlope)
# tomoIn, tomoOut = drawConnector("tomo_con", tomoIn, tomoOut, BondLength, BondSlope)
# qubitIn, qubitOut = drawConnector("qubit_con", qubitIn, qubitOut, BondLength, BondSlope)
# fluxIn, fluxOut = drawConnector("flux_con", fluxIn, fluxOut, BondLength, BondSlope)
# memoryIn, memoryOut = drawConnector("memory_con", memoryIn, memoryOut, BondLength, BondSlope)



# ######## Drawing of the Transmon #######
# # comment: this transmon is non-standard, it resembles in Xmon, and is
# # formed by two T junctions

# crossGap = LitExp("$cross_gap", 0.02)
# crossSize = LitExp("$cross_size", 0.3)

# # ports
# cross1In = [[memoryIn[POS][0]-crossGap-track/2, chipIn[POS][1]], [1, 0], track, crossGap]

# cross1OutRight, cross1OutLeft = [None, None, track, crossGap], [None, None, track, crossGap]
# cross1OutUp, cross1OutRight, cross1OutLeft = drawTriJunction("cross1", cross1In, cross1OutRight, cross1OutLeft)

# crossJSJuncIn, crossJSJuncOut = [cross1OutRight[POS], [0, 1], track, crossGap], [None, None, track, crossGap]
# crossJSJuncIn, crossJSJuncOut = drawJSJunc("memory_jsjunc", crossJSJuncIn, crossJSJuncOut, 0.02, 0.01, 0.02, 8.2)

# cross2In = [crossJSJuncOut[POS]+Vector(crossGap+track/2, crossGap+track/2), [-1, 0], track, crossGap]
# cross2OutLeft, cross2OutRight = [None, None, track, crossGap], [None, None, track, crossGap]
# cross2OutDown, cross2OutLeft, cross2OutRight = drawTriJunction("cross2", cross2In, cross2OutLeft, cross2OutRight)

# ####### BUFFER RESONATOR ##############
# # Dessin de la piste du buffer (buffer)

# chipCapaPos, chipWaveLength = LitExp("$chip_capa_pos", 0.6), LitExp("$chip_wave_length", 0.52)
# chipCapaLength, chipCapaSize = LitExp("$chip_capa_length", 0.13), LitExp("$chip_capa_size", 0.05)
# chipCapa2Length, chipCapa2Size = LitExp("$chip_capa2_length", 0.16), LitExp("$chip_capa2_size", 0.05)

# chipCapa1In, chipCapa1Out = [chipOut[POS]+Vector(chipCapaPos, 0), [1, 0], track, gap], [None, None, track, gap]
# chipCapa1In, chipCapa1Out = drawCapa("chip_capa1", chipCapa1In, chipCapa1Out, chipCapaLength, track/2, chipCapaSize)

# chipCapa2In, chipCapa2Out = [cross1OutUp[POS]+Vector(-crossSize, 0), [-1, 0], track, crossGap], [None, None, track, gap]
# chipCapa2In, chipCapa2Out = drawCapa("chip_capa2", chipCapa2In, chipCapa2Out, chipCapa2Length, track/2, chipCapa2Size)


# chipSquideIn, chipSquideOut = [chipCapa1Out[POS]*0.5+chipCapa2Out[POS]*0.5+Vector(-0.052, 0), [1, 0], track, gap], [None, None, track, gap]
# chipSquideIn, chipSquideOut, fluxFluxIn, _ = drawSquide("chip_squide", chipSquideIn, chipSquideOut, 0.02, 0.01, 0.02, 0.02, 2.62)


# chipWave1In, chipWave1Out = [chipCapa1Out[POS], [1, 0], track, gap], [None, None, track, gap]
# chipWave1In, chipWave1Out = drawWave("chip1_wave", chipWave1In, chipWave1Out, 1, chipWaveLength, chipWaveLength, 0.4)

# chipWave2In, chipWave2Out = [chipCapa2Out[POS], [-1, 0], track, gap], [None, None, track, gap]
# chipWave2In, chipWave2Out = drawWave("chip2_wave", chipWave2In, chipWave2Out, 1, chipWaveLength, chipWaveLength, 0.4)

# drawCable("chip_cable1", chipOut, chipCapa1In)
# drawCable("chip_cable2", chipWave1Out, chipSquideIn)
# drawCable("chip_cable3", chipSquideOut, chipWave2Out)
# drawCable("chip_cable4", chipCapa2In, cross1OutUp)

# ####### DIRECT INPUT LINE FOR THE TRANSMON #####
# #Dessin du la piste du transmon (transmon-out)

# qubitCapaLength, qubitCapaSize = LitExp("$qubit_capa_length", 0.11), LitExp("$qubit_capa_size", 0.31)

# qubitCapaIn, qubitCapaOut = [cross2OutDown[POS]+Vector(crossSize, 0), [1, 0], track, crossGap], [None, None, track, gap]
# qubitCapaIn, qubitCapaOut = drawCapa("qubit_capa", qubitCapaIn, qubitCapaOut, qubitCapaLength, track, qubitCapaSize)


# qubitDoubleJuncOut = [qubitCapaOut[POS]+Vector(0.5, 0.4), [0, 1], track, gap]
# qubitDoubleJuncIn = [qubitDoubleJuncOut[POS], [0, -1], track, gap]

# drawCable("qubit_cable1", qubitDoubleJuncOut, qubitOut)
# drawCable("qubit_cable2", qubitCapaOut, qubitDoubleJuncIn)
# drawCable("qubit_cable3", qubitCapaIn, cross2OutDown)

# ######## READOUT RESONATOR #############
# #Dessin de la piste de la tomographie (read-out)

# tomoDecal, tomoCapaPos = LitExp("$tomo_decal", 2.32), LitExp("$tomo_capa_pos", 0.2)
# tomoCapaLength, tomoCapaSize = LitExp("$tomo_capa_length", 0.13), LitExp("$tomo_capa_size", 0.05)
# tomoCapa2Length, tomoCapa2Size = LitExp("$tomo_capa2_length", 0.16), LitExp("$tomo_capa2_size", 0.05)

# tomoCapa1In, tomoCapa1Out = [tomoOut[POS]+Vector(0, -tomoCapaPos), [0, -1], track, gap], [None, None, track, gap]
# tomoCapa1In, tomoCapa1Out = drawCapa("tomo_capa1", tomoCapa1In, tomoCapa1Out, tomoCapaLength, track/2, tomoCapaSize)


# tomoCapa2In, tomoCapa2Out = [cross2OutRight[POS]+Vector(0, crossSize), [0, 1], track, crossGap], [None, None, track, gap]
# tomoCapa2In, tomoCapa2Out = drawCapa("tomo_capa2", tomoCapa2In, tomoCapa2Out, tomoCapa2Length, track/2, tomoCapa2Size)


# tomoDoubleJunc1In = [[tomoCapa1Out[POS][0], tomoCapa1Out[POS][1]*0.70+tomoCapa2Out[POS][1]*0.3]+Vector(-tomoDecal, 0), [0, 1], track, gap]
# tomoDoubleJunc1Out = [[tomoCapa1Out[POS][0], tomoCapa1Out[POS][1]*0.70+tomoCapa2Out[POS][1]*0.3]+Vector(-tomoDecal, 0), [0, -1], track, gap]

# tomoDoubleJunc2In = [[tomoCapa1Out[POS][0], tomoCapa1Out[POS][1]*0.5+tomoCapa2Out[POS][1]*0.5], [-1, 0], track, gap]
# tomoDoubleJunc2Out = [[tomoCapa1Out[POS][0], tomoCapa1Out[POS][1]*0.5+tomoCapa2Out[POS][1]*0.5], [1, 0], track, gap]

# drawCable("tomo_cable1", tomoOut, tomoCapa1In)
# drawCable("tomo_cable2", tomoCapa1Out, tomoDoubleJunc1In)
# drawCable("tomo_cable3", tomoDoubleJunc1Out, tomoDoubleJunc2In)
# drawCable("tomo_cable4", tomoDoubleJunc2Out, tomoCapa2Out)
# drawCable("tomo_cable5", tomoCapa2In, cross2OutRight)


# ######## MEMORY RESONATOR ###############
# #Dessin de la piste de la memoire (memory)

# memoryCapa1Pos, memoryWaveLength = LitExp("$memory_capa_pos", 0.2), LitExp("$memory_wave_length", 1.7)
# memoryCapaLength, memoryCapaSize = LitExp("$memory_capa_length", 0.05), LitExp("$memory_capa_size", 0.20)
# memoryCapa2Length, memoryCapa2Size = LitExp("$memory_capa2_length", 0.16), LitExp("$memory_capa2_size", 0.05)

# memoryCapa1In, memoryCapa1Out = [memoryOut[POS]+Vector(0, memoryCapa1Pos), [0, 1], track, gap], [None, None, track, gap]
# memoryCapa1In, memoryCapa1Out = drawCapa("memory_capa1", memoryCapa1In, memoryCapa1Out, memoryCapaLength, track/2, memoryCapaSize)

# memoryCapa2In, memoryCapa2Out = [cross1OutLeft[POS]+Vector(0, -crossSize), [0, -1], track, crossGap], [None, None, track, gap]
# memoryCapa2In, memoryCapa2Out = drawCapa("memory_capa2", memoryCapa2In, memoryCapa2Out, memoryCapa2Length, track/2, memoryCapa2Size)

# memoryWaveIn, memoryWaveOut = [memoryCapa1Out[POS], [0, 1], track, gap], [None, None, track, gap]
# memoryWaveIn, memoryWaveOut = drawWave("memory_wave", memoryWaveIn, memoryWaveOut, 2, memoryWaveLength, memoryWaveLength, 0.35)

# drawCable("memory_cable1", memoryOut, memoryCapa1In)
# drawCable("memory_cable2", memoryWaveOut, memoryCapa2Out)
# drawCable("memory_cable3", memoryCapa2In, cross1OutLeft)

# ######### FAST FLUX LINE ################
# #Dessin de la ligne d'application du flux (flux)

# fluxDecal = LitExp("$flux_decal", 0.0)

# fluxFluxIn[TRACK], fluxFluxIn[GAP] = 0.006, 0.022
# fluxFluxIn[POS] = fluxFluxIn[POS] + Vector(fluxDecal, -fluxFluxIn[TRACK])
# fluxFluxIn[ORI] = -fluxFluxIn[ORI]

# fluxDoubleJuncOut = [fluxFluxIn[POS]+Vector(0, -0.2), [0, 1], 0.006, 0.022]
# fluxDoubleJuncIn = [fluxFluxIn[POS]+Vector(0, -0.2), [0, -1], 0.006, 0.022]

# drawCable("flux_cable1", fluxOut, fluxDoubleJuncIn)
# drawCable("flux_cable2", fluxDoubleJuncOut, fluxFluxIn)

######## UNITE - SUBSTRACT - PerfE #####
#Finalisation globale du dessin



















#trackObject = unite(c.trackObjects)
#gapObject = unite(c.gapObjects)
#
#subtract(gapObject, "Ground_plane")
#assign_perfE([trackObject], "track")
#assign_perfE(c.bondwireObjects, "bondwire")


release()