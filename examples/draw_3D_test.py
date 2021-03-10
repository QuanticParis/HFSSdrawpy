# -*- coding: utf-8 -*-
"""
Created on Tue Oct 29 14:05:00 2019

@author: wcs
"""
import os

import HFSSdrawpy.libraries.example_elements as elt
from HFSSdrawpy import Body, Modeler
from HFSSdrawpy.parameters import GAP, TRACK

# import HFSSdrawpy.libraries.base_elements as base

pm = Modeler("hfss")

relative = pm.set_variable("1mm")

main = Body(pm, "main")

chip = Body(pm, "chip", rel_coor=[["1mm", "1mm", "1mm"], [1, 0, 0], [0, 0, 1]], ref_name="main")

chip1 = Body(pm, "chip1", rel_coor=[[0, 0, 0], [0, 1, 0], [1, 0, 0]], ref_name="chip")

chip2 = Body(pm, "chip2", rel_coor=[[0, 0, 0], [1, 0, 0], [0, 0, 1]], ref_name="chip")

track = pm.set_variable("20um")
gap = pm.set_variable("10um", name="gap")

track_big = pm.set_variable("25um")
gap_big = pm.set_variable("15um")

track_middle = pm.set_variable("22.5um")
gap_middle = pm.set_variable("12.5um")

offset = pm.set_variable("-50um")

# chip1

# default is the widths of track and gap
(port11,) = elt.create_port(chip1, [track, track + 2 * gap], name="port11")

with chip1(["2.0mm", "0.0mm"], [1, 0]):
    # default is the widths of track and gap
    (port12,) = elt.create_port(chip1, [track, track + 2 * gap], name="port12")

bond_length, bond_slope, pcb_track, pcb_gap = "200um", 0.5, "300um", "200um"

with chip1(["0.5mm", "0.5mm"], [0, 1]):
    (con_port1,) = elt.draw_connector(chip1, pcb_track, pcb_gap, bond_length, name="con_port1")

    with chip1(["1.5mm", "-1.0mm"], [0, 1]):
        (port13,) = elt.create_port(chip1, [track, track + 2 * gap], name="port13")

chip1.draw_cable(
    con_port1,
    port13,
    is_bond=True,
    fillet="100um",
    reverse_adaptor=False,
    to_meander=[0, 0, 0],
    meander_length=0,
    name="con_port1_port13",
)

ground_plane1 = chip1.rect([0, 0], ["3mm", "3mm"], layer=TRACK, name="gp1")

# chip2

# default is the widths of track and gap
(port21,) = elt.create_port(chip2, [track, track + 2 * gap], name="port21")

with chip2(["2.0mm", "0.0mm"], [1, 0]):
    # default is the widths of track and gap
    (port22,) = elt.create_port(chip2, [track, track + 2 * gap], name="port22")

bond_length, bond_slope, pcb_track, pcb_gap = "200um", 0.5, "300um", "200um"

with chip2(["0.5mm", "0.5mm"], [0, 1]):
    (con_port2,) = elt.draw_connector(chip2, pcb_track, pcb_gap, bond_length, name="con_port2")

    with chip2(["1.5mm", "-1.0mm"], [0, 1]):
        (port23,) = elt.create_port(chip2, [track, track + 2 * gap], name="port23")

chip2.draw_cable(
    con_port2,
    port23,
    is_bond=True,
    fillet="100um",
    reverse_adaptor=False,
    to_meander=[0, 0, 0],
    meander_length=0,
    name="con_port2_port23",
)
# # 3D
chip.box([0, 0, 0], ["3mm", "3mm", "3mm"], material="silicon")

ground_plane2 = chip2.rect([0, 0], ["3mm", "3mm"], layer=TRACK, name="gp2")

ground_plane1.subtract(chip1.entities[GAP])
ground_plane1.unite(chip1.entities[TRACK])
ground_plane1.assign_perfect_E()

ground_plane2.subtract(chip2.entities[GAP])
ground_plane2.unite(chip2.entities[TRACK])
ground_plane2.assign_perfect_E()

main.cylinder([0, 0, 0], "0.5mm", "0.7mm", "Z", name="tube")

# generate gds file
pm.generate_gds(os.path.join(os.getcwd(), "gds_files"), "cable_test")
