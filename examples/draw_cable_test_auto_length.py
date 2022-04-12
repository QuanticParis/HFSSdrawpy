import os

import HFSSdrawpy.libraries.example_elements as elt
from HFSSdrawpy import Body, Modeler
from HFSSdrawpy.path_finding.path_finder_auto import Meander
from HFSSdrawpy.parameters import GAP, TRACK

# import HFSSdrawpy.libraries.base_elements as base

pm = Modeler("gds")

relative = pm.set_variable("1mm")

chip = Body(pm, "chip")
track = pm.set_variable("20um")
gap = pm.set_variable("10um", name="gap")

track_big = pm.set_variable("25um")
gap_big = pm.set_variable("15um")

track_middle = pm.set_variable("22.5um")
gap_middle = pm.set_variable("12.5um")

offset = pm.set_variable("-50um")

fillet = pm.set_variable('100um')

# chip2
with chip(["0.5mm", "0.5mm"], [1, 0]):
    port0, = elt.create_port(
        chip, [track, track + 2 * gap], name="port0")  # default is the widths of track and gap

with chip(["2.2mm", "0.7mm"], [1, 0]):
    port1, = elt.create_port(chip, name="port1")  # default is the widths of track and gap

with chip(["2mm", "2mm"], [1, 0]):
    port2, = elt.create_port(
        chip, [track, track + 2 * gap], name="port2")  # default is the widths of track and gap

meander0 = Meander(2, None)
chip.draw_cable_auto(port0, port2, is_bond=True, fillet=fillet,
                     meanders=[meander0],
                     name='cable_1',
                     mid=0, target_length='5mm')


ground_plane = chip.rect([0, 0], ["3mm", "3mm"], layer=TRACK)

ground_plane.subtract(chip.entities[GAP])
ground_plane.unite(chip.entities[TRACK])

# generate gds file
pm.generate_gds(os.path.join(os.getcwd(), "gds_files"), "cable_test")
