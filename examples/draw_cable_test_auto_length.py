import os

import HFSSdrawpy.libraries.example_elements as elt
from HFSSdrawpy import Body, Modeler
from HFSSdrawpy.path_finding.path_finder_auto import Meander
from HFSSdrawpy.parameters import GAP, TRACK

pm = Modeler("gds")

chip = Body(pm, "chip")

track = pm.set_variable("40um")
gap = pm.set_variable("20um", name="gap")

fillet = pm.set_variable('200um')

with chip(["2.2mm", "0.1mm"], [1, 0]):
    portint, = elt.create_port(chip, name="portint")  # default is the widths of track and gap

oris = [[1, 0], [0, 1], [-1, 0], [0, -1]]
xs = [0, 2.5, 5, 7.5]

for ii, (ori, x) in enumerate(zip(oris, xs)):

    with chip([x*1e-3, '0mm'], [1, 0]):
        port0, = elt.create_port(chip, [track, track + 2 * gap], name="port0_%d"%ii)  # default is the widths of track and gap
    
    with chip([(x+1.5)*1e-3, "1.5mm"], ori):
        port1, = elt.create_port(chip, [track, track + 2 * gap], name="port1_%d"%ii)  # default is the widths of track and gap

    meander0 = Meander(0, length=None, n=2, stretch=False, pos=0)
    meander1 = Meander(1, length='0.5mm')
    chip.draw_cable_auto(port0, port1,
                         fillet='100um',
                         meanders=[meander0, meander1],
                         mid=1,
                         target_length='6.0mm',
                         name='cable_%d'%ii)
    
                         

ground_plane = chip.rect(['-1mm', '-1mm'], ["11mm", "3.5mm"], layer=TRACK)

ground_plane.subtract(chip.entities[GAP])
ground_plane.unite(chip.entities[TRACK])

# generate gds file
pm.generate_gds(os.path.join(os.getcwd(), "gds_files"), "cable_test")
