import numpy as np

from HFSSdrawpy import Modeler, Body

pm = Modeler("gds")
body = Body(pm, "chip")

dx = 1e1
dy = 1e-2

track = pm.set_variable("100mm")
gap = pm.set_variable("50mm")
fillet = track + 1.5 * gap

body.polyline([[-dx, -dy], [-dx, dy], [dx, dy], [dx, -dy]], layer=3)
mirror = body.mirror(np.pi / 2, 0)
with mirror:
    # body.polyline([[0, track], [track, 2 * track], [3 * track, 2 * track]], closed=True)
    with body([-10 * track, 10 * track], [1, 0]):
        (port_1,) = body.port(name="port_1", widths=[track, track + 2 * gap])
    with body([0, 2 * track], [0, -1]):
        (port_2,) = body.port(name="port_2", widths=[track, track + 2 * gap])
    body.draw_cable(port_1, port_2, fillet=fillet)

port_2_symmetric = mirror.get_symmetric_of(port_2)
body.draw_cable_auto(port_2, port_2_symmetric.r, fillet=fillet)

path = r"."
pm.generate_gds(path, "test_symmetry")
