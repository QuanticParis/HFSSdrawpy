import numpy as np

from HFSSdrawpy import Modeler, Body
from HFSSdrawpy.path_finding.path_finder import Path

pm = Modeler("gds")
body = Body(pm, "chip")

dx = 1e1
dy = 1e-2

track = pm.set_variable("100mm")
gap = pm.set_variable("50mm")
fillet = track + 1.5 * gap

with body([0, 0], [1, 0]):
    body.polyline([[-dx, -dy], [-dx, dy], [dx, dy], [dx, -dy]])
    mirror = body.mirror(np.pi / 2, 0)
    with mirror:
        body.polyline(
            [[0, track], [track, 2 * track], [3 * track, 2 * track]], closed=False
        )
        with body([0, 0.5], [1, 0]):
            with body([-1, 0.6], [1, 0]):
                (port_1,) = body.port(name="port_1", widths=[track, track + 2 * gap])
            with body([0, 0], [0, -1]):
                (port_2,) = body.port(name="port_2", widths=None)
            with body([1, 0.6], [0, 1]):
                (port_3,) = body.port(name="port_3", widths=[track, track + 2 * gap])
        path = Path("my_path", port_1, port_2, fillet)
        body.path(path.points, port_3.r, fillet=fillet)

# port_1_mirror = mirror.get_symmetric_of(port_1)
# body.draw_cable_auto(port_1, port_1_mirror.r, fillet="0mm")

path = r"."
pm.generate_gds(path, "test_symmetry")
