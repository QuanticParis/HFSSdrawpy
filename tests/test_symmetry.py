import numpy as np

from HFSSdrawpy import Modeler, Body

pm = Modeler("gds")
body = Body(pm, "chip")

dx = 1e1
dy = 1e-2

with body([0, 0], [1, 1]):
    body.polyline([[-dx, -dy], [-dx, dy], [dx, dy], [dx, -dy]])
    mirror = body.mirror(np.pi / 2, 0)
    with mirror:
        body.polyline([[0, 1], [-1, 2], [1, 2]])
        with body([0, 0.5], [1, 0]):
            with body([0, 0], [0, -1]):
                (port_1,) = body.port([0.1, 0.2], name="port_1")
            # with body([0.3, 0], [-1, 0]):
            #     (port_2,) = body.port([0.1, 0.2], name="port_2")
        # body.draw_cable(port_1, port_2)

    port_1_mirror = mirror.get_symmetric_of(port_1)
    # (port_2_mirror,) = mirror.get_symmetric_of(port_2)

    body.draw_cable(port_1, port_1_mirror.r)
    # body.draw_cable(port_2.r, port_2_mirror.r)

path = r"."
pm.generate_gds(path, "test_symmetry")
