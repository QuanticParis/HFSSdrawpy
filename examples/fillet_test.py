import os

import HFSSdrawpy.libraries.example_elements as elt
from HFSSdrawpy import Body, Modeler
from HFSSdrawpy.parameters import GAP, TRACK

# import HFSSdrawpy.libraries.base_elements as base

pm = Modeler("hfss")

chip1 = Body(pm, "chip1")


track = pm.set_variable("20um")
gap = pm.set_variable("10um")
radius1 = pm.set_variable("100um")
radius2 = pm.set_variable("400um")


rect1 = chip1.rect([0, 0], ["1mm", "1mm"], layer=TRACK)
rect2 = chip1.rect(["0.5mm", "0.5mm"], ["-1mm", "-1mm"], layer=GAP)

rect1.unite(rect2)

rect1.fillet([radius1, radius2], [[3, 1, 2, -1, -2, -3], [0, 4]])
# convention for fillet :
# if the geometry is a genuine base element, fillet indices are order in the
# natural way :
# - order or points for a polyline
# - origin then 'x' dimension point etc for a rectangle
# If the polygon result from a boolean operation, the fillets are order
# such as the 0th is the leftest among the lowest points. Indices increase
# in the trigonometric order.

# generate gds file
pm.generate_gds(os.path.join(os.getcwd(), "gds_files"), "fillet_test")
