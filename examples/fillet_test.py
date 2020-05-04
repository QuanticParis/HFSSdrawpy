# -*- coding: utf-8 -*-
"""
Created on Tue Oct 29 14:05:00 2019

@author: wcs
"""
import os

from HFSSdrawpy import Modeler, Body
from HFSSdrawpy.parameters import layer_TRACK, layer_GAP, layer_RLC
import HFSSdrawpy.libraries.example_elements as elt
# import HFSSdrawpy.libraries.base_elements as base


pm = Modeler('gds')

chip1 = Body(pm, 'chip1')


track = pm.set_variable('20um')
gap = pm.set_variable('10um')
radius = pm.set_variable('100um')

rect1 = chip1.rect([0, 0], ['1mm', '1mm'], layer=layer_TRACK)
rect2 = chip1.rect(['0.5mm', '0.5mm'], ['-1mm', '-1mm'], layer=layer_GAP)

rect1.unite(rect2)

# rect1.fillet(radius, [3, 1, 2, -1])

# generate gds file
# pm.generate_gds(os.getcwd(), 'test')
pm.generate_gds(os.path.join(os.getcwd(), 'gds_files'), 'fillet_test')
# pm.generate_gds('G:/Git', 'test')

# ModelEntity.print_instances()
