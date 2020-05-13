# -*- coding: utf-8 -*-
"""
Created on Tue Oct 29 14:05:00 2019

@author: wcs
"""
import os

from HFSSdrawpy import Modeler, Body
from HFSSdrawpy.parameters import TRACK, GAP
import HFSSdrawpy.libraries.example_elements as elt
# import HFSSdrawpy.libraries.base_elements as base

pm = Modeler('hfss')

relative = pm.set_variable('1mm')

chip = Body(pm, 'chip', rel_coor=[[0, 0, 0], [1, 0, 0], [0, 1, 0]])

chip1 = Body(pm, 'chip1', rel_coor=[['1mm', 0, 0], [1, 0, 0], [0, 1, 0]])

chip2 = Body(pm, 'chip2', rel_coor=[['1mm', 0, 0], [
             1, 0, 0], [0, 1, 0]],  ref_name="chip1")
