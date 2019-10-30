# -*- coding: utf-8 -*-
"""
Created on Wed Oct 30 15:39:44 2019

@author: antho
"""

import hfss
from designer import Object, Chip
from KeyElement import KeyElt
object = Object("Model","vacuum")

elt = KeyElt("Josephson")
chip = Chip()
chip.set_variable("A",2)