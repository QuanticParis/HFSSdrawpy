 # -*- coding: utf-8 -*-
"""
Created on Thu Oct 25/04/2020

@author: Ulysse
"""

##################################################################
# This file is ment to regroup the constant parameters of DrawPy #
# so that they are not defined more than one time, and we avoid  #
# circular inclusion at the same time.                           #
##################################################################


# PARAMETERS FOR THE GDS OUTPUT AND FOR FILLETS
eps = 1e-7
layer_TRACK = 1
layer_GAP = 0
layer_RLC = 2
layer_MESH = 3
layer_MASK = 4
layer_PORT = 10
layer_Default = 10