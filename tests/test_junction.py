# -*- coding: utf-8 -*-
"""
Created on Wed Mar 11 09:07:00 2020

@author: wcs
"""

'''
testing script for consolidation of junction functions in designer

changes:
    - 
'''  
  
from scripts.hfss import get_active_project, release, parse_entry
from scripts.designer import Circuit, KeyElt, Vector

project = get_active_project()
design = project.get_active_design()
modeler = design.modeler
modeler.set_units('mm')
modeler.delete_all_objects()

c = Circuit(design, modeler)

KeyElt.is_litho = False
KeyElt.is_hfss = True

### Drawing

c.set_variable('xx', '100um') # vertical spacing
c.set_variable('yy', '50um') # horizontal spacing

def pos(i, j):
    return [(i + 0.5)*c.xx, (j + 0.1)*c.yy]

# Finger junctions


c.key_elt('temp', pos(0, 0), [1, 0])
c.temp.draw_test_jcts(['10um', '10um'],
                '60um', '350nm', '400nm')

c.key_elt('temp', pos(0,1), [1,0])
c.temp.draw_test_jcts(['10um', '10um'], 
                       '60um', '350nm', '400nm', n_bridge=1, spacing_bridge='1.4um')


c.key_elt('temp', pos(0,2), [1,0])
c.temp.draw_test_jcts(['10um', '10um'], 
                 '60um', '350nm', '200nm')

c.key_elt('temp', pos(0,3), [1,0])
c.temp.draw_test_jcts(['10um', '10um'], 
                 '60um', '350nm', '200nm', overlap=2e-6)

c.key_elt('temp', pos(0,4), [1,0])
c.temp.draw_test_jcts(['10um', '10um'], 
                 '60um', '350nm', '200nm', rotspace='10um')

# Slab junctions

c.key_elt('temp', pos(1,0), [1,0])
c.temp.draw_test_jcts(['10um', '10um'], 
                 '60um', '350nm', '200nm', Width='500nm')

c.key_elt('temp', pos(1,1), [1,0])
c.temp.draw_test_jcts(['10um', '10um'], 
                 '60um', '350nm', '200nm', Width='500nm')

c.key_elt('temp', pos(1,2), [1,0])
c.temp.draw_test_jcts(['10um', '10um'], 
                 '60um', '350nm', '200nm', Width='500nm', overlap=2e-6)

c.key_elt('temp', pos(1,3), [1,0])
c.temp.draw_test_jcts(['10um', '10um'], 
                 '60um', '350nm', '200nm', Width='500nm', rotspace='10um')

# Cross junctions

c.key_elt('temp', pos(2,0), [1,0])
c.temp.draw_test_jct_cross(['10um', '10um'], 
                       '60um', '350nm', '200nm')

c.key_elt('temp', pos(2,1), [1,0])
c.temp.draw_test_jct_cross(['10um', '10um'], 
                       '60um', '350nm', '200nm', way=-1)

c.key_elt('temp', pos(2,2), [1,0])
c.temp.draw_test_jct_cross(['10um', '10um'], 
                       '60um', '350nm', '200nm')

c.key_elt('temp', pos(2,3), [1,0])
c.temp.draw_test_jct_cross(['10um', '10um'], 
                       '60um', '350nm', '200nm', overlap=2e-6)

c.key_elt('temp', pos(2,4), [1,0])
c.temp.draw_test_jct_cross(['10um', '10um'], 
                       '60um', '350nm', '200nm', rotspace='10um')

# Junction arrays

c.key_elt('temp', pos(3,0), [1,0])
c.temp.draw_test_jcts(['10um', '10um'], 
                       '60um', '400nm', '1.2um', n_bridge=10, spacing_bridge='1.4um')

c.key_elt('temp', pos(3,1), [1,0])
c.temp.draw_test_jcts(['10um', '10um'], 
                       '60um', '400nm', '1.2um', n_bridge=10, spacing_bridge='1.4um')

c.key_elt('temp', pos(3,2), [1,0])
c.temp.draw_test_jcts(['10um', '10um'], 
                       '60um', '400nm', '1.2um', n_bridge=10, spacing_bridge='1.4um', overlap=2e-6)

c.key_elt('temp', pos(3,3), [1,0])
c.temp.draw_test_jcts(['10um', '10um'], 
                       '60um', '400nm', '1.2um', n_bridge=10, spacing_bridge='1.4um', rotspace='20um')

### Other stuff

release()