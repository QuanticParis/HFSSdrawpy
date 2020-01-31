# -*- coding: utf-8 -*-
"""
Created on Mon Jan 27 13:53:38 2020

@author: Zaki
"""

import gdspy
# Create the geometry: a single rectangle.
rect = gdspy.Rectangle((0, 0), (2, 1))
cell = gdspy.Cell('FIRST')
cell.add(rect)
# Save all created cells in file 'first.gds'.
gdspy.write_gds('first.gds')
# Optionally, display all cells using the internal viewer.
#gdspy.LayoutViewer(gdspy.current_library)