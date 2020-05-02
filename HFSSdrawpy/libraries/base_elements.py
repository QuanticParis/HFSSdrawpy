# -*- coding: utf-8 -*-
"""
Created on Mon Oct 28 16:18:36 2019

@author: Ulysse
"""

#Binding the 2D elements of python_modeler

def port(self, *args, **kwargs):

    return self.port(*args, **kwargs)

def draw_cable(self, *args, **kwargs):

    return self.draw_cable(*args, **kwargs)

def disk_2D(self, pos, radius, axis, **kwargs):

    return self.disk_2D(pos, radius, axis, **kwargs)

def polyline_2D(self, points, closed=True, **kwargs):

    return self.polyline_2D(points, closed, **kwargs)

def path_2D(self, points, port, fillet, **kwargs):

    return self.path_2D(points, port, fillet, **kwargs)

def rect_corner_2D(self, pos, size, **kwargs):

    return self.rect_corner_2D(pos, size, **kwargs)

def rect_center_2D(self, pos, size, **kwargs):

    return self.rect_center_2D(pos, size, **kwargs)

def wirebond_2D(self, pos, ori, ymax, ymin, **kwargs):

    return self.wirebond_2D(pos, ori, ymax, ymin, **kwargs)

#Binding the 3D elements of python_modeler

def box_corner_3D(self, pos, size, **kwargs):

    return self.box_corner_3D(pos, size, **kwargs)

def box_center_3D(self, pos, size, **kwargs):

    return self.box_center_3D(pos, size, **kwargs)

def cylinder_3D(self, pos, radius, height, axis, **kwargs):

    return self.cylinder_3D(pos, radius, height, axis, **kwargs)