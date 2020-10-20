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

def disk(self, pos, radius, axis, **kwargs):

    return self.disk(pos, radius, axis, **kwargs)

def polyline(self, points, closed=True, **kwargs):

    return self.polyline(points, closed, **kwargs)

def path(self, points, port, fillet, **kwargs):

    return self.path(points, port, fillet, **kwargs)

def rect(self, pos, size, **kwargs):

    return self.rect(pos, size, **kwargs)

def rect_center(self, pos, size, **kwargs):

    return self.rect_center(pos, size, **kwargs)

def wirebond(self, pos, ori, ymax, ymin, **kwargs):

    return self.wirebond(pos, ori, ymax, ymin, **kwargs)

#Binding the 3D elements of python_modeler

def box(self, pos, size, **kwargs):

    return self.box(pos, size, **kwargs)

def box_center(self, pos, size, **kwargs):

    return self.box_center(pos, size, **kwargs)

def cylinder(self, pos, radius, height, axis, **kwargs):

    return self.cylinder(pos, radius, height, axis, **kwargs)

def sphere(self, pos, radius, **kwargs):

    return self.sphere(pos, radius, **kwargs)

def torus(self, pos, majorradius, minorradius, axis, **kwargs):
    
    return self.torus(pos, majorradius, minorradius, axis, **kwargs)