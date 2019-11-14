# -*- coding: utf-8 -*-
"""
Created on Thu Oct 31 14:14:51 2019

@author: antho
"""

from sympy.parsing import sympy_parser
from hfss import extract_value_unit, \
                 extract_value_dim, \
                 parse_entry, \
                 get_active_project, \
                 ModelEntity, \
                 VariableString 
from pint import UnitRegistry
from gds_modeler import GdsModeler
from Vector import Vector
from KeyElement import KeyElt

ureg = UnitRegistry()
Q = ureg.Quantity

LENGTH = '[length]'
INDUCTANCE = '[length] ** 2 * [mass] / [current] ** 2 / [time] ** 2'
CAPACITANCE = '[current] ** 2 * [time] ** 4 / [length] ** 2 / [mass]'
RESISTANCE = '[length] ** 2 * [mass] / [current] ** 2 / [time] ** 3'

LENGTH_UNIT = 'meter'
INDUCTANCE_UNIT = 'nH'
CAPACITANCE_UNIT = 'fF'
RESISTANCE_UNIT = 'ohm'


class PythonModeler():

    def __init__(self, name_interface, connector): #"Hfss" or "Gds"
        self.ports = {}
        self.variables = {}
        self.connector = connector
        
        if name_interface=="hfss":
            project = get_active_project()
            design = project.get_active_design()
            self.design = design
            self.modeler = design.modeler
            self.modeler.set_units('mm')
            self.modeler.delete_all_objects()
            self.interface = self.modeler

        if name_interface=="gds":
            self.interface = GdsModeler()
            
        self.mode = name_interface
    
    def body(self, body_name, coor_name='Global', coor_sys=None):
        if coor_name != 'Global':
            if not(coor_sys is None):
                coor_sys = parse_entry(coor_sys)
                self.interface.create_coor_sys(coor_name, coor_sys)
        return Body(self.interface, coor_name, body_name, self.mode, self.variables, self.connector)

    def set_units(self, units='m'):
        if (self.modeler != None):
            self.modeler.set_units('mm')
            self.interface = self.modeler
    
    def generate_gds(self, name_file):
        assert type(self.interface)==GdsModeler
        self.interface.generate_gds(name_file)
    
    def set_variable(self, name, value):
        """
        name (str): name of the variable in HFSS e.g. 'chip_length'
        value (str, VarStr, float): value of the variable if str will try to analyse the unit
        """
        self.store_variable(name, value) # for later parsing
        self.__dict__[name] = VariableString(name) # for handy use by user
        if self.mode=='hfss':
            self.design.set_variable(name, value) # for HFSS
    
    def store_variable(self, name, value): # put value in SI
        if not isinstance(value, VariableString):
            if LENGTH == extract_value_dim(value):
                self.variables[name]=extract_value_unit(value, LENGTH_UNIT)
            if INDUCTANCE == extract_value_dim(value):
                self.variables[name]=extract_value_unit(value, INDUCTANCE_UNIT)
            if CAPACITANCE == extract_value_dim(value):
                self.variables[name]=extract_value_unit(value, CAPACITANCE_UNIT)
            if RESISTANCE == extract_value_dim(value):
                self.variables[name]=extract_value_unit(value, RESISTANCE_UNIT)
        else:
            self.variables[name]=value
            
    def eval_var_str(self, name, unit=None): # return numerical value of a given expression
                                             # using the values stored in self.variables
        # TODO: parse several times
        # can only parse 2 times for now

        try:
            _val = float(eval(str(sympy_parser.parse_expr(str(sympy_parser.parse_expr(str(name), self.variables)), self.variables))))
        except Exception:
            msg = ('Parsed expression contains a string which does '
                             'not correspond to any design variable')
            raise ValueError(msg)
        if unit is not None:
            _val = str(_val)+' '+unit
        return _val
    
    def val(self, name): # to use if you want to compare two litt expressions
                         # Compute the numerical value using for eval_var_str
        if isinstance(name, list):
            name_list = []
            for elt in name:
                name_list.append(self.eval_var_str(elt))
            if isinstance(name, Vector):
                return Vector(name_list)
            else:
                return name_list
        else:
            return self.eval_var_str(name)

    def box_corner(self, pos, size, layer, **kwargs):
        name = self.interface.draw_box_corner(pos, size, **kwargs)
        return ModelEntity(name, 3, self.coor_sys, layer=layer)
    
    def box_center(self, pos, size, layer, **kwargs):
        name = self.interface.draw_box_center(pos, size, **kwargs)
        return ModelEntity(name, 3, self.coor_sys, layer=layer)
  
#    def polyline_3D(self, points3D, closed=True, **kwargs): # among kwargs, name should be given
#        if self.mode=='hfss':
#            name = self.interface.draw_polyline(points3D, closed=closed, **kwargs)
#        elif self.mode=='gds':
#            points2D = []
#            for point in points3D:
#                points2D.append(points3D[0:2])
#            name = self.interface.draw_polyline(points2D, [0,0], closed=closed, **kwargs)
#
#        dim = closed + 1 # 2D when closed, 1D when open
#        return ModelEntity(name, dim, self.coor_sys, layer=kwargs['layer'])
#
#    def rect_corner_3D(self, pos3D, size3D, **kwargs):
#        if self.mode=='hffs': 
#            name = self.interface.draw_rect_corner(pos3D, size3D, **kwargs)
#        elif self.mode=='gds':
#            pos2D = pos3D[0:2]
#            size2D = size3D[0:2]
#            name = self.interface.draw_rect_corner(pos2D, size2D, **kwargs)
#        return ModelEntity(name, 2, self.coor_sys, layer=kwargs['layer'])
#    
#    def rect_center_3D(self, pos3D, size3D, **kwargs):
#        if self.mode=='hffs': 
#            name = self.interface.draw_rect_corner(pos3D, size3D, **kwargs)
#        elif self.mode=='gds':
#            pos2D = pos3D[0:2]
#            size2D = size3D[0:2]
#            name = self.interface.draw_rect_corner(pos2D, size2D, **kwargs)
#        return ModelEntity(name, 2, self.coor_sys, layer=kwargs['layer'])
    
#    def polyline_2D(self, points2D, z=0, closed=True, **kwargs): # among kwargs, name should be given
#        points3D = []
#        for point in points2D:
#            points3D.append(point+(z,))
#        return self.polyline_3D(points3D, closed, **kwargs)
#
#    def rect_corner_2D(self, pos2D, size2D, z=0, **kwargs):
#        name = self.interface.draw_rect_corner(pos2D, size2D, **kwargs)
#        return ModelEntity(name, 2, self.coor_sys, z, layer=kwargs['layer'])
#    
#    def rect_center_2D(self, pos2D, size2D, z=0, **kwargs):
#        name = self.interface.draw_rect_center(pos2D, size2D, **kwargs)
#        return ModelEntity(name, 2, self.coor_sys, z, layer=kwargs['layer'])



    def polyline_2D(self, points2D, z=0, closed=True, **kwargs): # among kwargs, name should be given
        name = self.interface.draw_polyline(points2D, [0,0], closed=closed, **kwargs)
        dim = closed + 1
        return ModelEntity(name, dim, self.coor_sys, layer=kwargs['layer'])

    def rect_corner_2D(self, pos2D, size2D, z=0, **kwargs):
        name = self.interface.draw_rect_corner(pos2D, size2D, **kwargs)
        return ModelEntity(name, 2, self.coor_sys, z, layer=kwargs['layer'])
    
    def rect_center_2D(self, pos2D, size2D, z=0, **kwargs):
        name = self.interface.draw_rect_center(pos2D, size2D, **kwargs)
        return ModelEntity(name, 2, self.coor_sys, z, layer=kwargs['layer'])
    
    
    
    
    def rot(self, x, y=0):
        return Vector(x, y).rot(self.ori)

    def cylinder(self, name, pos, radius, height, axis, layer, **kwargs):
        name = self.interface.draw_cylinder(pos, radius, height, axis, **kwargs)
        return ModelEntity(name, 3, self.coor_sys, layer=kwargs['layer'])
 
    def disk(self, name, pos, radius, axis, layer, **kwargs):
        name = self.interface.draw_cylinder(pos, radius, axis, **kwargs)
        return ModelEntity(name, 2, self.coor_sys, layer=kwargs['layer'])
    
    def wirebond(self, name, pos, ori, width, layer, **kwargs):
        name = self.interface.draw_cylinder(pos, ori, width, **kwargs)
        return ModelEntity(name, 2, self.coor_sys, layer=kwargs['layer'])
    
    def connect_faces(self, name, entity1, entity2):
        assert entity1.dimension == entity2.dimension
        assert entity1.dimension == 2
        return ModelEntity(name, 3, entity1.coor_sys, entity1.model)
        
    def delete(self, entity):
        del entity
        
    def rename_entity(self, entity, name):
        entity.rename_entity(name)
        
    def unite(self, entities, name=None, keep_originals=False):
        loc_entities = entities.copy()
        dim_Union = 0;
        for entity in loc_entities:
            if isinstance(entity, tuple) or isinstance(entity, list):
                for element in entity:
                    loc_entities.append(element)
            if isinstance(entity, ModelEntity):
                if entity.dimension>dim_Union:
                    dim_Union = entity.dimension
                
                if not(keep_originals):
                    self.delete(entity)
                 
        if name is None:
            union = entities[0].copy(dim_Union)
        else:
            union = entities[0].copy(dim_Union)
            union.rename_entity(name)
        
        
        return union
        
    def intersect(self, entities, keep_originals = False):
        dim_Intersection = 3;
        for entity in entities:
            if entity.dimension<dim_Intersection:
                dim_Intersection = entity.dimension
        if dim_Intersection == 0:
            pass
        elif dim_Intersection == 1:
            Intersection = ModelEntity(entities[0].name, 1, entities[0].coor_sys, entities[0].model)
        elif dim_Intersection == 2:
            Intersection = ModelEntity(entities[0].name, 2, entities[0].coor_sys, entities[0].model, entities[0].boundaries)
        elif dim_Intersection == 3:
            Intersection = ModelEntity(entities[0].name, 3, entities[0].coor_sys, entities[0].model, entities[0].boundaries)
   
        #A factoriser
        if not(keep_originals):
            for entity in entities:
                self.delete(entity)
        
        return Intersection
    
    def subtract(self, blank_entity, tool_entities):
        self.interface.subtract(blank_entity, tool_entities)
        for i in tool_entities:
            self.delete(i)
        return blank_entity

    def translate(self, entities, vector=[0,0,0]):
        self.interface.translate(entities, vector)
    
    def rotate(self, entities, angle=0):
        self.interface.rotate(entities, angle)
    
    def separate_bodies(self,name):
        #This looks hard
        pass
    
    def _fillet(self, radius, vertex_index, entity):
        self.interface._fillet(radius, vertex_index, entity.name)
#        print(vertices)
#        print(radius)
    
    def assign_perfect_E(self, entity):
        entity.boundaries.append('PerfE')
    
    def mirrorZ(self, obj):
        pass
    
    def copy(self, entity):
        if entity.dimension ==0:
            pass
        elif entity.dimension ==1:
            return ModelEntity(entity.name, entity.coor_sys, entity.model)
        elif entity.dimension ==2:
            return ModelEntity(entity.name, entity.coor_sys, entity.model, entity.boundaries)
        elif entity.dimension ==3:
            return ModelEntity(entity.name, entity.coor_sys, entity.model, entity.boundaries)
        
    def duplicate_along_line(self, entity, vec):
        #Create clones
        pass
    
    def assign_lumped_RLC(self, entity, ori, parameters):
#        entity.boundaries.append('lumpedRLC')
# TODO
        pass
    
    def make_material(self, entity, material):
        # Problème pour les unions de matériaux si on ajoute un attribut
        pass
    
    def delete_all_objects(self, entities):
        for entity in entities:
            self.delete(entity)
    
    def mesh_zone(self, zone, mesh_length):
        mesh_length = parse_entry(mesh_length)
        self.interface.assign_mesh_length(zone, mesh_length)

        
        
    
#    def rectangle(pos,size):
#        interface.rectangle(pos,size)
#        return ModelEntity(name)
#    
#    def chip():
#        return chip(repere,modeler)
#    

#        
#class GdsModeler():
#    def rectangle():
#        plt.plot(pos,size)
#        
class Body(PythonModeler, KeyElt):

    pos_elt = [0,0]
    ori_elt = 0
    
    def __init__(self, interface, coor_sys, name, mode, variables, connector):
        self.interface = interface
        self.coor_sys = coor_sys
        self.name = name
        self.maskObjects = []
        self.trackObjects = []
        self.gapObjects = []
        self.mode = mode # 'hfss' or 'gds'
        self.variables = variables
        self.connector = connector


    def append_points(self, coor_list):
        # coor_list can be [()()] or [[][]]
        points = [(coor_list[0][0],coor_list[0][1])]

        for coor in coor_list[1:]:
            points.append((points[-1][0] + coor[0],points[-1][1] + coor[1]))
        return points

    def refx_points(self, coor_list, offset=0, absolute=False):
        points=[]
        for ii, coor in enumerate(coor_list):
            if ii>0 and not absolute:
                offset=0
            points.append(Vector(*coor).refx(offset))
        return points

    def refy_points(self, coor_list, offset=0, absolute=False):
        points=[]
        for ii, coor in enumerate(coor_list):
            if ii>0 and not absolute:
                offset=0
            points.append(Vector(*coor).refy(offset))
        return points
    
    def move_points(self, coor_list, move, absolute=False):
        points=[]
        for ii, coor in enumerate(coor_list):
            if ii>0 and not absolute:
                move=[0,0]
            points.append(Vector(*coor)+Vector(*move))
        return points


    def coor(self, vec): # Change of coordinate for a point
        return self.rot(*vec)+self.pos

    def coor_vec(self, vec): # Change of coordinate for a vector
        return self.rot(*vec)
    