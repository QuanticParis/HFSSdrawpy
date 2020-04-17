 # -*- coding: utf-8 -*-
"""
Created on Thu Oct 31 14:14:51 2019

@author: antho
"""

from sympy.parsing import sympy_parser
from pint import UnitRegistry
import numpy as np
# import sys
from functools import wraps
from inspect import currentframe, getfile

from . import Lib
from . import KeyElement
from . import CustomElement
from .variable_string import VariableString, \
                             extract_value_unit, \
                             extract_value_dim, \
                             parse_entry, \
                             _val, val, \
                             way, \
                             Vector
from .path_finder import Path

# ModelEntity should be defined here probably

# PARAMETERS FOR THE GDS OUTPUT AND FOR FILLETS
eps = 1e-7
layer_TRACK = 1
layer_GAP = 0
layer_RLC = 2
layer_MESH = 3
layer_MASK = 4
layer_PORT = 10
layer_Default = 10

##IMPORT KEY / CUSTOM Elements




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

def equal_float(float1, float2):
    if float1!=0:
        rel_diff = abs((float1-float2)/float1)
        if rel_diff<1e-5:
            return True
        else:
            return False

    elif float2!=0:
        rel_diff = abs((float1-float2)/float2)
        if rel_diff<1e-5:
            return True
        else:
            return False
    else:
        return True

class PythonModeler():
    """
    Modeler which defines basic operations and methods to perform on ModelEntity and on the chosen interface.
    To create a new interface, one needs to copy an existing one in a new file, and adapt all methods to the formalism of the new interface.

    The PythonMdlr class is called at the beginning of a script but then it is the body that is always used.
    The syntax is the following:
        import scripts
        PM = PythonMdlr("hfss")
        chip1 = PM.body(name, coordinate system)

    Inputs:
    -------
    name_interface: string in "gds" or "hfss"
    """
    is_overdev = False
    is_litho = False
    is_mask = False
    gap_mask = parse_entry('20um')
    overdev = parse_entry('0um')

    def __init__(self, name_interface):
        """
        Creates a PythonMdlr object based on the chosen interface.
        For now the interface cannot be changed during an execution, only at the beginning
        """
        self.variables = {}
        self.name_interface = name_interface
        if name_interface=="hfss":
            from .hfss import get_active_project
            project = get_active_project()
            design = project.get_active_design()
            self.design = design
            self.modeler = design.modeler
            self.modeler.set_units('mm')
            self.modeler.delete_all_objects()
            self.interface = self.modeler

        elif name_interface=="gds":
            from . import gds_modeler
            self.interface = gds_modeler.GdsModeler()

        else:
            print('Interface should be either hfss or gds')

        self.mode = name_interface

    def set_active_coor_system(func):
        """
        Defines a wrapper/decorator which allows the user to always work in the coordinate system of the chosen chip.
        """
        @wraps(func)
        def updated(*args, **kwargs):
            args[0].interface.set_coor_sys(args[0].coor_sys)
            return func(*args, **kwargs)
        return updated

    def append_lists(self):
        """
        We use a tree-like architecture to store the entities and port to be moved at the right place.
        """
        self.modelentities_to_move().append([])
        self.ports_to_move().append([])

    @classmethod
    def append_points(self, coor_list):
        """
        Almost depreciated.
        It allows the user to defines a polygon using its vertices.

        Inputs:
        -------
        coor_list: [[x1,y1],[x2,y2],...] or [(x1,y1),(x2,y2),...]

        Outputs:
        -------
        new_coor_list: new_coor_list[k] = Somme des k premiers vecteurs de coor_list
        """
        new_coor_list = [(coor_list[0][0],coor_list[0][1])]

        for coor in coor_list[1:]:
            new_coor_list.append((new_coor_list[-1][0] + coor[0],new_coor_list[-1][1] + coor[1]))
        return new_coor_list

    def assign_perfect_E(self, entities, name='perfE'):
        """
        Change the property of the HFFS objects which correspond to the entites.
        """
        if isinstance(entities, list):
            self.interface.assign_perfect_E(entities, entities[0].name+name)
        else:
            self.interface.assign_perfect_E(entities, entities.name+name)

    def body(self, body_name, coor_name='Global', coor_sys=None):
        """
        Creates a Body object which inherits from the current PythonModeler object.
        The body is associated with a coordinate system of choice.
        """
        if coor_name != 'Global':
            if not(coor_sys is None):
                coor_sys = parse_entry(coor_sys)
                self.interface.create_coor_sys(coor_name, coor_sys)
        N = Network(body_name, coor_name, self.interface, self.variables)
        B = Body(self.interface, coor_name, body_name, N, self.mode, self.variables)
        return B

    @set_active_coor_system
    def box_corner_3D(self, pos, size, **kwargs):
        """
        Draws a 3D box based on the coordinates of its corner.

        Inputs:
        -------
        pos: [x,y,z] the coordinates of the corner of the rectangle in the euclidian basis.
        size: [lx,ly,lz] the dimensions of the rectangle

        **kwargs include layer and name
        Outputs:
        -------
        box: Corresponding 3D Model Entity
        """

        layer = kwargs['layer'] if layer!=None else layer_Default
        name = self.interface.box_corner_3D(pos, size, **kwargs)
        box = ModelEntity(name, 3, self.coor_sys, layer=layer)
        return box

    @set_active_coor_system
    def box_center(self, pos, size, layer, **kwargs):
        """
        Draws a 3D box based on the coordinates of its center.

        Inputs:
        -------
        pos: [x,y,z] the coordinates of the center of the rectangle in the euclidian basis.
        size: [lx,ly,lz] the dimensions of the rectangle

        **kwargs include layer and name
        Outputs:
        -------
        box: Corresponding 3D Model Entity
        """
        layer = kwargs['layer'] if layer!=None else layer_Default
        name = self.interface.box_center_3D(pos, size, **kwargs)
        return ModelEntity(name, 3, self.coor_sys, layer=layer)

    @set_active_coor_system
    def cylinder_3D(self, pos, radius, height, axis, **kwargs):
        kwargs['coor_sys']=self.coor_sys
        name = self.interface.cylinder_3D(pos, radius, height, axis, **kwargs)
        return ModelEntity(name, 3, self.coor_sys, layer=kwargs['layer'])

    def connect_faces(self, name, entity1, entity2):
        assert entity1.dimension == entity2.dimension
        assert entity1.dimension == 2
        return ModelEntity(name, 3, entity1.coor_sys, entity1.model)

    def copy(self, entity):
        if entity.dimension ==0:
            pass
        elif entity.dimension ==1:
            return ModelEntity(entity.name, entity.coor_sys, entity.model)
        elif entity.dimension ==2:
            return ModelEntity(entity.name, entity.coor_sys, entity.model, entity.boundaries)
        elif entity.dimension ==3:
            return ModelEntity(entity.name, entity.coor_sys, entity.model, entity.boundaries)

    def delete(self, entity):
        if (entity.name in entity.__class__.dict_instances):
            entity.__class__.find_last_list(entity.__class__.instances_to_move).remove(entity)
            a = entity.__class__.dict_instances.pop(entity.name, None)
        del entity

    def delete_all_objects(self, entities):
        for entity in entities:
            self.delete(entity)

    @set_active_coor_system
    def disk_2D(self, pos, radius, axis, **kwargs):
        kwargs['coor_sys']=self.coor_sys
        if self.mode=='gds':
            pos = val(pos)
            radius = val(radius)
        name = self.interface.disk_2D(pos, radius, axis, **kwargs)
        return ModelEntity(name, 2, self.coor_sys, layer=kwargs['layer'])

    def duplicate_along_line(self, entity, vec):
        #Create clones
        pass

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

    def _fillet(self, radius, vertex_index, entity):
        if self.mode=='gds':
            radius = val(radius)
        try:
            self.interface._fillet(radius, vertex_index, entity)
        except Exception:
            print("Fillet operation resulted in an error")

    def _fillets(self, radius, entity, kind='open'):
        print(self)
        print(type(entity))
        vertices = self.interface.get_vertex_ids(entity)
        if entity.dimension==1:
            self.interface._fillets(radius, vertices[1:-1], entity)
        elif entity.dimension==2:
            self.interface._fillets(radius, entity)

    def generate_gds(self, name_file):
        '''Only for gds modeler'''
        self.interface.generate_gds(name_file)

    def intersect(self, entities, keep_originals = False):
        dim_Intersection = 3;
        for entity in entities:
            if entity.dimension<dim_Intersection:
                dim_Intersection = entity.dimension

        Intersection = ModelEntity(entities[0].name, dim_Intersection, entities[0].coor_sys, entities[0].model)

        if not(keep_originals):
            for entity in entities:
                self.delete(entity)
        return Intersection

    def make_rlc_boundary(self, corner, size, axis, r, l, c, name="LumpRLC"):
#        self.interface.make_rlc_boundary(corner, size, axis, r, l, c, name)
        pass

    def make_material(self, entity, material):
        self.interface.make_material(entity, material)
        pass

    def mesh_zone(self, entity, mesh_length):
        mesh_length = parse_entry(mesh_length)
#        print("mesh_length",mesh_length)
        self.interface.assign_mesh_length(entity, mesh_length)

    def mirrorZ(self, obj):
        pass

    def modelentities_to_move(self):
        inter1 = ModelEntity.instances_to_move
        inter2= ModelEntity.find_last_list(inter1)

        return inter2

    def move_points(self, coor_list, move, absolute=False):
        points=[]
        for ii, coor in enumerate(coor_list):
            if ii>0 and not absolute:
                move=[0,0]
            points.append(Vector(*coor)+Vector(*move))
        return points

    @set_active_coor_system
    def polyline_2D(self, points, closed=True, **kwargs): # among kwargs, name should be given
#        try:
#            print(kwargs['name'])
#        except:
#            pass
#
        kwargs['coor_sys']=self.coor_sys
        i = 0
        while i < len(points[:-1]):
            if np.array_equal(points[i], points[i+1]):
                points.pop(i)
            else:
                i+=1
        if self.mode=='gds':
            points = val(points)
        name = self.interface.polyline_2D(points, closed, **kwargs)
        dim = closed + 1
        return ModelEntity(name, dim, self.coor_sys, layer=kwargs['layer'])

    def path(self, points, port, fillet, **kwargs):
        if self.mode=='gds':
            points = val(points)
            fillet = val(fillet)
            _port = port.val()
            elements = self.interface.path(points, _port, fillet, **kwargs)
            return tuple(elements)
            # return ModelEntity(name, dim, self.coor_sys, layer=kwargs['layer'])
        elif self.mode=='hfss':
            # check that port is at the BEGINNING of the path
            self.polyline_2D(points, closed=True, **kwargs)
            ori = port.ori
            pos = port.pos
            entities = []
            for ii in port.N:
                offset = port.offsets[ii]
                width = port.widths[ii]
                subname = port.subnames[ii]
                layer = port.layers[ii]
                points_starter = [Vector(0, offset+width/2).rot(ori)+pos,
                                  Vector(0, offset-width/2).rot(ori)+pos]
                entity = self.polyline_2D(points_starter, closed=False,
                                          name=kwargs['name']+'_'+subname,
                                          layer=layer)
                path = self.polyline_2D(points, closed=False,
                                name=kwargs['name']+'_'+subname+'_path',
                                layer=layer)
                self._fillets(fillet, path, kind='open')

                entity = self._sweep_along_path(entity, path, name=None)
                entities.append(entity)

            return tuple(entities)

    def ports_to_move(self):
        inter1 = Port.instances_to_move
        inter2= Port.find_last_list(inter1)
        return inter2

    @set_active_coor_system
    def rect_corner_2D(self, pos, size, **kwargs):
        kwargs['coor_sys']=self.coor_sys
        if self.mode=='gds':
            pos = val(pos)
            size = val(size)
        name = self.interface.rect_corner_2D(pos, size, **kwargs)
        return ModelEntity(name, 2, self.coor_sys, layer=kwargs['layer'])

    @set_active_coor_system
    def rect_center_2D(self, pos, size, **kwargs):
        kwargs['coor_sys']=self.coor_sys
        if self.mode=='gds':
            pos = val(pos)
            size = val(size)
        name = self.interface.rect_center_2D(pos, size, **kwargs)
        return ModelEntity(name, 2, self.coor_sys, layer=kwargs['layer'])

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

    def rename_entity_1(self, entity, name):
        if entity.name in ModelEntity.dict_instances:
            ModelEntity.dict_instances.pop(entity.name,None)
            ModelEntity.dict_instances[name]=entity

        self.interface.rename_entity(entity,name)
        entity.rename_entity(name)


    @staticmethod
    def reset(self):
        Port.instances_to_move = []
        ModelEntity.instances_to_move = []

    def rot(self, x, y=0):
        return Vector(x, y).rot(self.ori)

    def rotate(self, entities, angle=0):
        if isinstance(angle, list) or isinstance(angle, np.ndarray):
            if len(angle)==2:
                new_angle= np.math.atan2(np.linalg.det([[1,0],angle]),np.dot([1,0],angle))
                new_angle= new_angle/np.pi*180
            else:
                raise Exception("angle should be either a float or a 2-dim array")
        elif isinstance(angle, float) or isinstance(angle, int):
            new_angle = angle
        else:
            raise Exception("angle should be either a float or a 2-dim array")

        self.interface.rotate(entities, new_angle)

    def separate_bodies(self, name):
        # This looks hard
        pass

    def set_current_coor(self, pos, ori):
        self.current_pos, self.current_ori = parse_entry(pos, ori)

    def set_variable(self, value, name=None):
        """
        name (str): name of the variable in HFSS e.g. 'chip_length'
        value (str, VarStr, float): value of the variable
                                    if str will try to analyse the unit
        """
        if name is None:
            # this auto-parsing is clearly a hack and not robust
            # but I find it convenient
            f = currentframe().f_back#.f_back
            filename = getfile(f)
            code_line = open(filename).readlines()[f.f_lineno - 1]
            name = code_line.split("=")[0].strip()

        # self.store_variable(name, value)  # for later parsing
        # self.__dict__[name] = VariableString(name)
        # for handy use by user / deprecated
        if self.mode == 'hfss':
            self.design.set_variable(name, value)  # for HFSS
        if not name in VariableString.variables.keys():
            return VariableString(name, value=value)
        else:
            VariableString.variables[name]=value
            return VariableString.instances[name]

    # def store_variable(self, name, value):  # put value in SI
    #     if not isinstance(value, VariableString):
    #         if LENGTH == extract_value_dim(value):
    #             self.variables[name] = extract_value_unit(value,
    #                                                       LENGTH_UNIT)
    #         if INDUCTANCE == extract_value_dim(value):
    #             self.variables[name] = extract_value_unit(value,
    #                                                       INDUCTANCE_UNIT)
    #         if CAPACITANCE == extract_value_dim(value):
    #             self.variables[name] = extract_value_unit(value,
    #                                                       CAPACITANCE_UNIT)
    #         if RESISTANCE == extract_value_dim(value):
    #             self.variables[name] = extract_value_unit(value,
    #                                                       RESISTANCE_UNIT)
    #     else:
    #         self.variables[name] = value

    def subtract(self, blank_entity, tool_entities, keep_originals=False):
        name = self.interface.subtract(blank_entity, tool_entities,
                                       keep_originals)
        if not(keep_originals):
            for tool_entity in tool_entities:
                self.delete(tool_entity)
#        return ModelEntity(name, 2, self.coor_sys, layer=blank_entity.layer)

    def _sweep_along_path(self, entity_to_sweep, path_entity, name=None):
        new_name = self.interface._sweep_along_path(entity_to_sweep,
                                                    path_entity)
#        self.delete(path_entity)
        entity_to_sweep.modify_dimension(2)
        path_entity.rename_entity(new_name)
        return path_entity

    def translate(self, entities, vector=[0, 0, 0]):
        if self.mode == 'gds':
            vector = val(vector)
        self.interface.translate(entities, vector)

    def unite(self, entities, name=None, keep_originals=False):
        loc_entities = []
        dim_Union = 0
        union = None
        for entity in entities:
            if entity is not None and isinstance(entity, ModelEntity):
                union = entity

            if union is not None:
                loc_entities.append(union)
                if union.dimension > dim_Union:
                    dim_Union = union.dimension
                if not(keep_originals):
                    self.delete(union)

        if len(loc_entities)>1:
            new_name = self.interface.unite(loc_entities)
            if name is None:
                union = entities[0].copy(dim_Union)
                self.rename_entity_1(union, new_name)

            else:
                union = entities[0].copy(dim_Union)
                self.rename_entity_1(union, name)

        return union

    # def val(self, name): # to use if you want to compare two litt expressions
    #                      # Compute the numerical value using for eval_var_str
    #     if isinstance(name, list) or isinstance(name, tuple):
    #         print('list')
    #         name_list = []
    #         for elt in name:
    #             if isinstance(elt, VariableString):
    #                 name_list.append(elt.value())
    #             else:
    #                 name_list.append(elt)
    #         if isinstance(name, Vector):
    #             return Vector(name_list)
    #         else:
    #             return name_list
    #     else:
    #         if isinstance(name, VariableString):
    #             return name.value()
    #         else:
    #             return

    @set_active_coor_system
    def wirebond_2D(self, pos, ori, ymax, ymin, **kwargs):
        kwargs['coor_sys']=self.coor_sys
        if self.mode=='gds':
            pos, ori, ymax, ymin = val(pos, ori, ymax, ymin)
            name_a, name_b = self.interface.wirebond_2D(pos, ori, ymax, ymin, **kwargs)
            return ModelEntity(name_a, 2, self.coor_sys, layer=kwargs['layer']), \
                    ModelEntity(name_b, 2, self.coor_sys, layer=kwargs['layer'])
        else:
            name = self.interface.wirebond_2D(pos, ori, ymax, ymin, **kwargs)
            return ModelEntity(name, 3, self.coor_sys, layer=kwargs['layer'])

class ModelEntity():
    instances_layered = {}
    dict_instances = {}
    instances_to_move = []
    def __init__(self, name, dimension, coor_sys, model = 'True', layer=layer_Default):# model,
        self.name = name
        self.dimension = dimension
        self.coor_sys = coor_sys
        self.history = []
        self.model = model
        self.layer = layer

        ModelEntity.dict_instances[name] = self
        ModelEntity.find_last_list(ModelEntity.instances_to_move).append(self)
        if layer in self.instances_layered.keys():
            ModelEntity.instances_layered[layer].append(self)
        else:
            ModelEntity.instances_layered[layer]=[self]

    def __str__(self):
        return self.name

    def check_name(self, name):
        i = 0
        new_name = name
        while(new_name in self.dict_instances.keys()):
            i+=1
            new_name = name+'_'+str(i)
        return new_name

    @classmethod
    def printInstances(cls):
        for instance in cls.instances:
            print(instance)

    @staticmethod
    def find_last_list(list_entities=instances_to_move):
#        print(list_entities)
#        if isinstance(list_entities, Port) :
#            return None
        if isinstance(list_entities, list):
            if len(list_entities)==0:
                return list_entities
            else:
                if isinstance(list_entities[-1], list):
                    return ModelEntity.find_last_list(list_entities[-1])
                else:
                    return list_entities
        else:
            return list_entities


#    @staticmethod
#    def merge(n):
#        key_list = []
#        for idx, key in enumerate(ModelEntity.instances_to_move):
#            if idx>=n:
#                key_list.append(key)
#        for key in key_list:
#            entity = ModelEntity.instances_to_move.pop(key, None)
#            ModelEntity.instances_moved[key]=entity

    @staticmethod
    def reset():
        ModelEntity.instances_layered = {}
        ModelEntity.dict_instances = {}
        ModelEntity.instances_to_move = []

    def append_history(self, entity):
        self.history.append(entity)

    def delete(self):
        del self

    def copy(self, dimension=None):
        #TODO complete
        if dimension==None:
            dimension = self.dimension
        return ModelEntity(self.name, dimension, self.coor_sys)

    def rename_entity(self, name):
        self.name = name

    def modify_dimension(self, dimension):
        self.dimension = dimension

    def _sweep_along_path(self, path_obj):
        self.dimension = 2
        self.append_history(path_obj)
        return self

    def sweep_along_vector(self, vect):
        "Duplicate the line and create a surface with the copy"
        self.dimension = 2
        return self

    def thicken_sheet(self, thickness, bothsides=False):
        self.dimension = 3
        return self

class Port():
    instances_to_move = []
    dict_instances  = {}

    def __init__(self, name, pos, ori, widths, subnames, layers, offsets, constraint_port, key='name'):
        new_name = name
        if not isinstance(key, Port):
            new_name = self.check_name(name)
        self.name = new_name
        self.pos = Vector(pos)
        self.ori = Vector(ori)
        self.constraint_port = constraint_port
        self.save = None
        if not constraint_port:
            self.widths = parse_entry(widths)
            self.subnames = subnames
            self.layers = layers
            self.offsets = parse_entry(offsets)
            self.N = len(widths)
        else:
            self.widths = widths
            self.subnames = subnames
            self.layers = layers
            self.offsets = offsets
            self.N = 0

        Port.find_last_list(Port.instances_to_move).append(self)
        if key=='name':  # normal initialisation
            self.dict_instances[name] = self

            # create a reversed version of the port that can be called by either
            # port.r or 'port_name.r'
            reversed_ori = -self.ori
            reversed_offsets = None
            if self.offsets is not None:
                reversed_offsets = []
                for ii in range(self.N):
                    reversed_offsets.append(-self.offsets[ii])
            self.r = Port(self.name+'_r', self.pos, reversed_ori, self.widths,
                          self.subnames, self.layers, reversed_offsets,
                          self.constraint_port, key=self)

        elif isinstance(key, Port):  # reverse initialisation, key is the previous port
            self.dict_instances[name] = self
            self.r = key
        else:
            pass  # when the port is only a float eval do not add it in dict



    @staticmethod
    def find_last_list(list_entities=instances_to_move):
        if isinstance(list_entities, list):
            if len(list_entities)==0:
                return list_entities
            else:
                if isinstance(list_entities[-1], list):
                    return ModelEntity.find_last_list(list_entities[-1])
                else:
                    return list_entities
        else:
            return list_entities

    @staticmethod
    def reset():
        Port.instances_to_move = []
        Port.dict_instances  = {}

    def check_name(self, name):
        i = 0
        new_name = name
        while(new_name in self.dict_instances.keys()):
            i+=1
            new_name = name+'_'+str(i)
        return new_name

    # def reverse(self):
    #     self.ori = -self.ori
    #     if self.offsets is not None:
    #         offsets = []
    #         for ii in range(self.N):
    #             offsets.append(-self.offsets[ii])
    #         self.offsets = offsets

    def compare(self, other, pm):
        points = []

        adapt_dist = pm.set_variable(1e-5, name=self.name+'_adapt')
        max_diff = 0
        for ii in range(self.N):
            if self.layers[ii]!=other.layers[ii]:
                raise ValueError('Tried to connect ports form different \
                                 layers: %s != %s'%(self.layers[ii],
                                                    other.layers[ii]))
            width1 = self.widths[ii]
            width2 = other.widths[ii]

            offset1 = self.offsets[ii]
            offset2 = -other.offsets[ii]

            if width1!=width2 or offset1!=offset2:
                # need adaptor
                points.append([Vector(0, offset1+width1/2).rot(self.ori)+self.pos,
                               Vector(adapt_dist, offset2+width2/2).rot(self.ori)+self.pos,
                               Vector(adapt_dist, offset2-width2/2).rot(self.ori)+self.pos,
                               Vector(0, offset1-width1/2).rot(self.ori)+self.pos])
            max_diff = max(max_diff, abs(_val(offset1+width1/2-(offset2+width2/2))),
                           abs(_val(offset2-width2/2-(offset1-width1/2))))

        adapt_dist = pm.set_variable(2*max_diff, name=self.name+'_adapt')

        if len(points) != 0:
            self.save = {'pos':self.pos, 'widths':self.widths,
                         'offsets':self.offsets}
            self.pos = self.pos + Vector(2*max_diff, 0).rot(self.ori)
            self.widths = other.widths
            self.offsets = [-offset for offset in other.offsets]

            self.r.save = {'pos':self.r.pos, 'widths':self.r.widths,
                         'offsets':self.r.offsets}
            self.r.pos = self.pos
            self.r.widths = self.widths
            self.r.offsets = other.offsets
        return points, 2*max_diff

    def val(self):
        _widths = []
        _offsets = []
        for ii in range(self.N):
            width = self.widths[ii]
            offset = self.offsets[ii]
            _widths.append(_val(width))
            _offsets.append(_val(offset))

        _pos = []
        for coor in self.pos:
            _pos.append(_val(coor))
        _pos = Vector(_pos)

        _ori = []
        for coor in self.ori:
            _ori.append(_val(coor))
        _ori = Vector(_ori)

        return Port(self.name, _pos, _ori, _widths, self.subnames,
                    self.layers, _offsets, self.constraint_port, key=None)

    def revert(self):
        if self.save is not None:
            self.pos = self.save['pos']
            self.widths = self.save['widths']
            self.offsets = self.save['offsets']

            self.r.pos = self.save['pos']
            self.r.widths = self.save['widths']
            reversed_offsets = []
            for ii in range(self.N):
                reversed_offsets.append(-self.offsets[ii])
            self.r.offsets = reversed_offsets

    def bond_params(self):
        y_max = -np.infty
        y_max_val = -np.infty
        y_min = np.infty
        y_min_val = np.infty
        for ii in range(self.N):
            # widths should not be negative
            _y_max_val = _val(self.offsets[ii]+self.widths[ii]/2)
            _y_min_val = _val(self.offsets[ii]-self.widths[ii]/2)
            if _y_max_val > y_max_val:
                y_max = self.offsets[ii]+self.widths[ii]/2
                y_max_val = _y_max_val
            if _y_min_val < y_min_val:
                y_min = self.offsets[ii]-self.widths[ii]/2
                y_min_val = _y_min_val
        return y_max, y_min

class Network(PythonModeler):
    variables= None
    to_bond=[]

    def __init__(self, name, coor_sys, interface, variables):
        self.interface = interface
        self.coor_sys = coor_sys
        self.name = name
        self.variables = variables
#        self.__class__.variables = variables


    def update(self, coor_sys):
        self.coor_sys = coor_sys

    def port(self, name, pos, ori, widths, subnames, layers, offsets, constraint_port):
        return Port(name, pos, ori, widths, subnames, layers, offsets, constraint_port)

    def draw_capa(self, name, iInPort, iOutPort, iLength, iWidth, iSize):
        '''
        Inputs:
        -------
        name: string name of object
        iIn: (position, direction, track, gap) defines the input port
        iOut: (position, direction, track, gap) defines the output port
               position and direction are None: this is calculated from
               other parameters
        iLength: (float) length of pads
        iWidth: (float) width of pads
        iSize: (float) spacing between pads (see drawing)

        Outputs:
        --------
        retIn: same as iIn, with flipped vector
        retOut: calculated output port to match all input dimensions

                 iSize
              +--+  +--+
              |  |  |  |
            +-+  |  |  +-+
        iIn |    |  |    | iOut
            +-+  |  |  +-+
              |  |  |  |
              +--+  +--+
        '''



        iLength, iWidth,iSize = parse_entry(iLength, iWidth, iSize)
        retIn = [[0,0,0], 0, iInPort.track, iInPort.gap]
        retOut = [[iInPort.gap+iOutPort.gap+iSize+2*iWidth,0],0 , iOutPort.track, iOutPort.gap]

        points1 = self.append_points([(iInPort.gap+iWidth, 0),
                                     (0, -iLength/2),
                                     (-iWidth, 0),
                                     (0, iLength/2-iInPort.track/2),
                                     (-iInPort.gap, 0),
                                     (0, iInPort.track),
                                     (iInPort.gap, 0),
                                     (0, iLength/2-iInPort.track/2),
                                   (iWidth, 0)])


        trackIn = self.polyline_2D(points1, name=name+'_track1', layer=layer_TRACK)
#        self.chip.trackObjects.append(trackIn)

        points2 = self.append_points([(iInPort.gap+iWidth+iSize, 0),
                                     (0, -iLength/2),
                                     (+iWidth, 0),
                                     (0, iLength/2-iOutPort.track/2),
                                     (+iOutPort.gap, 0),
                                     (0, iOutPort.track),
                                     (-iOutPort.gap, 0),
                                     (0, iLength/2-iOutPort.track/2),
                                     (-iWidth, 0)])
        trackOut = self.polyline_2D(points2, name=name+'_track2', layer=layer_TRACK)
#        self.chip.trackObjects.append(trackOut)

        points3 = self.append_points([(0, 0),
                                     (0, iLength/2+iInPort.gap),
                                     (iInPort.gap+iWidth+iSize/2, 0),
                                     (0, iOutPort.gap-iInPort.gap),
                                     (iOutPort.gap+iWidth+iSize/2, 0),
                                     (0, -iLength-2*iOutPort.gap),
                                     (-(iOutPort.gap+iWidth+iSize/2),0),
                                     (0, iOutPort.gap-iInPort.gap),
                                     (-(iInPort.gap+iWidth+iSize/2),0)
                                     ])
        gap1 = self.polyline_2D(points3, name=name+'_gap1', layer=layer_GAP)
#        self.chip.gapObjects.append(gap1)
#
#



#        TODO !!
#        if not self.is_litho:
#            self.draw(self.name+"_mesh", points)
#            self.modeler.assign_mesh_length(self.name+"_mesh",1/2*iLength)

        self.iIn = retIn
        self.iOut = retOut
#        return [retIn, retOut]


    def find_slanted_path(self, name, iInPort, iOutPort):

        iIn_pos = Vector(iInPort.pos)
        iIn_ori = Vector(iInPort.ori)
        iOut_pos = Vector(iOutPort.pos)
        iOut_ori = Vector(iOutPort.ori)

        if iIn_ori.dot(iOut_ori)!=-1:
            raise ValueError('Cannot find slanted path: ports are not oriented correctly')
        else:
            dist = (iOut_pos-iIn_pos).dot(iIn_ori)

        pointA = iIn_pos+iIn_ori*dist/3
        pointB = iOut_pos+iOut_ori*dist/3

        # TODO
        self.to_bond.append([iIn_pos, pointA])
        self.to_bond.append([pointB, iOut_pos])
        return [iIn_pos, pointA, pointB, iOut_pos], dist/3






@Lib.add_methods_from(KeyElement, CustomElement)
class Body(PythonModeler):

    def __init__(self, interface, coor_sys, name, network, mode, variables):
        self.interface = interface
        self.coor_sys = coor_sys
        self.name = name
        self.maskObjects = []
        self.trackObjects = []
        self.gapObjects = []
        self.mode = mode # 'hfss' or 'gds'
        self.variables = variables
        network.update(coor_sys)
        self.network = network



#    @staticmethod
#    def find_last_list(list_entities):
##        print(list_entities)
##        if isinstance(list_entities, Port) :
##            return None
#        if isinstance(list_entities, list):
#            if len(list_entities)==0:
#                return list_entities
#            else:
#                if isinstance(list_entities[-1], list):
#                    return Body.find_last_list(list_entities[-1])
#                else:
#                    return list_entities
#        else:
#            return list_entities
    def move_port(func):
        @wraps(func)
        def moved(*args, **kwargs):
            new_args = [args[0], args[1]]
            #  args[0] = PM
            #  args[1] = name
            compteur = 0
            for i, argument in enumerate(args[2:]):
                if isinstance(argument, str) and (argument in Port.dict_instances):
                    #  if argument is the sting representation of the port
                    new_args.append(Port.dict_instances[argument])
                    compteur+=1
                elif isinstance(argument, Port):
                    #  it the argument is the port itself
                    new_args.append(argument)
                    compteur+=1
                else:
                    new_args.append(argument)
                # else:
                #     error = '%s arg should be a port'%str(argument)
                #     raise Exception(error)
#            print("compteur",compteur)
            if compteur==0:
                raise Exception("Please indicate more than 0 port")

            previous_pos = args[0].current_pos
            previous_ori = args[0].current_ori
            #TODO
            if func.__name__=='draw_cable':
                #TODO It depends of a parameter of drawCable
                args[0].set_current_coor([0,0],[1,0])
            elif func.__name__=='find_path':
                args[0].set_current_coor([0,0],[1,0])

            #  the following is not robust
            elif compteur==1:
#                print(new_args[2].pos)
                args[0].set_current_coor(new_args[2].pos, new_args[2].ori)
            elif compteur==2:
                args[0].set_current_coor(1/2*(new_args[2].pos+new_args[3].pos), new_args[2].ori)
            new_args = tuple(new_args)
            return KeyElement._moved(func, previous_pos, previous_ori, *new_args, **kwargs)
        return moved

    def port(self, name, pos, ori, widths, subnames, layers, offsets, constraint_port):
        if constraint_port:
            pos, ori = parse_entry(pos, ori)
            offset=0
            width=50e-6  # 50um
            points = [(0, offset+width/2),
                      (width/3, offset),
                      (0, offset-width/2)]
            self.polyline_2D(points, name='_'+name, layer=layer_PORT)
        else:
            pos, ori, widths, offsets = parse_entry(pos, ori, widths, offsets)
            for ii in range(len(widths)):
                width = widths[ii]
                offset = offsets[ii]
                points = [(0, offset+width/2),
                          (width/3, offset),
                          (0, offset-width/2)]
                self.polyline_2D(points, name='_'+name+'_'+subnames[ii], layer=layer_PORT)

        return self.network.port(name, pos, ori, widths, subnames, layers, offsets, constraint_port)

    def double_port(self, name, pos, ori, track, gap):
        port1 = self.port(name+'_front', pos, ori, track, gap)
        port2 = self.port(name+'_back', pos, -Vector(ori), track, gap)
        return port1, port2

    @staticmethod
    def number_ports():
        return len(Port.instances_to_move)

    @staticmethod
    def number_modelentities():
        return len(ModelEntity.instances_to_move)

    def translate_port(self, ports, vector):
        for port in ports:

            port.pos = port.pos+Vector(vector)

    def rotate_port(self, ports, angle):
        if isinstance(angle, list):
            if len(angle)==2:
                new_angle= np.math.atan2(np.linalg.det([[1,0],angle]),np.dot([1,0],angle))
                new_angle= new_angle/np.pi*180
            else:
                raise Exception("angle should be either a float or a 2-dim array")
        else :
            new_angle=angle
        rad = new_angle/180*np.pi
        rotate_matrix = np.array([[np.cos(rad) ,np.sin(-rad)],[np.sin(rad) ,np.cos(rad)]])
        for port in ports:
            port.ori = rotate_matrix.dot(port.ori)
            posx = port.pos[0]*np.cos(rad)+port.pos[1]*np.sin(-rad)
            posy = port.pos[0]*np.sin(rad)+port.pos[1]*np.cos(rad)
            port.pos = Vector([posx, posy])
#    @staticmethod
#    def ports_to_move(n):
#        lastP = OrderedDict()
#        for idx, key in enumerate(Port.instances_to_move):
#            if n <= idx:
#                lastP[key] = Port.instances_to_move[key]
#        return lastP
#
#    @staticmethod
#    def modelentities_to_move(n):
#        lastM = OrderedDict()
#        for idx, key in enumerate(ModelEntity.instances_to_move):
#            if n <= idx:
#                lastM[key] = ModelEntity.instances_to_move[key]
#        return lastM


#    @staticmethod
#    def ports_to_move():
#        return Port.instances_to_move
#
#    @staticmethod
#    def modelentities_to_move():
#        return ModelEntity.instances_to_move
#
    @staticmethod
    def ports_were_moved(n):
        Port.merge(n)
    @staticmethod
    def modelentities_were_moved(n):
        ModelEntity.merge(n)

    @move_port
    def _connect_JJ(self, name, iInPort, iOutPort, iTrackJ, iLengthJ=None, iInduct='1nH', fillet=None):
        '''
        Draws a Joseph's Son Junction.

        Draws a rectangle, here called "junction",
        with Bondary condition :lumped RLC, C=R=0, L=iInduct in nH
        Draws needed adaptors on each side

        Inputs:
        -------
        name:
        iIn: (tuple) input port
        iOut: (tuple) output port
        iSize: (float) length of junction
        iWidth: (float) width of junction
        iLength: (float) distance between iIn and iOut, including
                 the adaptor length
        iInduct: (float in nH)

        Outputs:
        --------

        '''
        iLength = (iOutPort.pos-iInPort.pos).norm()
        if iLengthJ is None:
            iLengthJ = iTrackJ
        iTrackJ = parse_entry(iTrackJ)

        if 0:
            induc_H = val(iInduct)*1e-9
            print('induc'+str(induc_H))
            w_plasma = 2*np.pi*24*1e9
            capa_fF = 1/(induc_H*w_plasma**2)*10**15
            capa_plasma = self.set_variable(name+'_capa_plasma', str(capa_fF)+'fF')
            print(capa_fF)
            print(capa_plasma)
        else:
            capa_plasma = 0

        # No parsing needed, should not be called from outside
        iTrack1, iTrack2 = parse_entry(iInPort.track, iOutPort.track)

        adaptDist1 = iTrack1/2-iTrackJ/2
        adaptDist2 = iTrack2/2-iTrackJ/2
        if val(adaptDist1)>val(iLength/2-iTrackJ/2) or val(adaptDist2)>val(iLength/2-iTrackJ/2):
            raise ValueError('Increase iTrackJ %s' % name)

        raw_points = [(iLengthJ/2, iTrackJ/2),
                      ((iLength/2-iLengthJ/2-adaptDist1), 0),
                      (adaptDist1, (iTrack1-iTrackJ)/2),
                      (0, -iTrack1),
                      (-adaptDist1, (iTrack1-iTrackJ)/2),
                      (-(iLength/2-iLengthJ/2-adaptDist1), 0)]
        points = self.append_points(raw_points)
        right_pad = self.polyline_2D(points, name =name+"_pad1", layer=layer_TRACK)

        raw_points = [(iLengthJ/2, iTrackJ/2),
                      ((iLength/2-iLengthJ/2-adaptDist2), 0),
                      (adaptDist2, (iTrack2-iTrackJ)/2),
                      (0, -iTrack2),
                      (-adaptDist2, (iTrack2-iTrackJ)/2),
                      (-(iLength/2-iLengthJ/2-adaptDist2), 0)]
        points = self.append_points(self.refy_points(raw_points))
        left_pad = self.polyline_2D(points, name=name+"_pad2", layer=layer_TRACK)

        pads = self.unite([right_pad, left_pad], name=name+'_pads')

        if not self.is_litho:
            if val(iTrack1) > val(iTrack2):
                mesh = self.rect_center_2D([0,0], [iLength, iTrack1], name=name+'_mesh', layer=layer_MESH)
            else:
                mesh = self.rect_center_2D([0,0], [iLength, iTrack2], name=name+'_mesh', layer=layer_MESH)
            self.mesh_zone(mesh, iTrackJ/2)

            points = self.append_points([(iTrackJ/2,0),(-iTrackJ,0)])
            flow_line = self.polyline_2D(points, closed = False, name=name+'_flowline', layer=layer_MESH)
            JJ = self.rect_center_2D([0,0], [iLengthJ, iTrackJ], name=name+'_lumped', layer=layer_RLC)
            #TODO
            #self.assign_lumped_RLC(JJ,0, iInduct, capa_plasma, flow_line)

        return pads

    @move_port
    def draw_cable(self, name, *ports, fillet="0.3mm", is_bond=False, to_meander=[[]], meander_length=0, meander_offset=0, is_mesh=False, reverse_adaptor=False):
        '''
        Draws a CPW transmission line along the ports

        Be careful: if the ports are facing eachother and offset by a small distance, sometimes the cable cannot be drawn.

        Inputs:
        -------
        name: (string) base-name of object, draws 'name_adaptor' etc
        iIn: (tuple) input port
        iOut: (tuple) output port
        iMaxfillet: (float), maximum fillet radius
        reverseZ: performs a mirror operation along Z --> useful only when the thickening operation goes in the wrong direction

        '''
        #TODO change the format of the arguments
        meander_length, meander_offset, fillet = parse_entry(meander_length, meander_offset, fillet)

        # to_meander should be a list of list
        # meander_length, meander_offset should be lists
        if not isinstance(to_meander[0], list):
            to_meander = [to_meander]
        if not isinstance(meander_length, list):
            meander_length = [meander_length]*len(to_meander)
        if not isinstance(meander_offset, list):
            meander_offset = [meander_offset]*len(to_meander)

        ports = list(ports)

        # asserts neither in nor out port are constraint_ports
        if ports[0].constraint_port and ports[-1].constraint_port:
            raise ValueError('At least the first (%s) or last port (%s) \
                             should define the port parameters'%(ports[0].name,
                                                             ports[-1].name))
        elif ports[0].constraint_port:
            ports[0].widths = ports[-1].widths
            ports[0].offsets = ports[-1].offsets
            ports[0].layers = ports[-1].layers
            ports[0].subnames = ports[-1].subnames
            ports[0].N = ports[-1].N
            ports[0].constraint_port = False  # ports[0] is now defined
        elif ports[-1].constraint_port:
            ports[-1].widths = ports[0].widths
            ports[-1].offsets = ports[0].offsets
            ports[-1].layers = ports[0].layers
            ports[-1].subnames = ports[0].subnames
            ports[-1].N = ports[0].N
            ports[-1].constraint_port = False  # ports[-1] is now defined
        else:
            pass

        # recursive approach if there are intermediate non_constraint ports

        cable_portion = 0
        _ports = [ports[0]]
        for port in ports[1:-1]:
            _ports.append(port)
            if not port.constraint_port:
                print(to_meander[cable_portion])
                print(meander_length[cable_portion])
                self.draw_cable(name+'_%d'%cable_portion, *_ports,
                                fillet=fillet, is_bond=is_bond,
                                to_meander=[to_meander[cable_portion]],
                                meander_length=[meander_length[cable_portion]],
                                meander_offset=[meander_offset[cable_portion]],
                                is_mesh=is_mesh,
                                reverse_adaptor=reverse_adaptor)
                cable_portion += 1
                _ports = [port.r]
        if cable_portion != 0:
            name = name+'_%d'%cable_portion
            to_meander = [to_meander[cable_portion]]
            meander_length = [meander_length[cable_portion]]
            meander_offset = [meander_offset[cable_portion]]
        _ports.append(ports[-1])


        # at this stage first and last port are not constraint_port and all
        # intermediate port should be constraint_port

        ports = _ports

        # find and plot adaptor geometry
        if reverse_adaptor:
            points, length_adaptor = ports[-1].compare(ports[0], self)
            index_modified = -1
        else:
            points, length_adaptor = ports[0].compare(ports[-1], self)
            index_modified = 0

        # plot adaptors
        for jj, pts in enumerate(points):
            self.polyline_2D(pts, name=ports[index_modified].name+'_'+ports[index_modified].subnames[jj]+'_adapt', layer=ports[index_modified].layers[jj])

        # define the constraint_port parameters
        for port in ports[1:-1]:
            port.widths = ports[0].widths
            port.offsets = ports[0].offsets
            port.layers = ports[0].layers
            port.subnames = ports[0].subnames
            port.N = ports[0].N

        # find all intermediate paths
        total_path = None
        for ii in range(len(ports)-1):
            path = Path(name, ports[ii], ports[ii+1], fillet)
            if total_path is None:
                total_path = path
            else:
                total_path += path
            ports[ii+1] = ports[ii+1].r # reverse the last port

        total_path.clean()

        # do meandering
        total_path.meander(to_meander[0], meander_length[0], meander_offset[0])

        total_path.clean()

        # plot cable
        self.path(total_path.points, total_path.port_in, total_path.fillet,
                  name=name)

        # if bond plot bonds
        if is_bond:
            self.draw_bond(total_path.to_bond(), *ports[0].bond_params(), name=name+'_wb_%d'%ii)

        ports[0].revert()
        ports[-1].revert()

        length = total_path.length() + length_adaptor
        print('Cable "%s" length = %.3f mm'%(name, length*1000))
        return length

    def draw_bond(self, to_bond, ymax, ymin, name='wb', min_dist='0.5mm'):
        # to_bond list of segments
        ymax, ymin, min_dist = parse_entry(ymax, ymin, min_dist)

        min_dist = val(min_dist)
        n_segments = len(to_bond)
        jj=0
        while jj<n_segments:
            elt = to_bond[jj]
            A = elt[0]
            B = elt[1]
            if jj+1 < n_segments and to_bond[jj+1][0] == B:
                B = to_bond[jj+1][1]
                jj+=1
            val_BA = val(B-A)
            ori = way(val_BA)
            length = Vector(val_BA).norm()
            n_bond = int(length/_val(min_dist))+1
            spacing = (B-A).norm()/n_bond
            pos = A+ori*spacing/2
            self.wirebond_2D(pos, ori, ymax, ymin, layer=layer_Default, name=name+'_0')
            for ii in range(n_bond-1):
                pos = pos + ori*spacing
                self.wirebond_2D(pos, ori, ymax, ymin, layer=layer_Default, name=name+'_%d'%ii)
            jj+=1

    @move_port
    def _connect_snails2(self, name, iInPort, iOutPort, squid_size, width_top, n_top, width_bot, n_bot, N, width_bridge, spacing_bridge, litho='opt'):

        width_track = iInPort.track # assume both are equal
        spacing = (iOutPort.pos-iInPort.pos).norm()

        width_snail = squid_size[0]+6*width_track #ZL

        tot_width = width_snail*N
        if np.abs(val(tot_width))>np.abs(val(spacing)):
            raise ValueError("cannot put all snails in given space")
        offset = 5.0e-6-0.630e-6
        overlap=0
        if litho=='elec':
            overlap = 5e-6
        snails = []
        snails.append(self.rect_corner_2D([-tot_width/2, -width_track/2], [-(spacing-tot_width)/2-overlap,width_track], name=self.name+'_pad_left', layer=layer_TRACK))
        snails.append(self.rect_corner_2D([tot_width/2, -width_track/2], [(spacing-tot_width)/2+overlap,width_track], name=self.name+'_pad_right', layer=layer_TRACK))
        if N%2==1:
            x_pos=-(N//2)*width_snail
        else:
            x_pos=-(N//2-1/2)*width_snail

        for jj in range(int(N)):
            snail=[]
            snail.append(self.rect_corner_2D([x_pos-squid_size[0]/2, -0.1e-6+(-squid_size[1]/2-width_bot+0.1e-6)-((-squid_size[1]/2-width_bot+0.1e-6)+width_track/2)-offset], [-3*width_track/2, squid_size[1]+width_top+width_bot+2*0.1e-6], name=self.name+'_left', layer=layer_TRACK))
            snail.append(self.rect_corner_2D([x_pos+squid_size[0]/2, -squid_size[1]/2-((-squid_size[1]/2)+width_track/2)-offset], [3*width_track/2, squid_size[1]+width_top+width_bot+0.1e-6], name=self.name+'_right', layer=layer_TRACK))
            for width, n, way in [[width_top, n_top, 1], [width_bot, n_bot, -1]]:
                if n==1:
                    snail.append(self.rect_corner_2D([x_pos-squid_size[0]/2, way*squid_size[1]/2-((-squid_size[1]/2-width_bot-0.1e-6)+width_track/2)-offset], [(squid_size[0]-width_bridge)/2, way*width-2*0.1e-6], name=name+'_islandtop_left', layer=layer_TRACK))
                    snail.append(self.rect_corner_2D([x_pos+squid_size[0]/2, way*squid_size[1]/2-((-squid_size[1]/2-width_bot)+width_track/2)-offset], [-(squid_size[0]-width_bridge)/2, way*width], name=name+'_islandtop_right', layer=layer_TRACK))
                if n==3: #TODO
                    snail.append(self.rect_corner_2D([x_pos-squid_size[0]/2, way*squid_size[1]/2-((-squid_size[1]/2-width_bot+0.1e-6)+width_track/2)-offset], [(squid_size[0]-2*width_bridge-spacing_bridge)/2, way*width+2*0.1e-6], name=name+'_islandtop_left', layer=layer_TRACK))
                    snail.append(self.rect_corner_2D([x_pos+squid_size[0]/2, way*squid_size[1]/2-((-squid_size[1]/2-width_bot+0.1e-6)+width_track/2)-offset], [-(squid_size[0]-2*width_bridge-spacing_bridge)/2, way*width+2*0.1e-6], name=name+'_islandtop_right', layer=layer_TRACK))
                    snail.append(self.rect_corner_2D([x_pos-spacing_bridge/2, way*squid_size[1]/2-((-squid_size[1]/2-width_bot)+width_track/2)-offset], [spacing_bridge, way*width], name=name+'_island_middle_', layer=layer_TRACK) )
                if n==4:
                    snail.append(self.rect_corner_2D([x_pos-squid_size[0]/2, way*squid_size[1]/2-((-squid_size[1]/2-width_bot+0.1e-6)+width_track/2)-offset], [(squid_size[0]-3*width_bridge-2*spacing_bridge)/2, way*width+2*0.1e-6], name=name+'_islandtop_left', layer=layer_TRACK))
                    snail.append(self.rect_corner_2D([x_pos+squid_size[0]/2, way*squid_size[1]/2-((-squid_size[1]/2-width_bot)+width_track/2)-offset], [-(squid_size[0]-3*width_bridge-2*spacing_bridge)/2, way*width], name=name+'_islandtop_right', layer=layer_TRACK))
                    snail.append(self.rect_corner_2D([x_pos-spacing_bridge-width_bridge/2, way*squid_size[1]/2-((-squid_size[1]/2-width_bot)+width_track/2)-offset], [spacing_bridge, way*width], name=name+'_island_middle_left', layer=layer_TRACK) )
                    snail.append(self.rect_corner_2D([x_pos+spacing_bridge+width_bridge/2, way*squid_size[1]/2-((-squid_size[1]/2-width_bot+0.1e-6)+width_track/2)-offset], [-spacing_bridge, way*width+2*0.1e-6], name=name+'_island_middle_right', layer=layer_TRACK) )
            snail.append(self.rect_corner_2D([x_pos-squid_size[0]/2-3*width_track/2, -width_track/2], [-1.5*width_track, width_track], name=name+'_connect_left', layer=layer_TRACK)) #ZL
            snail.append(self.rect_corner_2D([x_pos+squid_size[0]/2+3*width_track/2, -width_track/2], [1.5*width_track, width_track], name=name+'_connect_right', layer=layer_TRACK))

            self.unite(snail, name=self.name+'_snail_'+str(jj))
            x_pos = x_pos+width_snail

    @move_port
    def _connect_jct(self, name, iInPort, iOutPort, width_bridge, n=1, spacing_bridge=0, assymetry=0.1e-6, overlap=5e-6, width_jct=None, thin=False): #opt assymetry=0.25e-6
        limit_dose = 8e-6
        width = iInPort.track # assume both are equal
        spacing = (iOutPort.pos-iInPort.pos).norm()
#        self.pos = (self.pos+self.posOut)/2
        n = int(n)

        tot_width = n*width_bridge+(n-1)*spacing_bridge

        if width_jct is not None:
            margin = 1e-6

        self.rect_corner_2D([-tot_width/2-limit_dose,-width/2-assymetry], [-(spacing-tot_width)/2-overlap+limit_dose, width+2*assymetry], name=name+'_left', layer=layer_TRACK)
        if thin and width_jct is not None:
            self.rect_corner_2D([-tot_width/2-margin,-width/2-assymetry], [-limit_dose+margin, width+2*assymetry], name=name+'_left2', layer=layer_TRACK)
            self.rect_corner_2D([-tot_width/2,-width_jct/2],[-margin, width_jct], name=name+'_left3', layer=layer_TRACK)
        else:
            self.rect_corner_2D([-tot_width/2,-width/2-assymetry], [-limit_dose, width+2*assymetry], name=name+'_left2', layer=layer_TRACK)
        #print(n)
        if n%2==0:
            _width_right = width+2*assymetry
        else:
            _width_right = width
        if width_jct is not None:
            self.rect_corner_2D([tot_width/2+margin,-_width_right/2], [limit_dose-margin, _width_right], name=name+'right2', layer=layer_TRACK)
            self.rect_corner_2D([tot_width/2,-width_jct/2], [margin, width_jct], name=name+'_right3', layer=layer_TRACK)
        else:
            self.rect_corner_2D([tot_width/2,-_width_right/2], [limit_dose, _width_right], name=name+'_right2', layer=layer_TRACK)

        self.rect_corner_2D([tot_width/2+limit_dose,-_width_right/2], [(spacing-tot_width)/2+overlap-limit_dose, _width_right], name=name+'_right', layer=layer_TRACK)

        x_pos = -(tot_width)/2+width_bridge

        for ii in range(n-1):
            if ii%2==1:
                _width_loc = width+2*assymetry
            else:
                _width_loc = width
            self.rect_corner_2D([x_pos,-_width_loc/2], [spacing_bridge, _width_loc], name=name+'_middle'+str(ii), layer=layer_TRACK)
            x_pos = x_pos+spacing_bridge+width_bridge

    @move_port
    def _connect_array(self, name, iInPort, iOutPort, width_bridge, spacing_bridge, n=1):
        width = iInPort.track # assume both are equal
        spacing = (iOutPort.pos-iInPort.pos).norm()
        n = int(n)

        tot_width = n*width_bridge+(n-1)*spacing_bridge

        self.rect_corner_2D([-tot_width/2,-width/2], [-(spacing-tot_width)/2, width], name=name+'_left', layer=layer_TRACK)
        self.rect_corner_2D([tot_width/2,-width/2], [(spacing-tot_width)/2, width], name=name+'_right', layer=layer_TRACK)

        x_pos = -(tot_width)/2+width_bridge

        for ii in range(n-1):
            self.rect_corner_2D([x_pos,-width/2], [spacing_bridge, width], name=name+'_middle'+str(ii), layer=layer_TRACK)
            x_pos = x_pos+spacing_bridge+width_bridge


    @move_port
    def cable_starter(self, name, iInPort, width = 'track', index=None, border=parse_entry('15um')): # width can also be 'gap'
        if width=='track' or width=='Track':
            points = [(0,iInPort.track/2), (0,-iInPort.track/2)]
        elif width=='gap' or width=='Gap':
            points = [(0,iInPort.gap+iInPort.track/2),(0,-iInPort.gap-iInPort.track/2)]
        elif width=='mask' or width=='Mask':
            points = [(0,iInPort.gap+iInPort.track/2+self.gap_mask),(0,-iInPort.gap-iInPort.track/2-self.gap_mask)]
        if index == None:
            track = self.polyline_2D(points, closed=False, name=name+'_starter_'+width, layer=layer_TRACK)
        else:
            track = self.polyline_2D(points, closed=False, name=name+'_starter_'+width, layer=layer_TRACK)

         #TODO Find unknown variables
    #        elif width=='dc_track':
    #            points = self.append_absolute_points([(0, self.rel_posIn[index]-self.widIn[index]/2),\
    #                                                  (0, self.rel_posIn[index]+self.widIn[index]/2)])
    #        elif width=='dc_gap':
    #            points = self.append_absolute_points([(0, self.rel_posIn[index]-self.widIn[index]/2-border),\
    #                                                  (0, self.rel_posIn[index]+self.widIn[index]/2+border)])
    #        elif width=='dc_cutout':
    #            points = self.append_absolute_points([(0, -self.cutIn/2),(0, self.cutIn/2)])

        return track

    @move_port
    def _connect_jct_corrected(self, width_bridge, iInduct='0nH', n=1, spacing_bridge=0, assymetry=0.1e-6, overlap=None, width_jct=None, thin=False, cross=False, override=False, dose=False, way=1): #opt assymetry=0.25e-6
        '''
        This function is very similar to _connect_jct
        It is called in draw_dose_test_jct
        In case of single junction not crossed, it does not draw rectangle

        Before in _connect_jct
        --------------
                |-----|-------
                |-----|-------
        -------------

        After in _connect_jct_corrected
        --------------
                |-------
                |-------
        -------------

        The rectangle of the finger of the junction does not go on the big rectangle
        The rectangle at the end of the big rectangle has been deleted
        Each aera for the junction is just exposed once.
        '''

        limit_dose = 1e-6
        width = self.inTrack # assume both are equal
        spacing = (self.posOut-self.pos).norm()
        #spacing = (iInPort.pos-iOutPort.pos).norm()

        self.pos = (self.pos+self.posOut)/2
        n = int(n)

        if width_jct is not None:
            margin = 2e-6
        if overlap is None:
            overlap = 0.0

        if cross and n==1:
            tot_width = 1.5*margin + width_bridge + width_jct
        else:
            tot_width = n*width_bridge+(n-1)*spacing_bridge

        if (self.is_litho and not override) or dose:
            if cross:
                self.draw_rect(self.name+'_left', self.coor([-tot_width/2,-width/2-assymetry]), self.coor_vec([-(spacing-tot_width)/2-overlap, width+2*assymetry]))
                self.draw_rect(self.name+'_detour', self.coor([-tot_width/2,-way*(margin-width_bridge+0.5*width_jct+width/2)-assymetry]), self.coor_vec([-width,way*(margin-width_bridge+0.5*width_jct)]))
                self.draw_rect(self.name+'_detour1', self.coor([-tot_width/2,-way*(margin-width_bridge+0.5*width_jct+width/2)-assymetry]), self.coor_vec([margin,way*width_jct]))
            else:
                self.draw_rect(self.name+'_left', self.coor([-tot_width/2-limit_dose,-width/2-assymetry]), self.coor_vec([-(spacing-tot_width)/2-overlap+limit_dose, width+2*assymetry]))
            if thin and width_jct is not None:
        #                self.draw_rect(self.name+'_left2', self.coor([-tot_width/2-margin,-width/2-assymetry]), self.coor_vec([-limit_dose+margin, width+2*assymetry]))
                self.draw_rect(self.name+'_left3', self.coor([-tot_width/2,-width_jct/2]), self.coor_vec([-margin+limit_dose, width_jct]))
            elif not cross:
                self.draw_rect(self.name+'_left2', self.coor([-tot_width/2,-width/2-assymetry]), self.coor_vec([-limit_dose, width+2*assymetry]))

            if n%2==0:
                _width_right = width+2*assymetry
            else:
                _width_right = width
            if width_jct is not None and not cross:
        #                self.draw_rect(self.name+'_right2', self.coor([tot_width/2+margin,-_width_right/2]), self.coor_vec([limit_dose-margin, _width_right]))
                self.draw_rect(self.name+'_right3', self.coor([tot_width/2,-width_jct/2]), self.coor_vec([margin-limit_dose, width_jct]))
            elif not cross:
                self.draw_rect(self.name+'_right2', self.coor([tot_width/2,-_width_right/2]), self.coor_vec([limit_dose, _width_right]))

            if cross:
                self.draw_rect(self.name+'_right', self.coor([tot_width/2-0.5*margin,-_width_right/2]), self.coor_vec([(spacing-tot_width)/2+overlap+0.5*margin, _width_right]))
                self.draw_rect(self.name+'_detour2', self.coor([tot_width/2-0.5*margin-width_jct,way*(_width_right/2)]), self.coor_vec([width_jct,-way*(margin+width)]))
            else:
                self.draw_rect(self.name+'_right', self.coor([tot_width/2+limit_dose,-_width_right/2]), self.coor_vec([(spacing-tot_width)/2+overlap-limit_dose, _width_right]))

            x_pos = -(tot_width)/2+width_bridge

            for ii in range(n-1):
                if ii%2==1:
                    _width_loc = width+2*assymetry
                else:
                    _width_loc = width
                self.draw_rect(self.name+'_middle', self.coor([x_pos,-_width_loc/2]), self.coor_vec([spacing_bridge, _width_loc]))
                x_pos = x_pos+spacing_bridge+width_bridge

        elif not override:
            mesh = self.draw_rect_center(self.name+'_mesh', self.coor([0,0]), self.coor_vec([spacing, width]))
            self.modeler.assign_mesh_length(mesh, 0.5*width)

            JJ = self.draw_rect_center(self.name, self.coor([0,0]), self.coor_vec([spacing, width]))
            self.assign_lumped_RLC(JJ, self.ori, (0, iInduct, 0))