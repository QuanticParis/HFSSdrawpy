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
                             var, \
                             Vector

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

def way(vec):
    if vec[1] != 0:
        if abs(vec[0]/vec[1])<1e-2:
            if vec[1]>0:
                return Vector(0,1)
            elif vec[1]<0:
                return Vector(0,-1)
    if vec[0] != 0 :
        if abs(vec[1]/vec[0])<1e-2:
            if vec[0]>0:
                return Vector(1,0)
            elif vec[0]<0:
                return Vector(-1,0)

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

# useful function to find cable path
def next_point(point1, point2, vec):
    choice1 = point1+vec*var((point2-point1).dot(vec))
    choice2 = point1+vec.orth()*var((point2-point1).dot(vec.orth()))
    return [[point1, choice1, point2], [point1, choice2, point2]]

def check(points):
    # points : series of points that represent the cable path
    # returns : the cost (number of turns) and make sure each point is at a
    # turn ie that there are no consecutive segments in the same direction.
    length = 0
    prev_point = Vector(points[0])
    _points = [Vector(points[0])]
    vecs = []
    for point in points[1:]:
        if not ( equal_float(point[0], prev_point[0]) and equal_float(point[1], prev_point[1])):
           # TODO vec should already be a Vector
           vec = point-prev_point
           length += vec.norm()
           vecs.append(way(vec))

           prev_point = point
           _points.append(point)
    cost = 0

    points = _points.copy()

    new_points = [points[0]]
    prev_vec = vecs[0]
    for ii, vec in enumerate(vecs[1:]):
        curr_vec = vec
        if curr_vec.dot(prev_vec)==0:
            new_points.append(points[ii+1])
        added_cost = cost_f(prev_vec.dot(curr_vec))
        cost += added_cost
        prev_vec = curr_vec
    new_points.append(points[-1])

    return cost, new_points, length

def cost_f(x):
    if x==1:
        return 0
    elif x==0:
        return 1
    else:
        return 100

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
            raise NotImplementedError()


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

    def set_variable(self, value):
        """
        name (str): name of the variable in HFSS e.g. 'chip_length'
        value (str, VarStr, float): value of the variable
                                    if str will try to analyse the unit
        """
        f = currentframe().f_back#.f_back
        filename = getfile(f)
        code_line = open(filename).readlines()[f.f_lineno - 1]
        name = code_line.split("=")[0].strip()

        # self.store_variable(name, value)  # for later parsing
        # self.__dict__[name] = VariableString(name)
        # for handy use by user / deprecated
        if self.mode == 'hfss':
            self.design.set_variable(name, value)  # for HFSS
        return VariableString(name, value=value)

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

    def compare(self, other):
        points = []
        adapt_dist = VariableString(self.name+'_adapt', value=1e-5)
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
        VariableString.variables[self.name+'_adapt'] = 2*max_diff
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
        return points

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
    def draw_adaptor(self, name, iInPort, track, gap, iSlope=0.33):
        '''
        Draws an adaptor between two ports.
        Given input port iIn, and slope of line iSlope, calculates iOut, and draws adpator.

        Inputs:
        -------
        name:
        iIn: tuple, input port
        iOut: tuple, output port, usually Nones
        iSlope: slope of line to connect iIn and iOut. If iSlope=1, 45 degrees.

        Returns:
        --------
        reversed iIn and calculated iOut
        '''
#        if not self.isOut:
#            # calculate the output
#            # do not forget to add the new port to dict
        track, gap = parse_entry(track, gap)
        adaptDist = abs(track/2-iInPort.track/2)/iSlope
#            outPort = [self.pos+self.ori*adaptDist, self.ori, self.outPort['track'], self.outPort['gap']]
#            self.ports[self.iIn+'_bis'] = outPort
#            self.__init__(self.name, self.iIn, self.iIn+'_bis')
#        else:

#        adaptDist = (iInPort.pos-iOutPort.pos).norm()

        points = self.append_points([(0, iInPort.track/2),
                                     (adaptDist, track/2-iInPort.track/2),
                                     (0, -track),
                                     (-adaptDist, track/2-iInPort.track/2)])
        track_points = self.polyline_2D(points, name=name+"_track", layer = layer_TRACK)
#        self.trackObjects.append(track)
        points = self.append_points([(0, iInPort.gap+iInPort.track/2),
                                     (adaptDist, (gap-iInPort.gap)+(track-iInPort.track)/2),
                                     (0, -2*gap-track),
                                     (-adaptDist, (gap-iInPort.gap)+(track-iInPort.track)/2)])
        gap_points = self.polyline_2D(points, name=name+"_GAP", layer = layer_GAP)
#        self.gapObjects.append(gap)

        mask = None
        if self.is_mask:
            points = self.append_points([(0, self.gap_mask+iInPort.gap+iInPort.track/2),
                             (adaptDist, (gap-iInPort.gap)+(track-iInPort.track)/2),
                             (0, -2*self.gap_mask-2*gap-track),
                             (-adaptDist, (gap-iInPort.gap)+(track-iInPort.track)/2)])
            mask = self.polyline_2D(points, name=name+"_MASK", layer = layer_MASK)
#            self.maskObjects.append(mask)
        #TODO create port out
        return iInPort.name+'_bis', adaptDist, track_points, gap_points, mask

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
    def draw_cable(self, name, *ports, fillet="0.3mm", is_bond=False, to_meanders=None, meander_length=0, meander_offset=0, is_mesh=False, reverse_adaptor=False):
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
        meander_length, meander_offset, fillet =parse_entry(meander_length, meander_offset, fillet)

        to_bond=[]

        ports = list(ports)

        # TODO check port compatibility
        # should return a set of trapeze points

        # inTrack = self.variables[ports[0].track]
        # outTrack = self.variables[ports[-1].track]
        # inGap = self.variables[ports[0].gap]
        # outGap = self.variables[ports[-1].gap]

        # if (not equal_float(inTrack, outTrack)) or (not equal_float(inGap, outGap)):
        #     if reverse_adaptor:
        #         if (inTrack+inGap) > (outTrack+outGap):
        #             iIn, adaptor_length, track_adaptor, gap_adaptor, mask_adaptor = self.draw_adaptor('reverse',ports[0].name,ports[-1].track, ports[-1].gap)
        #         else:
        #             iOut, adaptor_length, track_adaptor, gap_adaptor, mask_adaptor = self.draw_adaptor('reverse',ports[0].name,ports[-1].track, ports[-1].gap)
        #     else:
        #         if (inTrack+inGap) > (outTrack+outGap):
        #             iOut, adaptor_length, track_adaptor, gap_adaptor, mask_adaptor = self.draw_adaptor('reverse',ports[0].name,ports[-1].track, ports[-1].gap)
        #         else:
        #             iIn, adaptor_length, track_adaptor, gap_adaptor, mask_adaptor = self.draw_adaptor('reverse',ports[0].name,ports[-1].track, ports[-1].gap)

#            all_constrains = []
#            for constrain in constrains:
#                all_constrains.append(Port.dict_instances[constrain])
#                all_constrains.append(Port.dict_instances[constrain])

#                all_constrains.append(constrain)#[self.ports[constrain][POS], self.ports[constrain][ORI], self.ports[constrain][TRACK], self.ports[constrain][GAP]])
#                previous modification to tackle the case where a cable is drawn between two ports defined on the fly
        nb_constraints = len(ports)-1
        if 0:
            if not isinstance(to_meanders[0], list):
                to_meanders = [to_meanders for ii in range(nb_constraints+1)]
#        print(to_meanders)


        cable_length = []
        tracks = []
        gaps = []
        masks = []
#            port_names = [iInPort]+all_constrains+[iOutPort] # List of ports
#            print(meander_length)

        all_points = []
        index_prev = 0
        index_next = 0
        for ii in range(nb_constraints):

            if not ports[ii].constraint_port:
                index_next = ii+1
                while ports[index_next].constraint_port:
                    index_next += 1
                #index_next is the index of the next normal port

                # find and plot adaptor geometry
                # or reverse
                if reverse_adaptor:
                    points = ports[index_next].compare(ports[ii])
                    index_modified = index_next
                else:
                    points = ports[ii].compare(ports[index_next])
                    index_modified = ii

                for jj, pts in enumerate(points):
                    self.polyline_2D(pts, name=ports[index_modified].name+'_'+ports[index_modified].subnames[jj]+'_adapt', layer=ports[index_modified].layers[jj])

                for jj in range(ii+1, index_next):
                    ports[jj].widths = ports[index_next].widths
                    ports[jj].offsets = ports[index_next].offsets
                    ports[jj].layers = ports[index_next].layers
                    ports[jj].subnames = ports[index_next].subnames
                    ports[jj].N = ports[index_next].N

            # to_meander = to_meanders[ii]
            to_meander = []
            if isinstance(meander_length, (list, np.ndarray)):
                m_length = meander_length[ii]
            else:
                m_length = meander_length

            points = self.find_path('points', ports[ii].name, ports[ii+1].name, fillet, to_bond, to_meander, m_length, meander_offset)
            all_points.append(points)
            ports[ii+1] = ports[ii+1].r

            # plotting cable in between two normal ports
            if ii+1==index_next:
                points = []
                for jj, p in enumerate(all_points):
                    if jj==0:
                        points += p
                    else:
                        points += p[1:]
                self.path(points, ports[index_next], fillet, name=name, layer=layer_TRACK)
                if is_bond:
                    self.draw_bond(to_bond, *ports[index_next].bond_params(), name=name+'_wb_%d'%ii)
                index_prev = index_next
                all_points = []
                to_bond = []

            ports[ii].revert()
            ports[ii+1].revert()

#         for ii in range(nb_constraints):
#             to_meander = to_meanders[ii]
#             if isinstance(meander_length, (list, np.ndarray)):
#                 m_length = meander_length[ii]
#             else:
#                 m_length = meander_length
#             if nb_constraints!=0:
#                 to_add = '_'+str(ii)
#             else:
#                 to_add = ''
# #                print(port_names[2*ii:2*ii+2])
# #                self.__init__(self.name, *port_names[2*ii:2*ii+2]) CONNECT ELEMENT USELESS
#             points = self.find_path('points', ports[ii].name, ports[ii+1].name, fillet, is_meander, to_meander, m_length, meander_offset)
#             all_points.append(points)
#             ports[ii+1].reverse()

        if 0:
            points = []
            for ii, p in enumerate(all_points):
                if ii==0:
                    points += p
                else:
                    points += p[1:]
            if 1:
                self.path(points, ports[0], fillet, name=name, layer=layer_TRACK)
            else:
                connection_track = self.polyline_2D(points, name=name+'path_track'+to_add, closed=False ,layer = layer_Default)
    #            print('length_adaptor = %.3f'%(val(adaptor_length)*1000))
                cable_length.append(self.length(points, 0, len(points)-1, fillet)+val(adaptor_length))
                self._fillets(fillet-eps, connection_track)

                connection_gap = self.polyline_2D(points, name=name+'path_gap'+to_add, closed=False ,layer = layer_Default)
                self._fillets(fillet-eps, connection_gap)
    #                self.set_current_coor([0,0], [0, -1])
                track_starter = self.polyline_2D(points[0]+[[-ports[2*ii].track/2,0],[ports[2*ii].track/2,0]], closed=False, name = name+'start_gap', layer = layer_Default)
                gap_starter = self.polyline_2D([[-ports[2*ii].gap,0],[ports[2*ii].gap,0]], closed=False, name = name+'start_track', layer = layer_Default)

                if self.is_mask:
                    connection_mask = self.polyline_2D(points, name=name+'_mask'+to_add, closed=False ,layer = layer_MASK)
                    self.set_current_coor([0,0], [0, -1])
                    mask_starter = self.cable_starter(name, ports[2*ii].name)
                    mask = self._sweep_along_path(mask_starter, connection_mask)
    #                    masks.append(connection_mask.sweep_along_path(mask_starter))

                track = self._sweep_along_path(track_starter, connection_track)
                tracks.append(track)
                gap = self._sweep_along_path(gap_starter, connection_gap)
                gaps.append(gap)

            if 0:
                if is_bond:
                    self.draw_bond(to_bond, (ports[0].track+ports[0].gap*2)*1.5)

                if track_adaptor is not None:
                    tracks = [*tracks, track_adaptor]
                    gaps = [*gaps, gap_adaptor]
                    if self.is_mask:
        #                    self.maskObjects.pop()
                        masks = [*masks, mask_adaptor]

                if len(tracks)>1:
                    print([t.name for t in tracks])
                    track = self.unite(tracks, name=name+'union_track')
                    gap = self.unite(gaps, name=name+'union_track')
        #                if layer is None:
        #                    self.trackObjects.append(track)
        #                    self.gapObjects.append(gap)
        #                else:
        #                    self.layers[layer]['trackObjects'].append(track)
        #                    self.layers[layer]['gapObjects'].append(gap)
        #                if self.is_mask:
        #                    mask = self.unite(masks, names[2])
        #                    self.maskObjects.append(mask)
                else:
                    track = tracks[0]
                    gap = gaps[0]
        #                if layer is None:
        #                    self.trackObjects.append(track)
        #                    self.gapObjects.append(gap)
        #                else:
        #                    self.layers[layer]['trackObjects'].append(track)
        #                    self.layers[layer]['gapObjects'].append(gap)
        #                if self.is_mask:
        #                    self.maskObjects.append(*masks)
        #            print(track)
        #            if is_mesh is True:
        #                if not self.is_litho:
        #                    self.mesh_zone(track,2*iInPort.track)
        #
                for length in cable_length:
                    print('{0}_length = {1:.3f} mm'.format(name, length*1000))
                print('sum = %.3f mm'%(1000*np.sum(cable_length)))

    @move_port
    def find_simple_path(self, name, iInPort, iOutPort, fillet):

        iIn_pos = Vector(iInPort.pos)
        iIn_ori = Vector(iInPort.ori)
        iOut_pos = Vector(iOutPort.pos)
        iOut_ori = Vector(iOutPort.ori)
        room_bonding = 0*100e-6 #SMPD MANU BOND SPACE

        dist_y = (iOut_pos-iIn_pos).dot(iIn_ori.orth())
        slanted=False
        if iIn_ori.dot(iOut_ori)==-1 and abs(val(dist_y))<val(2*fillet):
            slanted=True  # need to do a slanted_path

        pointA = val(iIn_pos+iIn_ori*room_bonding)
        pointB = val(iOut_pos+iOut_ori*room_bonding)


        point1 = val(iIn_pos+iIn_ori*(1.1*fillet+room_bonding))
        point2 = val(iOut_pos+iOut_ori*(1.1*fillet+room_bonding))

        points_choices = Vector([])
        if iIn_ori.dot(iOut_ori)==-1:
            middle_point = (point1 + point2)/2

            choice_in = next_point(point1, middle_point, iIn_ori) #bon sens
            choice_out = next_point(point2, middle_point, iOut_ori) #à inverser
            for c_in in choice_in:
                for c_out in choice_out:
                    points_choices.append([pointA, *c_in, *c_out[:-1][::-1], pointB])
        else:
            choice_in = next_point(point1, point2, iIn_ori)
            for c_in in choice_in:
                points_choices.append([pointA, *c_in, pointB])

        final_choice= None
        cost=np.inf
        for ii, choice in enumerate(points_choices):
            new_cost, new_choice, new_length = check(choice)

            if new_cost<cost:
                final_choice = new_choice
                cost = new_cost
                length = new_length

        def working_points(points, min_dist, to_meander):
            min_dist = min_dist*1.1
            working_p_start = []
            working_p_end = []
            left_p_start=[points[0]]
            left_p_end=[points[-1]]
            success=False
            index_start = 0
            for ii, point in enumerate(points[1:]):
                A = left_p_start[-1]
                B = point
                AB = B-A
                vec = way(val(B-A))
                if val(AB).norm() > val(min_dist):
                    working_p_start.append(A+vec*min_dist/2)
                    success = True
                    index_start = ii+1
                    break
                else:
                    left_p_start.append(B)
                    to_meander.pop(0)


            if not success:
                print('Warning: Could not find points to elongate cable %s' %self.name)
                left_p = left_p_start+left_p_end[::-1]
                return [], left_p, 0
            else:
                success=False
                index_end = 0
                for ii, point in enumerate(points[::-1][1:]):
                    A = left_p_end[-1]
                    B = point
                    AB = B-A
                    vec = way(val(B-A))
                    if val(AB).norm() > val(min_dist):
                        working_p_end.append(A+vec*min_dist/2)
                        success = True
                        index_end = ii+1
                        break
                    else:
                        left_p_end.append(B)
                        to_meander.pop(-1)
                if not success:
                    print('Warning: Could not find points to elongate cable %s' %self.name)
                    left_p = left_p_start+left_p_end[::-1]
                    return [], left_p, 0


            working_p = working_p_start+points[index_start:-index_end]+working_p_end
            index_insertion = len(left_p_start)
            left_p = left_p_start+left_p_end[::-1]

            return working_p, left_p, index_insertion

        def right_left(points):
            vecs = []
            A = points[0]
            for B in points[1:]:
                vecs.append(way(val(B-A)))
                A=B
        #    print(points)
        #    print(vecs)
            vecA = vecs[0]
            r_l = [0]
            for vecB in vecs[1:]:
                r_l.append(vecA.cross(vecB))
                vecA=vecB
            r_l.append(0)
            return r_l

        def add_points(points, rl, min_dist, n_meander):
            min_dist = min_dist*1.1
            n_points = len(points)

            A = points[0]
            new_points =[]
            if n_points==2:
                new_points.append(A)
                B=points[-1]
                vec = way(val(B-A))
                AB = (B-A).norm()
                n_add = int(val(AB/min_dist))
                if rl[0]*rl[1]==1 and n_add%2==0:
                    n_add-=1
                if rl[0]*rl[1]==-1 and n_add%2==1:
                    n_add-=1
                if n_meander==-1 or n_meander>=n_add:
                    dist = AB/n_add
                    ignore=False
                elif n_meander<n_add:
                    n_add=n_meander
                    centerAB=(A+B)/2
                    addedA=centerAB-vec*n_add/2*min_dist
                    addedB=centerAB+vec*n_add/2*min_dist
                    dist=min_dist
                    A=addedA
                    new_points.append(addedA)
                    ignore=True
                new_points+=[A+vec*dist/2]
                new_points+=[A+vec*dist*(jj+1+1/2) for jj in range(n_add-1)]
                rl=None
                indices_corners=None
#                else:
#                    indices_corners=[]
#                    rl = right_left(points)
#
#                    for ii, B in enumerate(points[1:]):
#                        new_points.append(A)
#                        vec = way(val(B-A))
#                        AB = (B-A).norm()
#                        if ii==0 or ii==n_points-2:
#                            factor = 0.5
#                        else:
#                            factor = 1
#                        n_add = int(val(AB/min_dist)-factor)
#                        if not(ii==0 or ii==n_points-2):
#                            if rl[ii-1]*rl[ii]==-1 and n_add%2==1:
#                                n_add-=1
#                            if rl[ii-1]*rl[ii]==1 and n_add%2==0:
#                                n_add-=1
#
#                        dist = AB/(n_add+factor)
#                        if n_add>=1:
#                            if ii==0:
#                                new_points+=[A+vec*dist/2]
#                                new_points+=[A+vec*dist*(jj+1+1/2) for jj in range(n_add-1)]
#                            elif ii!=n_points-2:
#                                new_points+=[A+vec*dist*(jj+1) for jj in range(n_add)]
#                            else:
#                                new_points+=[A+vec*dist*(jj+1) for jj in range(n_add)]
#                        indices_corners.append(len(new_points))
#                        A=B
#                    indices_corners= indices_corners[:-1]

            if ignore:
                new_points.append(addedB)
            new_points.append(points[-1])
            return new_points, indices_corners, dist, ignore

        def displace(points, rl, min_dist, displacement=0, offset=0, n_meander=-1):
            if np.abs(val(displacement))<val(min_dist)*1.1:
                displacement = min_dist*1.1
            points, indices_corners, dist, ignore = add_points(points, rl, min_dist, n_meander=n_meander)
            new_points = [points[0]]
            parity = 1
            if indices_corners is not None:
                for ii, B in enumerate(points[1:-1]):
                    A = points[ii]
                    AB= B-A
                    vec = way(val(AB))
                    if ii==0:
                        parity = (-2*((indices_corners[0]-(ii+1))%2)+1)*(-rl[0])
                    else:
                        parity = -parity

                    if ii+1 not in indices_corners:
                        #regular point
                        new_points[ii+1] = points[ii+1]+vec.orth()*parity*min_dist
                    else:
                        new_points[ii+1] = points[ii+1]+(vec.orth()*parity+vec).unit()*min_dist
            else:
                if rl[0]!=0:
                    parity = -rl[0]
                else:
                    parity = (2*(len(points)%2)-1) * (-rl[1]*(rl[1]+1)+1)
                if ignore:
                    n_ignore=2
                    new_points.append(points[1])
                else:
                    n_ignore=1
                for ii, B in enumerate(points[n_ignore:-n_ignore]):
                    A=points[ii]
                    AB=B-A
                    vec=way(val(AB))
                    if parity==1:
                        new_points.append(points[ii+n_ignore]+vec.orth()*parity*(displacement+offset)-vec*dist/2)
                        new_points.append(points[ii+n_ignore]+vec.orth()*parity*(displacement+offset)+vec*dist/2)
                    else:
                        new_points.append(points[ii+n_ignore]+vec.orth()*parity*(displacement-offset)-vec*dist/2)
                        new_points.append(points[ii+n_ignore]+vec.orth()*parity*(displacement-offset)+vec*dist/2)
                    parity = -parity
                if ignore:
                    new_points.append(points[-2])
                new_points.append(points[-1])


            return new_points

        def meander(points, min_dist, to_meander, meander_length, meander_offset): # to_meander is list of segments to be meander
            n_points = len(points)
            n_to_meander = len(to_meander)
            if n_points-1>n_to_meander:
                to_meander = to_meander+[0 for ii in range(n_points-1-n_to_meander)]
            else:
                to_meander = to_meander[:n_points-1]
            working_p, left_p, index_insertion = working_points(points, min_dist, to_meander)

            if len(working_p) != 0:
                rl = right_left(working_p)

                working_ps = []
                for ii, isit in enumerate(to_meander):
                    if isit!=0:
                        new_working_p = displace(working_p[ii:ii+2], rl[ii:ii+2], min_dist, displacement=meander_length, offset=meander_offset, n_meander=isit) # n_meander=-1 -> auto
                    else:
                        new_working_p = working_p[ii:ii+2]
                    working_ps += new_working_p
#                        print(working_ps)

                left_p[index_insertion:index_insertion] = working_ps
            return  left_p#left_p#,

        if any(to_meander):
            min_dist = 2*fillet
            final_choice = meander(final_choice, min_dist, to_meander, meander_length, meander_offset)

        final_choice =  [val(iIn_pos)] + final_choice + [val(iOut_pos)]

# Needed to draw Manu bond
        def add_fillet_points(points, fillet):
            new_points = [points[0]]
            for ii, point in enumerate(points[1:-1]):
                index = ii+1
                p_vec = points[index-1]-point
                n_vec = points[index+1]-point
                new_points.append(point+way(val(p_vec))*fillet)
                new_points.append(point+way(val(n_vec))*fillet)
            new_points.append(points[-1])
            return new_points



#        new_points = add_fillet_points(final_choice, fillet)
#        for ii, point in enumerate(new_points[::2]):
#            self.draw('bef_test', [new_points[2*ii], new_points[2*ii+1]], closed=False)
#            self.to_bond.append([new_points[2*ii], new_points[2*ii+1]])



#        self.draw('test', new_points, closed=False)

# Needed for equidistant fillet
#
#
#
#

        def dist(points, A, B, fillet): # A and B are integer point indices
            if A<0 or A>=len(points):
                raise ValueError('First index should be within the point list')
            if B<0 or B>=len(points):
                raise ValueError('Second index should be within the point list')
            if A==B:
                return 0
            if abs(A-B)==1:
                if A<B:
                    if A%2==1:
                        return fillet*np.pi/2
                    else:
                        return (points[A]-points[B]).norm()
                else:
                    return dist(points, B, A, fillet)
            if abs(A-B)>1:
                if A<B:
                    return dist(points, A, B-1, fillet) + dist(points, B-1, B, fillet)
                else:
                    return dist(points, B, A, fillet)

        def where(points, length, fillet):
            n_points = len(points)
            for ii in range(n_points-1):
                distance = dist(points, ii, ii+1, fillet)
                if length <= distance:
                    if ii%2==0:
                        kind = 'normal'
                    else:
                        kind = 'fillet'
                    return [ii, ii+1], kind, length
                else:
                    length = length-distance
            raise ValueError('Length should be smaller than cable length')

        _, final_choice, _ = check(final_choice)

        if slanted:
            h = abs(dist_y/2)
            x = (fillet-h)/fillet
            d = h*x/(1-x**2)**0.5
            dist_x = (iOut_pos-iIn_pos).dot(iIn_ori)
            pointA = iIn_pos+iIn_ori*(dist_x/2-d)
            pointB = iOut_pos+iOut_ori*(dist_x/2-d)

            to_bond.append([iIn_pos, pointA])
            to_bond.append([pointB, iOut_pos])

            return [iIn_pos, pointA, pointB, iOut_pos]
        else:
            to_bond_points = add_fillet_points(final_choice, fillet)
            for ii, point in enumerate(to_bond_points[::2]):
                points = [to_bond_points[2*ii], to_bond_points[2*ii+1]]
    #            print("points", points)
                # self.polyline_2D(points , closed=False, name='bef_test'+str(ii), layer=layer_Default)
                to_bond.append(points)
            return final_choice

    def find_meander_path(points, fillet, to_bond, to_meander, meander_length, meander_offset))

    def length(self, points, A, B, fillet): # A and B are integer point indices
#        for point in points:
#            print(val(point[0]), val(point[1]))
        if A<0 or A>=len(points):
            raise ValueError('First index should be within the point list')
        if B<0 or B>=len(points):
            raise ValueError('Second index should be within the point list')
        if A==B:
            return 0
        if A<B:
            value = 0
            for ii in range(B-A):
                value+=val((points[A+ii+1]-points[A+ii]).norm())
            return value-(B-A-1)*val(fillet*(2-np.pi/2))
        else:
            return self.length(points, B, A, fillet)


    def draw_bond(self, to_bond, ymax, ymin, name='wb', min_dist='0.5mm'):
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