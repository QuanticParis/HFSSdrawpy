 # -*- coding: utf-8 -*-
"""
Created on Thu Oct 31 14:14:51 2019

@author: antho
"""

from sympy.parsing import sympy_parser
from pint import UnitRegistry
import numpy as np
import os
from functools import wraps
from inspect import currentframe, getfile

from . import Lib
from . import KeyElement
from . import CustomElement
from .utils import VariableString, \
                   extract_value_unit, \
                   extract_value_dim, \
                   parse_entry, \
                   _val, val, \
                   way, \
                   Vector

from .path_finder import Path

from .parameters import layer_TRACK, \
                        layer_GAP, \
                        layer_RLC, \
                        layer_MESH, \
                        layer_MASK, \
                        layer_Default, \
                        layer_PORT, \
                        eps

from .utils import LENGTH, \
                   INDUCTANCE, \
                   CAPACITANCE, \
                   RESISTANCE, \
                   LENGTH_UNIT, \
                   INDUCTANCE_UNIT, \
                   CAPACITANCE_UNIT, \
                   RESISTANCE_UNIT \

# ModelEntity should be defined here probably

# PARAMETERS FOR THE GDS OUTPUT AND FOR FILLETS
# NOTE: They are now defined in the parameter file.
                        
eps = eps
layer_TRACK = layer_TRACK
layer_GAP = layer_GAP
layer_RLC = layer_RLC
layer_MESH = layer_MESH
layer_MASK = layer_MASK
layer_PORT = layer_PORT
layer_Default = layer_Default

##IMPORT KEY / CUSTOM Elements




ureg = UnitRegistry()
Q = ureg.Quantity

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

def find_last_list(list_entities):
    # return the last list of a set of nested lists
    if isinstance(list_entities, list):
        if len(list_entities)==0:
            return list_entities
        else:
            if isinstance(list_entities[-1], list):
                return find_last_list(list_entities[-1])
            else:
                return list_entities
    else:
        raise TypeError('There are no list')

def add_to_corresponding_list(elt, nested_list, added_elt):
    # return the last list of a set of nested lists
    if isinstance(nested_list, list):
        if elt in nested_list:
            index = nested_list.index(elt)
            nested_list.insert(index+1, added_elt)
            return True
        else:
            for elt_list in nested_list:
                if isinstance(elt_list, list):
                    if add_to_corresponding_list(elt, elt_list, added_elt):
                        break
            else:
                return False
            return True
    else:
        raise TypeError('Argument is not a list')

def general_remove(elt, nested_list):
    # same as list.remove(elt) but for a nested list
    if isinstance(nested_list, list):
        if elt in nested_list:
            nested_list.remove(elt)
            return True
        else:
            for elt_list in nested_list:
                if isinstance(elt_list, list):
                    success = general_remove(elt, elt_list)
                    if success:
                        break
    else:
        raise TypeError('Argument is not a list')

def gen_name(name):
    # routine to mimic the default naming procedure of HFSS when object
    # already exists
    end = ''
    for ii in name[::-1]:
        if ii.isdigit():
            end+=ii
        else:
            break
    if end=='':
        return name+'1'
    number = int(end[::-1])
    if number==0:
        return name+'1'
    else:
        prefix = name[:-len(str(number))]
        suffix = str(number+1)
        return prefix+suffix

def entity_kwargs(kwargs, keys):
    entity_kwargs = {}
    for key in keys:
        if key in kwargs.keys():
            entity_kwargs[key] = kwargs[key]
    return entity_kwargs

def check_name(_class, name):
    i = 0
    new_name = name
    while(new_name in _class.dict_instances.keys()):
        new_name = name+'_'+str(i)
        i+=1
    if new_name != name:
        print("%s: changed '%s' name into '%s'"%(_class.__name__, name, new_name))
    return new_name

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
    mode: string in "gds" or "hfss"
    """
    is_overdev = False
    is_litho = False
    is_mask = False
    gap_mask = parse_entry('20um')
    overdev = parse_entry('0um')

    def __init__(self, mode):
        """
        Creates a PythonMdlr object based on the chosen interface.
        For now the interface cannot be changed during an execution, only at the beginning
        """
        self.mode = mode
        if mode == "hfss":
            from .hfss import get_active_project
            project = get_active_project()
            design = project.get_active_design()
            self.design = design
            self.modeler = design.modeler
            self.modeler.set_units('mm')
            self.modeler.delete_all_objects()
            self.interface = self.modeler
        elif mode=="gds":
            from . import gds_modeler
            self.interface = gds_modeler.GdsModeler()
        else:
            print('Mode should be either hfss or gds')

    ### Utils methods

    @staticmethod
    def reset(self):
        Port.instances_to_move = []
        ModelEntity.instances_to_move = []

    def set_active_coor_system(func):
        """
        Defines a wrapper/decorator which allows the user to always work in the coordinate system of the chosen chip.
        """
        @wraps(func)
        def updated(*args, **kwargs):
            args[0].interface.set_coor_sys(args[0].coor_sys)
            return func(*args, **kwargs)
        return updated

    def delete_all_objects(self, entities):
        for entity in entities:
            entity.delete()

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

        if self.mode == 'hfss':
            self.design.set_variable(name, value)  # for HFSS
        if not name in VariableString.variables.keys():
            return VariableString(name, value=value)
        else:
            VariableString.variables[name]=value
            return VariableString.instances[name]

    def body(self, body_name, rel_coor=None, ref_name='Global'):
        """
        Creates a Body object which inherits from the current PythonModeler object.
        The body is associated with a coordinate system of choice.
        """
        # Note: for now coordinate systems are not reactualized at each run
        if rel_coor is None:
            rel_coor = [[0, 0, 0],  # origin
                        [1, 0, 0],  # new_x
                        [0, 1, 0]]  # new_y
        else:
            rel_coor = parse_entry(rel_coor)
        self.interface.create_coor_sys(coor_sys=body_name, rel_coor=rel_coor,
                                       ref_name=ref_name)
        _body = Body(self, body_name, rel_coor, ref_name)
        return _body

    def generate_gds(self, folder, filename):
        file = os.path.join(folder, filename)
        if self.mode=='gds':
            self.interface.generate_gds(file)

    def make_material(self, material_params, name):

        raise NotImplementedError()
        #raise ImportWarning("make material is not yet implemented")

    ### Drawing methods

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
        name = self.interface.box_corner_3D(pos, size, **kwargs)
        kwargs = entity_kwargs(kwargs, ['layer', 'nonmodel'])
        return ModelEntity(name, 3, self, **kwargs)

    @set_active_coor_system
    def box_center_3D(self, pos, size, **kwargs):
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
        name = self.interface.box_center_3D(pos, size, **kwargs)
        kwargs = entity_kwargs(kwargs, ['layer', 'nonmodel'])
        return ModelEntity(name, 3, self, **kwargs)

    @set_active_coor_system
    def cylinder_3D(self, pos, radius, height, axis, **kwargs):
        name = self.interface.cylinder_3D(pos, radius, height, axis, **kwargs)
        kwargs = entity_kwargs(kwargs, ['layer', 'nonmodel'])
        return ModelEntity(name, 3, self, **kwargs)


    @set_active_coor_system
    def disk_2D(self, pos, radius, axis, **kwargs):
        if self.mode=='gds':
            pos = val(pos)
            radius = val(radius)
        name = self.interface.disk_2D(pos, radius, axis, **kwargs)
        kwargs = entity_kwargs(kwargs, ['layer', 'nonmodel'])
        return ModelEntity(name, 2, self, **kwargs)

    @set_active_coor_system
    def polyline_2D(self, points, closed=True, **kwargs): # among kwargs, name should be given
        i = 0
        while i < len(points[:-1]):
            points_equal = [equal_float(val(p0),val(p1))
                            for p0, p1 in zip(points[i], points[i+1])]
            if all(points_equal):
                points.pop(i)
                print('Warning: Delete two coinciding points on a polyline2D')
            else:
                i+=1
        if self.mode=='gds':
            points = val(points)
        name = self.interface.polyline_2D(points, closed, **kwargs)
        dim = closed + 1
        kwargs = entity_kwargs(kwargs, ['layer', 'nonmodel'])
        return ModelEntity(name, dim, self, **kwargs)

    @set_active_coor_system
    def path_2D(self, points, port, fillet, **kwargs):
        name = kwargs['name']
        model_entities = []
        if self.mode == 'gds':
            points = val(points)
            fillet = val(fillet)
            _port = port.val()
            names, layers = self.interface.path(points, _port, fillet, name=name)
            for name, layer in zip(names, layers):
                kwargs = {'layer':layer}  # model by default for now
                model_entities.append(ModelEntity(name, 2, self, **kwargs))
        elif self.mode == 'hfss':
            # check that port is at the BEGINNING of the path (hfss only)
            ori = port.ori
            pos = port.pos
            path_entity = self.polyline_2D(points, closed=False,
                                           name=name, layer=layer_Default)
            path_entity.fillets(fillet)

            for ii in range(port.N):
                offset = port.offsets[ii]
                width = port.widths[ii]
                subname = port.subnames[ii]
                layer = port.layers[ii]
                points_starter = [Vector(0, offset+width/2).rot(ori)+pos,
                                  Vector(0, offset-width/2).rot(ori)+pos]
                entity = self.polyline_2D(points_starter, closed=False,
                                          name=name+'_'+subname, layer=layer)
                path_name = name+'_'+subname+'_path'
                current_path_entity = path_entity.copy(new_name=path_name)
                self.interface._sweep_along_path(entity, current_path_entity)
                current_path_entity.delete()
                model_entities.append(entity)

            path_entity.delete()

        return model_entities

    @set_active_coor_system
    def rect_corner_2D(self, pos, size, **kwargs):
        kwargs['coor_sys']=self.coor_sys
        if self.mode=='gds':
            pos = val(pos)
            size = val(size)
        name = self.interface.rect_corner_2D(pos, size, **kwargs)
        kwargs = entity_kwargs(kwargs, ['layer', 'nonmodel'])
        return ModelEntity(name, 2, self, **kwargs)

    @set_active_coor_system
    def rect_center_2D(self, pos, size, **kwargs):
        kwargs['coor_sys']=self.coor_sys
        if self.mode=='gds':
            pos = val(pos)
            size = val(size)
        name = self.interface.rect_center_2D(pos, size, **kwargs)
        kwargs = entity_kwargs(kwargs, ['layer', 'nonmodel'])
        return ModelEntity(name, 2, self, **kwargs)





    @set_active_coor_system
    def wirebond_2D(self, pos, ori, ymax, ymin, **kwargs):
        kwargs['coor_sys']=self.coor_sys
        if self.mode=='gds':
            pos, ori, ymax, ymin = val(pos, ori, ymax, ymin)
            name_a, name_b = self.interface.wirebond_2D(pos, ori, ymax, ymin, **kwargs)
            kwargs = entity_kwargs(kwargs, ['layer', 'nonmodel'])
            return ModelEntity(name_a, 2, self, **kwargs), \
                    ModelEntity(name_b, 2, self, **kwargs)
        else:
            name = self.interface.wirebond_2D(pos, ori, ymax, ymin, **kwargs)
            kwargs = entity_kwargs(kwargs, ['layer', 'nonmodel'])
            return ModelEntity(name, 3, self, **kwargs)

    ### Methods acting on list of entities

    def intersect(self, entities, keep_originals = False):
        raise NotImplementedError()
        dim_Intersection = 3;
        for entity in entities:
            if entity.dimension<dim_Intersection:
                dim_Intersection = entity.dimension
        intersection = ModelEntity(entities[0].name, dim_Intersection, entities[0].coor_sys, entities[0].nonmodel)
        if not(keep_originals):
            for entity in entities:
                entity.delete()
        return intersection

    def unite(self, entities, main=None, new_name=None):
        # main: name or entity that should be returned/preserved/final union
        # if new_name (str) is provided, the original entities are kept and
        # the union is named new_name
        if not isinstance(entities, list):
            entities = [entities]
        entities = entities.copy()

        if new_name is None:
            keep_originals = False
        else:
            keep_originals = True

        if main is not None:
            if isinstance(main, str):
                main = ModelEntity.dict_instances[main]
            if main in entities:
                entities.remove(main)
            entities = [main] + entities

        name = entities[0].name

        if len(entities)!=1:
            if not all([entity.dimension == entities[0].dimension
                                                    for entity in entities]):
                raise TypeError('All united elements should have the \
                                same dimension')
            else:
                if keep_originals:
                    entities[0] = entities[0].copy()

                union_entity = self.interface.unite(entities, keep_originals=keep_originals)

                if not keep_originals:
                    for ii in range(len(entities)):
                        entities[ii].delete()
        else:
            union_entity = entities[0]

        if keep_originals:
            union_entity.rename(new_name)

        return union_entity

    def rotate(self, entities, angle=0):
        if isinstance(angle, (list, np.ndarray)):
            if len(angle)==2:
                angle = np.math.atan2(np.linalg.det([[1,0],angle]),np.dot([1,0],angle))
                angle = angle/np.pi*180
            else:
                raise Exception("angle should be either a float or a 2-dim array")
        elif not isinstance(angle, (float, int, VariableString)):
            raise Exception("angle should be either a float or a 2-dim array")
        if self.mode == 'gds':
            angle = val(angle)
        self.interface.rotate(entities, angle)  # angle in degrees

    def translate(self, entities, vector=[0, 0, 0]):
        vector = parse_entry(vector)
        if self.mode == 'gds':
            vector = val(vector)
        self.interface.translate(entities, vector)


class ModelEntity():
    # this should be the objects we are handling on the python interface
    # each method of this class should act in return in HFSS/GDS when possible
    instances_layered = {layer_TRACK:[], layer_GAP:[], layer_MASK:[],
                         layer_Default:[], layer_RLC:[]}
    dict_instances = {}
    instances_to_move = []
    def __init__(self, name, dimension, body, nonmodel=False,
                 layer=layer_Default, coor_sys=None, copy=None):
        name = check_name(self.__class__, name)
        self.name = name
        self.dimension = dimension
        self.body = body
        self.nonmodel = nonmodel
        self.layer = layer
        if coor_sys is None:
            self.coor_sys = body.coor_sys
        else:
            self.coor_sys = coor_sys

        ModelEntity.dict_instances[name] = self
        if layer in self.instances_layered.keys():
            ModelEntity.instances_layered[layer].append(self)
        else:
            ModelEntity.instances_layered[layer]=[self]

        if copy is None:
            find_last_list(ModelEntity.instances_to_move).append(self)
        else:
            # copy is indeed the original object
            # the new object should be put in the same list indent
            add_to_corresponding_list(copy, ModelEntity.instances_to_move,
                                      self)

    def __str__(self):
        return self.name

    def __repr__(self):
        return self.name

    @staticmethod
    def reset():
        ModelEntity.instances_layered = {}
        ModelEntity.dict_instances = {}
        ModelEntity.instances_to_move = []

    @classmethod
    def print_instances(cls):
        for instance_name in cls.dict_instances:
            print(instance_name)

    def delete(self):
        # deletes the modelentity and its occurences throughout the code
        self.body.interface.delete(self)
        self.dict_instances.pop(self.name)
        self.instances_layered[self.layer].remove(self)
        general_remove(self, self.instances_to_move)
        del self

    def copy(self, new_name=None):
        generated_name = gen_name(self.name)
        self.body.interface.copy(self)
        copied = ModelEntity(generated_name, self.dimension, self.body,
                             nonmodel=self.nonmodel, layer=self.layer,
                             coor_sys=self.coor_sys, copy=self)
        if new_name is not None:
            copied.rename(new_name)
        return copied

    def rename(self, new_name):
        self.dict_instances.pop(self.name)
        self.dict_instances[new_name] = self
        self.body.interface.rename(self, new_name)
        self.name = new_name

    def thicken_sheet(self, thickness, bothsides=False):
        raise NotImplementedError()

    def assign_perfect_E(self, suffix='perfE'):
        self.body.interface.assign_perfect_E(self, self.name+'_'+suffix)

    def connect_faces(self, name, entity1, entity2):
        raise NotImplementedError()

    def duplicate_along_line(self, vec):
        # not implemented with the HFSS function for handling the copy better
        # copy and translate the copy
        vec = Vector(vec)
        copy = entity.copy()
        copy.translate(vec)
        return copy

    def fillet(self, radius, vertex_indices):
        # fillet a subset of vertices
        # vertex_indices can be an int or a list of int
        if self.mode=='gds':
            raise NotImplementedError()
        else:
            self.body.interface.fillet(self, radius, vertex_indices)

    def fillets(self, radius):
        # fillet all corner of an entity
        self.body.interface.fillets(self, radius)



    # def make_rlc_boundary(self, corner, size, axis, r, l, c, name="LumpRLC"):
    #     raise NotImplementedError()
#        self.interface.make_rlc_boundary(corner, size, axis, r, l, c, name)

    def assign_material(self, material):
        self.body.interface.assign_material(self, material)

    def assign_mesh_length(self, mesh_length):
        mesh_length = parse_entry(mesh_length)
        self.body.interface.assign_mesh_length(self, mesh_length)

    def assign_lumped_RLC(self, points, rlc):
        points = parse_entry(points)
        # move the points coordinate in the global coordinate system
        if self.body.ref_name != 'Global':
            # TODO do recursive to handle this
            raise NotImplementedError('Do not handle 2nd order relative \
                                      coordinate system yet.')
        origin = self.body.rel_coor[0]
        new_x = self.body.rel_coor[1]
        new_y =self.body.rel_coor[2]
        point_0 = []
        point_1 = []
        for ii in range(3):
            point_0.append(origin[ii] + new_x[ii] * points[0][0] + new_y[ii] * points[0][1])
            point_1.append(origin[ii] + new_x[ii] * points[1][0] + new_y[ii] * points[1][1])

        r, l, c = rlc
        self.body.interface.assign_lumped_rlc(self, r, l, c, point_0,
                                              point_1, name="RLC")

    def mirrorZ(self):
        raise NotImplementedError()

    def rotate(self, angle):
        self.body.rotate(self, angle)

    def translate(self, vector):
        self.body.translate(self, vector)

    def subtract(self, tool_entities, keep_originals=False):
        """
        tool_entities: a list of ModelEntity or a ModelEntity
        keep_originals: Boolean, True : the tool entities still exist after
                        boolean operation
        """
        if not isinstance(tool_entities, list):
            tool_entities = [tool_entities]
        if len(tool_entities)==0:
            pass
        else:
            if not all([entity.dimension==self.dimension
                                                for entity in tool_entities]):
                raise TypeError('All subtracted elements should have the \
                                same dimension')
            else:
                self.body.interface.subtract(self, tool_entities,
                                           keep_originals=True)
            if not keep_originals:
                for ii in range(len(tool_entities)):
                    tool_entities[0].delete()

    def unite(self, tool_entities, new_name=None):
        """
        tool_entities: a list of ModelEntity or a ModelEntity
        if new_name (str) is provided, the tool_entities + self are kept and
        the union is named new_name
        """
        return self.body.unite(tool_entities, main=self, new_name=new_name)



class Port():
    instances_to_move = []
    dict_instances  = {}

    def __init__(self, name, pos, ori, widths, subnames, layers, offsets, constraint_port, key='name'):
        if not (isinstance(key, Port) or key is None):
            name = check_name(self.__class__, name)
        self.name = name
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

        find_last_list(Port.instances_to_move).append(self)
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

    def __str__(self):
        return self.name

    def __repr__(self):
        return self.name

    @staticmethod
    def reset():
        Port.instances_to_move = []
        Port.dict_instances  = {}



    @classmethod
    def print_instances(cls):
        for instance_name in cls.dict_instances:
            print(instance_name)#, cls.dict_instances[instance_name])

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

# class Network(PythonModeler):
#     variables= None
#     to_bond=[]

#     def __init__(self, name, coor_sys, interface, variables):
#         self.interface = interface
#         self.coor_sys = coor_sys
#         self.name = name
#         self.variables = variables
# #        self.__class__.variables = variables


#     def update(self, coor_sys):
#         self.coor_sys = coor_sys

#     def port(self, name, pos, ori, widths, subnames, layers, offsets, constraint_port):
#         return Port(name, pos, ori, widths, subnames, layers, offsets, constraint_port)

#     def draw_capa(self, name, iInPort, iOutPort, iLength, iWidth, iSize):
#         '''
#         Inputs:
#         -------
#         name: string name of object
#         iIn: (position, direction, track, gap) defines the input port
#         iOut: (position, direction, track, gap) defines the output port
#                position and direction are None: this is calculated from
#                other parameters
#         iLength: (float) length of pads
#         iWidth: (float) width of pads
#         iSize: (float) spacing between pads (see drawing)

#         Outputs:
#         --------
#         retIn: same as iIn, with flipped vector
#         retOut: calculated output port to match all input dimensions

#                  iSize
#               +--+  +--+
#               |  |  |  |
#             +-+  |  |  +-+
#         iIn |    |  |    | iOut
#             +-+  |  |  +-+
#               |  |  |  |
#               +--+  +--+
#         '''



#         iLength, iWidth,iSize = parse_entry(iLength, iWidth, iSize)
#         retIn = [[0,0,0], 0, iInPort.track, iInPort.gap]
#         retOut = [[iInPort.gap+iOutPort.gap+iSize+2*iWidth,0],0 , iOutPort.track, iOutPort.gap]

#         points1 = self.append_points([(iInPort.gap+iWidth, 0),
#                                      (0, -iLength/2),
#                                      (-iWidth, 0),
#                                      (0, iLength/2-iInPort.track/2),
#                                      (-iInPort.gap, 0),
#                                      (0, iInPort.track),
#                                      (iInPort.gap, 0),
#                                      (0, iLength/2-iInPort.track/2),
#                                    (iWidth, 0)])


#         trackIn = self.polyline_2D(points1, name=name+'_track1', layer=layer_TRACK)
# #        self.chip.trackObjects.append(trackIn)

#         points2 = self.append_points([(iInPort.gap+iWidth+iSize, 0),
#                                      (0, -iLength/2),
#                                      (+iWidth, 0),
#                                      (0, iLength/2-iOutPort.track/2),
#                                      (+iOutPort.gap, 0),
#                                      (0, iOutPort.track),
#                                      (-iOutPort.gap, 0),
#                                      (0, iLength/2-iOutPort.track/2),
#                                      (-iWidth, 0)])
#         trackOut = self.polyline_2D(points2, name=name+'_track2', layer=layer_TRACK)
# #        self.chip.trackObjects.append(trackOut)

#         points3 = self.append_points([(0, 0),
#                                      (0, iLength/2+iInPort.gap),
#                                      (iInPort.gap+iWidth+iSize/2, 0),
#                                      (0, iOutPort.gap-iInPort.gap),
#                                      (iOutPort.gap+iWidth+iSize/2, 0),
#                                      (0, -iLength-2*iOutPort.gap),
#                                      (-(iOutPort.gap+iWidth+iSize/2),0),
#                                      (0, iOutPort.gap-iInPort.gap),
#                                      (-(iInPort.gap+iWidth+iSize/2),0)
#                                      ])
#         gap1 = self.polyline_2D(points3, name=name+'_gap1', layer=layer_GAP)
# #        self.chip.gapObjects.append(gap1)
# #
# #



# #        TODO !!
# #        if not self.is_litho:
# #            self.draw(self.name+"_mesh", points)
# #            self.modeler.assign_mesh_length(self.name+"_mesh",1/2*iLength)

#         self.iIn = retIn
#         self.iOut = retOut
# #        return [retIn, retOut]


@Lib.add_methods_from(KeyElement, CustomElement)
class Body(PythonModeler):

    def __init__(self, pm, coor_sys, rel_coor, ref_name): #network
        self.pm = pm
        self.coor_sys = coor_sys # also called body_name
                                 # a body is equivalent to a coor_sys
        self.rel_coor = rel_coor
        self.ref_name = ref_name
        self.interface = pm.interface
        self.mode = pm.mode # 'hfss' or 'gds'
        self.current_pos = [0,0]
        self.current_ori = [1,0]

    def modelentities_to_move(self):
        inter1 = ModelEntity.instances_to_move
        inter2= find_last_list(inter1)
        return inter2

    def ports_to_move(self):
        inter1 = Port.instances_to_move
        inter2= find_last_list(inter1)
        return inter2

    def append_lists(self):
        """
        We use a tree-like architecture to store the entities and port to be moved at the right place.
        """
        self.modelentities_to_move().append([])
        self.ports_to_move().append([])

    def set_coor(self, pos, ori):
        self.current_pos, self.current_ori = parse_entry(pos, ori)

    def move_port(func):
        @wraps(func)
        def moved(*args, **kwargs):
            new_args = [args[0], args[1]]  # args[0] = chip, args[1] = name

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
            return func(*new_args, **kwargs)
                # else:
                #     error = '%s arg should be a port'%str(argument)
                #     raise Exception(error)
#            print("compteur",compteur)
#             if compteur==0:
#                 raise Exception("Please indicate more than 0 port")

#             previous_pos = args[0].current_pos
#             previous_ori = args[0].current_ori
#             #TODO
#             if func.__name__=='draw_cable':
#                 #TODO It depends of a parameter of drawCable
#                 args[0].set_coor([0,0],[1,0])
#             elif func.__name__=='find_path':
#                 args[0].set_coor([0,0],[1,0])

#             #  the following is not robust
#             elif compteur==1:
# #                print(new_args[2].pos)
#                 args[0].set_coor(new_args[2].pos, new_args[2].ori)
#             elif compteur==2:
#                 args[0].set_coor(1/2*(new_args[2].pos+new_args[3].pos), new_args[2].ori)
#             new_args = tuple(new_args)
#             return KeyElement._moved(func, previous_pos, previous_ori, *new_args, **kwargs)
        return moved

    def port(self, name, pos, ori, widths, subnames, layers, offsets, constraint_port):
        name = check_name(Port, name)
        if constraint_port:
            pos, ori = parse_entry(pos, ori)
            offset=0
            width=50e-6  # 50um
            points = [(0, offset+width/2),
                      (width/3, offset),
                      (0, offset-width/2)]
            self.polyline_2D(points, name='_'+name, layer=layer_PORT, nonmodel=True)
        else:
            pos, ori, widths, offsets = parse_entry(pos, ori, widths, offsets)
            for ii in range(len(widths)):
                width = widths[ii]
                offset = offsets[ii]
                points = [(0, offset+width/2),
                          (width/3, offset),
                          (0, offset-width/2)]
                self.polyline_2D(points, name='_'+name+'_'+subnames[ii], layer=layer_PORT, nonmodel=True)

        return Port(name, pos, ori, widths, subnames, layers, offsets, constraint_port)

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
            points, length_adaptor = ports[-1].compare(ports[0], self.pm)
            index_modified = -1
        else:
            points, length_adaptor = ports[0].compare(ports[-1], self.pm)
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
        self.path_2D(total_path.points, total_path.port_in, total_path.fillet,
                  name=name)

        # if bond plot bonds
        if is_bond:
            self.draw_bond(total_path.to_bond(), *ports[0].bond_params(), name=name+'_wb')

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
        bond_number = 0
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
            for ii in range(n_bond):
                self.wirebond_2D(pos, ori, ymax, ymin, layer=layer_Default, name=name+'_%d'%(bond_number))
                bond_number += 1
                pos = pos + ori*spacing
            jj+=1

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