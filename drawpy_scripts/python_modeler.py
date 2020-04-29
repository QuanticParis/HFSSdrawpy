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

from .model_entity import ModelEntity

from .utils import VariableString, \
                   extract_value_unit, \
                   extract_value_dim, \
                   parse_entry, \
                   _val, val, equal_float, \
                   way, \
                   Vector

from .path_finder import Path

from .parameters import layer_Default

# ModelEntity should be defined here probably

# PARAMETERS FOR THE GDS OUTPUT AND FOR FILLETS
# NOTE: They are now defined in the parameter file.

##IMPORT KEY / CUSTOM Elements

ureg = UnitRegistry()
Q = ureg.Quantity

def entity_kwargs(kwargs, keys):
    entity_kwargs = {}
    for key in keys:
        if key in kwargs.keys():
            entity_kwargs[key] = kwargs[key]
    return entity_kwargs

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
            from .hfss_modeler import get_active_project
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

    def generate_gds(self, folder, filename):
        file = os.path.join(folder, filename)
        if self.mode=='gds':
            self.interface.generate_gds(file)

    def make_material(self, material_params, name):
        raise NotImplementedError()

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
