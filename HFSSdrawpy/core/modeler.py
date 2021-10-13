import os
from inspect import currentframe, getfile

import numpy as np
import sympy
from pint import UnitRegistry

from ..utils import parse_entry, store_variable, val, variables
from .entity import Entity

sympy.init_printing(use_latex=False)


class Modeler:
    """
    Modeler which defines basic operations and methods to perform on Entity and on the chosen interface.
    To create a new interface, one needs to copy an existing one in a new file, and adapt all methods to the formalism of the new interface.

    The Modeler class is called at the beginning of a script but then it is the body that is always used.
    The syntax is the following:
        # TODO

    Inputs:
    -------
    mode: string in "gds" or "hfss"
    """

    is_overdev = False
    overdev = parse_entry("0um")

    def __init__(self, mode):
        """
        Creates a Modeler object based on the chosen interface.
        For now the interface cannot be changed during an execution, only at the beginning
        """
        self.mode = mode
        if mode == "hfss":
            from ..interfaces.hfss_modeler import get_desktop

            desktop = get_desktop()
            project = desktop.get_active_project()
            design = project.get_active_design()
            self.design = design
            self.modeler = design.modeler
            self.modeler.set_units("mm")
            self.modeler.delete_all_objects()
            desktop.clear_all_messages()
            self.interface = self.modeler
        elif mode == "gds":
            from ..interfaces import gds_modeler

            self.interface = gds_modeler.GdsModeler()
        else:
            print("Mode should be either hfss or gds")
            
        # default init for mask values (used to subtract holes from
        # critical areas)
        self.is_mask = False
        self.gap_mask = parse_entry("20um")
        self.is_litho = False

        # The list of bodies pointing to the current Modeler
        self.bodies = []

    ### Utils methods

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
            f = currentframe().f_back  # .f_back
            filename = getfile(f)
            code_line = open(filename).readlines()[f.f_lineno - 1]
            name = code_line.split("=")[0].strip()

        if self.mode == "hfss":
            self.design.set_variable(name, value)  # for HFSS
        symbol = sympy.symbols(name)
        store_variable(symbol, value)
        return symbol

    def generate_gds(self, folder, filename, max_points=0):
        file = os.path.join(folder, filename)
        if self.mode == "gds":
            self.interface.generate_gds(file, max_points)
            
    # def import_gds(self, folder, filename):
    #     file = os.path.join(folder, filename)
    #     if self.mode == "gds":
    #         self.interface.import_gds(file)

    def make_material(self, material_params, name):
        raise NotImplementedError()

    ### Methods acting on list of entities

    def intersect(self, entities, keep_originals=False):
        raise NotImplementedError()

    def unite(self, entities, main=None, keep_originals=False, new_name=None):
        # main: name or entity that should be returned/preserved/final union
        # if new_name (str) is provided, the original entities are kept and
        # the union is named new_name
        if not isinstance(entities, list):
            entities = [entities]
        entities = entities.copy()

        # if new_name is None:
        #     keep_originals = False
        # else:
        #     keep_originals = True

        if main is not None:
            if isinstance(main, str):
                main = Entity.dict_instances[main]
            if main in entities:
                entities.remove(main)
            entities = [main] + entities

        if len(entities) != 1:
            if not all([entity.dimension == entities[0].dimension for entity in entities]):
                raise TypeError(
                    "All united elements should have the \
                                same dimension"
                )
            else:
                if keep_originals:
                    entities[0] = entities[0].copy()

                union_entity = self.interface.unite(entities, keep_originals=keep_originals)
                union_entity.is_boolean = True
                list_fillet = [entity.is_fillet for entity in entities]
                union_entity.is_fillet = union_entity.is_fillet or any(list_fillet)

                if not keep_originals:
                    ents = entities.copy()
                    for entity in ents:
                        entity.delete()
        else:
            union_entity = entities[0]

        if new_name:
            union_entity.rename(new_name)

        return union_entity

    def subtract(self, blank_entities, tool_entities, keep_originals=False):
        """
        tool_entities: a list of Entity or a Entity
        keep_originals: Boolean, True : the tool entities still exist after
                        boolean operation
        """
        if not isinstance(blank_entities, list):
            blank_entities = [blank_entities]
        if not isinstance(tool_entities, list):
            tool_entities = [tool_entities]
        if len(blank_entities) == 0 or len(tool_entities) == 0:
            pass
        else:
            if not all(
                [entity.dimension == blank_entities[0].dimension for entity in blank_entities]
            ) or not all(
                [entity.dimension == tool_entities[0].dimension for entity in tool_entities]
            ):
                raise TypeError(
                    "All subtracted elements should have the \
                                same dimension"
                )
            else:
                self.interface.subtract(blank_entities, tool_entities, keep_originals=True)
                # actualize the properties of the blank_entities
                list_fillet_bool = any([entity.is_fillet for entity in tool_entities])
                for entity in blank_entities:
                    entity.is_boolean = True
                    entity.is_fillet = entity.is_fillet or list_fillet_bool
                    # this is not optimal fillet wise but hard to do better
            if not keep_originals:
                tools = tool_entities.copy()
                for tool_entity in tools:
                    tool_entity.delete()
    
    def delete_inside(self, polygon_set, mask, keep_originals=False):
        '''
        Test if the polygons within the polygon_set are in the mask object.
        If so, delete them.
        
        Not implemented in HFSS
        Parameters
        ----------
        polygon_set : Entity 
            Typically a hole array.
        mask : Entity
        keep_originals : bool, optional
            Shall we keep the mask element or not. The default is False.
        '''
        # TODO handle the case when the polygon_set is fully deleted
        # TODO handle the HFSS case by doing a substract
        self.interface.delete_inside(polygon_set, mask, keep_originals=keep_originals)

    def rotate(self, entities, angle=0):
        if isinstance(angle, (list, np.ndarray)):
            if len(angle) == 2:
                angle = np.math.atan2(np.linalg.det([[1, 0], angle]), np.dot([1, 0], angle))
                angle = angle / np.pi * 180
            else:
                raise Exception("angle should be either a float or a 2-dim array")
        elif not isinstance(angle, (float, int)):
            raise Exception("angle should be either a float or a 2-dim array")
        if self.mode == "gds":
            angle = val(angle)
        self.interface.rotate(entities, angle)  # angle in degrees

    def translate(self, entities, vector=[0, 0, 0]):
        vector = parse_entry(vector)
        if self.mode == "gds":
            vector = val(vector)
        self.interface.translate(entities, vector)
