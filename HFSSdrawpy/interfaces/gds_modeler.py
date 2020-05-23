# -*- coding: utf-8 -*-
"""
Created on Mon Nov  4 11:13:09 2019

@author: antho
"""

import numpy as np
import gdspy

from ..utils import parse_entry, var, Vector
from ..core.entity import gen_name

TOLERANCE = 1e-7 # for arcs
print("gdspy_version : ",gdspy.__version__)

class GdsModeler():
    gds_object_instances = {}
    gds_cells  = {}
    dict_units = {'km':1.0e3,'m':1.0,'cm':1.0e-2,'mm':1.0e-3}
    # coor_systems = {'Global':[[0,0,0],[1,0]]}
    # coor_system = coor_systems['Global']

    def __init__(self, unit=1.0e-6, precision=1.0e-9):
        self.unit = unit
        self.precision = precision
        gdspy.current_library = gdspy.GdsLibrary()

    @classmethod
    def print_instances(cls):
        for instance_name in cls.gds_object_instances:
            print(instance_name)

    def reset_cell(self):
        del self.cell

    def create_coor_sys(self, coor_sys='chip', rel_coor=None,
                        ref_name='Global'):
        # this creates a cell, should not care about the rel_coor
        if not (coor_sys in gdspy.current_library.cells.keys()):
            cell = gdspy.Cell(coor_sys)
            self.gds_cells[coor_sys] = cell
        else:
            cell = self.gds_cells[coor_sys]
        # active cell should be the new cell
        self.cell = cell

    def set_coor_sys(self, coor_sys):
        if coor_sys in self.gds_cells.keys():
            self.cell = self.gds_cells[coor_sys]
        else:
            raise ValueError('%s cell do not exist'%coor_sys)

    def copy(self, entity):
        new_polygon = gdspy.copy(self.gds_object_instances[entity.name], 0, 0)
        new_name = gen_name(entity.name)
        self.gds_object_instances[new_name] = new_polygon
        self.cell.add(new_polygon)

    def rename(self, entity, name):
        obj = self.gds_object_instances.pop(entity.name)
        self.gds_object_instances[name]=obj

    def generate_gds(self, file):
        for cell_name in self.gds_cells.keys():
            filename = file+'_%s.gds'%cell_name
            gdspy.write_gds(filename, cells=[cell_name],
                            unit=1.0, precision=1e-9)

    def get_vertices(self, entity):
        polygon = self.gds_object_instances[entity.name]
        return polygon.polygons[0]

    def set_units(self, units='m'):
        self.unit = self.dict_units[units]

    def box(self, pos, size, **kwargs):
        pass

    def box_center(self, pos, size, **kwargs):
        pass

    def polyline(self, points, closed, **kwargs):
        #TODO sace of open path
        #size is the thickness of the polyline for gds, must be a 2D-list with idential elements
        name = kwargs['name']
        layer = kwargs['layer']
        points = parse_entry(points)
        if closed:
            poly1 = gdspy.Polygon(points, layer=layer)
        else:
            poly1 = gdspy.FlexPath(points, 1e-9, layer=layer)

        self.gds_object_instances[name] = poly1
        self.cell.add(poly1)

    def rect(self, pos, size, **kwargs):
        pos, size = parse_entry(pos, size)
        name = kwargs['name']
        layer = kwargs['layer']
        #This function neglects the z coordinate
        points = [(pos[0],pos[1]), (pos[0]+size[0],pos[1]+0), (pos[0]+size[0],pos[1]+size[1]), (pos[0],pos[1]+size[1])]
        poly1 = gdspy.Polygon(points, layer)

        self.gds_object_instances[name] = poly1
        self.cell.add(poly1)

    def rect_center(self, pos, size, **kwargs):
        pos, size = parse_entry(pos, size)
        corner_pos = [var(p) - var(s)/2 for p, s in zip(pos, size)]
        self.rect(corner_pos, size, **kwargs)

    def cylinder(self, pos, radius, height, axis, **kwargs):
        pass

    def disk(self, pos, radius, axis, number_of_points=None, **kwargs):
        pos, radius = parse_entry(pos, radius)
        name = kwargs['name']
        layer = kwargs['layer']
        assert axis=='Z', "axis must be 'Z' for the gdsModeler"
        round1 = gdspy.Round((pos[0],pos[1]), radius, layer=layer, tolerance=TOLERANCE, number_of_points=number_of_points)
        self.gds_object_instances[name] = round1
        self.cell.add(round1)

    def wirebond(self, pos, ori, ymax, ymin, height='0.1mm', **kwargs): #ori should be normed
        bond_diam = '20um'
        pos, ori, ymax, ymin, heigth, bond_diam = parse_entry((pos, ori, ymax, ymin, height, bond_diam))
        bond1 = pos + ori.orth()*(ymax+2*bond_diam)
        bond2 = pos + ori.orth()*(ymin-2*bond_diam)
        self.disk(bond1, bond_diam/2, 'Z', layer=kwargs['layer'], name=kwargs['name']+'a', number_of_points=6)
        self.disk(bond2, bond_diam/2, 'Z', layer=kwargs['layer'], name=kwargs['name']+'b', number_of_points=6)

    def path(self, points, port, fillet, name=''):

        # use dummy layers to recover the right elements
        layers = [ii  for ii in range(len(port.widths))]
        cable = gdspy.FlexPath(points, port.widths, offset=port.offsets,
                               corners="circular bend",
                               bend_radius=fillet, gdsii_path=False,
                               tolerance=TOLERANCE, layer=layers, max_points=0) # tolerance (meter) is highly important here should be smaller than the smallest dim typ. 100nm

        polygons = cable.get_polygons()
        names = []
        layers = []
        for ii in range(len(polygons)):
            poly = gdspy.Polygon(polygons[ii])
            poly.layers = [port.layers[ii]]
            current_name = name+'_'+port.subnames[ii]
            names.append(current_name)
            layers.append(port.layers[ii])
            self.gds_object_instances[current_name] = poly
            self.cell.add(poly)
        return names, layers

    def connect_faces(self, entity1, entity2):
        pass

    def delete(self, entity):
        self.cell.polygons.remove(self.gds_object_instances[entity.name])
        self.gds_object_instances.pop(entity.name)

    def rename_entity(self, entity, name):
        polygon = self.gds_object_instances.pop(entity.name)
        self.gds_object_instances[name] = polygon

    def unite(self, entities, keep_originals=True):

        blank_entity = entities.pop(0)
        blank_polygon = self.gds_object_instances.pop(blank_entity.name)
        self.cell.polygons.remove(blank_polygon)

        tool_polygons = []
        for tool_entity in entities:
            tool_polygon = self.gds_object_instances[tool_entity.name]
            if isinstance(tool_polygon, gdspy.PolygonSet):
                for polygon in tool_polygon.polygons:
                    tool_polygons.append(polygon)
            else:
                tool_polygons.append(tool_polygon)

        #2 unite operation
        tool_polygon_set = gdspy.PolygonSet(tool_polygons, layer = blank_entity.layer)
        united = gdspy.boolean(blank_polygon, tool_polygon_set, 'or',
                               precision=TOLERANCE, max_points=0,
                               layer=blank_entity.layer)

        self.gds_object_instances[blank_entity.name] = united
        self.cell.add(united)

        return blank_entity

    def intersect(self, entities):
        raise NotImplementedError()

    def subtract(self, blank_entity, tool_entities, keep_originals=True):
        #1 We clear the cell of all elements and create lists to store the polygons
        blank_polygon = self.gds_object_instances.pop(blank_entity.name)
        self.cell.polygons.remove(blank_polygon)

        tool_polygons = []
        for tool_entity in tool_entities:
            tool_polygon = self.gds_object_instances[tool_entity.name]
            if isinstance(tool_polygon, gdspy.PolygonSet):
                for polygon in tool_polygon.polygons:
                    tool_polygons.append(polygon)
            else:
                tool_polygons.append(tool_polygon)


        #2 subtract operation
        tool_polygon_set = gdspy.PolygonSet(tool_polygons, layer = blank_entity.layer)
        subtracted = gdspy.boolean(blank_polygon, tool_polygon_set, 'not',
                                    precision=TOLERANCE, max_points=0,
                                    layer=blank_entity.layer)

        #3 At last we update the cell and the gds_object_instance
        self.gds_object_instances[blank_entity.name] = subtracted
        self.cell.add(subtracted)

    def assign_material(self, material):
        pass

    def assign_perfect_E(self, entity, name=None):
        pass

    def assign_perfect_E_faces(self, entity):
        pass

    def assign_mesh_length(self, entity, length):#, suff = '_mesh'):
        pass

    def assign_lumped_rlc(self, entity, r, l, c, start, end, name="RLC"):
        pass

    def create_object_from_face(self, name):
        pass

    def fillet(self, entity, radius, vertex_indices=None):
        polygon = self.gds_object_instances[entity.name]
        if vertex_indices is None:
            polygon.fillet(radius, max_points=0)
        else:
            vertices_number = len(polygon.polygons[0])
            radii = [0]*vertices_number
            for rad, indices in zip(radius, vertex_indices):
                for index in indices:
                    radii[index]=rad
            polygon.fillet([radii], max_points=0, precision=TOLERANCE)

    def get_vertex_ids(self, entity):
        return None



    def sweep_along_vector(self, names, vector):
        self._modeler.SweepAlongVector(self._selections_array(*names),
                                        	["NAME:VectorSweepParameters",
                                        		"DraftAngle:="		, "0deg",
                                        		"DraftType:="		, "Round",
                                        		"CheckFaceFaceIntersection:=", False,
                                        		"SweepVectorX:="	, vector[0],
                                        		"SweepVectorY:="	, vector[1],
                                        		"SweepVectorZ:="	, vector[2]
                                        	])


    def thicken_sheet(self, sheet, thickness, bothsides=False):
        self._modeler.ThickenSheet([
                                		"NAME:Selections", "Selections:=", sheet,
                                		"NewPartsModelFlag:="	, "Model"
                                	],
                                	[
                                		"NAME:SheetThickenParameters",
                                		"Thickness:="		, thickness,
                                		"BothSides:="		, bothsides
                                	])

    def mirrorZ(self, entity):
        pass


    def duplicate_along_line(self, obj, vec):
        self._modeler.DuplicateAlongLine(["NAME:Selections","Selections:=", obj,
                                          "NewPartsModelFlag:="	, "Model"],
                                        	["NAME:DuplicateToAlongLineParameters",
                                    		"CreateNewObjects:="	, True,
                                    		"XComponent:="		, vec[0],
                                    		"YComponent:="		, vec[1],
                                    		"ZComponent:="		, vec[2],
                                    		"NumClones:="		, "2"],
                                        	["NAME:Options",
                                        	"DuplicateAssignments:=", False],
                                        	["CreateGroupsForNewObjects:=", False	])

    def duplicate_along_line(self, obj, vec, n=2):
        self._modeler.DuplicateAlongLine(["NAME:Selections","Selections:=", obj,
                                          "NewPartsModelFlag:="	, "Model"],
                                        	["NAME:DuplicateToAlongLineParameters",
                                    		"CreateNewObjects:="	, True,
                                    		"XComponent:="		, vec[0],
                                    		"YComponent:="		, vec[1],
                                    		"ZComponent:="		, vec[2],
                                    		"NumClones:="		, str(n)],
                                        	["NAME:Options",
                                        	"DuplicateAssignments:=", False],
                                        	["CreateGroupsForNewObjects:=", False	])

    def _make_lumped_rlc(self, r, l, c, start, end, obj_arr, name="LumpLRC"):
        name = increment_name(name, self._boundaries.GetBoundaries())
        params = ["NAME:"+name]
        params += obj_arr
        params.append(["NAME:CurrentLine", "Start:=", start, "End:=", end])
        params += ["UseResist:=", r != 0, "Resistance:=", r,
                   "UseInduct:=", l != 0, "Inductance:=", l,
                   "UseCap:=", c != 0, "Capacitance:=", c]
        self._boundaries.AssignLumpedRLC(params)

    def _make_lumped_port(self, start, end, obj_arr, z0="50ohm", name="LumpPort"):
        name = increment_name(name, self._boundaries.GetBoundaries())
        params = ["NAME:"+name]
        params += obj_arr
        params += ["RenormalizeAllTerminals:=", True, "DoDeembed:=", False,
                   ["NAME:Modes", ["NAME:Mode1", "ModeNum:=", 1, "UseIntLine:=", True,
                                   ["NAME:IntLine", "Start:=", start, "End:=", end],
                                   "CharImp:=", "Zpi", "AlignmentGroup:=", 0, "RenormImp:=", "50ohm"]],
                   "ShowReporterFilter:=", False, "ReporterFilter:=", [True],
                   "FullResistance:=", "50ohm", "FullReactance:=", "0ohm"]

        self._boundaries.AssignLumpedPort(params)

    def make_material(self, obj, material):
        self._modeler.ChangeProperty(["NAME:AllTabs",
                                		["NAME:Geometry3DAttributeTab",
                                			["NAME:PropServers", obj],
                                			["NAME:ChangedProps",
                                				["NAME:Material","Value:=", material]
                                			]
                                		]
                                	])



    def eval_expr(self, expr, units="mm"):
        if not isinstance(expr, str):
            return expr
        return self.parent.eval_expr(expr, units)

    def eval_var_str(self, name, unit=None):
        if not isinstance(name, VariableString):
            if unit is not None:
                return str(name)+' '+unit
            else:
                return str(name)
        return self.parent.eval_var_str(name, unit=unit)

    def delete_all_objects(self):
        objects = []
        for ii in range(int(self._modeler.GetNumObjects())):
            objects.append(self._modeler.GetObjectName(str(ii)))
#        print(objects)
        self._modeler.Delete(self._selections_array(*objects))

    def get_matched_object_name(self, name):
        return self._modeler.GetMatchedObjectName(name+'*')

    def translate(self, entities, vector):
        '''vector is 3-dimentional but with a z=0 component'''
        if not isinstance(entities, list):
            entities = [entities]
        translation_vector = [vector[0], vector[1]]
        for entity in entities:
            # if entity!=None:
            gds_entity = self.gds_object_instances[entity.name]
            gds_entity.translate(*translation_vector)

    def rotate(self, entities, angle):
        if not isinstance(entities, list):
            entities = [entities]
        for entity in entities:
            # if entity!=None:
            gds_entity = self.gds_object_instances[entity.name]
            gds_entity.rotate(angle/360*2*np.pi)











