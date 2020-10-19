# -*- coding: utf-8 -*-
"""
Created on Mon Nov  4 11:13:09 2019

@author: antho
"""

import numpy as np
import gdspy

from ..utils import parse_entry, val, Vector
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

        #TODO, this is a dirty fixe cause of Vector3D

        points_2D = []
        for point in points:
            points_2D.append([point[0], point[1]])

        if closed:
            poly1 = gdspy.Polygon(points_2D, layer=layer)
        else:
            poly1 = gdspy.FlexPath(points_2D, 1e-9, layer=layer)

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
        corner_pos = [val(p) - val(s)/2 for p, s in zip(pos, size)]
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

    def path(self, points, port, fillet, name='', corner="circular bend"):

        #TODO, this is a dirty fixe cause of Vector3D

        points_2D = []
        for point in points:
            points_2D.append([point[0], point[1]])

        # use dummy layers to recover the right elements
        layers = [ii  for ii in range(len(port.widths))]
        cable = gdspy.FlexPath(points_2D, port.widths, offset=port.offsets,
                               corners=corner,
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
        self.cell = self.gds_cells[blank_entity.body.name]
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

    def subtract(self, blank_entities, tool_entities, keep_originals=True):
        if isinstance(blank_entities, list):
            for blank_entity in blank_entities:
                self.subtract(blank_entity, tool_entities,
                              keep_originals=keep_originals)
        else:
            blank_entity = blank_entities
            #1 We clear the cell of all elements and create lists to store the polygons
            blank_polygon = self.gds_object_instances.pop(blank_entity.name)
            self.cell = self.gds_cells[blank_entity.body.name] # assumes blank and tool are in same body
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
            if subtracted is not None:
                #3 At last we update the cell and the gds_object_instance
                self.gds_object_instances[blank_entity.name] = subtracted
                self.cell.add(subtracted)
            else:
                print('Warning: the entity %s was fully \
                      subtracted'%blank_entity.name)
                dummy = gdspy.Polygon([[0, 0]])
                self.gds_object_instances[blank_entity.name] = dummy
                self.cell.add(dummy)
                blank_entity.delete()

    def assign_material(self, material):
        pass

    def assign_perfect_E(self, entity, name=None):
        pass
    
    def assign_impedance(self, entities, ResistanceSq, ReactanceSq, name="impedance"):
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


    def rect_array(self, pos, size, columns,rows,spacing, origin=(0, 0), **kwargs):
        pos, size = parse_entry(pos, size)
        name = kwargs['name']
        layer = kwargs['layer']
        points = [(pos[0],pos[1]), (pos[0]+size[0],pos[1]+0), (pos[0]+size[0],pos[1]+size[1]), (pos[0],pos[1]+size[1])]
        poly1 = gdspy.Polygon(points, layer)

        self.gds_object_instances[name] = poly1


        cell_to_copy = gdspy.Cell('cell_to_copy')
        self.gds_cells['cell_to_copy'] = cell_to_copy
        cell_to_copy.add(poly1)


        spacing = parse_entry(spacing)

        cell_array=gdspy.CellArray(cell_to_copy, columns, rows, spacing, origin)
        polygon_list=cell_array.get_polygons()
        poly2 = gdspy.PolygonSet(polygon_list, layer)
        self.cell.add(poly2)

        self.gds_object_instances[name] = poly2