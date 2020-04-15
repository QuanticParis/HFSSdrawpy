# -*- coding: utf-8 -*-
"""
Created on Mon Nov  4 11:13:09 2019

@author: antho
"""

import numpy as np
import gdspy

from .variable_string import parse_entry, var, Vector

eps = 1e-7
TOLERANCE = 1e-7
print("gdspy_version : ",gdspy.__version__)

# Create the geometry: a single rectangle.
#class Cell( object ):
#    theWholeList= []
#    def __call__( self, *args, **kw ):
#         x= gdspy.Cell( *args, **kw )
#         self.theWholeList.append( x )
#         return x


class GdsModeler():
    gds_object_instances = {}
    gds_cells  = {'Global':[[0,0,0],[1,0,0],[0,1,0]]}
    dict_units = {'km':1.0e3,'m':1.0,'cm':1.0e-2,'mm':1.0e-3}
    coor_systems = {'Global':[[0,0,0],[1,0]]}
    coor_system = coor_systems['Global']

    def __init__(self, unit=1.0e-6, precision=1.0e-9):
        self.unit = unit
        self.precision = precision
        gdspy.current_library = gdspy.GdsLibrary()
        self.create_coor_sys()


    def reset_cell(self):
        del self.cell

    def create_coor_sys(self, coor_name='Global', coor_sys=[[0,0,0], [1,0]]):
        try:
            #Test if the cell already exists
            self.cell = gdspy.current_library.cells[coor_name]
            self.coor = GdsModeler.gds_cells[coor_name]
        except Exception:
            self.cell = gdspy.Cell(coor_name)
            self.coor_system = coor_sys
            self.coor_systems[coor_name] = coor_sys
            self.gds_cells[coor_name] = coor_sys

    def copy(self, polygon, name):
        new_polygon = gdspy.copy(polygon, 0, 0)
        self.gds_object_instances[name] = new_polygon

    def generate_gds(self,name_file):
        gdspy.write_gds(name_file, unit=1.0, precision=1e-9)
#        writer =
#        writer.write_cell(self.cell)

    def set_units(self, units='m'):
        self.unit = self.dict_units[units]

    def set_coor_sys(self, coor_name):
        self.coor_system = self.coor_systems[coor_name]

    def _attributes_array(self, name=None, nonmodel=False, color=None,
                          transparency=0.9, material=None, solve_inside = None):
        arr = ["NAME:Attributes", "PartCoordinateSystem:=", "Global"]
        if name is not None:
            arr.extend(["Name:=", name])
        if nonmodel:
            arr.extend(["Flags:=", "NonModel"])

        if color is not None:
            arr.extend(["Color:=", "(%d %d %d)" % color])
        if transparency is not None:
            arr.extend(["Transparency:=", transparency])
        if material is not None:
            arr.extend(["MaterialName:=", material])
        if solve_inside is not None:
            arr.extend(["SolveInside:=", solve_inside])
        return arr

    def _selections_array(self, *names):
        return ["NAME:Selections", "Selections:=", ",".join(names)]

    def box_corner_3D(self, pos, size, **kwargs):
        print("ERROR : The function --box_corner_3D-- cannot be used for GDSmodeler")
        pass
    def box_center_3D(self, pos, size, **kwargs):
        print("ERROR : The function --box_center_3D-- cannot be used for GDSmodeler")
        pass

    def polyline_2D(self, points, closed, **kwargs):
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
        return name

    def rect_corner_2D(self, pos, size, **kwargs):
        pos, size = parse_entry(pos, size)
        name = kwargs['name']
        layer = kwargs['layer']
        #This function neglects the z coordinate
        points = [(pos[0],pos[1]), (pos[0]+size[0],pos[1]+0), (pos[0]+size[0],pos[1]+size[1]), (pos[0],pos[1]+size[1])]
        poly1 = gdspy.Polygon(points, layer)

        self.gds_object_instances[name] = poly1
        self.cell.add(poly1)
        return name

    def rect_center_2D(self, pos, size, **kwargs):
        pos, size = parse_entry(pos, size)
        corner_pos = [var(p) - var(s)/2 for p, s in zip(pos, size)]
        return self.rect_corner_2D(corner_pos, size, **kwargs)

    def cylinder_3D(self, pos, radius, height, axis, **kwargs):
        print("ERROR : The function --cylinder_3D-- cannot be used for GDSmodeler")
        pass

    def cylinder_center_3D(self, pos, radius, height, axis, **kwargs):
        #Useless ?
        print("ERROR : The function --cylinder_center_3D-- cannot be used for GDSmodeler")
        pass

    def disk_2D(self, pos, radius, axis, number_of_points=0, **kwargs):
        pos, radius = parse_entry(pos, radius)
        name = kwargs['name']
        layer = kwargs['layer']
        assert axis=='Z', "axis must be 'Z' for the gdsModeler"
        round1 = gdspy.Round((pos[0],pos[1]), radius, layer=layer, tolerance=TOLERANCE, number_of_points=number_of_points)
        self.gds_object_instances[name] = round1
        self.cell.add(round1)
        return name

    def wirebond_2D(self, pos, ori, ymax, ymin, height='0.1mm', **kwargs): #ori should be normed
        bond_diam = '20um'
        pos, ori, ymax, ymin, heigth, bond_diam = parse_entry((pos, ori, ymax, ymin, height, bond_diam))
        bond1 = pos + ori.orth()*(ymax+2*bond_diam)
        bond2 = pos + ori.orth()*(ymin-2*bond_diam)
        name_a = self.disk_2D(bond1, bond_diam/2, 'Z', layer=kwargs['layer'], name=kwargs['name']+'a', number_of_points=6)
        name_b = self.disk_2D(bond2, bond_diam/2, 'Z', layer=kwargs['layer'], name=kwargs['name']+'b', number_of_points=6)
        return name_a, name_b

    def connect_faces(self, entity1, entity2):
        print("ERROR : The function --connect_faces-- cannot be used for GDSmodeler")
        pass

    def delete(self, entity):
        self.cell.remove_polygon_modified(self.gds_object_instances[entity.name])
        self.gds_object_instances.pop(entity.name, None)

    def rename_entity(self, entity, name):
        polygon = self.gds_object_instances.pop(entity.name, None)
        self.gds_object_instances[name] = polygon

    def polygonset_norm(self, polygonset):
        if len(polygonset.polygons)==1:
            return np.linalg.norm(polygonset.polygons)**2
        else:
            list_norm = []
            for polygon in polygonset.polygons:
                list_norm.append(np.linalg.norm(polygon)**2)
            return np.amin(np.array(list_norm))

    def gds_boolean(self, operand1, operand2, operation, precision=0.001, max_points=199, layer=0, datatype=0):
        ratio = 10**9
        operand1 = operand1.scale(ratio, ratio)
        operand2 = operand2.scale(ratio, ratio)
        result = gdspy.boolean(operand1, operand2, operation, precision, max_points, layer, datatype)
        result = result.scale(ratio**(-1), ratio**(-1))
        return result

    def unite(self, entities, name=None, keep_originals=False):
        blank_entity = entities.pop(0)
        blank_polygon = self.gds_object_instances[blank_entity.name]
        self.cell.polygons.remove(blank_polygon)
        tool_polygons = []
        for tool_entity in entities:
            tool_polygon = self.gds_object_instances[tool_entity.name]

            self.cell.polygons.remove(tool_polygon)
            if isinstance(tool_polygon, gdspy.PolygonSet):
                for polygon in tool_polygon.polygons:
                    tool_polygons.append(polygon)
            else:
                tool_polygons.append(tool_polygon)

        tool_polygon_set = gdspy.PolygonSet(tool_polygons, layer = blank_entity.layer)
        sub = self.gds_boolean(blank_polygon, tool_polygon_set, 'or')
        if name==None:
            self.gds_object_instances[blank_entity.name] = sub
        else:
            self.gds_object_instances[name] = sub

        self.cell.add(sub)
        return name

    def intersect(self, entities, name=None, keep_originals=False):
        if not(isinstance(entities, list)):
            raise Exception('Union takes a list of entities as an argument')
        if len(entities)==0:
            raise Exception('Union takes a non-empty list of entities as an argument')

        entity_0 = entities.pop(0)
        polygon_0 = self.gds_object_instances[entity_0.name]
        final_name = entities[0].name if name==None else name

        if len(entities)>=2:
            #We fuse all gds_entities on the first element of the list
            polygon_list = []
            for entity in entities:
                polygon_list.append(self.gds_object_instances[entity.name])
            polygon_set = gdspy.PolygonSet(polygon_list)
            fused_polygon = gdspy.fast_boolean(polygon_0, polygon_set , 'or')
            self.cell.add(fused_polygon, polygon_0.layer)

        if not(keep_originals):
            for entity in entities:
                polygon = self.gds_object_instances.pop(entity.name, None)
                self.cell.remove_polygons(polygon)
                self.cell.remove_labels(entity.name)

        self.gds_object_instances[final_name]=fused_polygon

        return final_name

    def subtract(self, blank_entity, tool_entities, keep_originals=True, name=None):
        if name!=None:
            final_name = name
        else:
            final_name = blank_entity.name

        #1 We clear the cell of all elements and create lists to store the polygons
        blank_polygon = self.gds_object_instances.pop(blank_entity.name, None)
        self.cell.polygons.remove(blank_polygon)

        tool_polygons = []
        for tool_entity in tool_entities:
            tool_polygon = self.gds_object_instances[tool_entity.name]
            tool_polygons.append(tool_polygon)
            if not(keep_originals):
                self.cell.polygons.remove(tool_polygon)
                self.gds_object_instances.pop(tool_entity.name, None)
#
#        #2 Then, we perform fast boolean
        sub = blank_polygon
        for tool_polygon in tool_polygons:
            try:
                sub = gdspy.Polygon(self.gds_boolean(sub, tool_polygon, 'not').polygons[0])

            except Exception:
                sub = blank_polygon

        #3 At last we update the cell and the gds_object_instance
        self.gds_object_instances[final_name] = sub
        self.cell.add(sub)


        return final_name

#    def make_perfect_E(self, *objects):
#        print(self._boundaries.GetBoundaries())
#        name = increment_name("PerfE", self._boundaries.GetBoundaries())
#        print(name)
#        self._boundaries.AssignPerfectE(["NAME:"+name, "Objects:=", objects, "InfGroundPlane:=", False])

    def assign_perfect_E(self, obj, name='PerfE'):
        '''useless'''
        pass

    def assign_perfect_E_faces(self, name):
        print("ERROR : The function --assign_perfect_E-- cannot be used for GDSmodeler")
        pass

    def assign_mesh_length(self, obj, length):#, suff = '_mesh'):
        print("ERROR : The function --assign_mesh_length-- cannot be used for GDSmodeler")
        pass

    def create_object_from_face(self, name):
        print("ERROR : The function --create_object_from_face-- cannot be used for GDSmodeler")
        pass

    def _fillet(self, radius, vertices, entity):
        #TODO : Correct the choice of vertices
        #1 We need to extract the associated polygon
        polygon = self.gds_object_instances[entity.name]
        points = polygon.polygons[0]
#        print("points = ", points)
#        #2 We adapt the format of the list of radius
#        vertices_number = len(points)
##        new radius = [(i in vertex_index) for i in range(vertices_number)]
        new_radius=[0]*len(points)
        for i in vertices:
            new_radius[i]=radius
#        #3 We apply fillet on the polygon
        polygon.fillet(new_radius)

    def _fillets(self, radius, vertices, entity):
        #Vertices = None in gds
        polygon = self.gds_object_instances[entity.name]
        points = gdspy.PolygonSet(polygon.get_polygons())
        self.cell.add(points)
        points.fillet(radius)

    def get_vertex_ids(self, entity):
        return None


    def path(self, points, port, fillet, **kwargs):

        name = kwargs['name']
        # use dummy layers to recover the right elements
        layers = [ii  for ii in range(len(port.widths))]
        cable = gdspy.FlexPath(points, port.widths, offset=port.offsets,
                               corners="circular bend",
                               bend_radius=fillet, gdsii_path=False,
                               tolerance=TOLERANCE, layer=layers, max_points=0) # tolerance (meter) is highly important here should be smaller than the smallest dim typ. 100nm

        polygons = cable.get_polygons()
        names = []
        # print('cable.max_points')
        # print(cable.max_points)
        # print(port.N)
        # print(len(polygons))
        # print('enumerate dict keys')
        # for key in cable._polygon_dict.keys():
        #     print(key)
        #     print(len(cable._polygon_dict[key]))
        #     for elt in cable._polygon_dict[key]:
        #         print(elt)

        for ii in range(len(polygons)):
            poly = gdspy.Polygon(polygons[ii])
            poly.layers = [port.layers[ii]]
            current_name = name+'_'+port.subnames[ii]
            names.append(current_name)
            self.gds_object_instances[current_name] = poly
            self.cell.add(poly)

        return names

    def _sweep_along_path(self, to_sweep, path_obj, name=None):
        try:
            length_to_sweep = (Vector(self.gds_object_instances[to_sweep.name].polygons[0][0])-Vector(self.gds_object_instances[to_sweep.name].polygons[0][1])).norm()
        except Exception:
            length_to_sweep = (Vector(self.gds_object_instances[to_sweep.name].points[0])-Vector(self.gds_object_instances[to_sweep.name].points[1])).norm()
        try:
            flexpath = gdspy.FlexPath(self.gds_object_instances[path_obj.name].polygons[0], length_to_sweep, layer = to_sweep.layer)
        except Exception:
            flexpath = gdspy.FlexPath(self.gds_object_instances[path_obj.name].points, length_to_sweep, layer = to_sweep.layer)


        self.cell.add(flexpath)

        if name==None:
            return to_sweep.name
        else:
            return name

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

    def mirrorZ(self, obj):
        self._modeler.Mirror([
                        		"NAME:Selections", "Selections:=", obj,
                        		"NewPartsModelFlag:="	, "Model"
                        	  ],
                        	  [
                        		"NAME:MirrorParameters",
                        		"MirrorBaseX:="		, "0mm",
                        		"MirrorBaseY:="		, "0mm",
                        		"MirrorBaseZ:="		, "0mm",
                        		"MirrorNormalX:="	, "0mm",
                        		"MirrorNormalY:="	, "0mm",
                        		"MirrorNormalZ:="	, "1mm"
                        	  ])
        return obj


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
        translation_vector = [vector[0], vector[1]]
        for entity in entities:
            if entity!=None:
                gds_entity = self.gds_object_instances[entity.name]
                gds_entity.translate(*translation_vector)


    def rotate(self, entities, angle):
        for entity in entities:
            if entity!=None:
                gds_entity = self.gds_object_instances[entity.name]
#                print(gds_entity.polygons)
                gds_entity.rotate(angle/360*2*np.pi)











