# -*- coding: utf-8 -*-
"""
Created on Mon Nov  4 11:13:09 2019

@author: antho
"""

import gdspy
from hfss import parse_entry, var
import numpy as np

# Create the geometry: a single rectangle.
#class Cell( object ):
#    theWholeList= []
#    def __call__( self, *args, **kw ):
#         x= gdspy.Cell( *args, **kw )
#         self.theWholeList.append( x )
#         return x


class GdsModeler():
    gds_object_instances = {}
    dict_units = {'km':1.0e3,'m':1.0,'cm':1.0e-2,'mm':1.0e-3}
    
    def __init__(self, unit=1.0e-6, precision=1.0e-9):
        self.unit = unit
        self.precision = precision
    

    def reset_cell(self):
        del self.cell
        
    def create_coor_sys(self, coor_name, coor_sys):
        print('initialisation cell')
        try:
            self.cell = gdspy.Cell(coor_name)
        except Exception:
            self.cell = gdspy.current_library.cell_dict[coor_name]
        self.coor_sys = coor_sys
        
        
    def generate_gds(self,name_file):
        writer = gdspy.GdsWriter(name_file)
        writer.write_cell(self.cell)
        writer.close()

    def set_units(self, units='m'):
        self.unit = self.dict_units[units]

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

    def draw_box_corner(self, pos, size, **kwargs):
        print("ERROR : The function --draw_box_corner-- cannot be used for GDSmodeler")
        pass
    def draw_box_center(self, pos, size, **kwargs):
        print("ERROR : The function --draw_box_center-- cannot be used for GDSmodeler")
        pass     
    
    def draw_polyline(self, points, size=[0,0], closed=True, **kwargs):
        #size is useless, I just needed a list as the second argument
        name = kwargs['name']
        layer = kwargs['layer']
        points = parse_entry(points)
        if (closed==True):
            poly1 = gdspy.Polygon(points, layer)
        else:
            poly1 = gdspy.PolyPath(points, layer)
        self.gds_object_instances[name] = poly1
        self.cell.add(poly1)
        return name
    
    def draw_rect_corner(self, pos, size, **kwargs):
        #The cas of an already existing name is not dealt with.
        pos = parse_entry(pos)
        size = parse_entry(size)
        name = kwargs['name']
        layer = kwargs['layer']
        #This function neglects the z coordinate
        pos = parse_entry(pos)
        size = parse_entry(size)
        points = [(pos[0],pos[1]), (pos[0]+size[0],pos[1]+0), (pos[0]+size[0],pos[1]+size[1]), (pos[0],pos[1]+size[1])]
        poly1 = gdspy.Polygon(points, layer)
        self.gds_object_instances[name] = poly1
        self.cell.add(poly1)
        return name
    
    def draw_rect_center(self, pos, size, **kwargs):
        pos = parse_entry(pos)
        size = parse_entry(size)
        corner_pos = [var(p) - var(s)/2 for p, s in zip(pos, size)]
        return self.draw_rect_corner(corner_pos, size, **kwargs)

    def draw_cylinder(self, pos, radius, height, axis, **kwargs):
        print("ERROR : The function --draw_cylinder-- cannot be used for GDSmodeler")
        pass
    
    def draw_cylinder_center(self, pos, radius, height, axis, **kwargs):
        #Useless ?
        print("ERROR : The function --draw_cylinder_center-- cannot be used for GDSmodeler")
        pass
    
    def draw_disk(self, pos, radius, axis, **kwargs):
        pos = parse_entry(pos)
        radius = parse_entry(radius)
        name = kwargs['name']
        layer = kwargs['layer']
        assert axis=='Z', "axis must be 'Z' for the gdsModeler"
        round1 = gdspy.Round((pos[0],pos[1]), radius, layer)
        
        self.gds_object_instances[name] = round1
        self.cell.add(round1)
        return name
        
    def draw_wirebond(self, pos, ori, width, height='0.1mm', **kwargs): #ori should be normed
        pass
    
    def connect_faces(self, entity1, entity2):
        print("ERROR : The function --connect_faces-- cannot be used for GDSmodeler")
        pass
    
    def delete(self, entity):
        self.cell.remove_polygon_modified(self.gds_object_instances[entity.name])
        self.gds_object_instances.pop(entity.name, None)
            
    def rename_entity(self, entity, name):
        polygon = self.gds_object_instances.pop(entity.name, None)
        self.gds_object_instances[name] = polygon

    def unite(self, entities, name=None, keep_originals=False):
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

    def subtract(self, blank_entity, tool_entities):
        final_name = blank_entity.name
        print("blank_entity", blank_entity.name)
        print("tool_entities", tool_entities[0].name)
        
        #1 We clear the cell of all elements and create lists to store the polygons
        blank_polygon = self.gds_object_instances[blank_entity.name]
        self.cell.polygons.remove(blank_polygon)
        
        tool_polygons = []
        for tool_entity in tool_entities:
            tool_polygon = self.gds_object_instances[tool_entity.name]
            tool_polygons.append(tool_polygon)
            self.cell.polygons.remove(tool_polygon)
            self.gds_object_instances.pop(tool_entity.name, None)
        
        #2 Then, we perform fast boolean
        for tool_polygon in tool_polygons:
            blank_polygon = gdspy.Polygon(gdspy.fast_boolean(blank_polygon, tool_polygon, 'not').polygons[0])
        
        
        #3 At last we update the cell and the gds_object_instance
        print(type(blank_polygon))
        self.gds_object_instances[final_name] = blank_polygon
        
        self.cell.add(blank_polygon)
        
#        for entity in tool_entities:
#            blank_polygon = self.gds_object_instances[final_name]
#            tool_polygon = self.gds_object_instances[entity.name]
#            self.cell.polygons.remove(self.gds_object_instances[final_name])
#            blank_polygon = gdspy.Polygon(gdspy.fast_boolean(blank_polygon, tool_polygon, 'not').polygons[0])
#            self.gds_object_instances[final_name] = blank_polygon
#            self.cell.polygons.remove(tool_polygon)
#            self.cell.add(blank_polygon)
#            print(polygon_0.polygons)
#            print(tool_polygon.polygons)
#            Not isn't working
#            polygon_0 = gdspy.Polygon(gdspy.fast_boolean(polygon_0, tool_polygon, 'not').polygons[0])
#            print("fused_polygon", polygon_0)
            
#            
#        self.cell.add(polygon_0)
#    
#        for entity in tool_entities:
#            polygon = self.gds_object_instances.pop(entity.name, None)
##            print("polygon",polygon.polygons[0])
##            self.cell.remove_polygon_modified(polygon.polygons[0])
#            self.cell.polygons.remove(polygon)
##            self.cell.remove_labels(lambda lbl : lbl.string==entity.name)
#        
#        self.gds_object_instances[final_name]=polygon_0
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

    def _fillet(self, radius, vertex_index, name_entity):     
        #1 We need to extract the associated polygon
        polygon = self.gds_object_instances[name_entity]
        points = polygon.polygons[0]
        print("points = ", points)
        #2 We adapt the format of the list of radius
        vertices_number = len(points)
#        new radius = [(i in vertex_index) for i in range(vertices_number)]
        new_radius=[0]*vertices_number
        for i in vertex_index:
            new_radius[i]=radius

        #3 We apply fillet on the polygon
        polygon.fillet(new_radius)
        

    def _sweep_along_path(self, to_sweep, path_obj):
        self.rename_obj(path_obj, str(path_obj)+'_path')
        new_name = self.rename_obj(to_sweep, path_obj)
        names = [path_obj, str(path_obj)+'_path']
        self._modeler.SweepAlongPath(self._selections_array(*names),
                                     ["NAME:PathSweepParameters",
                                		"DraftAngle:="		, "0deg",
                                		"DraftType:="		, "Round",
                                		"CheckFaceFaceIntersection:=", False,
                                		"TwistAngle:="		, "0deg"])
        return Polyline(new_name, self)
    
    
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


    def copy(self, obj):
        self._modeler.Copy(["NAME:Selections", "Selections:=", obj])
        new_obj = self._modeler.Paste()
        return new_obj[0]
    
    
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
                

    def get_face_ids(self, obj):
        return self._modeler.GetFaceIDs(obj)

    def get_vertex_ids(self, obj):
        return self._modeler.GetVertexIDsFromObject(obj)

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
#        print(angle)
        for entity in entities:
            if entity!=None:
                gds_entity = self.gds_object_instances[entity.name]
#                print(gds_entity.polygons)
                gds_entity.rotate(angle/360*2*np.pi)










