# -*- coding: utf-8 -*-
"""
Created on Mon Nov  4 11:13:09 2019

@author: antho
"""

import gdspy
from hfss import parse_entry, var
# Create the geometry: a single rectangle.


class GdsModeler():
    dict_units = {'km':1.0e3,'m':1.0,'cm':1.0e-2,'mm':1.0e-3}

    def __init__(self, body, unit=1.0e-6, precision=1.0e-9):
        self.unit = unit
        self.precision = precision
        self.body = body
        self.cell = gdspy.Cell(body)
        
    def change_body(self,body):
        self.body = body
        self.cell = gdspy.Cell(body)

    def generate_gds(self,name_file):
        gdspy.write_gds(name_file, unit= 1.0e-6, precision=10.e-9)

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


    def draw_polyline(self, points, layer, closed=True, **kwargs):
        points = parse_entry(points)
        if (closed==True):
            poly1 = gdspy.Polygon(points, layer)
        else:
            poly1 = gdspy.PolyPath(points, layer)
        self.cell.add(poly1)
        self.cell.add(gdspy.Label(kwargs['name'],(0, 0), 'nw', layer=10))

        return kwargs['name']

    def draw_rect_corner(self, pos, size=[0,0,0], **kwargs):
        pos = parse_entry(pos)
        size = parse_entry(size)
        points = [(pos[0],pos[1]), (pos[0]+size[0],pos[1]+0), (pos[0]+size[0],pos[1]+size[1]), (pos[0],pos[1]+size[1])]
        print(points)
        assert ('0' in size or 0 in size)
        poly1 = gdspy.Polygon(points, pos[2])
        self.cell.add(poly1)
        self.cell.add(gdspy.Label(kwargs['name'],(0, 0), 'nw', layer=10))

        return kwargs['name']

    def draw_rect_center(self, pos, size=[0,0,0], **kwargs):
        pos = parse_entry(pos)
        size = parse_entry(size)
        corner_pos = [var(p) - var(s)/2 for p, s in zip(pos, size)]
        return self.draw_rect_corner(corner_pos, size, **kwargs)

    def draw_cylinder(self, pos, radius, height, axis, **kwargs):
        pass
    
    def draw_cylinder_center(self, pos, radius, height, axis, **kwargs):
        pass
    
    def draw_disk(self, pos, radius, axis, layer, **kwargs):
        assert axis=='Z'
        self.cell.add(
              gdspy.Round(
                  (pos[0],pos[1]),
                  radius,
                  layer))
        self.cell.add(gdspy.Label(kwargs['name'],(0, 0), 'nw', layer=10))
        return kwargs['name']
           
    def draw_wirebond(self, pos, ori, width, height='0.1mm', **kwargs): #ori should be normed
        pass
    
    def connect_faces(self, entity1, entity2):
        pass
    
    def delete(self, entity):
        # a verifier
        pass
            
    def rename_entity(self, entity, name):
        pass

    def unite(self, entities, name=None, keep_originals=False):
        polygons = self.cell.get_polygons
        names_polygons = self.cell.get_labels
        
        if len(entities)>=2:
            polygon_0 = polygons[names_polygons.index(entities[0].name)]
            for i in range(len(entities)):
                assert entities[i].name in names_polygons
                if i>=1:
                    polygon_i = polygons[names_polygons.index(entities[i].name)]
                    self.cell.add(gdspy.fast_boolean(polygon_0, polygon_i , 'or'), entities[0].layer)

        if name is None:
            return entities[0].name
        else:
            return name

    def intersect(self, entities, name =None, keep_originals=False):
        polygons = self.cell.get_polygons
        names_polygons= self.cell.get_labels
        
        if len(entities)>=2:
            polygon_0 = polygons[names_polygons.index(entities[0].name)]
            for i in range(len(entities)):
                assert entities[i].name in names_polygons
                if i>=1:
                    polygon_i = polygons[names_polygons.index(entities[i].name)]
                    self.cell.add(gdspy.fast_boolean(polygon_0, polygon_i , 'or'), entities[0].layer)

        if name is None:
            return entities[0].name
        else:
            return name

    def subtract(self, blank_entity, tool_entities, keep_originals=False):
        polygons = self.cell.get_polygons
        names_polygons= self.cell.get_labels
        polygon_0 = polygons[names_polygons.index(blank_entity.name)]
        for i in range(len(tool_entities)):
            assert tool_entities[i].name in names_polygons
            polygon_i = polygons[names_polygons.index(tool_entities[i].name)]
            self.cell.add(gdspy.fast_boolean(polygon_0, polygon_i , 'or'), blank_entity.layer)


    def translate(self, name, vector):
        self._modeler.Move(
            self._selections_array(name),
            ["NAME:TranslateParameters",
        		"TranslateVectorX:="	, vector[0],
        		"TranslateVectorY:="	, vector[1],
        		"TranslateVectorZ:="	, vector[2]]
        )


    def separate_bodies(self, name):
        self._modeler.SeparateBody(["NAME:Selections",
                                		"Selections:=", name,
                                		"NewPartsModelFlag:="	, "Model"
                                	], 
                                	[
                                		"CreateGroupsForNewObjects:=", False
                                	])

#    def make_perfect_E(self, *objects):
#        print(self._boundaries.GetBoundaries())
#        name = increment_name("PerfE", self._boundaries.GetBoundaries())
#        print(name)
#        self._boundaries.AssignPerfectE(["NAME:"+name, "Objects:=", objects, "InfGroundPlane:=", False])

    def assign_perfect_E(self, obj, name='PerfE'):
#        print(obj)

        if isinstance(obj, list):
            self._boundaries.AssignPerfectE(["NAME:"+name, "Objects:=", obj, "InfGroundPlane:=", False])
        else:
            name = str(obj)+'_'+name
            self._boundaries.AssignPerfectE(["NAME:"+name, "Objects:=", [obj], "InfGroundPlane:=", False])
            
    def assign_perfect_E_faces(self, name):
        # this is very peculiar to cavity Si chips
        faces = list(self._modeler.GetFaceIDs(name))
        faces = [int(ii) for ii in faces]
        faces.sort()
        faces_perfE = [faces[1]]+faces[6:]
        faces_loss = faces[2:6]

        self._boundaries.AssignPerfectE(
        	[
        		"NAME:"+str(name)+'_PerfE',
        		"Faces:="		, faces_perfE,
        		"InfGroundPlane:="	, False
        	])
            
        self._boundaries.AssignFiniteCond([ "NAME:FiniteCond3",
                                        		"Faces:="		, faces_loss,
                                        		"UseMaterial:="		, True,
                                        		"Material:="		, "lossy conductor",
                                        		"UseThickness:="	, False,
                                        		"Roughness:="		, "0um",
                                        		"InfGroundPlane:="	, False,
                                        		"IsTwoSided:="		, False,
                                        		"IsInternal:="		, True ])

    def assign_mesh_length(self, obj, length):#, suff = '_mesh'):
        name = str(obj)
#        print(obj)
        params = ["NAME:"+name]
        params += ["RefineInside:=", False, "Enabled:=", True]
        ######## RefineInside Should be False for planar object
        params += ["Objects:=", [str(obj)]]
        params += ["RestrictElem:=", False,
			         "RestrictLength:=",  True,
			         "MaxLength:=", length]
        self._mesh.AssignLengthOp(params)
        
    def create_object_from_face(self, name):
        faces = list(self._modeler.GetFaceIDs(name))
        faces.sort()
        face = faces[0]
        self._modeler.CreateObjectFromFaces(["NAME:Selections",
                                            "Selections:="		, name,
                                            "NewPartsModelFlag:="	, "Model"],
            ["NAME:Parameters",["NAME:BodyFromFaceToParameters","FacesToDetach:="	, [int(face)]]], 
            ["CreateGroupsForNewObjects:=", False])
        return name+'_ObjectFromFace1'

    def _fillet(self, radius, vertex_index, obj):
        vertices = self._modeler.GetVertexIDsFromObject(obj)
        if isinstance(vertex_index, list):
            to_fillet = [int(vertices[v]) for v in vertex_index]
        else:
            to_fillet = [int(vertices[vertex_index])]
#        print(vertices)
#        print(radius)
        self._modeler.Fillet(["NAME:Selections", "Selections:=", obj],
                              ["NAME:Parameters",
                               ["NAME:FilletParameters",
                                "Edges:=", [],
                                "Vertices:=", to_fillet,
                                "Radius:=", radius,
                                "Setback:=", "0mm"]])
            
     
    def _fillet_edges(self, radius, edge_index, obj):
        edges = self._modeler.GetEdgeIDsFromObject(obj)
        print(edges)
        if isinstance(edge_index, list):
            to_fillet = [int(edges[e]) for e in edge_index]
        else:
            to_fillet = [int(edges[edge_index])]
#        print(vertices)
#        print(radius)
        self._modeler.Fillet(["NAME:Selections", "Selections:=", obj],
                              ["NAME:Parameters",
                               ["NAME:FilletParameters",
                                "Edges:=", to_fillet,
                                "Vertices:=", [],
                                "Radius:=", radius,
                                "Setback:=", "0mm"]])   
            

    def _fillets(self, radius, vertices, obj):
        self._modeler.Fillet(["NAME:Selections", "Selections:=", obj],
                              ["NAME:Parameters",
                               ["NAME:FilletParameters",
                                "Edges:=", [],
                                "Vertices:=", vertices,
                                "Radius:=", radius,
                                "Setback:=", "0mm"]])

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
    
    













