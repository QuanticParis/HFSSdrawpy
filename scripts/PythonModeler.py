# -*- coding: utf-8 -*-
"""
Created on Thu Oct 31 14:14:51 2019

@author: antho
"""

from sympy.parsing import sympy_parser
from pint import UnitRegistry
import numpy as np
from Vector import Vector
import sys
from functools import wraps
eps = 1e-7
layer_TRACK = 1
layer_GAP = 0
layer_RLC = 2
layer_MESH = 3
layer_MASK = 4
layer_Default = 10

#IMPORT GDS / HFSS Modelers
from hfss import extract_value_unit, \
                 extract_value_dim, \
                 parse_entry, \
                 get_active_project, \
                 ModelEntity, \
                 VariableString , \
                 var
import gds_modeler
##IMPORT KEY / CUSTOM Elements
from Lib import *
import KeyElement
import CustomElement


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

class PythonMdlr():
    is_overdev = False
    is_litho = False
    is_mask = False
    pcb_track = parse_entry('300um')
    pcb_gap = parse_entry('200um')
    gap_mask = parse_entry('20um')
    overdev = parse_entry('0um')

    def __init__(self, name_interface): #"Hfss" or "Gds"
        self.variables = {}     
        if name_interface=="hfss":
            project = get_active_project()
            design = project.get_active_design()
            self.design = design
            self.modeler = design.modeler
            self.modeler.set_units('mm')
            self.modeler.delete_all_objects()
            self.interface = self.modeler
        if name_interface=="gds":
            self.interface = gds_modeler.GdsModeler()
        self.mode = name_interface
        
    def set_active_coor_system(func):
        @wraps(func)
        def updated(*args, **kwargs):
            args[0].interface.set_coor_sys(args[0].coor_sys)
            return func(*args, **kwargs)
        return updated
    
    def append_lists(self):
        self.modelentities_to_move().append([])
        self.ports_to_move().append([])     
        
    @classmethod
    def append_points(self, coor_list):
        # coor_list can be [()()] or [[][]]
        points = [(coor_list[0][0],coor_list[0][1])]

        for coor in coor_list[1:]:
            points.append((points[-1][0] + coor[0],points[-1][1] + coor[1]))
        return points

    def assign_perfect_E(self, entities, name='perfE'):
        if isinstance(entities, list):
            self.interface.assign_perfect_E(entities, entities[0].name+name)
        else:
            self.interface.assign_perfect_E(entities, entities.name+name)

    def body(self, body_name, coor_name='Global', coor_sys=None):
        if coor_name != 'Global':
            if not(coor_sys is None):
                coor_sys = parse_entry(coor_sys)
                self.interface.create_coor_sys(coor_name, coor_sys)
        N = Network(body_name, coor_name, self.interface, self.variables)
        B = Body(self.interface, coor_name, body_name, N, self.mode, self.variables)
        return B
    
    @set_active_coor_system
    def box_corner(self, pos, size, layer=layer_Default, **kwargs):
        name = self.interface.draw_box_corner(pos, size, **kwargs)
        return ModelEntity(name, 3, self.coor_sys, layer=layer)
    
    @set_active_coor_system
    def box_center(self, pos, size, layer, **kwargs):
        name = self.interface.draw_box_center(pos, size, **kwargs)
        return ModelEntity(name, 3, self.coor_sys, layer=layer)
  
#    def polyline_3D(self, points3D, closed=True, **kwargs): # among kwargs, name should be given
#        if self.mode=='hfss':
#            name = self.interface.draw_polyline(points3D, closed=closed, **kwargs)
#        elif self.mode=='gds':
#            points2D = []
#            for point in points3D:
#                points2D.append(points3D[0:2])
#            name = self.interface.draw_polyline(points2D, [0,0], closed=closed, **kwargs)
#
#        dim = closed + 1 # 2D when closed, 1D when open
#        return ModelEntity(name, dim, self.coor_sys, layer=kwargs['layer'])
#
#    def rect_corner_3D(self, pos3D, size3D, **kwargs):
#        if self.mode=='hffs': 
#            name = self.interface.draw_rect_corner(pos3D, size3D, **kwargs)
#        elif self.mode=='gds':
#            pos2D = pos3D[0:2]
#            size2D = size3D[0:2]
#            name = self.interface.draw_rect_corner(pos2D, size2D, **kwargs)
#        return ModelEntity(name, 2, self.coor_sys, layer=kwargs['layer'])
#    
#    def rect_center_3D(self, pos3D, size3D, **kwargs):
#        if self.mode=='hffs': 
#            name = self.interface.draw_rect_corner(pos3D, size3D, **kwargs)
#        elif self.mode=='gds':
#            pos2D = pos3D[0:2]
#            size2D = size3D[0:2]
#            name = self.interface.draw_rect_corner(pos2D, size2D, **kwargs)
#        return ModelEntity(name, 2, self.coor_sys, layer=kwargs['layer'])
    
#    def polyline_2D(self, points2D, z=0, closed=True, **kwargs): # among kwargs, name should be given
#        points3D = []
#        for point in points2D:
#            points3D.append(point+(z,))
#        return self.polyline_3D(points3D, closed, **kwargs)
#
#    def rect_corner_2D(self, pos2D, size2D, z=0, **kwargs):
#        name = self.interface.draw_rect_corner(pos2D, size2D, **kwargs)
#        return ModelEntity(name, 2, self.coor_sys, z, layer=kwargs['layer'])
#    
#    def rect_center_2D(self, pos2D, size2D, z=0, **kwargs):
#        name = self.interface.draw_rect_center(pos2D, size2D, **kwargs)
#        return ModelEntity(name, 2, self.coor_sys, z, layer=kwargs['layer'])




    
    @set_active_coor_system
    def cylinder(self, pos, radius, height, axis, **kwargs):
        kwargs['coor_sys']=self.coor_sys
        name = self.interface.draw_cylinder(pos, radius, height, axis, **kwargs)
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
    def disk(self, pos, radius, axis, **kwargs):
        kwargs['coor_sys']=self.coor_sys
        if self.mode=='gds':
            pos = self.val(pos)
            radius = self.val(radius)
        name = self.interface.draw_cylinder(pos, radius, axis, **kwargs)
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
            radius = self.val(radius)
        try:
            self.interface._fillet(radius, vertex_index, entity)
        except Exception:
            print("Fillet operation resulted in an error")
    def _fillets(self, radius, entity):
        vertices = self.interface.get_vertex_ids(entity)
        self.interface._fillets(radius, vertices, entity)

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
        kwargs['coor_sys']=self.coor_sys
        i = 0
        while i < len(points[:-1]):
            if np.array_equal(points[i], points[i+1]):
                points.pop(i)
            else:
                i+=1
                
        if self.mode=='gds':
            points = self.val(points)

            
        name = self.interface.draw_polyline(points, closed=closed, **kwargs)
        dim = closed + 1
        return ModelEntity(name, dim, self.coor_sys, layer=kwargs['layer'])
    
    def ports_to_move(self):
        inter1 = Port.instances_to_move
        inter2= Port.find_last_list(inter1)
        return inter2
    
    @set_active_coor_system
    def rect_corner_2D(self, pos, size, **kwargs):
        kwargs['coor_sys']=self.coor_sys
        if self.mode=='gds':
            pos = self.val(pos)
            size = self.val(size)
        name = self.interface.draw_rect_corner(pos, size, **kwargs)
        return ModelEntity(name, 2, self.coor_sys, layer=kwargs['layer'])
    
    @set_active_coor_system
    def rect_center_2D(self, pos, size, **kwargs):
        kwargs['coor_sys']=self.coor_sys
        if self.mode=='gds':
            pos = self.val(pos)
            size = self.val(size)
        name = self.interface.draw_rect_center(pos, size, **kwargs)
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

            
    def separate_bodies(self,name):
        #This looks hard
        pass
    
    def set_current_coor(self, pos, ori):
        self.current_pos = parse_entry(pos)
        self.current_ori = parse_entry(ori)
        
    def set_variable(self, name, value):
        """
        name (str): name of the variable in HFSS e.g. 'chip_length'
        value (str, VarStr, float): value of the variable if str will try to analyse the unit
        """
        self.store_variable(name, value) # for later parsing
        self.__dict__[name] = VariableString(name) # for handy use by user
        if self.mode=='hfss':
            self.design.set_variable(name, value) # for HFSS
    
    def store_variable(self, name, value): # put value in SI
        if not isinstance(value, VariableString):
            if LENGTH == extract_value_dim(value):
                self.variables[name]=extract_value_unit(value, LENGTH_UNIT)
            if INDUCTANCE == extract_value_dim(value):
                self.variables[name]=extract_value_unit(value, INDUCTANCE_UNIT)
            if CAPACITANCE == extract_value_dim(value):
                self.variables[name]=extract_value_unit(value, CAPACITANCE_UNIT)
            if RESISTANCE == extract_value_dim(value):
                self.variables[name]=extract_value_unit(value, RESISTANCE_UNIT)
        else:
            self.variables[name]=value
            
    def subtract(self, blank_entity, tool_entities, keep_originals= False):
        name = self.interface.subtract(blank_entity, tool_entities, keep_originals)
        if not(keep_originals):
            for tool_entity in tool_entities:
                self.delete(tool_entity)
#        return ModelEntity(name, 2, self.coor_sys, layer=blank_entity.layer)
            
    def _sweep_along_path(self, entity_to_sweep, path_entity, name=None):
        new_name = self.interface._sweep_along_path(entity_to_sweep, path_entity, name)
#        self.delete(path_entity)
        entity_to_sweep.modify_dimension(2)
        path_entity.rename_entity(new_name)
        return path_entity
    
    def translate(self, entities, vector=[0,0,0]):
        if self.mode == 'gds':
            vector = self.val(vector)
        self.interface.translate(entities, vector)
        
    def unite(self, entities, name=None, keep_originals=False):
        loc_entities = []     
        dim_Union = 0
        union=None
        for entity in entities:
            if entity!=None and isinstance(entity, ModelEntity):
                union = entity

            if union!=None:    
                    loc_entities.append(union)
                    if union.dimension>dim_Union:
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

    def val(self, name): # to use if you want to compare two litt expressions
                         # Compute the numerical value using for eval_var_str
        if isinstance(name, list) or isinstance(name, tuple):
            name_list = []
            for elt in name:
                name_list.append(self.val(elt))
            if isinstance(name, Vector):
                return Vector(name_list)
            else:
                return name_list
        else:
            return self.eval_var_str(name)
        
    @set_active_coor_system
    def wirebond(self, pos, ori, width, **kwargs):
        kwargs['coor_sys']=self.coor_sys
        name = self.interface.draw_cylinder(pos, ori, width, **kwargs)
        return ModelEntity(name, 2, self.coor_sys, layer=kwargs['layer'])
    
class Port():
    instances_to_move = []
    dict_instances  = {}
    
    def __init__(self, name, pos, ori, track, gap):
        new_name = self.check_name(name)
        self.name = new_name
        self.pos = Vector(pos)
        self.ori = Vector(ori)
        self.track = parse_entry(track)
        self.gap = parse_entry(gap)
        Port.find_last_list(Port.instances_to_move).append(self)
        self.dict_instances[name]=self
        

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
            



class Network(PythonMdlr):
    variables= None
    to_bond=[]

    def __init__(self, name, coor_sys, interface, variables):
        self.interface = interface
        self.coor_sys = coor_sys
        self.name = name
        self.variables = variables
#        self.__class__.variables = variables
    
            
    def update(self,coor_sys):
        self.coor_sys = coor_sys

    


            
    def port(self, name, pos, ori, track, gap):
        return Port(name, pos, ori, track, gap)
    

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
        
        
        
        iLength, iWidth,iSize = parse_entry((iLength, iWidth, iSize))
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
#        print(points3)
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
    

    


    
@add_methods_from(KeyElement, CustomElement)
class Body(PythonMdlr):

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
            compteur = 0
            for i, argument in enumerate(args[2:]):
                if isinstance(argument, str) and (argument in Port.dict_instances):
                    new_args.append(Port.dict_instances[argument])
                    compteur+=1
                else:
                    new_args.append(argument)

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

                
            elif compteur==1:
#                print(new_args[2].pos)
                args[0].set_current_coor(new_args[2].pos, new_args[2].ori)
            elif compteur==2:
                args[0].set_current_coor(1/2*(new_args[2].pos+new_args[3].pos), new_args[2].ori)
            new_args = tuple(new_args)
            return KeyElement._moved(func,previous_pos, previous_ori,*new_args, **kwargs)
        return moved  
    
    def port(self, name, pos, ori, track, gap):
        pos, ori, track, gap = parse_entry((pos, ori, track, gap))
#        self.rect_center_2D(pos, [track*ori[0]/2+track*ori[1], track*ori[0]+track/2*ori[1]], name=name+'track', layer='PORT')
#        self.rect_center_2D(pos, [track*ori[0]/2+(2*gap+track)*ori[1], (2*gap+track)*ori[0]+track/2*ori[1]], name=name+'gap', layer='PORT')
        self.rect_corner_2D(pos, [(ori[0]+0.5)*(1e-5), (ori[1]+0.5)*(1e-5)],  name=name+'ori', layer=layer_Default)
#        
        return self.network.port(name, pos, ori, track, gap)
    
    def double_port(self, name, pos, ori, track, gap):
        port1 = self.port( name+'_front', pos, ori, track, gap)
        port2 = self.port( name+'_back', pos, -Vector(ori), track, gap)
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
        track, gap = parse_entry((track, gap))
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
            induc_H = self.val(iInduct)*1e-9
            print('induc'+str(induc_H))
            w_plasma = 2*np.pi*24*1e9
            capa_fF = 1/(induc_H*w_plasma**2)*10**15
            capa_plasma = self.set_variable(name+'_capa_plasma', str(capa_fF)+'fF')
            print(capa_fF)
            print(capa_plasma)
        else:
            capa_plasma = 0
        
        # No parsing needed, should not be called from outside
        iTrack1 = parse_entry(iInPort.track)
        iTrack2 = parse_entry(iOutPort.track)
        
        adaptDist1 = iTrack1/2-iTrackJ/2
        adaptDist2 = iTrack2/2-iTrackJ/2
        if self.val(adaptDist1)>self.val(iLength/2-iTrackJ/2) or self.val(adaptDist2)>self.val(iLength/2-iTrackJ/2):
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
            if self.val(iTrack1) > self.val(iTrack2):
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
    def draw_cable(self, name, *ports, fillet="0.3mm", is_bond=False, is_meander=False, to_meanders = [1,0,1,0,1,0,1,0,1,0], meander_length=0, meander_offset=0, is_mesh=False, reverse_adaptor=False, layer=layer_Default):
            '''
            Draws a CPW transmission line between iIn and iOut
    
            if iIn and iOut are facing eachother, and offset,
            draws a transmission line with two elbows half-way in between.
    
            if iIn and iOut are perpendicular,
            draws a transmission line with one elbow.
    
            if iIn and iOut do not have the same track/ gap size, this function calls
            drawAdaptor before iOut.
    
            N.B: do not separate the two ports by their track or gap size.
    
            Inputs:
            -------
            name: (string) base-name of object, draws 'name_adaptor' etc
            iIn: (tuple) input port
            iOut: (tuple) output port
            iMaxfillet: (float), maximum fillet radius
            reverseZ: performs a mirror operation along Z --> useful only when the thickening operation goes in the wrong direction
            
            '''
            #TODO change the format of the arguments
            fillet, meander_length, meander_offset=parse_entry((fillet, meander_length, meander_offset))

            self.to_bond=[]
            adaptor_length=0
            track_adaptor = None
            inTrack = self.variables[ports[0].track]
            outTrack = self.variables[ports[-1].track]
            inGap = self.variables[ports[0].gap]
            outGap = self.variables[ports[-1].gap]
            fillet = self.variables[fillet]
            
            if (not equal_float(inTrack, outTrack)) or (not equal_float(inGap, outGap)):
                if reverse_adaptor:
                    if (inTrack+inGap) > (outTrack+outGap):
                        iIn, adaptor_length, track_adaptor, gap_adaptor, mask_adaptor = self.draw_adaptor('reverse',ports[0].name,ports[-1].track, ports[-1].gap)
                    else:
                        iOut, adaptor_length, track_adaptor, gap_adaptor, mask_adaptor = self.draw_adaptor('reverse',ports[0].name,ports[-1].track, ports[-1].gap)
                else:
                    if (inTrack+inGap) > (outTrack+outGap):
                        iOut, adaptor_length, track_adaptor, gap_adaptor, mask_adaptor = self.draw_adaptor('reverse',ports[0].name,ports[-1].track, ports[-1].gap)
                    else:
                        iIn, adaptor_length, track_adaptor, gap_adaptor, mask_adaptor = self.draw_adaptor('reverse',ports[0].name,ports[-1].track, ports[-1].gap)
      
#            all_constrains = []
#            for constrain in constrains:
#                all_constrains.append(Port.dict_instances[constrain])                
#                all_constrains.append(Port.dict_instances[constrain])

#                all_constrains.append(constrain)#[self.ports[constrain][POS], self.ports[constrain][ORI], self.ports[constrain][TRACK], self.ports[constrain][GAP]])
#                previous modification to tackle the case where a cable is drawn between two ports defined on the fly
            nb_constraints = len(ports)//2-1
            if not isinstance(to_meanders[0], list):
                to_meanders = [to_meanders for ii in range(nb_constraints+1)]
    #        print(to_meanders)
                
            
            cable_length = []
            tracks = []
            gaps = []
            masks = []
#            port_names = [iInPort]+all_constrains+[iOutPort] # List of ports
#            print(meander_length)
            for ii in range(nb_constraints+1):
                to_meander = to_meanders[ii]
                if isinstance(meander_length, (list, np.ndarray)):
                    m_length = meander_length[ii]
                else:
                    m_length = meander_length
                if nb_constraints!=0:
                    to_add = '_'+str(ii)
                else:
                    to_add = ''
#                print(port_names[2*ii:2*ii+2])
#                self.__init__(self.name, *port_names[2*ii:2*ii+2]) CONNECT ELEMENT USELESS
                
                points = self.find_path('points', ports[2*ii].name, ports[2*ii+1].name, fillet, is_meander, to_meander, m_length, meander_offset)
                connection_track = self.polyline_2D(points, name=name+'path_track'+to_add, closed=False ,layer = layer_TRACK)
    #            print('length_adaptor = %.3f'%(self.val(adaptor_length)*1000))
                cable_length.append(self.length(points, 0, len(points)-1, fillet)+self.val(adaptor_length))
#                self._fillets(fillet-eps, connection_track)
        
                connection_gap = self.polyline_2D(points, name=name+'path_gap'+to_add, closed=False ,layer = layer_Default)
#                self._fillets(fillet-eps, connection_gap)
#                self.set_current_coor([0,0], [0, -1])
                track_starter = self.polyline_2D([[-ports[2*ii].track,0],[ports[2*ii].track,0]], name = name+'start_gap', layer = 0)
                gap_starter = self.polyline_2D([[-ports[2*ii].gap,0],[ports[2*ii].gap,0]], name = name+'start_track', layer = 1)
                
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
    
                
            if is_bond:
                self.draw_bond((ports[0].track+ports[0].gap*2)*1.5)
            
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
    def find_path(self, name, iInPort, iOutPort, fillet, is_meander, to_meander, meander_length, meander_offset):
        iIn_pos = Vector(iInPort.pos)
        iIn_ori = Vector(iInPort.ori)
        iOut_pos = Vector(iOutPort.pos)
        iOut_ori = Vector(iOutPort.ori)
        room_bonding = 0*100e-6 #SMPD MANU BOND SPACE
        print(str(1.1*fillet))
        
        pointA = self.val(iIn_pos+iIn_ori*room_bonding)
        pointB = self.val(iOut_pos+iOut_ori*room_bonding)
        
        
        point1 = self.val(iIn_pos+iIn_ori*(1.1*fillet+room_bonding))
        point2 = self.val(iOut_pos+iOut_ori*(1.1*fillet+room_bonding))
#        print(point1)
#        print(point2)
#        print(iIn_ori)
        
        def next_point(point1, point2, vec):
            choice1 = point1+vec*var((point2-point1).dot(vec))
            choice2 = point1+vec.orth()*var((point2-point1).dot(vec.orth()))
            return [[point1, choice1, point2], [point1, choice2, point2]]


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
            
        def cost_f(x):
            if x==1:
                return 0
            elif x==0:
                return 1
            else:
                return 100

        def check(points):
            length = 0
            prev_point = Vector(points[0])
            _points = Vector([points[0]])
            vecs = Vector([])
            for point in points[1:]:
                if not ( equal_float(point[0], prev_point[0]) and equal_float(point[1], prev_point[1])):
                   #☺ TODO vec should already be a Vector
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

        final_choice= None
        cost=np.inf
        for ii, choice in enumerate(points_choices):
            new_cost, new_choice, new_length = check(choice)
            if new_cost<cost:
                final_choice = new_choice
                cost = new_cost
                length = new_length


        length_fillet = length - cost*(2-np.pi/2)*fillet
        n_fillet = 10
        dist_fillet = length_fillet/n_fillet

        float_final_choice = []
        for point in final_choice:
            float_final_choice.append(self.val(point))

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
                vec = way(self.val(B-A))
                if self.val(AB).norm() > self.val(min_dist):
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
                    vec = way(self.val(B-A))
                    if self.val(AB).norm() > self.val(min_dist):
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
                vecs.append(way(self.val(B-A)))
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
                vec = way(self.val(B-A))
                AB = (B-A).norm()
                n_add = int(self.val(AB/min_dist))
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
#                        vec = way(self.val(B-A))
#                        AB = (B-A).norm()
#                        if ii==0 or ii==n_points-2:
#                            factor = 0.5
#                        else:
#                            factor = 1
#                        n_add = int(self.val(AB/min_dist)-factor)
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
            if np.abs(self.val(displacement))<self.val(min_dist)*1.1:
                displacement = min_dist*1.1
            points, indices_corners, dist, ignore = add_points(points, rl, min_dist, n_meander=n_meander)
            new_points = [points[0]]
            parity = 1
            if indices_corners is not None:
                for ii, B in enumerate(points[1:-1]):
                    A = points[ii]
                    AB= B-A
                    vec = way(self.val(AB))
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
                    vec=way(self.val(AB))
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

        if is_meander:
            min_dist = 2*fillet
            final_choice = meander(final_choice, min_dist, to_meander, meander_length, meander_offset)
            
        final_choice =  [self.val(iIn_pos)] + final_choice + [self.val(iOut_pos)]

# Needed to draw Manu bond
        def add_fillet_points(points, fillet):
            new_points = [points[0]]
            for ii, point in enumerate(points[1:-1]):
                index = ii+1
                p_vec = points[index-1]-point
                n_vec = points[index+1]-point
                new_points.append(point+way(self.val(p_vec))*fillet)
                new_points.append(point+way(self.val(n_vec))*fillet)
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


        def return_bonds(points, fillet, length_fillet, n_fillet): #lengh_fillet is the cable lenght with filleting taken into account
            # create bond at half dist_fillet
            prev_ori = way(self.val(points[1]-points[0]))
            unit_dist_fillet = length_fillet/n_fillet
            dist_fillet = unit_dist_fillet/2 #starting dist for the fillet
            for ii in range(n_fillet):
                indices, kind, remain = where(points, dist_fillet, fillet)
                A = points[indices[0]]
                B = points[indices[1]]
                if kind=='normal':
                    pos = A + remain*(B-A).unit()
                    ori = way(self.val(B-A))
                    width = 0.0004
                    self.draw_wirebond('wire', pos, ori, width)
                    prev_ori = ori
                else:
                    next_ori=way(self.val(points[indices[1]+1]-B)) #should be fine, if we have a fillet we have some straight portion after
                    print(f'kind={kind}')
                    ex = next_ori
                    ey = prev_ori
                    print(f'ex={ex}')
                    print(f'ey={ey}')
                    pos_center = A + ex*(B-A).dot(ex)
                    print(pos_center)
                    theta = remain/fillet
                    print(theta*180/np.pi)
                    pos = pos_center - ex*np.cos(theta)*fillet + ey * np.sin(theta)*fillet
                    print(f'pos={pos}')
                    ori = ey*np.cos(theta) + ex*np.sin(theta)
                    print(f'ori={ori}')
                    width = 0.0004
                    self.draw_wirebond('wire', pos, ori, width)
                dist_fillet += unit_dist_fillet

        _, final_choice, _ = check(final_choice)

        to_bond_points = add_fillet_points(final_choice, fillet)
        for ii, point in enumerate(to_bond_points[::2]):
            points = [to_bond_points[2*ii], to_bond_points[2*ii+1]]
#            print("points", points)
            self.polyline_2D(points , closed=False, name='bef_test'+str(ii), layer=layer_Default)
            self.to_bond.append(points)
        return final_choice     
        
    def length(self, points, A, B, fillet): # A and B are integer point indices
#        for point in points:
#            print(self.val(point[0]), self.val(point[1]))
        if A<0 or A>=len(points):
            raise ValueError('First index should be within the point list')
        if B<0 or B>=len(points):
            raise ValueError('Second index should be within the point list')
        if A==B:
            return 0
        if A<B:
            value = 0
            for ii in range(B-A):
                value+=self.val((points[A+ii+1]-points[A+ii]).norm())
            return value-(B-A-1)*self.val(fillet*(2-np.pi/2))
        else:
            return self.length(points, B, A, fillet)
    
    
    def draw_bond(self, width, min_dist='0.5mm'):
        width, min_dist = parse_entry((width, min_dist))

        min_dist = self.val(min_dist)
        for elt in self.to_bond:
            A = elt[0]
            B = elt[1]
            val_BA = self.val(B-A)
            ori = way(val_BA)
            length = Vector(val_BA).norm()
            n_bond = int(length/min_dist)+1
            spacing = (B-A).norm()/n_bond
            pos = A+ori*spacing/2
            self.draw_wirebond('wire', pos, ori, width)
            for ii in range(n_bond-1):
                pos = pos + ori*spacing
                self.draw_wirebond('wire', pos, ori, width)
              

    @move_port
    def _connect_snails2(self, name, iInPort, iOutPort, squid_size, width_top, n_top, width_bot, n_bot, N, width_bridge, spacing_bridge, litho='opt'):
        
        width_track = iInPort.track # assume both are equal
        spacing = (iOutPort.pos-iInPort.pos).norm()
        
        width_snail = squid_size[0]+6*width_track #ZL
        
        tot_width = width_snail*N
        if np.abs(self.val(tot_width))>np.abs(self.val(spacing)):
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