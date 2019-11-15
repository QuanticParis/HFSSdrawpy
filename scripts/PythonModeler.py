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



#IMPORT GDS / HFSS Modelers
from hfss import extract_value_unit, \
                 extract_value_dim, \
                 parse_entry, \
                 get_active_project, \
                 ModelEntity, \
                 VariableString 
import gds_modeler
#IMPORT KEY / CUSTOM Elements
import Lib
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


class PythonMdlr():
    is_overdev = False
    is_litho = False
    is_mask = False
    pcb_track = parse_entry('300um')
    pcb_gap = parse_entry('200um')
    gap_mask = parse_entry('20um')
    overdev = parse_entry('0um')

    def __init__(self, name_interface): #"Hfss" or "Gds"
        self.ports = {}
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
    

    def body(self, body_name, coor_name='Global', coor_sys=None):
        if coor_name != 'Global':
            if not(coor_sys is None):
                coor_sys = parse_entry(coor_sys)
                self.interface.create_coor_sys(coor_name, coor_sys)
        return Body(self.interface, coor_name, body_name, self.mode, self.variables), Network(self.interface, coor_name, body_name, self.mode, self.variables)
    
    @classmethod
    def append_points(self, coor_list):
        # coor_list can be [()()] or [[][]]
        points = [(coor_list[0][0],coor_list[0][1])]

        for coor in coor_list[1:]:
            points.append((points[-1][0] + coor[0],points[-1][1] + coor[1]))
        return points
    
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
    
    def move_points(self, coor_list, move, absolute=False):
        points=[]
        for ii, coor in enumerate(coor_list):
            if ii>0 and not absolute:
                move=[0,0]
            points.append(Vector(*coor)+Vector(*move))
        return points


    def coor(self, vec): # Change of coordinate for a point
        return self.rot(*vec)+self.pos

    def coor_vec(self, vec): # Change of coordinate for a vector
        return self.rot(*vec)
    def set_units(self, units='m'):
        if (self.modeler != None):
            self.modeler.set_units('mm')
            self.interface = self.modeler
    
    def generate_gds(self, name_file):
        assert type(self.interface)==GdsModeler
        self.interface.generate_gds(name_file)
    
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
            
    @classmethod        
    def eval_var_str(self, name, unit=None): # return numerical value of a given expression
                                             # using the values stored in self.variables
        # TODO: parse several times
        # can only parse 2 times for now
        print(self.variables)
        print(name)
        try:
            _val = float(eval(str(sympy_parser.parse_expr(str(sympy_parser.parse_expr(str(name), self.variables)), self.variables))))
        except Exception:
            msg = ('Parsed expression contains a string which does '
                             'not correspond to any design variable')
            raise ValueError(msg)
        if unit is not None:
            _val = str(_val)+' '+unit
        return _val
    
    def val(self, name): # to use if you want to compare two litt expressions
                         # Compute the numerical value using for eval_var_str
        if isinstance(name, list):
            name_list = []
            for elt in name:
                name_list.append(self.eval_var_str(elt))
            if isinstance(name, Vector):
                return Vector(name_list)
            else:
                return name_list
        else:
            return self.eval_var_str(name)

    def box_corner(self, pos, size, layer, **kwargs):
        name = self.interface.draw_box_corner(pos, size, **kwargs)
        return ModelEntity(name, 3, self.coor_sys, layer=layer)
    
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



    def polyline_2D(self, points2D, z=0, closed=True, **kwargs): # among kwargs, name should be given
        name = self.interface.draw_polyline(points2D, [0,0], closed=closed, **kwargs)
        dim = closed + 1
        return ModelEntity(name, dim, self.coor_sys, layer=kwargs['layer'])

    def rect_corner_2D(self, pos2D, size2D, z=0, **kwargs):
        name = self.interface.draw_rect_corner(pos2D, size2D, **kwargs)
        return ModelEntity(name, 2, self.coor_sys, z, layer=kwargs['layer'])
    
    def rect_center_2D(self, pos2D, size2D, z=0, **kwargs):
        name = self.interface.draw_rect_center(pos2D, size2D, **kwargs)
        return ModelEntity(name, 2, self.coor_sys, z, layer=kwargs['layer'])
    
    
    
    
    def rot(self, x, y=0):
        return Vector(x, y).rot(self.ori)

    def cylinder(self, name, pos, radius, height, axis, layer, **kwargs):
        name = self.interface.draw_cylinder(pos, radius, height, axis, **kwargs)
        return ModelEntity(name, 3, self.coor_sys, layer=kwargs['layer'])
 
    def disk(self, name, pos, radius, axis, layer, **kwargs):
        name = self.interface.draw_cylinder(pos, radius, axis, **kwargs)
        return ModelEntity(name, 2, self.coor_sys, layer=kwargs['layer'])
    
    def wirebond(self, name, pos, ori, width, layer, **kwargs):
        name = self.interface.draw_cylinder(pos, ori, width, **kwargs)
        return ModelEntity(name, 2, self.coor_sys, layer=kwargs['layer'])
    
    def connect_faces(self, name, entity1, entity2):
        assert entity1.dimension == entity2.dimension
        assert entity1.dimension == 2
        return ModelEntity(name, 3, entity1.coor_sys, entity1.model)
        
    def delete(self, entity):
        del entity
        
    def rename_entity(self, entity, name):
        entity.rename_entity(name)
        
    def unite(self, entities, name=None, keep_originals=False):
        loc_entities = entities.copy()
        dim_Union = 0;
        for entity in loc_entities:
            if isinstance(entity, tuple) or isinstance(entity, list):
                for element in entity:
                    loc_entities.append(element)
            if isinstance(entity, ModelEntity):
                if entity.dimension>dim_Union:
                    dim_Union = entity.dimension
                
                if not(keep_originals):
                    self.delete(entity)
                 
        if name is None:
            union = entities[0].copy(dim_Union)
        else:
            union = entities[0].copy(dim_Union)
            union.rename_entity(name)
        
        
        return union
        
    def intersect(self, entities, keep_originals = False):
        dim_Intersection = 3;
        for entity in entities:
            if entity.dimension<dim_Intersection:
                dim_Intersection = entity.dimension
        if dim_Intersection == 0:
            pass
        elif dim_Intersection == 1:
            Intersection = ModelEntity(entities[0].name, 1, entities[0].coor_sys, entities[0].model)
        elif dim_Intersection == 2:
            Intersection = ModelEntity(entities[0].name, 2, entities[0].coor_sys, entities[0].model, entities[0].boundaries)
        elif dim_Intersection == 3:
            Intersection = ModelEntity(entities[0].name, 3, entities[0].coor_sys, entities[0].model, entities[0].boundaries)
   
        #A factoriser
        if not(keep_originals):
            for entity in entities:
                self.delete(entity)
        
        return Intersection
    
    def subtract(self, blank_entity, tool_entities):
        self.interface.subtract(blank_entity, tool_entities)
        for i in tool_entities:
            self.delete(i)

    def translate(self, entities, vector=[0,0,0]):
        self.interface.translate(entities, vector)
    
    def rotate(self, entities, angle=0):
        if isinstance(angle, list):
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
    
    def _fillet(self, radius, vertex_index, entity):
        self.interface._fillet(radius, vertex_index, entity.name)
#        print(vertices)
#        print(radius)
    
    def assign_perfect_E(self, entity):
        entity.boundaries.append('PerfE')
    
    def mirrorZ(self, obj):
        pass
    
    def copy(self, entity):
        if entity.dimension ==0:
            pass
        elif entity.dimension ==1:
            return ModelEntity(entity.name, entity.coor_sys, entity.model)
        elif entity.dimension ==2:
            return ModelEntity(entity.name, entity.coor_sys, entity.model, entity.boundaries)
        elif entity.dimension ==3:
            return ModelEntity(entity.name, entity.coor_sys, entity.model, entity.boundaries)
        
    def duplicate_along_line(self, entity, vec):
        #Create clones
        pass
    
    def assign_lumped_RLC(self, entity, ori, parameters):
#        entity.boundaries.append('lumpedRLC')
# TODO
        pass
    
    def make_material(self, entity, material):
        # Problème pour les unions de matériaux si on ajoute un attribut
        pass
    
    def delete_all_objects(self, entities):
        for entity in entities:
            self.delete(entity)
    
    def mesh_zone(self, zone, mesh_length):
        mesh_length = parse_entry(mesh_length)
        self.interface.assign_mesh_length(zone, mesh_length)

        
        
    
#    def rectangle(pos,size):
#        interface.rectangle(pos,size)
#        return ModelEntity(name)
#    
#    def chip():
#        return chip(repere,modeler)
#    

#        
#class GdsModeler():
#    def rectangle():
#        plt.plot(pos,size)

class Port():
    instances = {}
    def __init__(self, name, pos, ori, track, gap):
        self.name = name
        self.pos = Vector(pos)
        self.ori = Vector(ori)
        self.track = track
        self.gap = gap
        self.__class__.instances[name] = self


class Network(PythonMdlr):
    variables= None
    def __init__(self, interface, coor_sys, name, mode, variables):
        self.interface = interface
        self.coor_sys = coor_sys
        self.name = name
        self.mode = mode # 'hfss' or 'gds'
        self.variables = variables
        self.__class__.variables = variables
        
    def decorator(func):
        @wraps(func)
        def decorated(*args, **kwargs):
            # parse the instructions of the user
            self1  = args[0]
            name = args[1]
            iIn = args[2]
            iOut = args[3]
            print('\n')
            print('Ports ', Port.instances)
            print('\n')
            if iIn in Port.instances:
                iInPort = Port.instances[iIn].__dict__
            else:
                raise ValueError('inPort %s does not exist' % iIn)
            if iOut in Port.instances:
                iOutPort = Port.instances[iOut].__dict__
            else:
                raise ValueError('outPort %s does not exist' % iOut)
            return func(*((self1, name, iInPort, iOutPort)+args[4:]))
        
        return decorated
    
    def port(self, name, pos, ori, track, gap):
        return Port(name, pos, ori, track, gap)
    @decorator
    def _connect_JJ(self, name, iInPort, iOutPort, iTrackJ, iInduct='1nH', fillet=None):
        '''
        Draws a Joseph's Son Junction.

        Draws a rectangle, here called "junction",
        with Bondary condition :lumped RLC, C=R=0, L=iInduct in nH
        Draws needed adaptors on each side

        Inputs:
        -------
        name:
        iIn: (tuple) input port
        iOut: (tuple) output port - None, ignored and recalculated
        iSize: (float) length of junction
        iWidth: (float) width of junction
        iLength: (float) distance between iIn and iOut, including
                 the adaptor length
        iInduct: (float in nH)

        Outputs:
        --------

        '''
        iLength = (iOutPort['pos']-iInPort['pos']).norm()
        
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
        iInPort['pos'] = (iInPort['pos']+iOutPort['pos'])/2
        iTrack = parse_entry(iInPort['track']) # assume both track are identical
        iTrackJ = parse_entry(iTrackJ)
        adaptDist = iTrack/2-iTrackJ/2

        if self.val(adaptDist)>self.val(iLength/2-iTrackJ/2):
            raise ValueError('Increase iTrackJ %s' % name)


        raw_points = [(iTrackJ/2, iTrackJ/2),
                      ((iLength/2-iTrackJ/2-adaptDist), 0),
                      (adaptDist, (iTrack-iTrackJ)/2),
                      (0, -iTrack),
                      (-adaptDist, (iTrack-iTrackJ)/2),
                      (-(iLength/2-iTrackJ/2-adaptDist), 0)]
        points = self.append_points(raw_points)
#        
#        fig, ax = plt.subplots()
#        patches = []
#        patches.append(Polygon(points))
#        p = PatchCollection(patches, alpha=0.4)
#        ax.add_collection(p)
        print(raw_points)
        print(points)
        right_pad = self.polyline_2D(points, name =name+"_pad1", layer="TRACK")

        points = self.append_points(self.refy_points(raw_points))
        left_pad = self.polyline_2D(points, name=name+"_pad2", layer="TRACK")

        pads = self.unite([right_pad, left_pad], name=name+'_pads')
        
        if not self.is_litho:
            mesh = self.rect_center_2D([0,0], [iLength, iTrack], name=name+'_mesh', layer='MESH')
            self.mesh_zone(mesh, iTrackJ/2)
    
            points = self.append_points([(iTrackJ/2,0),(-iTrackJ,0)])
            self.polyline_2D(points, closed = False, name=name+'_line', layer="MESH")
    
            JJ = self.rect_center_2D([0,0], [iTrackJ, iTrackJ], name=name, layer="RLC")
            #TODO
            #self.assign_lumped_RLC(JJ, self.ori, (0, iInduct, capa_plasma))

        return pads
    
    @decorator
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
        retIn = [[0,0], 0, iInPort['track'], iInPort['gap']]
        retOut = [[0,0]+[iInPort['gap']+iOutPort['gap']+iSize+2*iWidth,0],0 , iOutPort['track'], iOutPort['gap']]

        points1 = self.append_points([(iInPort['gap']+iWidth, 0),
                                     (0, -iLength/2),
                                     (-iWidth, 0),
                                     (0, iLength/2-iInPort['track']/2),
                                     (-iInPort['gap'], 0),
                                     (0, iInPort['track']),
                                     (iInPort['gap'], 0),
                                     (0, iLength/2-iInPort['track']/2),
                                   (iWidth, 0)])

    
        trackIn = self.polyline_2D(points1, name=name+'_track1', layer='TRACK')
        self.chip.trackObjects.append(trackIn)

        points2 = self.append_points([(iInPort['gap']+iWidth+iSize, 0),
                                     (0, -iLength/2),
                                     (+iWidth, 0),
                                     (0, iLength/2-iOutPort['track']/2),
                                     (+iOutPort['gap'], 0),
                                     (0, iOutPort['track']),
                                     (-iOutPort['gap'], 0),
                                     (0, iLength/2-iOutPort['track']/2),
                                     (-iWidth, 0)])
        trackOut = self.polyline_2D(points2, name=name+'_track2', layer='TRACK')
        self.chip.trackObjects.append(trackOut)

        points3 = self.append_points([(0, 0),
                                     (0, iLength/2+iInPort['gap']),
                                     (iInPort['gap']+iWidth+iSize/2, 0),
                                     (0, iOutPort['gap']-iInPort['gap']),
                                     (iOutPort['gap']+iWidth+iSize/2, 0),
                                     (0, -iLength-2*iOutPort['gap']),
                                     (-(iOutPort['gap']+iWidth+iSize/2),0),
                                     (0, iOutPort['gap']-iInPort['gap']),
                                     (-(iInPort['gap']+iWidth+iSize/2),0)
                                     ])
        print(points3)
        gap1 = self.polyline_2D(points3, name=name+'_gap1', layer='GAP')
        self.chip.gapObjects.append(gap1)
#    
#    

        
        
#        TODO !!
#        if not self.is_litho:
#            self.draw(self.name+"_mesh", points)
#            self.modeler.assign_mesh_length(self.name+"_mesh",1/2*iLength)

        self.iIn = retIn
        self.iOut = retOut
#        return [retIn, retOut]

    @decorator
    def find_slanted_path(self, name, iInPort, iOutPort):
        
        iIn_pos = Vector(iInPort['pos'])
        iIn_ori = Vector(iInPort['ori'])
        iOut_pos = Vector(iOutPort['pos'])
        iOut_ori = Vector(iOutPort['ori'])
        
        if iIn_ori.dot(iOut_ori)!=-1:
            raise ValueError('Cannot find slanted path: ports are not oriented correctly')
        else:
            dist = (iOut_pos-iIn_pos).dot(iIn_ori)
        
        pointA = iIn_pos+iIn_ori*dist/3
        pointB = iOut_pos+iOut_ori*dist/3
        
        self.to_bond.append([iIn_pos, pointA])
        self.to_bond.append([pointB, iOut_pos])
        return [iIn_pos, pointA, pointB, iOut_pos], dist/3      
    
    @decorator
    def find_path(self, name, iInPort, iOutPort, fillet, is_meander, to_meander, meander_length, meander_offset):
        
        iIn_pos = Vector(iInPort['pos'])
        iIn_ori = Vector(iInPort['ori'])
        iOut_pos = Vector(iOutPort['pos'])
        iOut_ori = Vector(iOutPort['ori'])
        
        room_bonding = 0*100e-6 #SMPD MANU BOND SPACE
        print(str(1.1*fillet))
        
        pointA = iIn_pos+iIn_ori*room_bonding
        pointB = iOut_pos+iOut_ori*room_bonding

        point1 = iIn_pos+iIn_ori*(1.1*fillet+room_bonding)
        point2 = iOut_pos+iOut_ori*(1.1*fillet+room_bonding)

        def next_point(point1, point2, vec):
            choice1 = point1+vec*np.dot(point2-point1,vec)
            choice2 = point1+vec.orth()*(point2-point1).dot(vec.orth())
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
                if not equal_float(self.val(point)[0], self.val(prev_point)[0]) or not equal_float(self.val(point)[1], self.val(prev_point)[1]):
                   #☺ TODO vec should already be a Vector


                    vec = self.val(point-prev_point)
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
            
        final_choice =  [iIn_pos] + final_choice + [iOut_pos]

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
#            self.draw('bef_test', [to_bond_points[2*ii], to_bond_points[2*ii+1]], closed=False)
            self.to_bond.append([to_bond_points[2*ii], to_bond_points[2*ii+1]])

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
        
    def cable_starter(self, width = 'track', index=None, border=parse_entry('15um')): # width can also be 'gap'
        if width=='track' or width=='Track':
            points = self.append_points([(0, self.inTrack/2),
                                         (0, -self.inTrack)])
        elif width=='gap' or width=='Gap':
            points = self.append_points([(0, self.inGap+self.inTrack/2),
                                         (0, -2*self.inGap-self.inTrack)])
        elif width=='mask' or width=='Mask':
            points = self.append_points([(0, self.inGap+self.inTrack/2+self.gap_mask),
                                         (0, -2*self.inGap-self.inTrack-2*self.gap_mask)])
        elif width=='dc_track':
            points = self.append_absolute_points([(0, self.rel_posIn[index]-self.widIn[index]/2),\
                                                  (0, self.rel_posIn[index]+self.widIn[index]/2)])
        elif width=='dc_gap':
            points = self.append_absolute_points([(0, self.rel_posIn[index]-self.widIn[index]/2-border),\
                                                  (0, self.rel_posIn[index]+self.widIn[index]/2+border)])
        elif width=='dc_cutout':
            points = self.append_absolute_points([(0, -self.cutIn/2),(0, self.cutIn/2)])
            
        return self.draw(self.name+'_width_'+width+'_'+str(index), points, closed=False) #used to be +'_width'

    def draw_cable(self, name, iInPort, iOutPort, fillet="0.3mm", is_bond=False, is_meander=False, to_meanders = [1,0,1,0,1,0,1,0,1,0], meander_length=0, meander_offset=0, is_mesh=False, constrains=[], reverse_adaptor=False, layer=None):
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
            
            fillet, meander_length, meander_offset=parse_entry((fillet, meander_length, meander_offset))
    #        inPos, inOri = Vector(iIn[POS]), Vector(iIn[ORI])
    #        outPos, outOri = Vector(iOut[POS]), Vector(iOut[ORI])
    #        _, _, track, gap = iIn
            self.to_bond=[]
            adaptor_length=0
            track_adaptor = None
            if (not equal_float(self.val(self.inTrack), self.val(self.outTrack))) or (not equal_float(self.val(self.inGap), self.val(self.outGap))):
                if reverse_adaptor:
                    if self.val(self.inTrack+self.inGap) > self.val(self.outTrack+self.outGap):
                        adaptor = ConnectElt(self.name+'_adaptor', self.iIn, [self.outTrack, self.outGap])
                        iIn, adaptor_length, track_adaptor, gap_adaptor, mask_adaptor = adaptor.draw_adaptor()
                        self.__init__(self.name, iIn, self.iOut)
                    else:
                        adaptor = ConnectElt(self.name+'_adaptor', self.iOut, [self.inTrack, self.inGap])
                        iOut, adaptor_length, track_adaptor, gap_adaptor, mask_adaptor = adaptor.draw_adaptor()
                        self.__init__(self.name, self.iIn, iOut)
                else:
                    if self.val(self.inTrack+self.inGap) > self.val(self.outTrack+self.outGap):
                        adaptor = ConnectElt(self.name+'_adaptor', self.iOut, [self.inTrack, self.inGap])
                        iOut, adaptor_length, track_adaptor, gap_adaptor, mask_adaptor = adaptor.draw_adaptor()
                        self.__init__(self.name, self.iIn, iOut)
                    else:
                        adaptor = ConnectElt(self.name+'_adaptor', self.iIn, [self.outTrack, self.outGap])
                        iIn, adaptor_length, track_adaptor, gap_adaptor, mask_adaptor = adaptor.draw_adaptor()
                        self.__init__(self.name, iIn, self.iOut)
      
            all_constrains = []
            for constrain in constrains:
                all_constrains.append([self.ports[constrain][POS], -self.ports[constrain][ORI], self.ports[constrain][TRACK], self.ports[constrain][GAP]])
                all_constrains.append(constrain)#[self.ports[constrain][POS], self.ports[constrain][ORI], self.ports[constrain][TRACK], self.ports[constrain][GAP]])
                # preivous modification to tackle the case where a cable is drawn between two ports defined on the fly
                
            if not isinstance(to_meanders[0], list):
                to_meanders = [to_meanders for ii in range(len(constrains)+1)]
    #        print(to_meanders)
                
            
            cable_length = []
            tracks = []
            gaps = []
            masks = []
            port_names = [self.iIn]+all_constrains+[self.iOut]
            print(port_names)
            for ii in range(len(constrains)+1):
                to_meander = to_meanders[ii]
                if isinstance(meander_length, (list, np.ndarray)):
                    m_length = meander_length[ii]
                else:
                    m_length = meander_length
                if len(constrains)!=0:
                    to_add = '_'+str(ii)
                else:
                    to_add = ''
                print(port_names[2*ii:2*ii+2])
                self.__init__(self.name, *port_names[2*ii:2*ii+2])
                
                points = self.find_path(fillet, is_meander, to_meander, m_length, meander_offset)
                connection = self.draw(self.name+'_track'+to_add, points, closed=False)
    #            print('length_adaptor = %.3f'%(self.val(adaptor_length)*1000))
                cable_length.append(self.length(points, 0, len(points)-1, fillet)+self.val(adaptor_length))
                connection.fillets(fillet-eps)
        
                connection_gap = connection.copy(self.name+"_gap"+to_add)
        
                track_starter = self.cable_starter('track')
                gap_starter = self.cable_starter('gap')
            
                if self.is_mask:
                    connection_mask = connection.copy(self.name+"_mask"+to_add)
                    mask_starter = self.cable_starter('mask')
                    masks.append(connection_mask.sweep_along_path(mask_starter))
        
                tracks.append(connection.sweep_along_path(track_starter))
                gaps.append(connection_gap.sweep_along_path(gap_starter))
    
                
            if is_bond:
                self.draw_bond((self.inTrack+self.inGap*2)*1.5)
            
            if track_adaptor is not None:
                self.trackObjects.pop()
                self.gapObjects.pop()
                tracks = [*tracks, track_adaptor]
                gaps = [*gaps, gap_adaptor]
                if self.is_mask:
                    self.maskObjects.pop()
                    masks = [*masks, mask_adaptor]
                    
            if len(tracks)>1:
                names = [self.name+'_track', self.name+'_gap', self.name+'_mask']
    #            if track_adaptor is not None:
    #                names = [self.name+'_track_1', self.name+'_gap_1', self.name+'_mask_1']
                track = self.unite(tracks, names[0])
                gap = self.unite(gaps, names[1])
                if layer is None:
                    self.trackObjects.append(track)
                    self.gapObjects.append(gap)
                else:
                    self.layers[layer]['trackObjects'].append(track)
                    self.layers[layer]['gapObjects'].append(gap)
                if self.is_mask:
                    mask = self.unite(masks, names[2])
                    self.maskObjects.append(mask)
            else:
                track = tracks[0]
                gap = gaps[0]
                if layer is None:
                    self.trackObjects.append(track)
                    self.gapObjects.append(gap)
                else:
                    self.layers[layer]['trackObjects'].append(track)
                    self.layers[layer]['gapObjects'].append(gap)
                if self.is_mask:
                    self.maskObjects.append(*masks)
            
            if is_mesh is True:
                if not self.is_litho:
                    self.modeler.assign_mesh_length(track,2*self.inTrack)
                    
            for length in cable_length:
                print('{0}_length = {1:.3f} mm'.format(self.name, length*1000))
            print('sum = %.3f mm'%(1000*np.sum(cable_length)))
            
@Lib.add_methods_from(KeyElement, CustomElement)
class Body(PythonMdlr):

    pos_elt = [0,0]
    ori_elt = 0

    
    def __init__(self, interface, coor_sys, name, mode, variables):
        self.interface = interface
        self.coor_sys = coor_sys
        self.name = name
        self.maskObjects = []
        self.trackObjects = []
        self.gapObjects = []
        self.mode = mode # 'hfss' or 'gds'
        self.variables = variables
        
    def new_connector(self):
        from ConnectElement2 import ConnectElt2
        self.connector = ConnectElt2(self)


    
    
