# -*- coding: utf-8 -*-
"""
Created on Mon Oct 28 16:18:36 2019

@author: Zaki
"""
#from designer import Vector, eps
#import numpy as np
#from hfss import parse_entry
#from CustomElement import CustomElt
#from CustomElement import Port
#
import Lib
from functools import wraps
from hfss import parse_entry
from Vector import Vector
import numpy as np
eps = 1e-7

# _DataStore

__methods__ = [] # self is a DataStore
register_method = Lib.register_method(__methods__)

#@register_method
#def big_method(self):
#    return self._b

#
#@register_method
#def huge_method(self):
#    return self._c


#class KeyElt(CustomElt):


#    @property
#    def pcb_track(self):
#        return self._pcb_track
#
#    @pcb_track.setter
#    def pcb_track(self, new_pcb_track):
#        self._pcb_track = parse_entry(new_pcb_track)
#
#    @property
#    def pcb_gap(self):
#        return self._pcb_gap
#
#    @pcb_gap.setter
#    def pcb_gap(self, new_pcb_gap):
#        self._pcb_gap = parse_entry(new_pcb_gap)


@register_method
def move(func):
    @wraps(func)
    # At the begining of the execution, we decide that all elements created go to instances_to_move
    def moved(*args, **kwargs):
        '''No movement along z axis'''
        #1 We need to keep track of the entities created during the execution of a function
#        nb_ports = args[0].number_ports()
#        nb_modelentities = args[0].number_modelentities()
        list_entities_old =args[0].modelentities_to_move()
        list_ports_old = args[0].ports_to_move()
        args[0].append_lists()      
        #2 We store the current position on the chip where the piece is to be drawn. 
        #2 It is modified during the execution of innner functions, which forces us to store them BEFORE the execution
        pos = args[0].current_pos
        angle = args[0].current_ori
        
        #3 We execute the actual function
        result = func(*args, **kwargs)

        #4 To avoid moving two times the same object, we splitted ModelEntity.instances in 'to_move' and 'moved'
        #4 But it does NOT mean that every 'to_move' object should be moved. We only move entities that were created during the execution
        list_entities_new = args[0].modelentities_to_move()
        list_ports_new = args[0].ports_to_move()
        print([i.name for i in list_entities_new])

        #5 We move the entities_to_move with the right operation
        if len(list_entities_new)>0:
            #TODO Move this in rotate and translate
            args[0].rotate(list_entities_new, angle=angle)
            args[0].translate(list_entities_new, vector=[pos[0], pos[1], 0])
            
        if len(list_ports_new)>0:            
            args[0].network.rotate(list_ports_new, angle)
            args[0].network.translate(list_ports_new, vector=[pos[0], pos[1], 0])
           
        #3 We empty a part of the 'to_move' dictionnaries in 'moved' dictionnaries
        if isinstance(list_entities_old[-1], list):
            a = list_entities_old.pop(-1)
            for entity in a:
                list_entities_old.append(entity)
        else:
            args[0].reset()
            
        if isinstance(list_ports_old[-1], list):
            a = list_ports_old.pop(-1)
            for entity in a:
                list_ports_old.append(entity)
        else:
            args[0].reset()
        
        return result
    return moved

@register_method
def create_port(self, name, iTrack=0, iGap=0):
    iTrack, iGap = parse_entry((iTrack, iGap))
    portOut = self.port(name, [0,0], [1,0], iTrack+2*self.overdev, iGap-2*self.overdev)
    return portOut
@register_method
def create_dc_port(self, name, layer, cut, rel_pos, wid):
    portOut = [layer+'_'+ name, [0,0], [1,0], cut, rel_pos, wid, len(rel_pos)]
    return portOut

#    draw_connector(name, pos, ori, iTrack, iGap, iBondLength, iSlope=1, pcbTrack=None, pcbGap=None, tr_line=True)
    
# a decorator wraps the function, pos and ori are not arguments at the
# the definition level, but should be used when using the function
@register_method
@move
def draw_connector(self, name, iTrack, iGap, iBondLength, iSlope=1, pcbTrack=None, pcbGap=None, tr_line=True):
    '''
    Draws a CPW connector for inputs and outputs.

    Inputs:
    -------

    iBondLength: (float) corresponds to dimension a in the drawing
    iSlope (float): 1 is 45 degrees to adpat lengths between Bonding pad and iout
    iLineTest (Bool): unclear, keep False

        ground plane
        +------+
        |       \
        |        \
        |   +-a-+\|
    iIn |   |.....| iOut
        |   +---+/|
        |        /
        |       /
        +------+

    Outputs:
    --------
    returns iIn and recalculated iOut
    '''
    
    iTrack, iGap = parse_entry((iTrack, iGap))
    iBondLength, iSlope = parse_entry((iBondLength, iSlope))
    
    if pcbGap is not None:
        pcbGap = parse_entry(pcbGap)
        self.pcb_gap = pcbGap
        
    if pcbTrack is not None:
        pcbTrack = parse_entry(pcbTrack)
        self.pcb_track = pcbTrack
        
    adaptDist = (self.pcb_track/2-iTrack/2)/iSlope

    portOut = self.port('iOut', [adaptDist+self.pcb_gap+iBondLength,0], [1,0], iTrack+2*self.overdev, iGap-2*self.overdev)
#        print(self.pos, self.ori)
#        print(adaptDist)
#        print(self.pos+self.ori*(adaptDist+iGap+iBondLength), self.ori)
    points = [(self.pcb_gap-self.overdev, self.pcb_track/2+self.overdev),
              (self.pcb_gap+iBondLength, self.pcb_track/2+self.overdev),
              (self.pcb_gap+iBondLength+adaptDist, self.overdev+iTrack/2),
              (self.pcb_gap+iBondLength+adaptDist, -iTrack/2-self.overdev),
              (self.pcb_gap+iBondLength, -self.pcb_track/2-self.overdev),
              (self.pcb_gap-self.overdev, -self.pcb_track/2-self.overdev)]

    track = self.polyline_2D(points, name=name+'_track', layer='TRACK')
#        self.trackObjects.append(self.draw(points, ))
   
    points = [(self.pcb_gap/2+self.overdev, self.pcb_gap+self.pcb_track/2-self.overdev),
             (self.pcb_gap+iBondLength, self.pcb_gap+self.pcb_track/2-self.overdev),
             (self.pcb_gap+iBondLength+adaptDist, iGap+iTrack/2-self.overdev),
             (self.pcb_gap+iBondLength+adaptDist, -iGap-iTrack/2+self.overdev),
             (self.pcb_gap+iBondLength, -self.pcb_gap-self.pcb_track/2+self.overdev),
             (self.pcb_gap/2+self.overdev, -self.pcb_gap-self.pcb_track/2+self.overdev)]

    gap = self.polyline_2D(points, name=name+'_gap', layer='GAP')
#        self.gapObjects.append(self.draw(self.name+"_gap", points))

    if self.is_mask:
        points =[(self.pcb_gap/2-self.gap_mask, self.pcb_gap+self.pcb_track/2+self.gap_mask),
                 (self.pcb_gap+iBondLength, self.pcb_gap+self.pcb_track/2+self.gap_mask),
                 (self.pcb_gap+iBondLength+adaptDist, iGap+iTrack/2+self.gap_mask),
                 (self.pcb_gap+iBondLength+adaptDist, -iGap-self.gap_mask),
                 (self.pcb_gap+iBondLength, (self.pcb_gap)+(iTrack-self.pcb_track)*0.5-self.gap_mask),
                 (self.pcb_gap/2-self.gap_mask, (self.pcb_gap)+(iTrack-self.pcb_track)*0.5-self.gap_mask)]
        
        mask = self.polyline_2D(points, name=self.name+"_mask", layer="MASK")

        self.maskObjects.append(mask)  
        
#        if not self.is_litho and tr_line:
#            points = self.append_points([(self.pcb_gap/2+self.overdev, self.pcb_track/2+self.overdev),
#                                         (self.pcb_gap/2-2*self.overdev, 0),
#                                         (0, -self.pcb_track-2*self.overdev),
#                                         (-self.pcb_gap/2+2*self.overdev, 0)])
#            ohm = self.draw(points, name=name+'_ohm')
#            self.assign_lumped_RLC(ohm, self.ori, ('50ohm',0,0))
#            self.modeler.assign_mesh_length(ohm, self.pcb_track/10)
##            self.trackObjects.append(ohm)
#            points = self.append_points([(self.pcb_gap/2+self.overdev,0),(self.pcb_gap/2-2*self.overdev,0)])
#            self.draw(self.name+'_line', points, closed=False)
        
    return [portOut], [track, gap, mask]
    #on renvoie track, gap, mask
    
@register_method
@move
def draw_quarter_circle(self, name, layer, fillet):
    # Usual case which is then rotated with ori and translated with pos
    temp = self.rect_corner_2D([0,0], 2*fillet*Vector([1,1]), name=name ,layer=layer)
    temp_fillet = self.rect_corner_2D([0,0], 2*(fillet)*Vector([1,1]), name=name+'_fillet', layer=layer)
    #Minus eps because hfss can't substract object of same size
    self._fillet(fillet-eps, [0], temp_fillet)
    self.subtract(temp, [temp_fillet])
    return [[],[temp]]

@register_method
@move        
def cutout(self, name, zone_size):
    zone_size = parse_entry(zone_size)
    cutout_rect = self.rect_center_2D([0, 0], zone_size, name=name, layer='mask')
    #self.maskObjects.append(cutout_rect)
    return [[], [cutout_rect]]


@register_method
@move
def draw_T(self, name, iTrack, iGap):
    if not self.is_overdev or self.val(self.overdev<0):
        cutout = self.rect_center_2D([0,self.overdev/2], [2*iGap+iTrack, 2*iGap+iTrack-self.overdev], name=self.name+'_cutout', layer='gap')
        self.gapObjects.append(cutout)
    else:
        points = [       (-(iGap+iTrack/2),-iTrack/2-iGap+self.overdev),
                         (-(iGap+iTrack/2), iGap-self.overdev),
                         ((iGap+iTrack/2), iGap-self.overdev),
                         ((iGap+iTrack/2), -(iGap+iTrack)+self.overdev),
                         ((iGap+iTrack/2)-self.overdev, -(iGap+iTrack)+self.overdev),
                         ((iGap+iTrack/2)-self.overdev, -(iGap+iTrack)), 
                         (-iTrack/2-iGap+self.overdev, -(iGap+iTrack)),
                         (-iTrack/2-iGap+self.overdev, self.overdev-(iGap+iTrack)) ]
        
        cutout = self.polyline_2D(points, name=name+'_cutout', layer='gap')
        self.gapObjects.append(cutout)
    
    if self.is_mask:
        mask = self.rect_center_2D([-iGap-iTrack/2,-iGap-iTrack/2], [2*iGap+iTrack, 2*iGap+iTrack+self.gap_mask], name = self.name+'_mask', layer='mask')
        self.maskObjects.append(mask)
        
    points = [                   (-(iGap+iTrack/2),-iTrack/2-self.overdev),
                                 (-(iGap+iTrack/2), iTrack/2+self.overdev),
                                 ((iGap+iTrack/2), iTrack/2+self.overdev),
                                 ((iGap+iTrack/2), -iTrack/2-self.overdev),
                                 (iTrack/2+self.overdev, -iTrack/2-self.overdev),
                                 (iTrack/2+self.overdev, -iGap-iTrack/2), 
                                 (-iTrack/2-self.overdev, -iGap-iTrack/2),
                                 (-iTrack/2-self.overdev, -iTrack/2-self.overdev)]
    track = self.polyline_2D(points, name = self.name+'_track', layer='track')
    if iGap<iTrack:
        #not sure if self.val(iGap)<self.val(iTrack) was correct
        fillet=iGap
    else:
        fillet=iTrack
    self._fillet(fillet-eps,[4,7], track)
    
    self.trackObjects.append(track)
    
    portOut1 = self.port('portOut1', [iTrack/2+iGap, 0], [1,0], iTrack+2*self.overdev, iGap-2*self.overdev)
    #ports[self.name+'_1'] = portOut1
    portOut2 = self.port('portOut2', [-(iTrack/2+iGap), 0], [-1,0], iTrack+2*self.overdev, iGap-2*self.overdev)
    #ports[self.name+'_2'] = portOut2
    portOut3 = self.port('portOut2', [0, -(iTrack/2+iGap)], [0,-1], iTrack+2*self.overdev, iGap-2*self.overdev)
    #ports[self.name+'_3'] = portOut3
    
    return [[portOut1,portOut2,portOut3], [track, cutout]]

@register_method
@move
def draw_end_cable(self, name, iTrack, iGap, typeEnd = 'open', fillet=None):
    iTrack, iGap = parse_entry((iTrack, iGap))
    track, mask = None, None # track is not created in every cases

    if typeEnd=='open' or typeEnd=='Open':
        cutout = self.rect_center_2D([iGap,-(iTrack+2*iGap)/2+self.overdev], [-iGap+self.overdev, iTrack+2*iGap-2*self.overdev], name=name+'_cutout', layer='gap')
        if fillet is not None:
            cutout.fillet(iGap-self.overdev-eps,[2,1])

        self.gapObjects.append(cutout)
        portOut = self.port('portOut'+name, [iGap, 0], [1,0], iTrack+2*self.overdev, iGap-2*self.overdev)
        
        if self.is_overdev:
            track = self.rect_center_2D([iGap,-iTrack/2-self.overdev], [-self.overdev, iTrack+2*self.overdev], name=name+'_track' ,layer='track')
            self.trackObjects.append(track)
        if self.is_mask:
            mask = self.rect_center_2D([-self.gap_mask,-(iTrack+2*iGap)/2-self.gap_mask],[iGap+self.gap_mask, iTrack+2*iGap+2*self.gap_mask],name=name+'_mask', layer='mask')
            if fillet is not None:
                mask.fillet(iGap+self.gap_mask-eps,[0,3])

            self.maskObjects.append(mask)
            
    elif typeEnd=='short' or typeEnd=='Short':
        cutout1 = self.rect_center_2D([iGap/2,-(iTrack/2+iGap)+self.overdev], [-iGap/2+self.overdev, iGap-2*self.overdev], name=name+'_cutout1', layer='gap')
        cutout2 = self.rect_center_2D([iGap/2,(iTrack/2+iGap)-self.overdev], [-iGap/2+self.overdev, -iGap+2*self.overdev], name=name+'_cutout2', layer='gap')
        if fillet is not None:
            cutout1.fillet(iGap/2-self.overdev-eps,[2,1])
            cutout2.fillet(iGap/2-self.overdev-eps,[2,1])

        cutout = self.unite([cutout1, cutout2], self.name+'_cutout')
        self.gapObjects.append(cutout)
        portOut = self.port([iGap/2, 0], [1,0], iTrack+2*self.overdev, iGap-2*self.overdev)
        
        if self.is_mask:
            mask = self.rect_center_2D([-self.gap_mask,-(iTrack+2*iGap)/2-self.gap_mask], [iGap/2+self.gap_mask, iTrack+2*iGap+2*self.gap_mask], name=name+'_mask', layer='mask')
            if fillet is not None:
                mask.fillet(iGap/2+self.gap_mask-eps,[0,3])

            self.maskObjects.append(mask)

    else:
        raise ValueError("typeEnd should be 'open' or 'short', given %s" % typeEnd)
        
    return [[portOut], [track, cutout, mask]]



@register_method
def size_dc_gap(self, length, positions, widths, border):
    # Ne pas s'en préoccuper
    length, positions, widths, border = parse_entry((length, positions, widths, border))
    
    low = np.argmin(positions)
    high = np.argmax(positions)
    pos_cutout = [-length/2, positions[low]-widths[low]/2-border]
    width_cutout = positions[high]-positions[low]+widths[high]/2+widths[low]/2+2*border
    vec_cutout = [length, width_cutout]
  
    return pos_cutout, width_cutout, vec_cutout

@register_method
@move
def cavity_3D(self, name, inner_radius, outer_radius, cylinder_height, antenna_radius, antenna_height):
    inner_cylinder = self.cylinder([0,0,0], inner_radius , cylinder_height, "Z", name=name+"_inner_cylinder", layer="l1")
    outer_cylinder = self.cylinder([0,0,0], outer_radius , cylinder_height, "Z", name=name+"_outer_cylinder", layer="l1")
    tube = self.subtract(outer_cylinder, [inner_cylinder])
    antenna = self.cylinder([0,0,0], antenna_radius ,antenna_height, "Z", name=self.name+"_antenna", layer="l1")
    return [[],[tube, antenna]] 

@register_method
@move
def cavity_3D_simple(self, name, radius, cylinder_height, antenna_radius, antenna_height, insert_radius, depth):
    radius, cylinder_height, antenna_radius, antenna_height, insert_radius, depth = parse_entry((radius, cylinder_height, antenna_radius, antenna_height, insert_radius, depth))
    cylinder = self.cylinder([0,0,0], radius , cylinder_height, "Z", name=name+"_inner_cylinder", layer="l1")
    self.assign_perfect_E(cylinder, name+'_PerfE')
    antenna = self.cylinder([0,0,0], antenna_radius ,antenna_height, "Z", name=name+"_antenna", layer="l1")
    self.make_material(antenna, "\"aluminum\"")
    cylinder_side = self.cylinder([0,radius-depth,antenna_height], insert_radius, radius/2, "Y", name=name+"_insert", layer="l1")
    self.unite([cylinder, cylinder_side], name="union_cylinder")

    return [[],[cylinder, antenna]] 
