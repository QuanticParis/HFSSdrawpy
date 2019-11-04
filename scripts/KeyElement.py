# -*- coding: utf-8 -*-
"""
Created on Mon Oct 28 16:18:36 2019

@author: Zaki
"""
from designer import Vector, way, equal_float, eps
import numpy as np
from hfss import parse_entry
from hfss import VariableString
from CustomElement import CustomElt


TOP = [0, 1]
DOWN = [0, -1]
RIGHT = [1, 0]
LEFT = [-1, 0]

POS = 0
ORI = 1
TRACK = 2
GAP = 3

class Port():
    def __init__(self, pos, ori, track, gap):
        self.pos = pos
        self.ori = ori
        self.track = track
        self.gap = gap
        
    # TODO several track within one gap for a port
    # TODO Be able to split a port which returns several ports
    # TODO Draw a port
    # TODO Adaptor should be defined here
    

class KeyElt(CustomElt):

    pcb_track = parse_entry('300um')
    pcb_gap = parse_entry('200um')
    is_mask = False
    gap_mask = parse_entry('20um')
    overdev = parse_entry('0um')
    is_overdev = False
    is_litho = False

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

    def move(func):
        def moved(*args, **kwargs):
            pos = args[2]
            angle = args[3]
            ports, entities = func(*(args[:2]+args[4:]), **kwargs)
            args[0].rotate(entities, angle=angle)
            args[0].translate(entities, vector=[pos[0], pos[1], 0])
            return ports, entities
        return moved
    
    def create_port(self, iTrack=0, iGap=0):
        iTrack, iGap = parse_entry((iTrack, iGap))
        portOut = [self.coor([0,0]), self.coor_vec([1,0]), iTrack+2*self.overdev, iGap-2*self.overdev]
        self.ports[self.name] = portOut
    
    def create_dc_port(self, layer, cut, rel_pos, wid):
        portOut = [self.coor([0,0]), self.coor_vec([1,0]), cut, rel_pos, wid, len(rel_pos)]
        self.ports_dc[layer+'_'+self.name] = portOut
    
#    draw_connector(name, pos, ori, iTrack, iGap, iBondLength, iSlope=1, pcbTrack=None, pcbGap=None, tr_line=True)
        
    # a decorator wraps the function, pos and ori are not arguments at the
    # the definition level, but should be used when using the function
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

        portOut = Port([adaptDist+self.pcb_gap+iBondLength,0], [1,0], iTrack+2*self.overdev, iGap-2*self.overdev)
#        print(self.pos, self.ori)
#        print(adaptDist)
#        print(self.pos+self.ori*(adaptDist+iGap+iBondLength), self.ori)
        points = [(self.pcb_gap-self.overdev, self.pcb_track/2+self.overdev),
                  (self.pcb_gap+iBondLength, self.pcb_track/2+self.overdev, 0),
                  (self.pcb_gap+iBondLength+adaptDist, self.overdev+iTrack/2),
                  (self.pcb_gap+iBondLength+adaptDist, -iTrack/2-self.overdev),
                  (self.pcb_gap+iBondLength, -self.pcb_track/2-self.overdev),
                  (self.pcb_gap-self.overdev, -self.pcb_track/2-self.overdev)]
        track = self.polyline_2D(points, name=name+'_track')
        
#        self.trackObjects.append(self.draw(points, ))
       
        points = [(self.pcb_gap/2+self.overdev, self.pcb_gap+self.pcb_track/2-self.overdev),
                 (self.pcb_gap+iBondLength, self.pcb_gap+self.pcb_track/2-self.overdev),
                 (self.pcb_gap+iBondLength+adaptDist, iGap+iTrack/2-self.overdev),
                 (self.pcb_gap+iBondLength+adaptDist, -iGap-iTrack/2+self.overdev),
                 (self.pcb_gap+iBondLength, -self.pcb_gap-self.pcb_track/2+self.overdev),
                 (self.pcb_gap/2+self.overdev, -self.pcb_gap-self.pcb_track/2+self.overdev)]
    
        gap = self.polyline_2D(points, name=name+'_gap')
#        self.gapObjects.append(self.draw(self.name+"_gap", points))

        if self.is_mask:
            points = self.append_points([(self.pcb_gap/2-self.gap_mask, self.pcb_gap+self.pcb_track/2+self.gap_mask),
                             (self.pcb_gap/2+iBondLength+self.gap_mask, 0),
                             (adaptDist, (iGap-self.pcb_gap)+(iTrack-self.pcb_track)*0.5),
                             (0, -2*iGap-iTrack-2*self.gap_mask),
                             (-adaptDist, (iGap-self.pcb_gap)+(iTrack-self.pcb_track)*0.5),
                             (-(self.pcb_gap/2+iBondLength)-self.gap_mask, 0)])
            self.maskObjects.append(self.draw(self.name+"_mask", points))


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
            
        return [portOut], [track, gap]


    def draw_quarter_circle(self, name, fillet, coor, ori=Vector([1,1])):
        ori=Vector(ori)
        temp = self.draw_rect(self.name+"_"+name, self.coor(coor), self.coor_vec(ori*2*fillet))
        temp_fillet = self.draw_rect(self.name+"_"+name+'f', self.coor(coor), self.coor_vec(ori*2*fillet))
        temp_fillet.fillet(fillet, 0)
        
        quarter = self.subtract(temp, [temp_fillet])
        return quarter
    
    def mesh_zone(self, zone_size, mesh_length):
        zone_size, mesh_length = parse_entry((zone_size, mesh_length))
        if not self.is_litho:
            zone = self.draw_rect_center(self.name, self.coor([0, 0]), self.coor_vec(zone_size))
            self.modeler.assign_mesh_length(zone, mesh_length)
            
    def cutout(self, zone_size):
        zone_size = parse_entry(zone_size)
        self.maskObjects.append(self.draw_rect_center(self.name, self.coor([0, 0]), self.coor_vec(zone_size)))
                    
    def draw_T(self, iTrack, iGap):
        
        if not self.is_overdev or self.val(self.overdev<0):
            cutout = self.draw_rect_center(self.name+'_cutout', self.coor([0,self.overdev/2]), self.coor_vec([2*iGap+iTrack, 2*iGap+iTrack-self.overdev]))
            self.gapObjects.append(cutout)
        else:
            points = self.append_points([(-(iGap+iTrack/2),-iTrack/2-iGap+self.overdev),
                             (0, 2*iGap+iTrack-2*self.overdev),
                             ((iGap+iTrack/2)*2, 0),
                             (0, -(2*iGap+iTrack)+2*self.overdev),
                             (-self.overdev, 0),
                             (0, -self.overdev), 
                             (-iTrack-2*iGap+2*self.overdev, 0),
                             (0, self.overdev)])
            cutout = self.draw(self.name+'_cutout', points)
            self.gapObjects.append(cutout)
        
        if self.is_mask:
            mask = self.draw_rect(self.name+'_mask', self.coor([-iGap-iTrack/2,-iGap-iTrack/2]), self.coor_vec([2*iGap+iTrack, 2*iGap+iTrack+self.gap_mask]))
            self.maskObjects.append(mask)
            
        points = self.append_points([(-(iGap+iTrack/2),-iTrack/2-self.overdev),
                                     (0, iTrack+2*self.overdev),
                                     ((iGap+iTrack/2)*2, 0),
                                     (0, -iTrack-2*self.overdev),
                                     (-iGap+self.overdev, 0),
                                     (0, -iGap+self.overdev), 
                                     (-iTrack-2*self.overdev, 0),
                                     (0, iGap-self.overdev)])
        track = self.draw(self.name+'_track', points)
        if self.val(iGap)<self.val(iTrack):
            fillet=iGap
        else:
            fillet=iTrack
        track.fillet(fillet-eps,[4,7])
        
        self.trackObjects.append(track)
        
        portOut1 = [self.coor([iTrack/2+iGap, 0]), self.coor_vec([1,0]), iTrack+2*self.overdev, iGap-2*self.overdev]
        self.ports[self.name+'_1'] = portOut1
        portOut2 = [self.coor([-(iTrack/2+iGap), 0]), self.coor_vec([-1,0]), iTrack+2*self.overdev, iGap-2*self.overdev]
        self.ports[self.name+'_2'] = portOut2
        portOut3 = [self.coor([0, -(iTrack/2+iGap)]), self.coor_vec([0,-1]), iTrack+2*self.overdev, iGap-2*self.overdev]
        self.ports[self.name+'_3'] = portOut3

    def draw_end_cable(self, iTrack, iGap, typeEnd = 'open', fillet=None):
        iTrack, iGap = parse_entry((iTrack, iGap))
        if typeEnd=='open' or typeEnd=='Open':
            cutout = self.draw_rect(self.name+'_cutout', self.coor([iGap,-(iTrack+2*iGap)/2+self.overdev]), self.coor_vec([-iGap+self.overdev, iTrack+2*iGap-2*self.overdev]))
            if fillet is not None:
                if abs(self.ori[0])==1:
                    cutout.fillet(iGap-self.overdev-eps,[2,1])
                else:
                    cutout.fillet(iGap-self.overdev-eps,[2,3])
            self.gapObjects.append(cutout)
            portOut = [self.coor([iGap, 0]), self.ori, iTrack+2*self.overdev, iGap-2*self.overdev]
            
            if self.is_overdev:
                track = self.draw_rect(self.name+'_track', self.coor([iGap,-iTrack/2-self.overdev]), self.coor_vec([-self.overdev, iTrack+2*self.overdev]))
                self.trackObjects.append(track)
            if self.is_mask:
                mask = self.draw_rect(self.name+'_mask', self.coor([-self.gap_mask,-(iTrack+2*iGap)/2-self.gap_mask]), self.coor_vec([iGap+self.gap_mask, iTrack+2*iGap+2*self.gap_mask]))
                if fillet is not None:
                    if abs(self.ori[0])==1:
                        mask.fillet(iGap+self.gap_mask-eps,[0,3])
                    else:
                        mask.fillet(iGap+self.gap_mask-eps,[0,1])
                self.maskObjects.append(mask)
                
        elif typeEnd=='short' or typeEnd=='Short':
            cutout1 = self.draw_rect(self.name+'_cutout1', self.coor([iGap/2,-(iTrack/2+iGap)+self.overdev]), self.coor_vec([-iGap/2+self.overdev, iGap-2*self.overdev]))
            cutout2 = self.draw_rect(self.name+'_cutout2', self.coor([iGap/2,(iTrack/2+iGap)-self.overdev]), self.coor_vec([-iGap/2+self.overdev, -iGap+2*self.overdev]))
            if fillet is not None:
                if abs(self.ori[0])==1:
                    cutout1.fillet(iGap/2-self.overdev-eps,[2,1])
                    cutout2.fillet(iGap/2-self.overdev-eps,[2,1])
                else:
                    cutout1.fillet(iGap/2-self.overdev-eps,[2,3])
                    cutout2.fillet(iGap/2-self.overdev-eps,[2,3])
            cutout = self.unite([cutout1, cutout2], self.name+'_cutout')
            self.gapObjects.append(cutout)
            portOut = [self.coor([iGap/2, 0]), self.ori, iTrack+2*self.overdev, iGap-2*self.overdev]
            
            if self.is_mask:
                mask = self.draw_rect(self.name+'_mask', self.coor([-self.gap_mask,-(iTrack+2*iGap)/2-self.gap_mask]), self.coor_vec([iGap/2+self.gap_mask, iTrack+2*iGap+2*self.gap_mask]))
                if fillet is not None:
                    if abs(self.ori[0])==1:
                        mask.fillet(iGap/2+self.gap_mask-eps,[0,3])
                    else:
                        mask.fillet(iGap/2+self.gap_mask-eps,[0,1])
                self.maskObjects.append(mask)
        else:
            raise ValueError("typeEnd should be 'open' or 'short', given %s" % typeEnd)
        self.ports[self.name] = portOut
        
    def size_dc_gap(self, length, positions, widths, border):
        
        length, positions, widths, border = parse_entry((length, positions, widths, border))
        
        low = np.argmin(positions)
        high = np.argmax(positions)
        pos_cutout = self.coor([-length/2, positions[low]-widths[low]/2-border])
        width_cutout = positions[high]-positions[low]+widths[high]/2+widths[low]/2+2*border
        vec_cutout = self.coor_vec([length, width_cutout])
      
        return pos_cutout, width_cutout, vec_cutout
