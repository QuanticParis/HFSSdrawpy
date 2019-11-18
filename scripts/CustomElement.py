# -*- coding: utf-8 -*-
"""
Created on Mon Oct 28 16:27:24 2019

@author: Zaki
"""

from designer import Vector, eps
import numpy as np
from hfss import parse_entry
#from PythonModeler import PythonModeler

from functools import wraps
import Lib

__methods__ = [] # self is a DataStore
register_method = Lib.register_method(__methods__)

@register_method
def draw_JJ(self, name, iTrack, iGap, iTrackJ, iLength, iInduct='1nH', fillet=None):
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
    iTrack, iGap, iTrackJ, iLength = parse_entry((iTrack, iGap, iTrackJ, iLength))

    portOut1 = self.port(name+'_1',[iLength/2,0], [1,0], iTrack, iGap)
    portOut2 = self.port(name+'_2',[-iLength/2,0], [-1,0], iTrack, iGap)

    junction = self.connect_elt(name, portOut1, portOut2)
    pads = junction._connect_JJ(iTrackJ, iInduct=iInduct, fillet=None)
    self.trackObjects.append(pads)
    
    self.gapObjects.append(self.draw_rect_center(self.name+"_cutout", self.coor([0,0]), self.coor_vec([iLength, iTrack+2*iGap])))

def draw_IBM_tansmon(self,
                   cutout_size,
                   pad_spacing,
                   pad_size,
                   Jwidth,
                   track,
                   gap,
                   Jinduc,
                   nport=1,
                   fillet=None):
    

    

    cutout_size, pad_spacing, pad_size, Jwidth, track, gap = parse_entry((cutout_size, pad_spacing, pad_size, Jwidth, track, gap))
    cutout_size = Vector(cutout_size)
    pad_size = Vector(pad_size)
    
    cutout = self.rect_center_2D([0,0], cutout_size-Vector([2*self.overdev, 2*self.overdev]), name=self.name+"_cutout", layer="GAP")

    if self.is_mask:
        mask = self.rect_center_2D([0,0], cutout_size+Vector([track,track]),name=self.name+"_mask",layer="MASK")
        self.maskObjects.append(mask)
    if not self.is_litho:
        mesh = self.rect_center_2D([0,0], cutout_size, name=self.name+"_mesh", layer="MESH")
        self.mesh_zone(mesh, 2*track)

    track_J=Jwidth*4.
    in_junction = self.port(self.name+'_in_jct', [-pad_spacing/2+self.overdev, 0], [1,0], track_J+2*self.overdev, 0)
    out_junction = self.port(self.name+'_out_jct', [pad_spacing/2-self.overdev, 0], [1,0], track_J+2*self.overdev, 0)
#        self.ports[self.name+'_in_jct'] = in_junction
#        self.ports[self.name+'_out_jct'] = out_junction
#        junction = self.connect_elt(self.name+'_junction', self.name+'_in_jct', self.name+'_out_jct')
    junction_pad = self.connector.f('JJ', self.name+'_in_jct', self.name+'_out_jct', Jwidth+2*self.overdev, iInduct=Jinduc, fillet=None)

    right_pad = self.rect_center_2D(Vector(pad_spacing+pad_size[0],0)/2, pad_size+Vector([2*self.overdev, 2*self.overdev]), name=self.name+"_pad1", layer="TRACK")
    left_pad = self.rect_center_2D(-Vector(pad_spacing+pad_size[0],0)/2, pad_size+Vector([2*self.overdev, 2*self.overdev]), name=self.name+"_pad2", layer="TRACK")
    

    pads = self.unite([right_pad, left_pad, junction_pads], name=self.name+'_pads')
    
    if fillet is not None:
        self._fillet(fillet+self.overdev,[19], pads)
        self._fillet(0.75*fillet-self.overdev,[18,17,14,13], pads)
        self._fillet(fillet+self.overdev,[12,11,10,9,8,7], pads)
        self._fillet(0.75*fillet-self.overdev,[6,5,2,1], pads)
        self._fillet(fillet+self.overdev,0, pads)
  
    self.trackObjects.append(pads)

    if nport==1:
#            self.ports[self.name+'_1'] = [self.coor(cutout_size.px()/2), self.ori, track, gap]
        self.port(self.name+'_1',cutout_size.px()/2, self.ori, track, gap)
    elif nport==2:
#            self.ports[self.name+'_1'] = [self.coor(cutout_size.px()/2), self.ori, track, gap]
#            self.ports[self.name+'_2'] = [self.coor(-cutout_size.px()/2), -self.ori, track, gap]
        self.port(self.name+'_1',cutout_size.px()/2, self.ori, track, gap)
        self.port(self.name+'_2',-cutout_size.px()/2, -self.ori, track, gap)
    elif nport==3:
#            self.ports[self.name+'_1a'] = [self.coor(cutout_size.px()/2+pad_size.py()/2), self.ori, track, gap]
#            self.ports[self.name+'_1b'] = [self.coor(cutout_size.px()/2-pad_size.py()/2), self.ori, track, gap]
#            self.ports[self.name+'_2'] = [self.coor(-cutout_size.px()/2), -self.ori, track, gap]
        self.port(self.name+'_1a',cutout_size.px()/2+pad_size.py()/2, self.ori, track, gap)
        self.port(self.name+'_1b',cutout_size.px()/2-pad_size.py()/2, self.ori, track, gap)

        self.port(self.name+'_2',-cutout_size.px()/2, -self.ori, track, gap)
    elif nport==4:
        pass
        # TODO
#            self.ports[self.name+'_1a'] = [self.coor(cutout_size.px()/2+pad_size.py()/2), self.ori, track, gap]
#            self.ports[self.name+'_1b'] = [self.coor(cutout_size.px()/2-pad_size.py()/2), self.ori, track, gap]
#            self.ports[self.name+'_2a'] = [self.coor(-cutout_size.px()/2+pad_size.py()/2), -self.ori, track, gap]
#            self.ports[self.name+'_2b'] = [self.coor(-cutout_size.px()/2-pad_size.py()/2), -self.ori, track, gap]
    elif nport==5:
        pass
#            self.ports[self.name+'_1'] = [self.coor(cutout_size.px()/2), self.ori, track+2*self.overdev, gap-2*self.overdev]
#            self.ports[self.name+'_2'] = [self.coor(-cutout_size.px()/2), -self.ori, track+2*self.overdev, gap-2*self.overdev]
#            self.ports[self.name+'_3a'] = [self.coor(Vector(pad_spacing/2+pad_size[0],-cutout_size[1]/2)), -self.ori.orth() ,track/2+2*self.overdev,gap/2-2*self.overdev]
#            if self.is_overdev:
#                sub_1 = self.draw_rect(self.name + '_sub_1', self.coor([cutout_size[0]/2, -track/2-gap+self.overdev]), self.coor_vec([-self.overdev, track+2*gap-2*self.overdev]))
#                sub_2 = self.draw_rect(self.name + '_sub_2', self.coor([-cutout_size[0]/2, -track/2-gap+self.overdev]), self.coor_vec([self.overdev, track+2*gap-2*self.overdev]))
#                sub_3a = self.draw_rect(self.name + '_sub_3a', self.coor([-track/2-gap+self.overdev+pad_spacing/2+pad_size[0],-cutout_size[1]/2]), self.coor_vec([track+2*gap-2*self.overdev, self.overdev]))

#        if self.is_overdev:
#            cutout = self.unite([cutout, sub_1, sub_2, sub_3a])
#            
    self.gapObjects.append(cutout)
    
#        self.draw_rect_center(self.name+"check1", self.coor(cutout_size.px()/2)+self.ori*pad_size[0]/20, self.rot(*pad_size)/10)
#        self.draw_rect_center(self.name+"check2", self.coor(-cutout_size.px()/2)-self.ori*pad_size[0]/20, self.rot(*pad_size)/10)
#        self.draw_rect_center(self.name+"check3", self.coor(Vector(pad_spacing/2+pad_size[0],cutout_size[1]/2))+self.ori.orth()*pad_size[1]/20, self.rot(*pad_size)/10)

@register_method
def draw_ZR_transmon(self,
                    cutout_size,
                    pad_spacing,
                    pad_size_right,
                    track_right,
                    gap_right,
                    length_right,
                    spacing_right,
                    short_right,
                    Jwidth,
                    Jinduc,
                    pad_size_left=None,
                    track_left=None,
                    gap_left=None,
                    length_left=None,
                    spacing_left=None,
                    short_left=None,
                    fillet=None):
    
    # Short should be 0 for no short

    parsed = parse_entry((cutout_size,
                          pad_spacing,
                          pad_size_right,
                          track_right,
                          gap_right,
                          length_right,
                          spacing_right,
                          short_right,
                          Jwidth,
                          Jinduc,
                          pad_size_left,
                          track_left,
                          gap_left,
                          length_left, 
                          spacing_left, 
                          short_left))
                    
    (cutout_size,
     pad_spacing,
     pad_size_right,
     track_right,
     gap_right,
     length_right,
     spacing_right,
     short_right,
     Jwidth,
     Jinduc,
     pad_size_left,
     track_left,
     gap_left,
     length_left,
     spacing_left,
     short_left) = parsed
     
    if pad_size_left is None:
        pad_size_left=pad_size_right
    if track_left is None:
        track_left=track_right
    if gap_left is None:
        gap_left=gap_right
    if length_left is None:
        length_left=length_right
    if spacing_left is None:
        spacing_left=spacing_right
    if short_left is None:
        short_left=short_right
    
    cutout_size = Vector(cutout_size)
    pad_size_left = Vector(pad_size_left)
    pad_size_right = Vector(pad_size_right)

    cutout = self.rect_center_2D([0,0], cutout_size-Vector([self.overdev*2, self.overdev*2]),name =self.name+"_cutout", layer = "GAP")
    if not self.is_litho:
        mesh = self.rect_center_2D([0,0], cutout_size, name=self.name+"_mesh", layer="MESH")
    if self.is_mask:
        mask = self.rect_center_2D([0,0], cutout_size+Vector([self.gap_mask,self.gap_mask])*2, name=self.name+"_mask", layer="MASK")
    
    track_J=Jwidth*4.
    in_junction = self.port(self.name+'_in_jct', [-pad_spacing/2+self.overdev, 0], [1,0], track_J+2*self.overdev, 0)
    out_junction = self.port(self.name+'_out_jct', [+pad_spacing/2-self.overdev, 0], [-1,0], track_J+2*self.overdev, 0)
    
#        self.ports[self.name+'_in_jct'] = in_junction
#        self.ports[self.name+'_out_jct'] = out_junction
#        junction = self.connect_elt(self.name+'_junction', self.name+'_in_jct', self.name+'_out_jct')
    junction_pads = self.network._connect_JJ('JJ', self.name+'_in_jct', self.name+'_out_jct', Jwidth+2*self.overdev, iInduct=Jinduc, fillet=None)
    
    raw_points = [(pad_spacing/2-self.overdev, -pad_size_right[1]/2-self.overdev),
                  (pad_size_right[0]+2*self.overdev, 0),
                  (0, pad_size_right[1]/2-spacing_right-short_right-gap_right-track_right/2+2*self.overdev),
                  (-pad_size_right[0]+length_right, 0),
                  (0, (spacing_right+short_right+gap_right+track_right/2)*2-2*self.overdev),
                  (pad_size_right[0]-length_right, 0),
                  (0, pad_size_right[1]/2-spacing_right-short_right-gap_right-track_right/2+2*self.overdev),
                  (-pad_size_right[0]-2*self.overdev, 0)]
    points = self.append_points(raw_points)
    right_pad = self.polyline_2D(points, name=self.name+"_pad1", layer="TRACK") 
        
    raw_points = [(-pad_spacing/2+self.overdev, -pad_size_left[1]/2-self.overdev),
                  (-pad_size_left[0]-2*self.overdev, 0),
                  (0, pad_size_left[1]/2-spacing_left-short_left-gap_left-track_left/2+2*self.overdev),
                  (pad_size_left[0]-length_left, 0),
                  (0, (spacing_left+short_left+gap_left+track_left/2)*2-2*self.overdev),
                  (-pad_size_left[0]+length_left, 0),
                  (0, pad_size_left[1]/2-spacing_left-short_left-gap_left-track_left/2+2*self.overdev),
                  (pad_size_left[0]+2*self.overdev, 0)]
    points = self.append_points(raw_points)
    left_pad = self.polyline_2D(points, name=self.name+"_pad2", layer="TRACK") 
    
    right_track_raw_points = [(cutout_size[0]/2-self.overdev,-track_right/2-self.overdev),
                             (-(cutout_size[0]/2-pad_spacing/2-length_right-gap_right-short_right-spacing_right), 0),
                             (0, track_right+2*self.overdev),
                             ((cutout_size[0]/2-pad_spacing/2-length_right-gap_right-short_right-spacing_right), 0)]
    right_track_points = self.append_points(right_track_raw_points)
    right_track = self.polyline_2D(right_track_points, name=self.name+"_track1", layer="TRACK") 
    
    left_track_raw_points = [(-cutout_size[0]/2+self.overdev,-track_left/2-self.overdev),
                             ((cutout_size[0]/2-pad_spacing/2-length_left-gap_left-short_left-spacing_left), 0),
                             (0, track_left+2*self.overdev),
                             (-(cutout_size[0]/2-pad_spacing/2-length_left-gap_left-short_left-spacing_left), 0)]
    left_track_points = self.append_points(left_track_raw_points)
    left_track = self.polyline_2D(left_track_points, name=self.name+"_track2", layer="TRACK") 
    
    if short_right!=0:
        raw_points = [(cutout_size[0]/2-self.overdev,-track_right/2-gap_right+self.overdev),
                      (-(cutout_size[0]/2-pad_spacing/2-length_right-short_right-spacing_right)+2*self.overdev, 0),
                      (0, 2*gap_right+track_right-2*self.overdev),
                      ((cutout_size[0]/2-pad_spacing/2-length_right-short_right-spacing_right)-2*self.overdev, 0),
                      (0, short_right+2*self.overdev),
                      (-(cutout_size[0]/2-pad_spacing/2-length_right-spacing_right), 0),
                      (0, -(2*gap_right+track_right+2*short_right)-2*self.overdev),
                      ((cutout_size[0]/2-pad_spacing/2-length_right-spacing_right), 0)]
        points = self.append_points(raw_points)
        right_short = self.polyline_2D(points, name=self.name+"_short1", layer="TRACK") 
        
    if short_left!=0:
        raw_points = [(-cutout_size[0]/2+self.overdev,-track_left/2-gap_left+self.overdev),
                      ((cutout_size[0]/2-pad_spacing/2-length_left-short_left-spacing_left)-2*self.overdev, 0),
                      (0, 2*gap_left+track_left-2*self.overdev),
                      (-(cutout_size[0]/2-pad_spacing/2-length_left-short_left-spacing_left)+2*self.overdev, 0),
                      (0, short_left+2*self.overdev),
                      ((cutout_size[0]/2-pad_spacing/2-length_left-spacing_left), 0),
                      (0, -(2*gap_left+track_left+2*short_left)-2*self.overdev),
                      (-(cutout_size[0]/2-pad_spacing/2-length_left-spacing_left), 0)]
        points = self.append_points(raw_points)
        left_short = self.polyline_2D(points, name=self.name+"_short2", layer="TRACK") 
        
    if fillet is not None:
        self._fillet(track_right/2-eps+self.overdev,[2,1], right_track)
        self._fillet(track_left/2-eps+self.overdev,[2,1], left_track)
        self._fillet(cutout_size[1]/6-self.overdev,[0,1,2,3], cutout)
        if not self.is_litho:
            self._fillet(cutout_size[1]/6,[0,1,2,3], mesh)
        
        if self.is_mask:
            self._fillet(cutout_size[1]/6+self.gap_mask,[0,1,2,3], mask)
        
        if short_right!=0:
            self._fillet(track_right/2+gap_right+short_right-eps+self.overdev,[6,5], right_short)
            self._fillet(track_right/2+gap_right-eps-self.overdev,[2,1], right_short)
            
            fillet_right1 = (pad_size_right[1]/2-spacing_right-short_right-gap_right-track_right/2)/4-self.overdev
            right_quarter_up1 = self.draw_quarter_circle('right_quarter_up1', [cutout_size[0]/2-self.overdev, track_right/2+gap_right+short_right+self.overdev], [1,0], 'TRACK', fillet_right1)
            right_quarter_down1 = self.draw_quarter_circle('right_quarter_down1', [cutout_size[0]/2-self.overdev, -(track_right/2+gap_right+short_right)-self.overdev], [0,-1], 'TRACK', fillet_right1)
            right_short = self.unite([right_short, right_quarter_up1, right_quarter_down1])
        else:
            fillet_right2 = (pad_size_right[1]/2-spacing_right-gap_right-track_right/2)/4+self.overdev
            right_quarter_up2 = self.draw_quarter_circle('right_quarter_up2', [cutout_size[0]/2-self.overdev, track_right/2+gap_right-self.overdev], [1,0], 'TRACK', fillet_right2)
            right_quarter_down2 = self.draw_quarter_circle('right_quarter_down2', [cutout_size[0]/2-self.overdev, -(track_right/2+gap_right)+self.overdev], [0,-1] ,'TRACK' , fillet_right2)
            cutout = self.unite([cutout, right_quarter_up2, right_quarter_down2])
            
        if short_left!=0:
            self._fillet(track_left/2+gap_left+short_left-eps+self.overdev,[6,5], left_short)
            self._fillet(track_left/2+gap_left-eps-self.overdev,[2,1], left_short)
            
            fillet_left1= (pad_size_left[1]/2-spacing_left-short_left-gap_left-track_left/2)/4-self.overdev
            left_quarter_up1 = self.draw_quarter_circle('left_quarter_up1', [-cutout_size[0]/2+self.overdev,track_left/2+gap_left+short_left+self.overdev], [-1,0], 'TRACK', fillet_left1)
            left_quarter_down1 = self.draw_quarter_circle('left_quarter_down1', [-cutout_size[0]/2+self.overdev,-(track_left/2+gap_left+short_left)-self.overdev], [1,0], 'TRACK', fillet_left1)
            left_short = self.unite([left_short, left_quarter_up1, left_quarter_down1])
        else:
            fillet_left2 = (pad_size_left[1]/2-spacing_left-gap_left-track_left/2)/4+self.overdev
            left_quarter_up2 = self.draw_quarter_circle('left_quarter_up2', [-cutout_size[0]/2+self.overdev,track_left/2+gap_left-self.overdev], [-1,0], 'TRACK', fillet_left2)
            left_quarter_down2 = self.draw_quarter_circle('left_quarter_down2', [-cutout_size[0]/2+self.overdev,-(track_left/2+gap_left)+self.overdev], [1,0], 'TRACK', fillet_left2)
            cutout = self.unite([cutout, left_quarter_up2, left_quarter_down2])
            
        self._fillet(pad_size_right[0]/4+self.overdev,[7,0], right_pad)
        self._fillet((pad_size_right[1]/2-spacing_right-short_right-gap_right-track_right/2)/4+self.overdev,[1,2,5,6], right_pad)
        self._fillet(track_right/2+gap_right+short_right+spacing_right-eps-self.overdev,[5,6], right_pad)

        self._fillet(pad_size_left[0]/4+self.overdev,[7,0], left_pad)
        self._fillet((pad_size_left[1]/2-spacing_left-short_left-gap_left-track_left/2)/4+self.overdev,[1,2,5,6], left_pad)
        self._fillet(track_left/2+gap_left+short_left+spacing_left-eps-self.overdev,[5,6], left_pad)
    
    if not self.is_litho:
        self.mesh_zone(mesh, 4*Jwidth)
    
    to_unite = [right_pad, left_pad, right_track, left_track, junction_pads]
    if short_right!=0:
        to_unite.append(right_short)
    if short_left!=0:
        to_unite.append(left_short)
    
    if self.is_overdev:
        added_gap_left = self.rect_center_2D([-cutout_size[0]/2, -(track_left+2*gap_left)/2+self.overdev], [self.overdev, (track_left+2*gap_left)-2*self.overdev], name=self.name+'_added_gap_left', layer='GAP')
        added_gap_right = self.rect_center_2D([cutout_size[0]/2, -(track_right+2*gap_right)/2+self.overdev], [-self.overdev, (track_right+2*gap_right)-2*self.overdev], name=self.name+'_added_gap_right', layer='GAP' )

        added_track_left = self.rect_center_2D([-cutout_size[0]/2, -(track_left)/2-self.overdev], [self.overdev, (track_left)+2*self.overdev], name=self.name+'_added_track_left', layer="TRACK")
        added_track_right = self.rect_center_2D([cutout_size[0]/2, -(track_right)/2-self.overdev], [-self.overdev, (track_right)+2*self.overdev], name=self.name+'_added_track_right', layer="TRACK")
        to_unite = to_unite+[added_track_left, added_track_right]
        
        cutout = self.unite([cutout, added_gap_left, added_gap_right])
        
    pads = self.unite(to_unite, name=self.name+'_pads')
    self.trackObjects.append(pads)
    
    self.gapObjects.append(cutout)
    
    if self.is_mask:
        self.maskObjects.append(mask)
    
    portOut1 = self.port('portOut1', [cutout_size[0]/2,0], [1,0], track_right+2*self.overdev, gap_right-2*self.overdev)
#        self.ports[self.name+'_1'] = portOut1
    portOut2 = self.port('portOut2', [-cutout_size[0]/2,0], [-1,0], track_left+2*self.overdev, gap_left-2*self.overdev)
#        self.ports[self.name+'_2'] = portOut2
    
    
def draw_capa(self, iTrack, iGap, pad_spacing, pad_size, half=False):
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


          +--+  +--+
          |  |  |  |
        +-+  |  |  +-+
    iIn |    |  |    | iOut
        +-+  |  |  +-+
          |  |  |  |
          +--+  +--+
    '''
    _pos = self.pos
    if half:
        self.pos = self.pos - self.ori*(pad_spacing/2+pad_size[0]+iGap)
        self.pos2 = self.pos
    iTrack, iGap, pad_spacing, pad_size = parse_entry((iTrack, iGap, pad_spacing, pad_size))
    pad_size = Vector(pad_size)

    portOut1 = [self.pos+self.ori*(pad_spacing/2+pad_size[0]+iGap), self.ori, iTrack+2*self.overdev, iGap-2*self.overdev]
    if not half:
        portOut2 = [self.pos-self.ori*(pad_spacing/2+pad_size[0]+iGap), -self.ori, iTrack+2*self.overdev, iGap-2*self.overdev]

    raw_points = [(pad_spacing/2-self.overdev, pad_size[1]/2+self.overdev),
                  (pad_size[0]+2*self.overdev, 0),
                  (0, -(pad_size[1]-iTrack)/2),
                  (iGap-self.overdev, 0),
                  (0, -iTrack-2*self.overdev),
                  (-iGap+self.overdev, 0),
                  (0, -(pad_size[1]-iTrack)/2),
                  (-pad_size[0]-2*self.overdev, 0)]
    points = self.append_points(raw_points)
    right_pad = self.draw(self.name+"_pad1", points)

    if not half:
        points = self.append_points(self.refy_points(raw_points))
        left_pad = self.draw(self.name+"_pad2", points)

    padlist = [right_pad]
    if not half:
        padlist.append(left_pad)
    pads = self.unite(padlist, name=self.name+'_pads')
    self.trackObjects.append(pads)

    jj = 2 if not half else 1
    pos_cutout = self.pos if not half else self.pos+self.ori*(pad_spacing/4+pad_size[0]/2+iGap/2)
    cutout = self.draw_rect_center(self.name+"_cutout", pos_cutout, self.coor_vec([pad_spacing*jj/2 + jj*pad_size[0]+jj*iGap-jj*self.overdev, pad_size[1] + 2*iGap-2*self.overdev]))
    
    if self.is_overdev:
        sub_1 = self.draw_rect(self.name + '_sub_1', self.coor([pad_spacing/2+pad_size[0]+iGap, -iTrack/2-iGap+self.overdev]), self.coor_vec([-self.overdev, iTrack+2*iGap-2*self.overdev]))
        if not half:
            sub_2 = self.draw_rect(self.name + '_sub_1', self.coor([-pad_spacing/2-pad_size[0]-iGap, -iTrack/2-iGap+self.overdev]), self.coor_vec([self.overdev, iTrack+2*iGap-2*self.overdev]))
        cutout_list = [cutout, sub_1]
        if not half:
            cutout_list.append(sub_2)
        cutout = self.unite(cutout_list)
    
    self.gapObjects.append(cutout)
    if self.is_mask:
        self.maskObjects.append(self.draw_rect_center(self.name+"_mask", self.coor([0,0]), self.coor_vec([pad_spacing + 2*pad_size[0]+4*iGap, pad_size[1] + 4*iGap])))
    if not self.is_litho:
        self.draw(self.name+"_mesh", points)
        self.modeler.assign_mesh_length(self.name+"_mesh",iTrack)

    self.ports[self.name+'_1'] = portOut1
    if not half:
        self.ports[self.name+'_2'] = portOut2

    self.pos = _pos
    
def draw_capa_inline(self, iTrack, iGap, capa_length, pad_spacing, n_pad=1, iTrack_capa=None, iGap_capa=None, premesh=True, tight=False): #iGap_capa is added gap
    
    iTrack, iGap, capa_length, pad_spacing, n_pad = parse_entry((iTrack, iGap, capa_length, pad_spacing, n_pad))
    if iTrack_capa is not None:
        iTrack_capa = parse_entry((iTrack_capa))
        usual=False
    else:
        iTrack_capa = iTrack
        usual=True
        
    if iGap_capa is not None:
        iGap_capa = parse_entry((iGap_capa))
    else:
        iGap_capa=0

    drawn_pads = []
    if n_pad==1:
        drawn_pads.append(self.draw_rect(self.name+'_pad_left', self.coor([-capa_length/2,-iTrack/2-self.overdev]), self.coor_vec([capa_length/2-pad_spacing/2+self.overdev,iTrack+2*self.overdev])))
        drawn_pads.append(self.draw_rect(self.name+'_pad_right', self.coor([capa_length/2,-iTrack/2-self.overdev]), self.coor_vec([-(capa_length/2-pad_spacing/2)-self.overdev,iTrack+2*self.overdev])))
    else:
        pad_width = (iTrack_capa-(n_pad-1)*pad_spacing)/n_pad
        if not usual:
#                pad_length = capa_length-pad_spacing-2*iTrack
            pad_length = capa_length-pad_spacing-2*pad_spacing
#                offset = iTrack
            offset = pad_spacing
        else:
            pad_length = capa_length-pad_spacing
            offset = 0
        curr_height = -iTrack_capa/2
        pad_size = Vector([pad_length, pad_width])
        for ii in range(int(n_pad/2)):
            drawn_pads.append(self.draw_rect(self.name+"_pad"+str(ii), self.coor([-capa_length/2+offset, curr_height-self.overdev]), self.coor_vec(pad_size+Vector([self.overdev, 2*self.overdev]))))
            drawn_pads.append(self.draw_rect(self.name+"_pad"+str(ii)+'b', self.coor([-capa_length/2+offset+pad_spacing-self.overdev, curr_height+pad_width+pad_spacing-self.overdev]), self.coor_vec(pad_size+Vector([self.overdev, 2*self.overdev]))))
            curr_height = curr_height+2*(pad_width+pad_spacing)
        if n_pad%2!=0:
            drawn_pads.append(self.draw_rect(self.name+"_pad"+str(ii+1), self.coor([-capa_length/2+offset, curr_height-self.overdev]), self.coor_vec(pad_size+Vector([self.overdev, 2*self.overdev]))))
        else:
            curr_height = curr_height-2*(pad_width+pad_spacing)
        if not usual:
#                drawn_pads.append(self.draw_rect(self.name+'_connect_left', self.coor([-capa_length/2-self.overdev, -iTrack_capa/2-self.overdev]), self.coor_vec([iTrack+2*self.overdev, curr_height+iTrack_capa/2+pad_width+2*self.overdev])))
            drawn_pads.append(self.draw_rect(self.name+'_connect_left', self.coor([-capa_length/2-self.overdev, -iTrack_capa/2-self.overdev]), self.coor_vec([offset+2*self.overdev, curr_height+iTrack_capa/2+pad_width+2*self.overdev])))

            if n_pad%2!=0:
#                    drawn_pads.append(self.draw_rect(self.name+'_connect_right', self.coor([capa_length/2+self.overdev, -iTrack_capa/2+pad_width+pad_spacing-self.overdev]), self.coor_vec([-iTrack-2*self.overdev, curr_height+(iTrack_capa/2-2*pad_width-2*pad_spacing)+pad_width+2*self.overdev])))
                drawn_pads.append(self.draw_rect(self.name+'_connect_right', self.coor([capa_length/2+self.overdev, -iTrack_capa/2+pad_width+pad_spacing-self.overdev]), self.coor_vec([-offset-2*self.overdev, curr_height+(iTrack_capa/2-2*pad_width-2*pad_spacing)+pad_width+2*self.overdev])))
            else:
#                    drawn_pads.append(self.draw_rect(self.name+'_connect_right', self.coor([capa_length/2+self.overdev, -iTrack_capa/2+pad_width+pad_spacing-self.overdev]), self.coor_vec([-iTrack-2*self.overdev, curr_height+iTrack_capa/2+pad_width+2*self.overdev])))
                drawn_pads.append(self.draw_rect(self.name+'_connect_right', self.coor([capa_length/2+self.overdev, -iTrack_capa/2+pad_width+pad_spacing-self.overdev]), self.coor_vec([-offset-2*self.overdev, curr_height+iTrack_capa/2+pad_width+2*self.overdev])))
#                    
#                    drawn_pads.append(self.draw_rect(self.name+'_connect_right_b', self.coor([capa_length/2+self.overdev, -iTrack/2-self.overdev]), self.coor_vec([-iTrack-2*self.overdev, (iTrack/2+iTrack_capa/2)+2*self.overdev])))
#                    drawn_pads.append(self.draw_rect(self.name+'_connect_left_b', self.coor([-capa_length/2-self.overdev, iTrack/2+self.overdev]), self.coor_vec([iTrack+2*self.overdev, -(iTrack/2+iTrack_capa/2)-2*self.overdev])))
                drawn_pads.append(self.draw_rect(self.name+'_connect_right_b', self.coor([capa_length/2+self.overdev, -iTrack/2-self.overdev]), self.coor_vec([-offset-2*self.overdev, (iTrack/2+iTrack_capa/2)+2*self.overdev])))
                drawn_pads.append(self.draw_rect(self.name+'_connect_left_b', self.coor([-capa_length/2-self.overdev, iTrack/2+self.overdev]), self.coor_vec([offset+2*self.overdev, -(iTrack/2+iTrack_capa/2)-2*self.overdev])))

    if tight:
        portOut1 = [self.pos+self.ori*capa_length/2, self.ori, iTrack_capa-2*pad_width-2*pad_spacing+2*self.overdev, iGap-(iTrack_capa-2*pad_width-2*pad_spacing-iTrack)/2-2*self.overdev]
    else:
        portOut1 = [self.pos+self.ori*capa_length/2, self.ori, iTrack+2*self.overdev, iGap-2*self.overdev]

#            portOut2 = [self.pos-self.ori*capa_length/2, -self.ori, iTrack_capa+2*self.overdev, iGap-(iTrack_capa-iTrack)/2-2*self.overdev]
#        else:
    portOut2 = [self.pos-self.ori*capa_length/2, -self.ori, iTrack+2*self.overdev, iGap-2*self.overdev]

    if self.is_overdev:
        drawn_pads.append(self.draw_rect(self.name+'_added_right', self.coor([capa_length/2, -iTrack/2]), self.coor_vec([-self.overdev, iTrack])))
        drawn_pads.append(self.draw_rect(self.name+'_added_left', self.coor([-capa_length/2, -iTrack/2]), self.coor_vec([self.overdev, iTrack])))

    pads = self.unite(drawn_pads, name=self.name+'_pads')
    
    self.trackObjects.append(pads)
    
    self.gapObjects.append(self.draw_rect_center(self.name+"_cutout", self.coor([0,0]), self.coor_vec([capa_length, iTrack + 2*iGap+2*iGap_capa-2*self.overdev])))
    if self.is_mask:
        self.maskObjects.append(self.draw_rect_center(self.name+"_mask", self.coor([0,0]), self.coor_vec([capa_length, iTrack + 2*iGap+2*iGap_capa +2*self.gap_mask])))
    
    if premesh:
        if not self.is_litho:
            self.draw_rect_center(self.name+"_mesh", self.coor([0,0]), self.coor_vec([capa_length, iTrack_capa+2*self.overdev]))
            self.modeler.assign_mesh_length(self.name+"_mesh",pad_spacing)

    self.ports[self.name+'_1'] = portOut1
    self.ports[self.name+'_2'] = portOut2
    
def draw_capa_interdigitated(self, iTrack, iGap, teeth_size,gap_size, N_period, fillet):
    '''
    '''
    iTrack, iGap = parse_entry((iTrack, iGap))
    fillet = parse_entry(fillet)
    teeth_size=parse_entry(teeth_size)
    gap_size=parse_entry(gap_size)
    teeth_size = Vector(teeth_size)
    portOut1 = [self.pos+self.ori*(teeth_size[0]+iTrack+iGap), self.ori, iTrack+2*self.overdev, iGap-2*self.overdev]
    portOut2 = [self.pos-self.ori*(teeth_size[0]+iTrack+iGap), -self.ori, iTrack+2*self.overdev, iGap-2*self.overdev]


    N_teeth=2*N_period+1    
    raw_points = [(teeth_size[0], -N_teeth*teeth_size[1]-self.overdev)]
    raw_points.append((teeth_size[0], (-N_teeth+1)*teeth_size[1]))
    for i in range(-N_teeth+1,N_teeth-1,4):
        raw_points.append((-teeth_size[0], i*teeth_size[1]))
        raw_points.append((-teeth_size[0], (i+2)*teeth_size[1]))
        raw_points.append((teeth_size[0], (i+2)*teeth_size[1]))
        raw_points.append((teeth_size[0], (i+4)*teeth_size[1]))
    raw_points.append((-teeth_size[0], (N_teeth-1)*teeth_size[1]))
    raw_points.append((-teeth_size[0], N_teeth*teeth_size[1]+self.overdev))

    points = self.append_absolute_points(raw_points)
    connection = self.draw(self.name+"_capagap", points, closed=False)
    
    
    connection.fillets(fillet)
    raw_points=[(-gap_size-teeth_size[0]+self.overdev,N_teeth*teeth_size[1]+self.overdev),(gap_size-teeth_size[0]-self.overdev,N_teeth*teeth_size[1]+self.overdev)]
    points=self.append_absolute_points(raw_points)
    capagap_starter = self.draw(self.name+'_width', points, closed=False)
    
    capagap = connection.sweep_along_path(capagap_starter)
    
   
    
    raw_points = [(-teeth_size[0]-iTrack-self.overdev, -N_teeth*teeth_size[1]-self.overdev),
                  (-teeth_size[0]-iTrack-self.overdev,-iTrack/2-self.overdev),
                  (-teeth_size[0]-iTrack-iGap,-iTrack/2-self.overdev),
                  (-teeth_size[0]-iTrack-iGap, iTrack/2+self.overdev),
                  (-teeth_size[0]-iTrack-self.overdev, iTrack/2+self.overdev),
                  (-teeth_size[0]-iTrack-self.overdev, N_teeth*teeth_size[1]+self.overdev),
                  (teeth_size[0]+iTrack+self.overdev,  N_teeth*teeth_size[1]+self.overdev),
                  (teeth_size[0]+iTrack+self.overdev,iTrack/2+self.overdev),
                  (teeth_size[0]+iTrack+iGap,iTrack/2+self.overdev),
                  (teeth_size[0]+iTrack+iGap,-iTrack/2-self.overdev),
                  (teeth_size[0]+iTrack+self.overdev, -iTrack/2-self.overdev),
                  (teeth_size[0]+iTrack+self.overdev, -N_teeth*teeth_size[1]-self.overdev)]
    points = self.append_absolute_points(raw_points)
    pads = self.draw(self.name+"_pads", points)
    #####Filets on edges of the capa
    pads.fillet(fillet+self.overdev,11)
    pads.fillet(fillet+self.overdev,6)
    pads.fillet(fillet+self.overdev,5)
    pads.fillet(fillet+self.overdev,0)
    
    pads_sub = self.subtract(pads, [capagap])
    #print(pads_sub.vertices())
    
    #####Filets on edge
    pads.fillet(fillet-self.overdev,73)
    pads.fillet(fillet-self.overdev,70)

    #####Filets on last teeth
    pads.fillet(0.5*fillet+self.overdev,67)
    pads.fillet(0.5*fillet+self.overdev,38)
    pads.fillet(0.5*fillet+self.overdev,10)

    #####Filets on edge
    pads.fillet(fillet-self.overdev,7)
    pads.fillet(fillet-self.overdev,4)

    #####Filets on last teeth
    pads.fillet(0.5*fillet+self.overdev,1)





    if not self.is_overdev:
        self.gapObjects.append(self.draw_rect_center(self.name+"_gap", self.coor([0,0]), self.coor_vec([2*teeth_size[0]+2*iTrack+2*iGap, 2*N_teeth*teeth_size[1]+2*iGap])))
    else:
        raw_points = [(-teeth_size[0]-iTrack-iGap, iTrack/2+iGap-self.overdev),
                      (-teeth_size[0]-iTrack-iGap+self.overdev, iTrack/2+iGap-self.overdev),
                      (-teeth_size[0]-iTrack-iGap+self.overdev, N_teeth*teeth_size[1]+iGap-self.overdev),
                      (teeth_size[0]+iTrack+iGap-self.overdev, N_teeth*teeth_size[1]+iGap-self.overdev),
                      (teeth_size[0]+iTrack+iGap-self.overdev, iTrack/2+iGap-self.overdev),
                      (teeth_size[0]+iTrack+iGap, iTrack/2+iGap-self.overdev),
                      (teeth_size[0]+iTrack+iGap, -iTrack/2-iGap+self.overdev),
                      (teeth_size[0]+iTrack+iGap-self.overdev, -iTrack/2-iGap+self.overdev),
                      (teeth_size[0]+iTrack+iGap-self.overdev, -N_teeth*teeth_size[1]-iGap+self.overdev),
                      (-teeth_size[0]-iTrack-iGap+self.overdev, -N_teeth*teeth_size[1]-iGap+self.overdev),
                      (-teeth_size[0]-iTrack-iGap+self.overdev, -iTrack/2-iGap+self.overdev),
                      (-teeth_size[0]-iTrack-iGap, -iTrack/2-iGap+self.overdev)]
        points = self.append_absolute_points(raw_points)
        self.gapObjects.append(self.draw(self.name+"_gap", points))
    if self.is_mask:
        self.maskObjects.append(self.draw_rect_center(self.name+"_mask", self.coor([0,0]), self.coor_vec([2*teeth_size[0]+2*iTrack+4*iGap, 2*N_teeth*teeth_size[1]+4*iGap])))

        
    if not self.is_litho:
        self.draw(self.name+"_mesh", points)
        self.modeler.assign_mesh_length(self.name+"_mesh",iTrack)

    self.trackObjects.append(pads_sub)



    self.ports[self.name+'_1'] = portOut1
    self.ports[self.name+'_2'] = portOut2
#    

def draw_squid(self, iTrack, iGap, squid_size, iTrackPump, iGapPump, iTrackSquid=None, iTrackJ=None, Lj_down='1nH', Lj_up=None,  typePump='down', doublePump=False, iSlope=1, iSlopePump=0.5, fillet=None): #for now assume left and right tracks are the same width
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
    if iTrackSquid is None:
        iTrackSquid = iTrack/4
    if iTrackJ is None:
        iTrackJ = iTrackSquid/2
    if Lj_up is None:
        Lj_up = Lj_down
    iTrack, iGap, squid_size, iTrackPump, iGapPump, iTrackSquid, iTrackJ = parse_entry((iTrack, iGap, squid_size, iTrackPump, iGapPump, iTrackSquid, iTrackJ))
    squid_size = Vector(squid_size)

    adapt_dist = squid_size[1]/2 #if slope ==1 !!!!

    #adapta
    raw_points_a = [(squid_size[0]/2+adapt_dist,0),
                    (0, iTrack/2),
                    (-adapt_dist, squid_size[1]/2+iTrackSquid-iTrack/2),
                    (-squid_size[0]/2+iTrackSquid/2, 0),
                    (0, -iTrackSquid),
                    (squid_size[0]/2-iTrackSquid/2, 0)]
    points_a = self.append_points(raw_points_a)
    track_a = self.draw(self.name+"_track_a", points_a)

    raw_points_b = self.refx_points(raw_points_a)
    points_b = self.append_points(raw_points_b)
    track_b = self.draw(self.name+"_track_b", points_b)

    raw_points_c = self.refy_points(raw_points_a)
    points_c = self.append_points(raw_points_c)
    track_c = self.draw(self.name+"_track_c", points_c)

    raw_points_d = self.refy_points(raw_points_b)
    points_d = self.append_points(raw_points_d)
    track_d = self.draw(self.name+"_track_d", points_d)

    #junction up
    print(self.ori)
    in_junction_up = [self.coor([-iTrackSquid/2,squid_size[1]/2+iTrackSquid/2]), self.coor_vec([1,0]), iTrackSquid, 0]
    out_junction_up = [self.coor([iTrackSquid/2,squid_size[1]/2+iTrackSquid/2]), self.coor_vec([-1,0]), iTrackSquid, 0]
    junction = self.connect_elt(self.name+'_junction_up', in_junction_up, out_junction_up)
    junction_pads_up = junction._connect_JJ(iTrackJ, iInduct=Lj_up, fillet=fillet)
    #junction down
    in_junction_down = [self.coor([-iTrackSquid/2,-squid_size[1]/2-iTrackSquid/2]), self.coor_vec([1,0]), iTrackSquid, 0]
    out_junction_down = [self.coor([iTrackSquid/2,-squid_size[1]/2-iTrackSquid/2]), self.coor_vec([-1,0]), iTrackSquid, 0]
    junction = self.connect_elt(self.name+'_junction_down', in_junction_down, out_junction_down)
    junction_pads_down = junction._connect_JJ(iTrackJ, iInduct=Lj_up, fillet=fillet)

    right_track = self.draw_rect_center(self.name+"_added_track1", self.coor([2*(squid_size[0]/2+adapt_dist),0]), self.coor_vec([squid_size[0]+2*adapt_dist, iTrack]))
    left_track = self.draw_rect_center(self.name+"_added_track2", self.coor([-2*(squid_size[0]/2+adapt_dist),0]), self.coor_vec([squid_size[0]+2*adapt_dist, iTrack]))

    squid = self.unite([right_track, left_track, track_a, track_b, track_c, track_d, junction_pads_down, junction_pads_up], name=self.name)
    self.trackObjects.append(squid)


    if fillet is not None:
        fillet=parse_entry(fillet)
        squid.fillet(fillet,[32,31,26,25,24,19,18,16,11,10,9,4,3,0])

    adapt_dist_pump = 4*iTrackPump#(4*iTrackPump - 2*iTrackSquid)/2/iSlopePump


    self.gapObjects.append(self.draw_rect_center(self.name+"_added_gap", self.coor([0,0]), self.coor_vec([3*(squid_size[0]+2*adapt_dist), iTrack+2*iGap])))
    if self.is_mask:
        self.maskObjects.append(self.draw_rect_center(self.name+"_added_gap", self.coor([0,0]), self.coor_vec([3*(squid_size[0]+2*adapt_dist), iTrack+3*iGap])))

    raw_points_adapt_pump_a = [(3/2*(squid_size[0]+2*adapt_dist),-iTrack/2-iGap-iTrackSquid),
                                     (-(squid_size[0]+2*adapt_dist)+iTrackSquid,0),
                                     (iTrackPump/2-iTrackSquid, -adapt_dist_pump),
                                     (iGapPump, 0)]
    if self.is_mask:
        raw_points_adapt_pump_mask = [(3/2*(squid_size[0]+2*adapt_dist)+iGapPump,-iTrack/2-iGap-iTrackSquid-iGapPump),
                                      (0,iTrackSquid+iGapPump),
                                         (-(squid_size[0]+2*adapt_dist)+iTrackSquid-iGapPump-iTrackPump/2,0),
                                         (iTrackPump/2-iTrackSquid, -adapt_dist_pump-iTrackSquid),
                                         (3*iGapPump+iTrackPump/2, 0)]
        self.maskObjects.append(self.draw(self.name+"_cutout_pump_a_mask", self.append_points(raw_points_adapt_pump_mask)))
        raw_points_adapt_pump_mask_b = self.refy_points(raw_points_adapt_pump_mask, offset = (squid_size[0]+2*adapt_dist)/2)
        self.maskObjects.append(self.draw(self.name+"_cutout_pump_b_mask", self.append_points(raw_points_adapt_pump_mask_b)))

    if typePump == 'up' or typePump == 'Up':
        raw_points_adapt_pump_a = self.refx_points(raw_points_adapt_pump_a)
        ori_pump = [0,1]
        pos_pump = [(squid_size[0]+2*adapt_dist)/2, iTrack/2+iGap+iTrackSquid+adapt_dist_pump]
    elif typePump =='down' or typePump == 'Down':
        ori_pump = [0,-1]
        pos_pump = [(squid_size[0]+2*adapt_dist)/2, -iTrack/2-iGap-iTrackSquid-adapt_dist_pump]
    else:
        raise ValueError("typePump should be 'up' or 'down', given %s" % typePump)
    points_adapt_pump_a = self.append_points(raw_points_adapt_pump_a)
    cutout_pump_a=self.draw(self.name+"_cutout_pump_a", points_adapt_pump_a)
    self.gapObjects.append(cutout_pump_a)
    
    raw_points_adapt_pump_b = self.refy_points(raw_points_adapt_pump_a, offset = (squid_size[0]+2*adapt_dist)/2)
    points_adapt_pump_b = self.append_points(raw_points_adapt_pump_b)
    cutout_pump_b=self.draw(self.name+"_cutout_pump_b", points_adapt_pump_b)
    self.gapObjects.append(cutout_pump_b)
    
    if fillet is not None:
        cutout_pump_a.fillet(fillet/2,1)
        cutout_pump_a.fillet(fillet/4,0)
        cutout_pump_b.fillet(fillet/2,1)
        cutout_pump_b.fillet(fillet/4,0)
    if doublePump:
        raw_points_adapt_pump_c = self.refx_points(raw_points_adapt_pump_a)
        points_adapt_pump_c = self.append_points(raw_points_adapt_pump_c)
        self.gapObjects.append(self.draw(self.name+"_cutout_pump_c", points_adapt_pump_c))

        raw_points_adapt_pump_d = self.refx_points(raw_points_adapt_pump_b)
        points_adapt_pump_d = self.append_points(raw_points_adapt_pump_d)
        self.gapObjects.append(self.draw(self.name+"_cutout_pump_d", points_adapt_pump_d))

    portOut1 = [self.pos+self.ori*(squid_size[0]/2+adapt_dist)*3, self.ori, iTrack, iGap]
    portOut2 = [self.pos-self.ori*(squid_size[0]/2+adapt_dist)*3, -self.ori, iTrack, iGap]
    self.ports[self.name+'_1'] = portOut1
    self.ports[self.name+'_2'] = portOut2

    if doublePump:
        portOutpump1 = [self.coor(pos_pump), self.coor_vec(ori_pump), iTrackPump, iGapPump]
        self.ports[self.name+'_pump1'] = portOutpump1
        portOutpump2 = [self.coor(Vector(pos_pump).refx()), -self.coor_vec(ori_pump), iTrackPump, iGapPump]
        self.ports[self.name+'_pump2'] = portOutpump2
    else:
        portOutpump1 = [self.coor(pos_pump), self.coor_vec(ori_pump), iTrackPump, iGapPump]
        self.ports[self.name+'_pump'] = portOutpump1


def draw_snails(self, iTrack, iGap, array_room, array_offset, iTrackPump, iGapPump, mode='litho', iTrackMinPump=None, iTrackSnail=None, fillet=None, N_snails=1, snail_dict={'loop_width':20e-6, 'loop_length':20e-6, 'length_big_junction':10e-6, 'length_small_junction':2e-6, 'bridge':1e-6, 'bridge_spacing':1e-6}, L_eq = '1nH'): #for now assume left and right tracks are the same width
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
    if iTrackSnail is None:
        iTrackSnail = iTrack/10
    if iTrackMinPump is None:
        iTrackMinPump=iTrackSnail

    iTrack, iGap, array_room, array_offset, iTrackPump, iGapPump, iTrackMinPump, iTrackSnail = parse_entry((iTrack, iGap, array_room, array_offset, iTrackPump, iGapPump, iTrackMinPump, iTrackSnail))

    adapt_dist = iTrack/2 #if slope ==1 !!!!
    
    if self.val(array_offset)>=0:
        raw_points_a = [(array_room/2+adapt_dist-iTrack,-iTrack/2),
                        (iTrack,0),
                        (0, iTrack),
                        (-adapt_dist, -iTrack/2+iTrackSnail/2+array_offset),
                        (0, -iTrackSnail),
                        (-iTrack, 0)]
    else: 
        array_offset = -array_offset
        raw_points_a = [(array_room/2+adapt_dist-iTrack,-iTrack/2),
            (iTrack,0),
            (0, iTrack),
            (-adapt_dist, -iTrack/2+iTrackSnail/2+array_offset),
            (0, -iTrackSnail),
            (-iTrack, 0)]
        raw_points_a = self.refx_points(raw_points_a)
        array_offset = -array_offset
        
    
    points_a = self.append_points(raw_points_a)
    track_a = self.draw(self.name+"_track_a", points_a)
    
    raw_points_c = self.refy_points(raw_points_a)
    points_c = self.append_points(raw_points_c)
    track_c = self.draw(self.name+"_track_c", points_c)

    #snail array
    if mode=='litho':
        in_array = [self.coor([array_room/2, array_offset]), -self.ori, iTrackSnail, 0]
        out_array = [self.coor([-array_room/2, array_offset]), self.ori, iTrackSnail, 0]      
        snail_array = self.connect_elt(self.name+'_array', out_array, in_array)
        snail_track = snail_array._connect_snails2([snail_dict['loop_width'], snail_dict['loop_length']], snail_dict['length_big_junction'], 4, snail_dict['length_small_junction'], 1, N_snails, snail_dict['bridge'], snail_dict['bridge_spacing'])
    
    if 0:
        in_array = [self.coor([array_room/2, array_offset]), -self.ori, iTrackSnail, 0]
        out_array = [self.coor([-array_room/2, array_offset]), self.ori, iTrackSnail, 0]      
        snail_array = self.connect_elt(self.name+'_array', in_array, out_array)
        snail_track = snail_array._connect_JJ(iTrack, iInduct=L_eq)
    
    if mode=='equivalent':
        connect_left = self.draw_rect(self.name+"_left", self.coor([array_room/2,array_offset-iTrackSnail/2]), self.coor_vec([-iTrackSnail, iTrackSnail]))
        connect_right = self.draw_rect(self.name+"_right", self.coor([-array_room/2,array_offset-iTrackSnail/2]), self.coor_vec([iTrackSnail, iTrackSnail]))
        array_eq = self.draw_rect_center(self.name+"_array_eq", self.coor([0,array_offset]), self.coor_vec([array_room-2*iTrackSnail, iTrackSnail]))
        self.assign_lumped_RLC(array_eq, self.ori, (0, L_eq, 0))
        points = self.append_points([(-array_room/2+iTrackSnail,array_offset),(array_room-2*iTrackSnail,0)])
        self.draw(self.name+'_array_eq_line', points, closed=False)
        
#        #junction up
#        print(self.ori)
#        in_junction_up = [self.coor([-iTrackSnail/2,squid_size[1]/2+iTrackSnail/2]), self.coor_vec([1,0]), iTrackSnail, 0]
#        out_junction_up = [self.coor([iTrackSnail/2,squid_size[1]/2+iTrackSnail/2]), self.coor_vec([-1,0]), iTrackSnail, 0]
#        junction = self.connect_elt(self.name+'_junction_up', in_junction_up, out_junction_up)
#        junction_pads_up = junction._connect_JJ(iTrackJ, iInduct=Lj_up, fillet=fillet)
#        
#        #junction down
#        in_junction_down = [self.coor([-iTrackSnail/2,-squid_size[1]/2-iTrackSnail/2]), self.coor_vec([1,0]), iTrackSnail, 0]
#        out_junction_down = [self.coor([iTrackSnail/2,-squid_size[1]/2-iTrackSnail/2]), self.coor_vec([-1,0]), iTrackSnail, 0]
#        junction = self.connect_elt(self.name+'_junction_down', in_junction_down, out_junction_down)
#        junction_pads_down = junction._connect_JJ(iTrackJ, iInduct=Lj_up, fillet=fillet)

    right_track = self.draw_rect_center(self.name+"_added_track1", self.coor([2*(array_room/2+adapt_dist),0]), self.coor_vec([array_room+2*adapt_dist, iTrack+2*self.overdev]))
    left_track = self.draw_rect_center(self.name+"_added_track2", self.coor([-2*(array_room/2+adapt_dist),0]), self.coor_vec([array_room+2*adapt_dist, iTrack+2*self.overdev]))
    if mode=='equivalent':
        squid = self.unite([right_track, left_track, track_a, track_c, connect_left, connect_right], name=self.name+'_temp')
    else:
        squid = self.unite([right_track, left_track, track_a, track_c], name=self.name+'_temp')
    if fillet is not None:
        squid.fillet(iTrack/2,[0, 3, 8, 12])
        
    if 0:
        squid = self.unite([squid, snail_track], name=self.name)
    self.trackObjects.append(squid)

    adapt_dist_pump = 4*iTrackPump
    
    gaps=[]
    masks=[]

    gaps.append(self.draw_rect_center(self.name+"_added_gap", self.coor([0,0]), self.coor_vec([3*(array_room+2*adapt_dist), iTrack+2*iGap-2*self.overdev])))
    if self.is_mask:
        masks.append(self.draw_rect_center(self.name+"_added_gap", self.coor([0,0]), self.coor_vec([3*(array_room+2*adapt_dist), iTrack+2*iGap+2*self.gap_mask])))

    raw_points_adapt_pump_a = [(3/2*(array_room+2*adapt_dist)-self.overdev,-iTrack/2-iGap-iTrackMinPump-self.overdev),
                               (-((array_room+2*adapt_dist))+iTrackMinPump+2*self.overdev,0),
                               (iTrackPump/2-iTrackMinPump, -adapt_dist_pump+self.overdev),
                               (iGapPump-2*self.overdev, 0)]
    if self.is_mask:
        raw_points_adapt_pump_a_mask = [(3/2*(array_room+2*adapt_dist)+self.gap_mask,-iTrack/2-iGap-iTrackMinPump+self.gap_mask),
                                        (-((array_room+2*adapt_dist))+iTrackMinPump-2*self.gap_mask,0),
                                        (iTrackPump/2-iTrackMinPump-(iTrackPump/2-self.gap_mask), -adapt_dist_pump-self.gap_mask),
                                        (iGapPump+2*self.gap_mask+(iTrackPump/2-self.gap_mask), 0)]
    typePump='down'
    doublePump=False

    if typePump == 'up' or typePump == 'Up':
        raw_points_adapt_pump_a = self.refx_points(raw_points_adapt_pump_a)
        if self.is_mask:
            raw_points_adapt_pump_a_mask = self.refx_points(raw_points_adapt_pump_a_mask)
        ori_pump = [0,1]
        pos_pump = [((array_room+2*adapt_dist))/2, iTrack/2+iGap+iTrackMinPump+adapt_dist_pump]
    elif typePump =='down' or typePump == 'Down':
        ori_pump = [0,-1]
        pos_pump = [((array_room+2*adapt_dist))/2, -iTrack/2-iGap-iTrackMinPump-adapt_dist_pump]
    else:
        raise ValueError("typePump should be 'up' or 'down', given %s" % typePump)
    points_adapt_pump_a = self.append_points(raw_points_adapt_pump_a)
    cutout_pump_a=self.draw(self.name+"_cutout_pump_a", points_adapt_pump_a)
    gaps.append(cutout_pump_a)
    if self.is_mask:
        points_adapt_pump_a_mask = self.append_points(raw_points_adapt_pump_a_mask)
        masks.append(self.draw(self.name+"_cutout_pump_a_mask", points_adapt_pump_a_mask))
        
    raw_points_adapt_pump_b = self.refy_points(raw_points_adapt_pump_a, offset = ((array_room+2*adapt_dist))/2)
    points_adapt_pump_b = self.append_points(raw_points_adapt_pump_b)
    cutout_pump_b=self.draw(self.name+"_cutout_pump_b", points_adapt_pump_b)
    gaps.append(cutout_pump_b)
    if self.is_mask:
        raw_points_adapt_pump_b_mask = self.refy_points(raw_points_adapt_pump_a_mask, offset = ((array_room+2*adapt_dist))/2)
        points_adapt_pump_b_mask = self.append_points(raw_points_adapt_pump_b_mask)
        masks.append(self.draw(self.name+"_cutout_pump_b_mask", points_adapt_pump_b_mask))
    

    if doublePump:
        raw_points_adapt_pump_c = self.refx_points(raw_points_adapt_pump_a)
        points_adapt_pump_c = self.append_points(raw_points_adapt_pump_c)
        gaps.append(self.draw(self.name+"_cutout_pump_c", points_adapt_pump_c))
        if self.is_mask:
            raw_points_adapt_pump_c_mask = self.refx_points(raw_points_adapt_pump_a_mask)
            points_adapt_pump_c_mask = self.append_points(raw_points_adapt_pump_c_mask)
            masks.append(self.draw(self.name+"_cutout_pump_c_mask", points_adapt_pump_c_mask))

        raw_points_adapt_pump_d = self.refx_points(raw_points_adapt_pump_b)
        points_adapt_pump_d = self.append_points(raw_points_adapt_pump_d)
        gaps.append(self.draw(self.name+"_cutout_pump_d", points_adapt_pump_d))
        if self.is_mask:
            raw_points_adapt_pump_d_mask = self.refx_points(raw_points_adapt_pump_b_mask)
            points_adapt_pump_d_mask = self.append_points(raw_points_adapt_pump_d_mask)
            masks.append(self.draw(self.name+"_cutout_pump_d_mask", points_adapt_pump_d_mask))

    gaps = self.unite(gaps, self.name+'_cutout')
    self.gapObjects.append(gaps)
    if self.is_mask:
        masks = self.unite(masks, self.name+'_mask')
        self.maskObjects.append(masks)
        
    portOut1 = [self.coor([(3*(array_room+2*adapt_dist))/2,0]), self.ori, iTrack+2*self.overdev, iGap-2*self.overdev]
    portOut2 = [self.coor([-(3*(array_room+2*adapt_dist))/2,0]), -self.ori, iTrack+2*self.overdev, iGap-2*self.overdev]
    self.ports[self.name+'_1'] = portOut1
    self.ports[self.name+'_2'] = portOut2

    if doublePump:
        portOutpump1 = [self.coor(pos_pump), self.coor_vec(ori_pump), iTrackPump+2*self.overdev, iGapPump-2*self.overdev]
        self.ports[self.name+'_pump1'] = portOutpump1
        portOutpump2 = [self.coor(Vector(pos_pump).refx()), -self.coor_vec(ori_pump), iTrackPump+2*self.overdev, iGapPump-2*self.overdev]
        self.ports[self.name+'_pump2'] = portOutpump2
    else:
        portOutpump1 = [self.coor(pos_pump), self.coor_vec(ori_pump), iTrackPump+2*self.overdev, iGapPump-2*self.overdev]
        self.ports[self.name+'_pump'] = portOutpump1


def draw_snails_Zaki(self, iTrack, iGap, array_room, array_offset, iTrackPump,
                iGapPump, iTrackSnail=None, fillet=None, N_snails=1,
                snail_dict={'loop_width':20e-6, 'loop_length':20e-6,
                            'length_big_junction':10e-6,
                            'length_small_junction':2e-6, 'bridge':1e-6,
                            'bridge_spacing':1e-6}, L_eq = '1nH'): #for now assume left and right tracks are the same width
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
    if iTrackSnail is None:
        iTrackSnail = iTrack/10

    iTrack, iGap, array_room, array_offset, iTrackPump, iGapPump, iTrackSnail = parse_entry((iTrack, iGap, array_room, array_offset, iTrackPump, iGapPump, iTrackSnail))

    adapt_dist = iTrack/2 #if slope ==1 !!!!

    #adapta
    raw_points_a = [(array_room/2+adapt_dist,-iTrack/2),
                    (0, iTrack),
                    (-adapt_dist, -iTrack/2+iTrackSnail/2+array_offset),
                    (0, -iTrackSnail)]
    
    points_a = self.append_points(raw_points_a)
    track_a = self.draw(self.name+"_track_a", points_a)
    
    raw_points_c = self.refy_points(raw_points_a)
    points_c = self.append_points(raw_points_c)
    track_c = self.draw(self.name+"_track_c", points_c)

    #snail array
    if L_eq is None:
        in_array = [self.coor([array_room/2, array_offset]), -self.ori, iTrackSnail, 0]
        out_array = [self.coor([-array_room/2, array_offset]), self.ori, iTrackSnail, 0]      
        snail_array = self.connect_elt(self.name+'_array', in_array, out_array)
        snail_track = snail_array._connect_snails_Zaki([snail_dict['loop_width'], snail_dict['loop_length']], snail_dict['length_big_junction'], 3, snail_dict['length_small_junction'], 1, N_snails, snail_dict['bridge'], snail_dict['bridge_spacing'])
    
    if L_eq is not None:
        array_eq = self.draw_rect_center(self.name+"_array_eq", self.coor([0,array_offset]), self.coor_vec([array_room, iTrackSnail]))
        self.assign_lumped_RLC(array_eq, self.ori, (0, L_eq, 0))
        points = self.append_points([(-array_room/2,array_offset),(array_room,0)])
        self.draw(self.name+'_array_eq_line', points, closed=False)
        
#        #junction up
#        print(self.ori)
#        in_junction_up = [self.coor([-iTrackSnail/2,squid_size[1]/2+iTrackSnail/2]), self.coor_vec([1,0]), iTrackSnail, 0]
#        out_junction_up = [self.coor([iTrackSnail/2,squid_size[1]/2+iTrackSnail/2]), self.coor_vec([-1,0]), iTrackSnail, 0]
#        junction = self.connect_elt(self.name+'_junction_up', in_junction_up, out_junction_up)
#        junction_pads_up = junction._connect_JJ(iTrackJ, iInduct=Lj_up, fillet=fillet)
#        
#        #junction down
#        in_junction_down = [self.coor([-iTrackSnail/2,-squid_size[1]/2-iTrackSnail/2]), self.coor_vec([1,0]), iTrackSnail, 0]
#        out_junction_down = [self.coor([iTrackSnail/2,-squid_size[1]/2-iTrackSnail/2]), self.coor_vec([-1,0]), iTrackSnail, 0]
#        junction = self.connect_elt(self.name+'_junction_down', in_junction_down, out_junction_down)
#        junction_pads_down = junction._connect_JJ(iTrackJ, iInduct=Lj_up, fillet=fillet)

    right_track = self.draw_rect_center(self.name+"_added_track1", self.coor([2*(array_room/2+adapt_dist),0]), self.coor_vec([array_room+2*adapt_dist, iTrack]))
    left_track = self.draw_rect_center(self.name+"_added_track2", self.coor([-2*(array_room/2+adapt_dist),0]), self.coor_vec([array_room+2*adapt_dist, iTrack]))

    squid = self.unite([right_track, left_track, track_a, track_c], name=self.name)
    self.trackObjects.append(squid)


    if fillet is not None:
        squid.fillet(iTrack/2,[0, 3, 7, 10])

    adapt_dist_pump = 4*iTrackPump#(4*iTrackPump - 2*iTrackSnail)/2/iSlopePump


    self.gapObjects.append(self.draw_rect_center(self.name+"_added_gap", self.coor([0,0]), self.coor_vec([3*(array_room+2*adapt_dist), iTrack+2*iGap])))
    if self.is_mask:
        self.maskObjects.append(self.draw_rect_center(self.name+"_added_gap", self.coor([0,0]), self.coor_vec([3*(array_room+2*adapt_dist), iTrack+3*iGap])))

    raw_points_adapt_pump_a = [(3/2*(array_room+2*adapt_dist),-iTrack/2-iGap-iTrackSnail),
                                     (-(array_room+2*adapt_dist)+iTrackSnail,0),
                                     (iTrackPump/2-iTrackSnail, -adapt_dist_pump),
                                     (iGapPump, 0)]
    if self.is_mask:
        raw_points_adapt_pump_mask = [(3/2*(array_room+2*adapt_dist)+iGapPump,-iTrack/2-iGap-iTrackSnail-iGapPump),
                                      (0,iTrackSnail+iGapPump),
                                         (-(array_room+2*adapt_dist)+iTrackSnail-iGapPump-iTrackPump/2,0),
                                         (iTrackPump/2-iTrackSnail, -adapt_dist_pump-iTrackSnail),
                                         (3*iGapPump+iTrackPump/2, 0)]
        self.maskObjects.append(self.draw(self.name+"_cutout_pump_a_mask", self.append_points(raw_points_adapt_pump_mask)))
        raw_points_adapt_pump_mask_b = self.refy_points(raw_points_adapt_pump_mask, offset = (array_room+2*adapt_dist)/2)
        self.maskObjects.append(self.draw(self.name+"_cutout_pump_b_mask", self.append_points(raw_points_adapt_pump_mask_b)))

    typePump='up'
    doublePump=False

    if typePump == 'up' or typePump == 'Up':
        raw_points_adapt_pump_a = self.refx_points(raw_points_adapt_pump_a)
        ori_pump = [0,1]
        pos_pump = [(array_room+2*adapt_dist)/2, iTrack/2+iGap+iTrackSnail+adapt_dist_pump]
    elif typePump =='down' or typePump == 'Down':
        ori_pump = [0,-1]
        pos_pump = [(array_room+2*adapt_dist)/2, -iTrack/2-iGap-iTrackSnail-adapt_dist_pump]
    else:
        raise ValueError("typePump should be 'up' or 'down', given %s" % typePump)
    points_adapt_pump_a = self.append_points(raw_points_adapt_pump_a)
    cutout_pump_a=self.draw(self.name+"_cutout_pump_a", points_adapt_pump_a)
    self.gapObjects.append(cutout_pump_a)
    
    raw_points_adapt_pump_b = self.refy_points(raw_points_adapt_pump_a, offset = (array_room+2*adapt_dist)/2)
    points_adapt_pump_b = self.append_points(raw_points_adapt_pump_b)
    cutout_pump_b=self.draw(self.name+"_cutout_pump_b", points_adapt_pump_b)
    self.gapObjects.append(cutout_pump_b)
    
    if fillet is not None and False:    
        cutout_pump_a.fillet(fillet/2,1)
        cutout_pump_a.fillet(fillet/4,0)
        cutout_pump_b.fillet(fillet/2,1)
        cutout_pump_b.fillet(fillet/4,0)
    if doublePump:
        raw_points_adapt_pump_c = self.refx_points(raw_points_adapt_pump_a)
        points_adapt_pump_c = self.append_points(raw_points_adapt_pump_c)
        self.gapObjects.append(self.draw(self.name+"_cutout_pump_c", points_adapt_pump_c))

        raw_points_adapt_pump_d = self.refx_points(raw_points_adapt_pump_b)
        points_adapt_pump_d = self.append_points(raw_points_adapt_pump_d)
        self.gapObjects.append(self.draw(self.name+"_cutout_pump_d", points_adapt_pump_d))

    portOut1 = [self.pos+self.ori*(array_room/2+adapt_dist)*3, self.ori, iTrack, iGap]
    portOut2 = [self.pos-self.ori*(array_room/2+adapt_dist)*3, -self.ori, iTrack, iGap]
    self.ports[self.name+'_1'] = portOut1
    self.ports[self.name+'_2'] = portOut2

    if doublePump:
        portOutpump1 = [self.coor(pos_pump), self.coor_vec(ori_pump), iTrackPump, iGapPump]
        self.ports[self.name+'_pump1'] = portOutpump1
        portOutpump2 = [self.coor(Vector(pos_pump).refx()), -self.coor_vec(ori_pump), iTrackPump, iGapPump]
        self.ports[self.name+'_pump2'] = portOutpump2
    else:
        portOutpump1 = [self.coor(pos_pump), self.coor_vec(ori_pump), iTrackPump, iGapPump]
        self.ports[self.name+'_pump'] = portOutpump1
			

def draw_snails_sym(self, iTrack, iGap, array_room, iTrackPump, iGapPump, approach='20um', mode='litho', iTrackMinPump=None, iTrackSnail=None, fillet=None, N_snails=1, snail_dict={'loop_width':20e-6, 'loop_length':20e-6, 'length_big_junction':10e-6, 'length_small_junction':2e-6, 'bridge':1e-6, 'bridge_spacing':1e-6}, L_eq = '1nH'): #for now assume left and right tracks are the same width
    '''
    --------

    '''
    if iTrackSnail is None:
        iTrackSnail = iTrack/10
    if iTrackMinPump is None:
        iTrackMinPump=iTrackSnail

    iTrack, iGap, array_room, iTrackPump, iGapPump, iTrackMinPump, iTrackSnail, approach = parse_entry((iTrack, iGap, array_room, iTrackPump, iGapPump, iTrackMinPump, iTrackSnail, approach))

    adapt_dist = iTrack/2 #if slope ==1 !!!!
    
    track_a = self.draw_rect(self.name+"_track_a", self.coor([array_room/2,-iTrackSnail/2]), self.coor_vec([adapt_dist, iTrackSnail]))
    track_c = self.draw_rect(self.name+"_track_c", self.coor([-array_room/2,-iTrackSnail/2]), self.coor_vec([-adapt_dist, iTrackSnail]))
    

    #snail array
    if mode=='litho':
        in_array = [self.coor([array_room/2, 0]), -self.ori, iTrackSnail, 0]
        out_array = [self.coor([-array_room/2, 0]), self.ori, iTrackSnail, 0]      
        snail_array = self.connect_elt(self.name+'_array', out_array, in_array)
        snail_track = snail_array._connect_snails2([snail_dict['loop_width'], snail_dict['loop_length']], snail_dict['length_big_junction'], 4, snail_dict['length_small_junction'], 1, N_snails, snail_dict['bridge'], snail_dict['bridge_spacing'])
    
    if 0:
        in_array = [self.coor([array_room/2, 0]), -self.ori, iTrackSnail, 0]
        out_array = [self.coor([-array_room/2, 0]), self.ori, iTrackSnail, 0]      
        snail_array = self.connect_elt(self.name+'_array', in_array, out_array)
        snail_track = snail_array._connect_JJ(iTrack, iInduct=L_eq)
    
    if mode=='equivalent':
        connect_left = self.draw_rect(self.name+"_left", self.coor([array_room/2,-iTrackSnail/2]), self.coor_vec([-iTrackSnail, iTrackSnail]))
        connect_right = self.draw_rect(self.name+"_right", self.coor([-array_room/2,-iTrackSnail/2]), self.coor_vec([iTrackSnail, iTrackSnail]))
        if not self.is_litho:
            array_eq = self.draw_rect_center(self.name+"_array_eq", self.coor([0,0]), self.coor_vec([array_room-2*iTrackSnail, iTrackSnail]))
            self.assign_lumped_RLC(array_eq, self.ori, (0, L_eq, 0))
            points = self.append_points([(-array_room/2+iTrackSnail,0),(array_room-2*iTrackSnail,0)])
            self.draw(self.name+'_array_eq_line', points, closed=False)
        
#        #junction up
#        print(self.ori)
#        in_junction_up = [self.coor([-iTrackSnail/2,squid_size[1]/2+iTrackSnail/2]), self.coor_vec([1,0]), iTrackSnail, 0]
#        out_junction_up = [self.coor([iTrackSnail/2,squid_size[1]/2+iTrackSnail/2]), self.coor_vec([-1,0]), iTrackSnail, 0]
#        junction = self.connect_elt(self.name+'_junction_up', in_junction_up, out_junction_up)
#        junction_pads_up = junction._connect_JJ(iTrackJ, iInduct=Lj_up, fillet=fillet)
#        
#        #junction down
#        in_junction_down = [self.coor([-iTrackSnail/2,-squid_size[1]/2-iTrackSnail/2]), self.coor_vec([1,0]), iTrackSnail, 0]
#        out_junction_down = [self.coor([iTrackSnail/2,-squid_size[1]/2-iTrackSnail/2]), self.coor_vec([-1,0]), iTrackSnail, 0]
#        junction = self.connect_elt(self.name+'_junction_down', in_junction_down, out_junction_down)
#        junction_pads_down = junction._connect_JJ(iTrackJ, iInduct=Lj_up, fillet=fillet)

    right_track = self.draw_rect_center(self.name+"_added_track1", self.coor([2*(array_room/2+adapt_dist)+iTrackMinPump/2,0]), self.coor_vec([array_room+2*adapt_dist+iTrackMinPump, iTrack+2*self.overdev]))
    left_track = self.draw_rect_center(self.name+"_added_track2", self.coor([-2*(array_room/2+adapt_dist)-iTrackMinPump/2,0]), self.coor_vec([array_room+2*adapt_dist+iTrackMinPump, iTrack+2*self.overdev]))
    if mode=='equivalent':
        squid = self.unite([right_track, left_track, track_a, track_c, connect_left, connect_right], name=self.name+'_temp')
    else:
        squid = self.unite([right_track, left_track, track_a, track_c], name=self.name+'_temp')
        
    self.trackObjects.append(squid)

    adapt_dist_pump = 4*iTrackPump
    
    gaps=[]
    masks=[]
    
    raw_points_gap = [(-(3*(array_room+2*adapt_dist)+2*iTrackMinPump)/2,iTrack/2+iGap-self.overdev),
                      (iTrackMinPump-self.overdev, 0),
                      (0, -approach),
                      (array_room+2*adapt_dist+2*adapt_dist+array_room+2*self.overdev,0),
                      (0,approach),
                      (array_room+2*adapt_dist+iTrackMinPump-self.overdev,0),
                      (0, -(iTrack/2+iGap-self.overdev)),
                      (-(3*(array_room+2*adapt_dist)+2*iTrackMinPump),0)]
    points_gap = self.append_points(raw_points_gap)
    gap1 = self.draw(self.name+"_gap_1", points_gap)
    gap1.fillet(2*iTrackMinPump,[2,3])
    gaps.append(gap1)
    points_gap = self.append_points(self.refy_points(self.refx_points(raw_points_gap)))
    gap2 = self.draw(self.name+"_gap_2", points_gap)
    gap2.fillet(2*iTrackMinPump,[2,3])
    gaps.append(gap2)
    
#        gaps.append(self.draw_rect_center(self.name+"_added_gap", self.coor([0,0]), self.coor_vec([3*(array_room+2*adapt_dist)+2*iTrackMinPump, iTrack+2*iGap-2*self.overdev])))
    if self.is_mask:
        masks.append(self.draw_rect_center(self.name+"_added_gap", self.coor([0,0]), self.coor_vec([3*(array_room+2*adapt_dist), iTrack+2*iGap+2*self.gap_mask])))

    raw_points_adapt_pump_a = [(3/2*(array_room+2*adapt_dist)-self.overdev-iTrackMinPump,-iTrack/2-iGap-iTrackMinPump-self.overdev+approach),
                               (-((array_room+2*adapt_dist))+2*iTrackMinPump+2*self.overdev,0),
                               (iTrackPump/2-iTrackMinPump, -adapt_dist_pump+self.overdev),
                               (iGapPump-2*self.overdev, 0)]
    if self.is_mask:
        raw_points_adapt_pump_a_mask = [(3/2*(array_room+2*adapt_dist)+self.gap_mask,-iTrack/2-iGap-iTrackMinPump+self.gap_mask),
                                        (-((array_room+2*adapt_dist))+iTrackMinPump-2*self.gap_mask,0),
                                        (iTrackPump/2-iTrackMinPump-(iTrackPump/2-self.gap_mask), -adapt_dist_pump-self.gap_mask),
                                        (iGapPump+2*self.gap_mask+(iTrackPump/2-self.gap_mask), 0)]



    

    points_adapt_pump_a = self.append_points(raw_points_adapt_pump_a)
    cutout_pump_a=self.draw(self.name+"_cutout_pump_a", points_adapt_pump_a)
    cutout_pump_a.fillet(iTrackMinPump,[0,1])
    gaps.append(cutout_pump_a)
    if self.is_mask:
        points_adapt_pump_a_mask = self.append_points(raw_points_adapt_pump_a_mask)
        masks.append(self.draw(self.name+"_cutout_pump_a_mask", points_adapt_pump_a_mask))
        
    raw_points_adapt_pump_b = self.refy_points(raw_points_adapt_pump_a, offset = ((array_room+2*adapt_dist))/2)
    points_adapt_pump_b = self.append_points(raw_points_adapt_pump_b)
    cutout_pump_b=self.draw(self.name+"_cutout_pump_b", points_adapt_pump_b)
    cutout_pump_b.fillet(iTrackMinPump,[0,1])
    gaps.append(cutout_pump_b)
    if self.is_mask:
        raw_points_adapt_pump_b_mask = self.refy_points(raw_points_adapt_pump_a_mask, offset = ((array_room+2*adapt_dist))/2)
        points_adapt_pump_b_mask = self.append_points(raw_points_adapt_pump_b_mask)
        masks.append(self.draw(self.name+"_cutout_pump_b_mask", points_adapt_pump_b_mask))
    
    ori_pump = [0,-1]
    pos_pump = [((array_room+2*adapt_dist))/2, -iTrack/2-iGap-iTrackMinPump-adapt_dist_pump+approach]
    
    raw_points_adapt_pump_c = self.refy_points(self.refx_points(raw_points_adapt_pump_a))
    points_adapt_pump_c = self.append_points(raw_points_adapt_pump_c)
    cutout_pump_c = self.draw(self.name+"_cutout_pump_c", points_adapt_pump_c)
    cutout_pump_c.fillet(iTrackMinPump,[0,1])
    gaps.append(cutout_pump_c)
    if self.is_mask:
        raw_points_adapt_pump_c_mask = self.refy_points(self.refx_points(raw_points_adapt_pump_a_mask))
        points_adapt_pump_c_mask = self.append_points(raw_points_adapt_pump_c_mask)
        masks.append(self.draw(self.name+"_cutout_pump_c_mask", points_adapt_pump_c_mask))

    raw_points_adapt_pump_d = self.refy_points(self.refx_points(raw_points_adapt_pump_b))
    points_adapt_pump_d = self.append_points(raw_points_adapt_pump_d)
    cutout_pump_d = self.draw(self.name+"_cutout_pump_d", points_adapt_pump_d)
    cutout_pump_d.fillet(iTrackMinPump,[0,1])
    gaps.append(cutout_pump_d)
    if self.is_mask:
        raw_points_adapt_pump_d_mask = self.refy_points(self.refx_points(raw_points_adapt_pump_b_mask))
        points_adapt_pump_d_mask = self.append_points(raw_points_adapt_pump_d_mask)
        masks.append(self.draw(self.name+"_cutout_pump_d_mask", points_adapt_pump_d_mask))
    
#        if fillet is not None:
#            cutout_pump_a.fillet(fillet/2,1)
#            cutout_pump_a.fillet(fillet/4,0)
#            cutout_pump_b.fillet(fillet/2,1)
#            cutout_pump_b.fillet(fillet/4,0)

    gaps = self.unite(gaps, self.name+'_cutout')
    self.gapObjects.append(gaps)
    if self.is_mask:
        masks = self.unite(masks, self.name+'_mask')
        self.maskObjects.append(masks)
        
    portOut1 = [self.coor([(3*(array_room+2*adapt_dist))/2+iTrackMinPump,0]), self.ori, iTrack+2*self.overdev, iGap-2*self.overdev]
    portOut2 = [self.coor([-(3*(array_room+2*adapt_dist))/2-iTrackMinPump,0]), -self.ori, iTrack+2*self.overdev, iGap-2*self.overdev]
    self.ports[self.name+'_1'] = portOut1
    self.ports[self.name+'_2'] = portOut2


    portOutpump1 = [self.coor(pos_pump), self.coor_vec(ori_pump), iTrackPump+2*self.overdev, iGapPump-2*self.overdev]
    self.ports[self.name+'_pump_1'] = portOutpump1
    portOutpump2 = [self.coor(Vector(pos_pump).refx().refy()), -self.coor_vec(ori_pump), iTrackPump+2*self.overdev, iGapPump-2*self.overdev]
    self.ports[self.name+'_pump_2'] = portOutpump2
#        else:
#            portOutpump1 = [self.coor(pos_pump), self.coor_vec(ori_pump), iTrackPump+2*self.overdev, iGapPump-2*self.overdev]
#            self.ports[self.name+'_pump'] = portOutpump1
        
def draw_squid_protect(self, iTrack, iGap, squid_size, iTrackPump, iGapPump, iTrackSquid=None, iTrackJ=None, Lj_down='1nH', Lj_up=None,  typePump='down', doublePump=False, iSlope=1, iSlopePump=0.5, fillet=None): #for now assume left and right tracks are the same width
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
        
    if iTrackSquid is None:
        iTrackSquid = iTrack/4
    if iTrackJ is None:
        iTrackJ = iTrackSquid/2
    if Lj_up is None:
        Lj_up = Lj_down
    iTrack, iGap, squid_size, iTrackPump, iGapPump, iTrackSquid, iTrackJ = parse_entry((iTrack, iGap, squid_size, iTrackPump, iGapPump, iTrackSquid, iTrackJ))
    squid_size = Vector(squid_size)

    adapt_dist = squid_size[1]/2 #if slope ==1 !!!!

    #adapta
    raw_points_out = [(squid_size[0]/2+iTrackSquid,squid_size[1]/2+iTrackSquid),
                    (-(squid_size[0]+2*iTrackSquid), 0),
                    (0, -(squid_size[1]/2+iTrackSquid)+2*iTrackSquid),
                    (-squid_size[0]*2/3, 0),
                    (0, -4*iTrackSquid),
                    (squid_size[0]*2/3, 0),
                    (0, -(squid_size[1]/2+iTrackSquid)+2*iTrackSquid),
                    ((squid_size[0]+2*iTrackSquid), 0),]
    points_out = self.append_points(raw_points_out)
    track_out = self.draw(self.name+"_track_out", points_out)
    
    raw_points_in = [(squid_size[0]/2,squid_size[1]/2),
                    (-(squid_size[0]), 0),
                    (0, -(squid_size[1]/2)+iTrackSquid),
                    (-squid_size[0]/3, 0),
                    (0, iTrackSquid),
                    (-iTrackSquid,0),
                    (0,-iTrackSquid),
                    (-squid_size[0]/3+iTrackSquid,0),
                    (0, -2*iTrackSquid),
                    (squid_size[0]/3-iTrackSquid, 0),
                    (0, -iTrackSquid),
                    (iTrackSquid,0),
                    (0,iTrackSquid),
                    (squid_size[0]/3, 0),
                    (0, -(squid_size[1]/2)+iTrackSquid),
                    ((squid_size[0]), 0),]
    points_in = self.append_points(raw_points_in)
    track_in = self.draw(self.name+"_track_in", points_in)
    
    track_out=self.subtract(track_out, [track_in])

    #junction up
    print(self.ori)
    in_junction_up = [self.coor([-squid_size[0]/2-squid_size[0]/3-iTrackSquid,3/2*iTrackSquid]), self.coor_vec([1,0]), iTrackSquid, 0]
    out_junction_up = [self.coor([-squid_size[0]/2-squid_size[0]/3,3/2*iTrackSquid]), self.coor_vec([-1,0]), iTrackSquid, 0]
    junction = self.connect_elt(self.name+'_junction_up', in_junction_up, out_junction_up)
    junction_pads_up = junction._connect_JJ(iTrackJ, iInduct=Lj_up, fillet=fillet)
    #junction down
    in_junction_down = [self.coor([-squid_size[0]/2-squid_size[0]/3-iTrackSquid,-3/2*iTrackSquid]), self.coor_vec([1,0]), iTrackSquid, 0]
    out_junction_down = [self.coor([-squid_size[0]/2-squid_size[0]/3,-3/2*iTrackSquid]), self.coor_vec([-1,0]), iTrackSquid, 0]
    junction = self.connect_elt(self.name+'_junction_down', in_junction_down, out_junction_down)
    junction_pads_down = junction._connect_JJ(iTrackJ, iInduct=Lj_up, fillet=fillet)

    right_track = self.draw_rect(self.name+"_added_track1", self.coor([squid_size[0]/2+iTrackSquid,-iTrack/2-self.overdev]), self.coor_vec([squid_size[0]+iTrackSquid, iTrack+2*self.overdev]))
    left_track = self.draw_rect(self.name+"_added_track2", self.coor([-squid_size[0]/2-iTrackSquid-squid_size[0]*2/3,-iTrack/2-self.overdev]), self.coor_vec([-squid_size[0]-iTrackSquid+squid_size[0]*2/3, iTrack+2*self.overdev]))

    squid = self.unite([right_track, left_track, track_out, junction_pads_down, junction_pads_up], name=self.name)
    self.trackObjects.append(squid)


    adapt_dist_pump = 4*iTrackPump#(4*iTrackPump - 2*iTrackSquid)/2/iSlopePump


    gaps=[]
    masks=[]
    
    gaps.append(self.draw_rect_center(self.name+"_added_gap", self.coor([0,0]), self.coor_vec([3*squid_size[0]+4*iTrackSquid, iTrack+2*iGap-2*self.overdev])))
    if self.is_mask:
        masks.append(self.draw_rect_center(self.name+"_mask", self.coor([0,0]), self.coor_vec([3*squid_size[0]+4*iTrackSquid, iTrack+2*iGap+2*self.gap_mask])))

    raw_points_adapt_pump_a = [(3/2*squid_size[0]+2*iTrackSquid-self.overdev,-iTrack/2-iGap-iTrackSquid-self.overdev),
                               (-(squid_size[0])+2*self.overdev,0),
                               (iTrackPump/2-iTrackSquid, -adapt_dist_pump+self.overdev),
                               (iGapPump-2*self.overdev, 0)]
    if self.is_mask:
        raw_points_adapt_pump_a_mask = [(3/2*squid_size[0]+2*iTrackSquid+self.gap_mask,-iTrack/2-iGap-iTrackSquid+self.gap_mask),
                                        (-(squid_size[0])-2*self.gap_mask,0),
                                        (iTrackPump/2-iTrackSquid-(iTrackPump/2-self.gap_mask), -adapt_dist_pump-self.gap_mask),
                                        (iGapPump+2*self.gap_mask+(iTrackPump/2-self.gap_mask), 0)]

    if typePump == 'up' or typePump == 'Up':
        raw_points_adapt_pump_a = self.refx_points(raw_points_adapt_pump_a)
        if self.is_mask:
            raw_points_adapt_pump_a_mask = self.refx_points(raw_points_adapt_pump_a_mask)
        ori_pump = [0,1]
        pos_pump = [(squid_size[0])/2+iTrackSquid, iTrack/2+iGap+iTrackSquid+adapt_dist_pump]
    elif typePump =='down' or typePump == 'Down':
        ori_pump = [0,-1]
        pos_pump = [(squid_size[0])/2+iTrackSquid, -iTrack/2-iGap-iTrackSquid-adapt_dist_pump]
    else:
        raise ValueError("typePump should be 'up' or 'down', given %s" % typePump)
    points_adapt_pump_a = self.append_points(raw_points_adapt_pump_a)
    cutout_pump_a=self.draw(self.name+"_cutout_pump_a", points_adapt_pump_a)
    gaps.append(cutout_pump_a)
    if self.is_mask:
        points_adapt_pump_a_mask = self.append_points(raw_points_adapt_pump_a_mask)
        masks.append(self.draw(self.name+"_cutout_pump_a_mask", points_adapt_pump_a_mask))
        
    raw_points_adapt_pump_b = self.refy_points(raw_points_adapt_pump_a, offset = (squid_size[0])/2+iTrackSquid)
    points_adapt_pump_b = self.append_points(raw_points_adapt_pump_b)
    cutout_pump_b=self.draw(self.name+"_cutout_pump_b", points_adapt_pump_b)
    gaps.append(cutout_pump_b)
    if self.is_mask:
        raw_points_adapt_pump_b_mask = self.refy_points(raw_points_adapt_pump_a_mask, offset = (squid_size[0])/2+iTrackSquid)
        points_adapt_pump_b_mask = self.append_points(raw_points_adapt_pump_b_mask)
        masks.append(self.draw(self.name+"_cutout_pump_b_mask", points_adapt_pump_b_mask))
    
    if fillet is not None:
        cutout_pump_a.fillet(fillet/2,1)
        cutout_pump_a.fillet(fillet/4,0)
        cutout_pump_b.fillet(fillet/2,1)
        cutout_pump_b.fillet(fillet/4,0)
    if doublePump:
        raw_points_adapt_pump_c = self.refx_points(raw_points_adapt_pump_a)
        points_adapt_pump_c = self.append_points(raw_points_adapt_pump_c)
        gaps.append(self.draw(self.name+"_cutout_pump_c", points_adapt_pump_c))
        if self.is_mask:
            raw_points_adapt_pump_c_mask = self.refx_points(raw_points_adapt_pump_a_mask)
            points_adapt_pump_c_mask = self.append_points(raw_points_adapt_pump_c_mask)
            masks.append(self.draw(self.name+"_cutout_pump_c_mask", points_adapt_pump_c_mask))

        raw_points_adapt_pump_d = self.refx_points(raw_points_adapt_pump_b)
        points_adapt_pump_d = self.append_points(raw_points_adapt_pump_d)
        gaps.append(self.draw(self.name+"_cutout_pump_d", points_adapt_pump_d))
        if self.is_mask:
            raw_points_adapt_pump_d_mask = self.refx_points(raw_points_adapt_pump_b_mask)
            points_adapt_pump_d_mask = self.append_points(raw_points_adapt_pump_d_mask)
            masks.append(self.draw(self.name+"_cutout_pump_d_mask", points_adapt_pump_d_mask))

    gaps = self.unite(gaps, self.name+'_cutout')
    self.gapObjects.append(gaps)
    if self.is_mask:
        masks = self.unite(masks, self.name+'_mask')
        self.maskObjects.append(masks)
        
    portOut1 = [self.coor([(3*squid_size[0]+4*iTrackSquid)/2,0]), self.ori, iTrack+2*self.overdev, iGap-2*self.overdev]
    portOut2 = [self.coor([-(3*squid_size[0]+4*iTrackSquid)/2,0]), -self.ori, iTrack+2*self.overdev, iGap-2*self.overdev]
    self.ports[self.name+'_1'] = portOut1
    self.ports[self.name+'_2'] = portOut2

    if doublePump:
        portOutpump1 = [self.coor(pos_pump), self.coor_vec(ori_pump), iTrackPump+2*self.overdev, iGapPump-2*self.overdev]
        self.ports[self.name+'_pump1'] = portOutpump1
        portOutpump2 = [self.coor(Vector(pos_pump).refx()), -self.coor_vec(ori_pump), iTrackPump+2*self.overdev, iGapPump-2*self.overdev]
        self.ports[self.name+'_pump2'] = portOutpump2
    else:
        portOutpump1 = [self.coor(pos_pump), self.coor_vec(ori_pump), iTrackPump, iGapPump]
        self.ports[self.name+'_pump'] = portOutpump1

def draw_inductance(self, iTrack, iGap, all_length, ind_length, mode='litho', L_eq = '1nH'): #for now assume left and right tracks are the same width
    '''
    --------

    '''

    iTrack, iGap, all_length, ind_length = parse_entry((iTrack, iGap, all_length, ind_length))

    #snail array
#        if mode=='litho':
#            in_array = [self.coor([all_length/2, 0]), -self.ori, iTrackSnail, 0]
#            out_array = [self.coor([-all_length/2, 0]), self.ori, iTrackSnail, 0]      
#            snail_array = self.connect_elt(self.name+'_array', out_array, in_array)
#            snail_track = snail_array._connect_snails2([snail_dict['loop_width'], snail_dict['loop_length']], snail_dict['length_big_junction'], 4, snail_dict['length_small_junction'], 1, N_snails, snail_dict['bridge'], snail_dict['bridge_spacing'])
#        
#        if 0:
#            in_array = [self.coor([array_room/2, 0]), -self.ori, iTrackSnail, 0]
#            out_array = [self.coor([-array_room/2, 0]), self.ori, iTrackSnail, 0]      
#            snail_array = self.connect_elt(self.name+'_array', in_array, out_array)
#            snail_track = snail_array._connect_JJ(iTrack, iInduct=L_eq)
    
    if mode=='equivalent':
        connect_left = self.draw_rect(self.name+"_left", self.coor([ind_length/2,-iTrack/2]), self.coor_vec([all_length/2-ind_length/2, iTrack]))
        connect_right = self.draw_rect(self.name+"_right", self.coor([-ind_length/2,-iTrack/2]), self.coor_vec([-(all_length/2-ind_length/2), iTrack]))
        if not self.is_litho:
            array_eq = self.draw_rect_center(self.name+"_eq", self.coor([0,0]), self.coor_vec([ind_length, iTrack]))
            self.assign_lumped_RLC(array_eq, self.ori, (0, L_eq, 0))
            points = self.append_points([(-ind_length/2,0),(ind_length,0)])
            self.draw(self.name+'_eq_line', points, closed=False)
 

    connect = self.unite([connect_left, connect_right])
    self.trackObjects.append(connect)
    
    gap = self.draw_rect(self.name+"_cutout", self.coor([-all_length/2,-iTrack/2-iGap]), self.coor_vec([all_length, iTrack+2*iGap]))
    self.gapObjects.append(gap)
    if self.is_mask:
        mask = self.draw_rect(self.name+"_mask", self.coor([-all_length/2,-iTrack/2-iGap-self.gap_mask]), self.coor_vec([all_length, iTrack+2*iGap+2*self.gap_mask]))
        self.maskObjects.append(mask)
        
    portOut1 = [self.coor([all_length/2,0]), self.ori, iTrack+2*self.overdev, iGap-2*self.overdev]
    portOut2 = [self.coor([-all_length/2,0]), -self.ori, iTrack+2*self.overdev, iGap-2*self.overdev]
    self.ports[self.name+'_1'] = portOut1
    self.ports[self.name+'_2'] = portOut2
    
    
def draw_halfJRM(self, iTrack, iGap, all_length, mode='litho', LS = '0.1nH',
                 Lind ='1nH', LJ='5nH'): #for now assume left and right tracks are the same width
    '''
    --------

    '''

    iTrack, iGap, all_length = parse_entry((iTrack, iGap, all_length))
    T = iTrack # Track
    A = all_length # legnth of halfJRM element top to bottom
    C = iGap # Gap. Width of halfJRM element is T+2*C
    L = C/3 # length of lumped inductor
    IJ = (C-L)/2 # Island of metal linking junctions to ground C = L+2*IJ
    CLR = (A-3*T)/4
    I = (A-2*CLR-3*T-2*L)/4 # Island of metal linking shunt junctions

    #snail array
#        if mode=='litho':
#            in_array = [self.coor([all_length/2, 0]), -self.ori, iTrackSnail, 0]
#            out_array = [self.coor([-all_length/2, 0]), self.ori, iTrackSnail, 0]      
#            snail_array = self.connect_elt(self.name+'_array', out_array, in_array)
#            snail_track = snail_array._connect_snails2([snail_dict['loop_width'], snail_dict['loop_length']], snail_dict['length_big_junction'], 4, snail_dict['length_small_junction'], 1, N_snails, snail_dict['bridge'], snail_dict['bridge_spacing'])
#        
#        if 0:
#            in_array = [self.coor([array_room/2, 0]), -self.ori, iTrackSnail, 0]
#            out_array = [self.coor([-array_room/2, 0]), self.ori, iTrackSnail, 0]      
#            snail_array = self.connect_elt(self.name+'_array', in_array, out_array)
#            snail_track = snail_array._connect_JJ(iTrack, iInduct=L_eq)

    raw_points_center = [(I+T/2, T/2),(-I-T/2, T/2),(-I-T/2, -T/2),(I+T/2, -T/2)]
    raw_points_right = [(3*I/2+L+T+I/2+T/2, T/2), (3*I/2+L+T-I/2-T/2, T/2),
                        (3*I/2+L+T-I/2-T/2, -T/2), (3*I/2+L+T+I/2+T/2, -T/2)]
    raw_points_left = self.refy_points(raw_points_right, absolute=True)
    raw_points_connector_top = [(T/2, T/2), (T/2, T/2+C), (-T/2, T/2+C),(-T/2, T/2)]
    raw_points_center1 =[(T/2, -T/2),(-T/2, -T/2),(-T/2, -T/2-IJ),(T/2, -T/2-IJ)]
    raw_points_center2 =self.move_points(raw_points_center1, [0, -IJ-L], absolute=True)
    raw_points_left1 =self.move_points(raw_points_center1, [-L-T-2*I, 0], absolute=True)
    raw_points_left2 =self.move_points(raw_points_center2, [-L-T-2*I, 0], absolute=True)
    raw_points_right1 =self.move_points(raw_points_center1, [L+T+2*I, 0], absolute=True)
    raw_points_right2 =self.move_points(raw_points_center2, [L+T+2*I, 0], absolute=True)

    x0 = -T/2-I-L-I-T
    raw_points_connector_left = [(x0, T/2),(x0-CLR, T/2),(x0-CLR, -T/2),(x0,-T/2)]
    raw_points_connector_right = self.refy_points(raw_points_connector_left, absolute=True)
    
    island_center = self.draw(self.name+"_center", self.append_absolute_points(raw_points_center))
    island_right = self.draw(self.name+"_right",self.append_absolute_points( raw_points_right))
    island_left = self.draw(self.name+"_left", self.append_absolute_points(raw_points_left))
    island_center1 = self.draw(self.name+"_center1", self.append_absolute_points(raw_points_center1))
    island_center2 = self.draw(self.name+"_center2", self.append_absolute_points(raw_points_center2))
    island_left1 = self.draw(self.name+"_left1", self.append_absolute_points(raw_points_left1))
    island_left2 = self.draw(self.name+"_left2", self.append_absolute_points(raw_points_left2))
    island_right1 = self.draw(self.name+"_right1", self.append_absolute_points(raw_points_right1))
    island_right2 = self.draw(self.name+"_right2", self.append_absolute_points(raw_points_right2))
    connector_top = self.draw(self.name+"_conn_top", self.append_absolute_points(raw_points_connector_top))
    connector_left = self.draw(self.name+"_conn_left", self.append_absolute_points(raw_points_connector_left))
    connector_right = self.draw(self.name+"_conn_right", self.append_absolute_points(raw_points_connector_right))
    
    if not self.is_litho:
        array_eqS1 = self.draw_rect(self.name+"_LS1", self.coor([T/2+I, -T/2]), self.coor_vec([L, T]))
        array_eqS2 = self.draw_rect(self.name+"_LS2", self.coor([-T/2-I-L, -T/2]), self.coor_vec([L, T]))
        array_eqL = self.draw_rect(self.name+"_L", self.coor([-T/2, -T/2-IJ]), self.coor_vec([T, -L]))
        array_eqJ1 = self.draw_rect(self.name+"_LJ1", self.coor([2*I+L+T-T/2, -T/2-IJ]), self.coor_vec([T, -L]))
        array_eqJ2 = self.draw_rect(self.name+"_LJ2", self.coor([-2*I-L-T-T/2, -T/2-IJ]), self.coor_vec([T, -L]))
        self.assign_lumped_RLC(array_eqS1, self.ori, (0, LS, 0))
        self.assign_lumped_RLC(array_eqS2, self.ori, (0, LS, 0))
        self.assign_lumped_RLC(array_eqL, self.ori, (0, Lind, 0))
        self.assign_lumped_RLC(array_eqJ1, self.ori, (0, LJ, 0))
        self.assign_lumped_RLC(array_eqJ2, self.ori, (0, LJ, 0))
        pointsS1 = [(I+T/2,0),(I+T/2+L,0)]
        self.draw(self.name+'_eq_lineS1', self.append_absolute_points(pointsS1), closed=False)
        pointsS2 = self.refy_points(pointsS1, absolute=True)
        self.draw(self.name+'_eq_lineS2', self.append_absolute_points(pointsS2), closed=False)
        pointsL = [(0,-T/2-IJ),(0,-T/2-IJ-L)]
        self.draw(self.name+'_eq_lineL', self.append_absolute_points(pointsL), closed=False)
        pointsJ1 = self.move_points(pointsL, [T+2*I+L, 0], absolute=True)
        self.draw(self.name+'_eq_lineJ1', self.append_absolute_points(pointsJ1), closed=False)
        pointsJ2 = self.move_points(pointsL, [-T-2*I-L, 0], absolute=True)
        self.draw(self.name+'_eq_lineJ2', self.append_absolute_points(pointsJ2), closed=False)
 

    connect = self.unite([island_center, island_right, island_left,
                          connector_top, island_center1, island_center2,
                          island_left1,island_left2,
                          island_right1,island_right2,
                          connector_left, connector_right])
    self.trackObjects.append(connect)
    

    gap = self.draw_rect(self.name+"_cutout", self.coor([-A/2,-T/2-C]), self.coor_vec([A, T+2*C]))
    self.gapObjects.append(gap)
#        if self.is_mask:
#            mask = self.draw_rect(self.name+"_mask", self.coor([-all_length/2,-iTrack/2-iGap-self.gap_mask]), self.coor_vec([all_length, iTrack+2*iGap+2*self.gap_mask]))
#            self.maskObjects.append(mask)
#            
    portT = [self.coor([0, T/2+C]), self.coor_vec([0,1]), iTrack+2*self.overdev, iGap-2*self.overdev]
    portL = [self.coor([-A/2,0]), self.coor_vec([-1,0]), iTrack+2*self.overdev, iGap-2*self.overdev]
    portR = [self.coor([A/2,0]), self.coor_vec([1,0]), iTrack+2*self.overdev, iGap-2*self.overdev]
    portB = [self.coor([0,-T/2-IJ-L-IJ]), self.coor_vec([0,-1]), iTrack+2*self.overdev, iGap-2*self.overdev]
    self.ports[self.name+'_T'] = portT
    self.ports[self.name+'_L'] = portL
    self.ports[self.name+'_R'] = portR
    self.ports[self.name+'_B'] = portB


def draw_T_join(self, iTrack, iGap):
    
    iTrack, iGap = parse_entry((iTrack, iGap))
    
    if not self.is_overdev or self.val(self.overdev<0):
        cutout = self.draw_rect_center(self.name+'_cutout', self.coor([0,self.overdev/2]), self.coor_vec([4*iGap+2*iTrack, 2*iGap+iTrack-self.overdev]))
        self.gapObjects.append(cutout)
    else:
        points = self.append_points([(-(iGap*2+iTrack),-iTrack/2-iGap+self.overdev),
                         (0, 2*iGap+iTrack-2*self.overdev),
                         ((iGap*2+iTrack)*2, 0),
                         (0, -(2*iGap+iTrack)+2*self.overdev),
                         (-self.overdev, 0),
                         (0, -self.overdev), 
                         (-iTrack*2-4*iGap+2*self.overdev, 0),
                         (0, self.overdev)])
        cutout = self.draw(self.name+'_cutout', points)
        self.gapObjects.append(cutout)
    
    if self.is_mask:
        mask = self.draw_rect(self.name+'_mask', self.coor([-iGap*2-iTrack,-iGap-iTrack/2]), self.coor_vec([4*iGap+2*iTrack, 2*iGap+iTrack+self.gap_mask]))
        self.maskObjects.append(mask)
        
    points = self.append_points([(-(iGap+iTrack/2)*2,-iTrack/2-self.overdev),
                                 (0, iTrack+2*self.overdev),
                                 ((iGap+iTrack/2)*4, 0),
                                 (0, -iTrack-2*self.overdev),
                                 (-iGap*2+self.overdev, 0),
                                 (0, -iGap+self.overdev), 
                                 (-iTrack*2-2*self.overdev, 0),
                                 (0, iGap-self.overdev)])
    track = self.draw(self.name+'_track', points)
    if self.val(iGap)<self.val(iTrack):
        fillet=iGap
    else:
        fillet=iTrack
    track.fillet(fillet-eps,[4,7])
    
    self.trackObjects.append(track)
    
    portOut1 = [self.coor([iTrack+iGap*2, 0]), self.coor_vec([1,0]), iTrack+2*self.overdev, iGap-2*self.overdev]
    self.ports[self.name+'_1'] = portOut1
    portOut2 = [self.coor([-(iTrack+iGap*2), 0]), self.coor_vec([-1,0]), iTrack+2*self.overdev, iGap-2*self.overdev]
    self.ports[self.name+'_2'] = portOut2
    portOut3 = [self.coor([0, -(iTrack/2+iGap)]), self.coor_vec([0,-1]), iTrack*2+2*self.overdev, iGap*2-2*self.overdev]
    self.ports[self.name+'_3'] = portOut3

    
def draw_fluxline(self, iTrack, iGap, length, track_flux, slope=0.5, sym='center', return_spacing=0, return_depth=0, opposite=False):
    is_fillet=True
    iTrack, iGap, length, track_flux, return_spacing, return_depth = parse_entry((iTrack, iGap, length, track_flux, return_spacing, return_depth))
    if sym == 'center':
        offset_gap_down = 0
        offset_track_down = 0
        offset_gap_up = 0
        offset_track_up = 0
        adapt_length = (iTrack/2-track_flux)/slope
    else:
        adapt_length = (iTrack/2-track_flux+length/2)/slope
        if sym=='down':
            offset_gap_down = iGap+2*track_flux
            offset_track_down = length/2+track_flux
            offset_gap_up = 0
            offset_track_up = 0
        elif sym=='up':
            offset_gap_down = 0
            offset_track_down = 0
            offset_gap_up = iGap+2*track_flux
            offset_track_up = length/2+track_flux
        else:
            raise ValueError("sym should be either True, 'down' or 'up'.")
    
    points = self.append_absolute_points([(adapt_length, iGap+iTrack/2),
                                          (adapt_length-track_flux, iGap+iTrack/2),
                                          (track_flux/2, length/2+offset_gap_up),
                                          (track_flux/2, -length/2-offset_gap_down),
                                          (adapt_length-track_flux, -(iGap+iTrack/2)),
                                          (adapt_length, -(iGap+iTrack/2))])
            

    if not sym == 'center':
        if return_spacing != 0:
            if not opposite:
                test1 = self.draw_rect(self.name+'_test1', self.coor([return_depth, length/2+offset_gap_up]), self.coor_vec([-(return_spacing+track_flux/2+return_depth), return_spacing+track_flux]))
            else:
                test1 = self.draw_rect(self.name+'_test1', self.coor([return_depth, -(length/2+offset_gap_up)]), self.coor_vec([-(return_spacing+track_flux/2+return_depth), -(return_spacing+track_flux)]))
            test1.fillet(return_spacing+track_flux+eps, 2)
            self.gapObjects.append(test1)                    
                
        
        
    gap = self.draw(self.name+'_gap', points)
    
    
    track = self.draw(self.name+'_track_guide', points, closed=False)
    gapext = self.draw(self.name+'_gapext_guide', points, closed=False)
    fillet = eps
    fillet2 = track_flux+eps
    if is_fillet:
        track.fillet(fillet2, [1,4])
        gapext.fillet(fillet2, [1,4])
        gap.fillet(fillet2,[1,4])
        
        track.fillet(fillet, [3,4])
        gapext.fillet(fillet, [3,4])
        gap.fillet(fillet,[3,4])

    points_starter = self.append_absolute_points([(adapt_length, iGap+iTrack/2),
                                          (adapt_length, iGap+iTrack/2+track_flux)])
    track_starter = self.draw(self.name+'_track', points_starter, closed=False)
    gapext_starter = self.draw(self.name+'_gapext', points_starter, closed=False)
    
    
    track = track.sweep_along_path(track_starter)
    gapext = gapext.sweep_along_path(gapext_starter)
    
    gap = self.unite([gap, gapext])
    
    if self.is_mask:
        mask = self.draw(self.name+'_mask', points)
        maskext = self.draw(self.name+'_mask_ext_guide', points, closed=False)
        points_starter_mask = self.append_absolute_points([(adapt_length, iGap+iTrack/2),
                                          (adapt_length, iGap+iTrack/2+self.gap_mask)])
        maskext_starter = self.draw(self.name+'_maskext', points_starter_mask, closed=False)
        maskext = maskext.sweep_along_path(maskext_starter)
        
        mask = self.unite([mask, maskext])
        self.maskObjects.append(mask)

    points = self.append_absolute_points([(adapt_length, iTrack/2),
                                          (adapt_length-track_flux, iTrack/2),
                                          (track_flux/2, track_flux-offset_track_down+offset_track_up),
                                          (track_flux/2, -track_flux-offset_track_down+offset_track_up),
                                          (adapt_length-track_flux, -iTrack/2),
                                          (adapt_length, -iTrack/2)])
              
    track_adapt = self.draw(self.name+'_track_adapt', points)
    if is_fillet:
        track_adapt.fillet(fillet2, [1,4])
    
    track = self.unite([track, track_adapt])
#            track.fillet(fillet, [15,20])
#            track.fillet(fillet2, [17,20])
    
    self.gapObjects.append(gap)
    self.trackObjects.append(track)

    portOut = [self.coor([adapt_length, 0]), self.coor_vec([1,0]), iTrack+2*self.overdev, iGap-2*self.overdev]
    self.ports[self.name] = portOut
    

def draw_double_fluxline(self, iTrack, iGap, length, track_flux, slope=0.7, sym='center', return_spacing=0, return_depth=0, opposite=False):
    is_fillet = True
    iTrack, iGap, length, track_flux, return_spacing, return_depth = parse_entry((iTrack, iGap, length, track_flux, return_spacing, return_depth))

    offset_gap_down = 0
    offset_track_down = 0
    offset_gap_up = 0
    offset_track_up = 0
    adapt_length = (iTrack/2-track_flux/2)/slope
    small_gap = length/4-3*track_flux/2/2
    big_length = 3*iTrack + 4*iGap
    
    points = self.append_absolute_points([(adapt_length, big_length/2),
                                          (adapt_length-track_flux,big_length/2),
                                          (track_flux/2, length/2+offset_gap_up),
                                          (track_flux/2, -length/2-offset_gap_down),
                                          (adapt_length-track_flux, -big_length/2),
                                          (adapt_length, -big_length/2)])
                  
                
    gap = self.draw(self.name+'_gap', points)
    
    track = self.draw(self.name+'_track_guide', points, closed=False)
    
    gapext = self.draw(self.name+'_gapext_guide', points, closed=False)


    points_starter = self.append_absolute_points([(adapt_length, big_length/2),
                                          (adapt_length, big_length/2+track_flux)])
    track_starter = self.draw(self.name+'_track', points_starter, closed=False)
    gapext_starter = self.draw(self.name+'_gapext', points_starter, closed=False)
    
    
    track = track.sweep_along_path(track_starter)
    gapext = gapext.sweep_along_path(gapext_starter)
    
    gap = self.unite([gap, gapext])
    
    if self.is_mask:
        mask = self.draw(self.name+'_mask', points)
        maskext = self.draw(self.name+'_mask_ext_guide', points, closed=False)
        points_starter_mask = self.append_absolute_points([(adapt_length, iGap+iTrack/2),
                                          (adapt_length, iGap+iTrack/2+self.gap_mask)])
        maskext_starter = self.draw(self.name+'_maskext', points_starter_mask, closed=False)
        maskext = maskext.sweep_along_path(maskext_starter)
        
        mask = self.unite([mask, maskext])
        self.maskObjects.append(mask)
        
    raw_points_top = [(adapt_length,big_length/2-iGap),
                      (adapt_length-track_flux, big_length/2-iGap),
                      (track_flux/2, small_gap+3*track_flux/2),
                      (track_flux/2, small_gap+3*track_flux/2-2/2*track_flux),
                      (adapt_length-track_flux, big_length/2-iGap-iTrack),
                      (adapt_length, big_length/2-iGap-iTrack)]
        
    raw_points_middle = [(adapt_length, iTrack/2),
                      (adapt_length-track_flux, iTrack/2),
                      (track_flux/2, track_flux/2-offset_track_down+offset_track_up),
                      (track_flux/2, -track_flux/2-offset_track_down+offset_track_up),
                      (adapt_length-track_flux, -iTrack/2),
                      (adapt_length, -iTrack/2)]
    
    raw_points_bottom = self.refx_points(raw_points_top)

    points_top = self.append_absolute_points(raw_points_top)
    points_middle = self.append_absolute_points(raw_points_middle)
    points_bottom = self.append_absolute_points(raw_points_bottom)
              
    track_adapt_top = self.draw(self.name+'_track_adapt_top', points_top)
    track_adapt_middle = self.draw(self.name+'_track_adapt_middle',
                                   points_middle)
    track_adapt_bottom = self.draw(self.name+'_track_adapt_bottom',
                                   points_bottom)
    

    
    track = self.unite([track, track_adapt_top,
                        track_adapt_middle,
                        track_adapt_bottom])
#            track.fillet(fillet, [15,20])
#            track.fillet(fillet2, [17,20])

    fillet2 = 2*track_flux/2+eps
    if is_fillet:
        track.fillet(fillet2, [1, 2, 3, 4, 7,8,9,10,13,14,15,16,19,20,21,22,25,26,27,28])
        gap.fillet(fillet2, [2,3,4, 5])
    
    self.gapObjects.append(gap)
    self.trackObjects.append(track)

    portOut_top = [self.coor([adapt_length, iTrack/2+iGap+iTrack/2]), self.coor_vec([1,0]), iTrack+2*self.overdev, iGap-2*self.overdev]
    portOut_bottom = [self.coor([adapt_length, -iTrack/2-iGap-iTrack/2]), self.coor_vec([1,0]), iTrack+2*self.overdev, iGap-2*self.overdev]
    self.ports[self.name+'_top'] = portOut_top
    self.ports[self.name+'_bottom'] = portOut_bottom
 


def draw_alignement_mark(self, iSize, iXdist, iYdist):
    iXdist, iYdist, iSize=parse_entry((iXdist, iYdist, iSize))
    raw_points = [(iXdist,iYdist),
                  (iSize/2,iSize/2),
                  (-iSize,0)]
    mark1=self.draw(self.name+"_mark_a", self.append_points(raw_points))
    
    raw_points = [(iXdist,iYdist),
                          (-iSize/2,-iSize/2),
                          (iSize,0)]
    mark2=self.draw(self.name+"_mark_b", self.append_points(raw_points))
    self.gapObjects.append(self.unite([mark1,mark2]))
    
def draw_alignement_mark_r(self, size, disp, suff=''):
    size, disp = parse_entry((size, disp))
    raw_points = [(*disp,),
                  (size/2,size/2),
                  (0,-size)]        
    marks = []
    marks.append(self.draw(self.name+'_'+suff+"a", self.append_points(raw_points)))
    marks.append(self.draw(self.name+'_'+suff+"b", self.append_points(self.refy_points(raw_points, disp[0]))))
    self.gapObjects.append(self.unite([*marks]))
    if self.is_mask:
        raw_points_mask = [(disp[0]-self.gap_mask,disp[1]),
                           (size/2+2*self.gap_mask,size/2+self.gap_mask*2),
                           (0,-size-self.gap_mask*4)]
        marks_mask = []
        marks_mask.append(self.draw(self.name+suff+"a_mask", self.append_points(raw_points_mask)))
        marks_mask.append(self.draw(self.name+suff+"b_mask", self.append_points(self.refy_points(raw_points_mask, disp[0]))))
        self.maskObjects.append(self.unite(marks_mask))
        
def draw_alignement_marks(self, size, disp, dict_except=None):
    size, disp = parse_entry((size, disp))
    disp = Vector(disp)
    moves = []
    directions_name = ['NE', 'NO', 'SE', 'SO']        
    directions_exp = [[1,1], [-1,1], [1,-1], [-1,-1]]
    if dict_except is not None:
        for direction_name in directions_name:
            try:
                moves.append(parse_entry(dict_except[direction_name]))
            except Exception:
                moves.append([0,0])
    else:
        moves = [[0,0] for ii in range(4)]
    
    for move, direction_exp, direction_name in zip(moves, directions_exp, directions_name):
        self.draw_alignement_mark_r(size, disp*direction_exp+move, suff=direction_name)      
    
def draw_dose_test_Nb(self, pad_size, pad_spacing, array):
    pad_size, pad_spacing = parse_entry((pad_size, pad_spacing))
    pad_size = Vector(pad_size)
    
    pads = []
    pos_x = pad_spacing
    pos_y = pad_spacing
    for jj in range(array[1]):
        for ii in range(array[0]):
            pads.append(self.draw_rect(self.name+'_%d_%d'%(ii, jj), self.coor([pos_x, pos_y]), self.coor_vec(pad_size)))
            pos_x += pad_size[0]+pad_spacing
            pads.append(self.draw_rect(self.name+'_%d_%d'%(ii, jj), self.coor([pos_x, pos_y]), self.coor_vec(pad_size)))
            pos_x += pad_size[0]+2*pad_spacing
        pos_y += pad_size[1]+pad_spacing
        pos_x = pad_spacing
            
    cutout = self.draw_rect(self.name+'_cutout', self.coor([0, 0]), self.coor_vec([pad_size[0]*2*array[0]+pad_spacing*3*array[0],pad_size[1]*array[1]+pad_spacing*(array[1]+1)]))
    if self.is_mask:
        mask = self.draw_rect(self.name+'_cutout', self.coor([-self.gap_mask, -self.gap_mask]), self.coor_vec([pad_size[0]*2*array[0]+pad_spacing*3*array[0]+2*self.gap_mask,pad_size[1]*array[1]+pad_spacing*(array[1]+1)+2*self.gap_mask]))
        self.maskObjects.append(mask)
    cutout = self.subtract(cutout, pads)
    self.gapObjects.append(cutout)

def draw_dose_test(self, pad_size, pad_spacing, iTrack, bridge, N, bridge_spacing, length_big_junction, length_small_junction):
    pad_size, pad_spacing, iTrack, bridge, bridge_spacing, length_big_junction, length_small_junction = parse_entry((pad_size, pad_spacing, iTrack, bridge, bridge_spacing, length_big_junction, length_small_junction))
    pad_size = Vector(pad_size)
    
    self.draw_rect(self.name+'_left', self.coor([-pad_spacing/2, -pad_size[1]/2]), self.coor_vec([-pad_size[0],pad_size[1]]))
    self.draw_rect(self.name+'_right', self.coor([pad_spacing/2, pad_size[1]/2]), self.coor_vec([pad_size[0],-pad_size[1]]))
    
    portOut1 = [self.coor([pad_spacing/2, 0]), -self.ori, iTrack, 0]
    portOut2 = [self.coor([-pad_spacing/2, 0]), self.ori, iTrack, 0]
    self.ports[self.name+'_1'] = portOut1
    self.ports[self.name+'_2'] = portOut2
    
    in_array = portOut2
    out_array = portOut1
    snail_array = self.connect_elt(self.name+'_junction', in_array, out_array)
    snail_track = snail_array._connect_snails2([20e-6,20e-6], length_big_junction, 3, length_small_junction, 1, N, bridge, bridge_spacing)#(squid_size, width_top, n_top, width_bot, n_bot, N, width_bridge)

def draw_dose_test_snail_Zaki(self, pad_size, pad_spacing, iTrack, N_snails=1,
                              snail_dict={'loop_width':20e-6, 'loop_length':20e-6,
                                            'length_big_junction':10e-6,
                                            'length_small_junction':2e-6, 'bridge':1e-6,
                                            'bridge_spacing':1e-6}):
    pad_size, pad_spacing, iTrack = parse_entry((pad_size, pad_spacing, iTrack))
    pad_size = Vector(pad_size)
    
    self.draw_rect(self.name+'_left', self.coor([-pad_spacing/2, -pad_size[1]/2]), self.coor_vec([-pad_size[0],pad_size[1]]))
    self.draw_rect(self.name+'_right', self.coor([pad_spacing/2, pad_size[1]/2]), self.coor_vec([pad_size[0],-pad_size[1]]))
    
    portOut1 = [self.coor([pad_spacing/2, 0]), -self.ori, iTrack, 0]
    portOut2 = [self.coor([-pad_spacing/2, 0]), self.ori, iTrack, 0]
    self.ports[self.name+'_1'] = portOut1
    self.ports[self.name+'_2'] = portOut2
    
    in_array = portOut2
    out_array = portOut1
    snail_array = self.connect_elt(self.name+'_junction', in_array, out_array)
    snail_track = snail_array._connect_snails_Zaki([snail_dict['loop_width'], snail_dict['loop_length']], snail_dict['length_big_junction'], 3, snail_dict['length_small_junction'], 1, N_snails, snail_dict['bridge'], snail_dict['bridge_spacing'])#(squid_size, width_top, n_top, width_bot, n_bot, N, width_bridge)
    
def draw_dose_test_junction(self, pad_size, pad_spacing, width, width_bridge, n_bridge=1, spacing_bridge=0, alternate_width=True):
    pad_size, pad_spacing,width, spacing_bridge, width_bridge = parse_entry((pad_size, pad_spacing, width, spacing_bridge, width_bridge))
    pad_size = Vector(pad_size)
    if self.val(width)<1.5e-6 and n_bridge==1:
        width_jct = width
        width = 1.5e-6
    else: 
        width_jct=None
    
    self.draw_rect(self.name+'_left', self.coor([-pad_spacing/2, -pad_size[1]/2]), self.coor_vec([-pad_size[0],pad_size[1]]))
    self.draw_rect(self.name+'_right', self.coor([pad_spacing/2, pad_size[1]/2]), self.coor_vec([pad_size[0],-pad_size[1]]))
    
    portOut1 = [self.coor([pad_spacing/2, 0]), -self.ori, width, 0]
    portOut2 = [self.coor([-pad_spacing/2, 0]), self.ori, width, 0]
    self.ports[self.name+'_1'] = portOut1
    self.ports[self.name+'_2'] = portOut2
    
    in_array = portOut2
    out_array = portOut1
    jcts = self.connect_elt(self.name+'_junction', in_array, out_array)
    print(n_bridge)
    if alternate_width:
        jcts._connect_jct(width_bridge, n=n_bridge, spacing_bridge=spacing_bridge, width_jct=width_jct)
    else:
        jcts._connect_jct(width_bridge, n=n_bridge, spacing_bridge=spacing_bridge, assymetry=0, width_jct=width_jct)

def draw_cav_Si(self, size, chip_height, tunnel_length, size_plot=None, height_plot=None, n=8, track='42um', gap='25um', ports_senw = ['500um',0,0,0], margin_litho='10um'):
    # can only make square cavity
    # size : intra size of the cavity : 4mm->15GHz
    # chip_heigth : silicon height
    # tunnel length : length of the silicon tunnel between 2 via
    # size_plot : base size of the central plot
    # height_plot : thickness of the capacitance at the top of the plot
    # n : nb of vias
    # track, gap : if we have an access port, track and gap of this port 
    #   currently all port have same gap and track
    # ports_senw = list giving the length of the cable coming from each port ; south, east, north, west
    #   0 means there should be no port
    # ports will be numeroted from 1 to n in the trigo way
    # margin_litho : width of the Si3N4 mask between the vias
    
    ports_senw_bool = [bool(ii) for ii in ports_senw]
    size, chip_height, tunnel_length, track, gap, margin_litho, ports_senw = parse_entry((size, chip_height, tunnel_length, track, gap, margin_litho, ports_senw))
#        self.draw_rect(self.name+'_rect', self.coor([-size/2, -size/2]), self.coor_vec([size,size]))
    access_width = track+2*gap+2*margin_litho
    if size_plot is None:
        size_plot=size/2
    else:
        size_plot = parse_entry((size_plot))
    
    if height_plot is None:
        height_plot=chip_height/2
    else:
        height_plot = parse_entry((height_plot))
        
    # bottom row
    width_trapeze = (size-(n+1)*margin_litho)/n
    track = '42um'
    gap = '25um'
    name_s, name_e, name_n, name_w = self.name+'_s', self.name+'_e', self.name+'_n', self.name+'_w'
    port_nb = 0
    for ii, name in enumerate([name_s, name_e, name_n, name_w]):
        if ports_senw_bool[ii]:
            port_nb += 1 
            _width_trapeze = width_trapeze-access_width/n
            x_pos = -size/2+margin_litho+_width_trapeze/2
            _n = n//2
            self.draw_trapeze(name, self.coor([x_pos, -size/2-tunnel_length/2]), 0, self.coor_vec([_width_trapeze, tunnel_length]), -chip_height/2)
            self.draw_trapeze(name+'_'+str('bis'), self.coor([x_pos, -size/2-tunnel_length/2]), -chip_height, self.coor_vec([_width_trapeze, tunnel_length]), chip_height/2)
            self.unite([name, name+'_'+str('bis')])

            self.duplicate_along_line(name, self.coor_vec([_width_trapeze+margin_litho, 0]), n=_n)
            unite_list = [name]+[name+'_'+str(ii+1) for ii in range(_n-1)]
            self.unite(unite_list)
            
            self.duplicate_along_line(name, self.coor_vec([_n*(_width_trapeze+margin_litho)+access_width, 0]), n=2)
            unite_list = [name, name+'_'+str(_n)]
            self.unite(unite_list)
            
            portOut = [self.coor([0, -size/2-tunnel_length]), self.coor_vec([0,-1]), track, gap]
            self.ports[self.name+'_port_'+str(port_nb)] = portOut
            
            end = self.key_elt(self.name+'_end_'+str(ii), self.coor([0, -size/2-tunnel_length+ports_senw[ii]]), self.coor_vec([0,-1]))
            end.draw_end_cable(track, gap, fillet=gap)
            cable = self.connect_elt(self.name+'_cable_'+str(ii), self.name+'_end_'+str(ii), self.name+'_port_'+str(port_nb))
            cable.draw_cable(is_bond=False)
        else:
            x_pos = -size/2+margin_litho+width_trapeze/2
    
            self.draw_trapeze(name, self.coor([x_pos, -size/2-tunnel_length/2]), 0, self.coor_vec([width_trapeze, tunnel_length]), -chip_height/2)
            self.draw_trapeze(name+'_'+str('bis'), self.coor([x_pos, -size/2-tunnel_length/2]), -chip_height, self.coor_vec([width_trapeze, tunnel_length]), chip_height/2)
            self.unite([name, name+'_'+str('bis')])
            
            self.duplicate_along_line(name, self.coor_vec([width_trapeze+margin_litho, 0]), n=n)
            unite_list = [name]+[name+'_'+str(ii+1) for ii in range(n-1)]
            self.unite(unite_list)
            
        self.ori = self.ori.rot(Vector([0,1]))

    name_c = self.name+'_c'
    self.draw_trapeze(name_c, self.coor([-size/2-tunnel_length/2, -size/2-tunnel_length/2]), 0, self.coor_vec([tunnel_length, tunnel_length]), -chip_height/2)
    self.draw_trapeze(name_c+'_bis', self.coor([-size/2-tunnel_length/2, -size/2-tunnel_length/2]), -chip_height, self.coor_vec([tunnel_length, tunnel_length]), chip_height/2)
    self.unite([name_c, name_c+'_bis'])
    
    self.duplicate_along_line(name_c, self.coor_vec([size+tunnel_length, 0]), n=2)
    unite_list = [name_c, name_c+'_1']
    self.unite(unite_list)
    
    self.duplicate_along_line(name_c, self.coor_vec([0, size+tunnel_length]), n=2)
    unite_list = [name_c, name_c+'_2']
    self.unite(unite_list)
    
    name_p = name+'_p'
    self.draw_trapeze(name_p, self.coor([0, 0]), -chip_height, self.coor_vec([size_plot, size_plot]), chip_height-height_plot)


    unite_list = [name_s, name_e, name_n, name_w, name_c, name_p]
    self.unite(unite_list)
    self.name = name_s
    
        
def draw_dc_Nrect(self, layer_name, length, rel_pos, widths, border='10um'):
    '''
    REL_POS HAS TO BE A LIST
    '''

    if isinstance(widths, list):
        if len(widths) is not len(rel_pos):
            raise ValueError('widths and rel_pos do not have same length')
    elif isinstance(widths, str):
        widths = [widths for ii in range(len(rel_pos))]
    else:
        raise ValueError('widths should be a string or a list of strings')
    mult = len(rel_pos)
    
    rel_pos, length, widths, border = parse_entry((rel_pos, length, widths, border))
    
    rect = []
    for ii in range(mult):
        rect.append(self.draw_rect(layer_name+'_'+self.name+'_track_'+str(ii), \
                                    self.coor([-length/2, rel_pos[ii]-widths[ii]/2]), \
                                    self.coor_vec([length, widths[ii]])))
    if len(rect) > 1:
        rect = self.unite(rect, name=layer_name+'_'+self.name+'_track')
    else:
        rect = self.rename(rect[0], layer_name+'_'+self.name+'_track')
    self.layers[layer_name]['trackObjects'].append(rect)
    
    pos_cutout, width_cutout, vec_cutout = self.size_dc_gap(length, rel_pos, widths, border)
    cutout = self.draw_rect(layer_name+'_'+self.name+'_cutout', pos_cutout, vec_cutout)
    self.layers[layer_name]['gapObjects'].append(cutout)
    
    portOut1 = [self.coor([length/2, 0]), self.ori, rel_pos, widths, width_cutout, mult] #portOut 
    portOut2 = [self.coor([-length/2, 0]), -self.ori, rel_pos[::-1], widths, width_cutout, mult] #portIn
    self.ports_dc[self.name+'_1'] = portOut1
    self.ports_dc[self.name+'_2'] = portOut2
    
def draw_dc_test_Nrect(self, layer_name, N, length, rel_pos, widths, border='10um'):
    '''
<<<<<<< HEAD
    Poorly coded: rel_pos is the center of one edge but we use draw_rect instead of draw_rect_center
=======
>>>>>>> simplified routines draw rect
    REL_POS HAS TO BE A LIST
    '''

    if isinstance(widths, list):
        if len(widths) is not len(rel_pos):
            raise ValueError('widths and rel_pos do not have same length')
    elif isinstance(widths, str):
        widths = [widths for ii in range(len(rel_pos))]
    else:
        raise ValueError('widths should be a string or a list of strings')
    mult = len(rel_pos)
    
    rel_pos, length, widths, border = parse_entry((rel_pos, length, widths, border))
    
    rect = []
    for kk in range(N):
        y0 = parse_entry((str(kk*100)+'um'))
        for ii in range(mult):
            rect.append(self.draw_rect(layer_name+'_'+self.name+'_track_'+str(ii), \
                                        self.coor([-length/2, y0+rel_pos[ii]-widths[ii]/2]), \
                                        self.coor_vec([length, widths[ii]])))
    if len(rect) > 1:
        rect = self.unite(rect, name=layer_name+'_'+self.name+'_track')
    else:
        rect = self.rename(rect[0], layer_name+'_'+self.name+'_track')
    self.layers[layer_name]['trackObjects'].append(rect)
    
#        pos_cutout, width_cutout, vec_cutout = self.size_dc_gap(length, rel_pos, widths, border)
#        cutout = self.draw_rect(layer_name+'_'+self.name+'_cutout', pos_cutout, vec_cutout)
#        self.layers[layer_name]['gapObjects'].append(cutout)
    
#        portOut1 = [self.coor([length/2, 0]), self.ori, rel_pos, widths, width_cutout, mult] #portOut 
#        portOut2 = [self.coor([-length/2, 0]), -self.ori, list(-np.array(rel_pos)), widths, width_cutout, mult] #portIn
#        self.ports_dc[self.name+'_1'] = portOut1
#        self.ports_dc[self.name+'_2'] = portOut2
    
    
def draw_dc_pad(self, layer_name, iTrack, iGap, xlength='250um', ylength='250um'):
    
    iTrack, iGap, xlength, ylength = parse_entry((iTrack, iGap, xlength, ylength))
    
    pad = self.draw_rect('pad',\
                         self.coor([-xlength/2, -ylength/2]),\
                         self.coor_vec([xlength, ylength]))
    padout = self.draw_rect('padout',\
                            self.coor([xlength/2, -iTrack/2]), \
                            self.coor_vec([2*iGap, iTrack]))
    pad = self.unite([pad, padout], name=layer_name+'_'+self.name+'_track')
    pad_gap = self.draw_rect(layer_name+'_'+self.name+'_gap',\
                             self.coor([-xlength/2-2*iGap, -ylength/2-2*iGap]), \
                             self.coor_vec([xlength+4*iGap, ylength+4*iGap]))

    self.layers[layer_name]['trackObjects'].append(pad)
    self.layers[layer_name]['gapObjects'].append(pad_gap)

    portOut = [self.coor([xlength/2+2*iGap, 0]), self.ori, iTrack, iGap]
    self.ports[self.name] = portOut
    
    
def draw_dc_pads(self, layer_name, rel_pos, widths, gaps, xlength='250um', ylength='250um'):

    if isinstance(widths, list):
        if len(widths) is not len(rel_pos):
            raise ValueError('width and rel_pos lists do not have same length')
    elif isinstance(widths, str):
        widths = [widths for ii in range(len(rel_pos))]
    else:
        raise ValueError('width should be a string or a list of strings')
        
    if isinstance(gaps, list):
        if len(gaps) is not len(rel_pos):
            raise ValueError('gap and rel_pos lists do not have same length')
    elif isinstance(gaps, str):
        gaps = [gaps for ii in range(len(rel_pos))]
    else:
        raise ValueError('gap should be a string or a list of strings')
        
    mult = len(rel_pos)
    
    rel_pos, widths, gaps, xlength, ylength = parse_entry((rel_pos, widths, gaps, xlength, ylength))
    
    pads = []
    pad_gaps = []
    for ii in range(mult):
        pad = self.draw_rect('pad'+str(ii), \
                             self.coor([-xlength/2, rel_pos[ii]-ylength/2]), \
                             self.coor_vec([xlength, ylength]))
        padout = self.draw_rect('pad_out'+str(ii), \
                                self.coor([xlength/2, rel_pos[ii]-widths[ii]/2]), \
                                self.coor_vec([4*gaps[ii], widths[ii]+2*gaps[ii]]))
        pads.append(self.unite([pad, padout]))
        pad_gap = self.draw_rect('gap'+str(ii), \
                                 self.coor([-xlength/2-3*gaps[ii], rel_pos[ii]-ylength/2-3*gaps[ii]]), \
                                 self.coor_vec([xlength+6*gaps[ii], ylength+6*gaps[ii]]))
        pad_gapout = self.draw_rect('gap_out'+str(ii), \
                                    self.coor([xlength/2+3*gaps[ii], rel_pos[ii]-widths[ii]/2-gaps[ii]]), \
                                    self.coor_vec([gaps[ii], widths[ii]+2*gaps[ii]]))
        pad_gaps.append(self.unite([pad_gap, pad_gapout]))
    if len(pads) > 1:
        self.unite(pads, name=layer_name+'_'+self.name+'_track')
        self.unite(pad_gaps, name=layer_name+'_'+self.name+'_gap')
    else:
        self.rename(pads[0], layer_name+'_'+self.name+'_track')
        self.rename(pad_gaps[0], layer_name+'_'+self.name+'_gap')
    self.layers[layer_name]['trackObjects'].append(pads[0])
    self.layers[layer_name]['gapObjects'].append(pad_gaps[0])
    
    cutout_list = [widths[ii]+2*gaps[ii] for ii in range(mult)]
    track_list = [widths[ii] for ii in range(mult)]

    portOut = [self.coor([xlength/2+4*gaps[ii], 0]), self.ori, cutout_list, rel_pos, track_list, mult]
    self.ports_dc[self.name] = portOut
  
def draw_alignment_marks(self, layer_name, isLayer63=True):
    
    w_thin, l_thin, w_large, l_large, square = '0.1um', '4um', '1um', '3um', '12um'
    w_thin, l_thin, w_large, l_large, square = parse_entry((w_thin, l_thin, w_large, l_large, square))
    marks = []
    if isLayer63 is True:
        self.new_layer('layer63')
        squares = []
    count = -1
    for x0 in ['-42um', '42um']:
        count += 1
        for y0 in ['-42um', '42um']:
            count += 1
            mark = []
            x0, y0 = parse_entry((x0, y0))
            mark.append(self.draw_rect('thin_x', self.coor([x0-l_thin/2, y0-w_thin/2]), self.coor_vec([l_thin, w_thin])))
            mark.append(self.draw_rect('thin_y', self.coor([x0-w_thin/2, y0-l_thin/2]), self.coor_vec([w_thin, l_thin])))
            mark.append(self.draw_rect('large_right', self.coor([x0+l_thin/2, y0-w_large/2]), self.coor_vec([l_large, w_large])))
            mark.append(self.draw_rect('large_left', self.coor([x0-l_thin/2, y0-w_large/2]), self.coor_vec([-l_large, w_large])))
            mark.append(self.draw_rect('large_top', self.coor([x0-w_large/2, y0+l_thin/2]), self.coor_vec([w_large, l_large])))
            mark.append(self.draw_rect('large_bot', self.coor([x0-w_large/2, y0-l_thin/2]), self.coor_vec([w_large, -l_large])))
            mark = self.unite(mark, name=layer_name+'_'+self.name+'_alignement_mark_'+str(count))
            marks.append(mark)
            if isLayer63 is True:
                squares.append(self.draw_rect_center('layer63', self.coor([x0, y0]), self.coor_vec([square, square])))
                
    marks = self.unite(marks, name=layer_name+'_'+self.name+'_alignement_mark')
    self.layers[layer_name]['trackObjects'].append(marks)
    if isLayer63 is True:
        squares = self.unite(squares, name='layer63_'+self.name)
        
    
def draw_pits(self, rel_pos, length, width):
    '''
    rel_pos gives the position of the fine structure
    length is the length along the CNT region
    width is in the transverse direction, spans the whole comb
    '''
    self.new_layer('pits')
    rel_pos, length, width = parse_entry((rel_pos, length, width))
    
    print(rel_pos)
    dist2CNT = max([abs(pos) for pos in rel_pos])
    dist2CNT += parse_entry('10um')
    space = parse_entry('40um')
    close_pit_left = self.draw_rect('close_pit_left', \
                                    self.coor([-length/2, dist2CNT]), \
                                    self.coor_vec([length, space]))
    close_pit_right = self.draw_rect('close_pit_right', \
                                     self.coor([-length/2, -dist2CNT]), \
                                     self.coor_vec([length, -space]))
    pit_left = self.draw_rect('pit_left', \
                              self.coor([-length, dist2CNT+space]), \
                              self.coor_vec([2*length, width]))
    pit_right = self.draw_rect('pit_right', \
                                self.coor([-length, -dist2CNT-space]), \
                                self.coor_vec([2*length, -width]))
    pit_left = self.unite([pit_left, close_pit_left], name='pits_'+self.name+'_left')
    pit_right = self.unite([pit_right, close_pit_right], name='pits_'+self.name+'_right')
    self.layers['pits']['trackObjects'].append(pit_left)
    self.layers['pits']['trackObjects'].append(pit_right)
    
 