# -*- coding: utf-8 -*-
"""
Created on Mon Oct 28 16:27:24 2019

@author: Zaki
"""

from designer import Vector, eps
import numpy as np
from hfss import parse_entry
#from PythonModeler import PythonModeler
from KeyElement import move
from functools import wraps
import Lib

__methods__ = []
register_method = Lib.register_method(__methods__)


@register_method
def draw_JJ(self, name, iInPort, iOutPort, iInduct='1nH', fillet=None):
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
#    iTrack, iGap, iTrackJ, iLength = parse_entry((iTrack, iGap, iTrackJ, iLength))
#    portOut1 = self.port(name+'_1',[iLength/2,0], [1,0], iTrack, iGap)
#    portOut2 = self.port(name+'_2',[-iLength/2,0], [-1,0], iTrack, iGap)
    
    pads = self.network._connect_JJ(name+'JJ', name+'_1', name+'_2', iInduct)
#    self.trackObjects.append(pads)
    gap = self.rect_center_2D([0,0], [iLength, iTrack+2*iGap], name=name+'_gap', layer='GAP')
    return gap, portOut1, portOut2, pads

@register_method
@move
def draw_IBM_transmon(self,
                   name, 
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
    
    cutout = self.rect_center_2D([0,0], cutout_size-Vector([2*self.overdev, 2*self.overdev]), name=name+"_cutout", layer="GAP")

    if self.is_mask:
        mask = self.rect_center_2D([0,0], cutout_size+Vector([track,track]),name=name+"_mask",layer="MASK")
        self.maskObjects.append(mask)
    if not self.is_litho:
        mesh = self.rect_center_2D([0,0], cutout_size, name=name+"_mesh", layer="MESH")
        self.mesh_zone(mesh, 2*track)

    right_pad = self.rect_center_2D(Vector(pad_spacing+pad_size[0],0)/2, pad_size+Vector([2*self.overdev, 2*self.overdev]), name=name+"_pad1", layer="TRACK")
    left_pad = self.rect_center_2D(-Vector(pad_spacing+pad_size[0],0)/2, pad_size+Vector([2*self.overdev, 2*self.overdev]), name=name+"_pad2", layer="TRACK")
    

    track_J=Jwidth*4.
    in_junction = self.port(name+'_in_jct', [-pad_spacing/2+self.overdev, 0], [1,0], track_J+2*self.overdev, 0)
    out_junction = self.port(name+'_out_jct', [pad_spacing/2-self.overdev, 0], [1,0], track_J+2*self.overdev, 0)
#        self.ports[name+'_in_jct'] = in_junction
#        self.ports[name+'_out_jct'] = out_junction
#        junction = self.connect_elt(name+'_junction', name+'_in_jct', name+'_out_jct')
    junction_pad = self.network._connect_JJ(name+'JJ', name+'_in_jct', name+'_out_jct', Jwidth+2*self.overdev, iInduct=Jinduc, fillet=None)
    pads = self.unite([right_pad, left_pad, junction_pad], name=name+'_pads')
    return pads, cutout, in_junction, out_junction

    
    if fillet is not None:
        self._fillet(fillet+self.overdev,[19], pads)
        self._fillet(0.75*fillet-self.overdev,[18,17,14,13], pads)
        self._fillet(fillet+self.overdev,[12,11,10,9,8,7], pads)
        self._fillet(0.75*fillet-self.overdev,[6,5,2,1], pads)
        self._fillet(fillet+self.overdev,0, pads)
    if nport==1:
        self.port(name+'_1',cutout_size.px()/2, self.ori, track, gap)
        
@register_method
@move
def draw_ZR_transmon(self,
                    name,
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

    cutout = self.rect_center_2D([0,0], cutout_size-Vector([self.overdev*2, self.overdev*2]),name =name+"_cutout", layer = "GAP")
    if not self.is_litho:
        mesh = self.rect_center_2D([0,0], cutout_size, name=name+"_mesh", layer="MESH")
    if self.is_mask:
        mask = self.rect_center_2D([0,0], cutout_size+Vector([self.gap_mask,self.gap_mask])*2, name=name+"_mask", layer="MASK")
    
    track_J=Jwidth*4.
    in_junction = self.port(name+'_in_jct', [-pad_spacing/2+self.overdev, 0], [1,0], track_J+2*self.overdev, 0)
    out_junction = self.port(name+'_out_jct', [+pad_spacing/2-self.overdev, 0], [-1,0], track_J+2*self.overdev, 0)
    
#        self.ports[name+'_in_jct'] = in_junction
#        self.ports[name+'_out_jct'] = out_junction
#        junction = self.connect_elt(name+'_junction', name+'_in_jct', name+'_out_jct')
    junction_pads = self.network._connect_JJ(name+'JJ', name+'_in_jct', name+'_out_jct', Jwidth+2*self.overdev, iInduct=Jinduc, fillet=None)
    
    raw_points = [(pad_spacing/2-self.overdev, -pad_size_right[1]/2-self.overdev),
                  (pad_size_right[0]+2*self.overdev, 0),
                  (0, pad_size_right[1]/2-spacing_right-short_right-gap_right-track_right/2+2*self.overdev),
                  (-pad_size_right[0]+length_right, 0),
                  (0, (spacing_right+short_right+gap_right+track_right/2)*2-2*self.overdev),
                  (pad_size_right[0]-length_right, 0),
                  (0, pad_size_right[1]/2-spacing_right-short_right-gap_right-track_right/2+2*self.overdev),
                  (-pad_size_right[0]-2*self.overdev, 0)]
    points = self.append_points(raw_points)
    right_pad = self.polyline_2D(points, name=name+"_pad1", layer="TRACK") 
        
    raw_points = [(-pad_spacing/2+self.overdev, -pad_size_left[1]/2-self.overdev),
                  (-pad_size_left[0]-2*self.overdev, 0),
                  (0, pad_size_left[1]/2-spacing_left-short_left-gap_left-track_left/2+2*self.overdev),
                  (pad_size_left[0]-length_left, 0),
                  (0, (spacing_left+short_left+gap_left+track_left/2)*2-2*self.overdev),
                  (-pad_size_left[0]+length_left, 0),
                  (0, pad_size_left[1]/2-spacing_left-short_left-gap_left-track_left/2+2*self.overdev),
                  (pad_size_left[0]+2*self.overdev, 0)]
    points = self.append_points(raw_points)
    left_pad = self.polyline_2D(points, name=name+"_pad2", layer="TRACK") 
    
    right_track_raw_points = [(cutout_size[0]/2-self.overdev,-track_right/2-self.overdev),
                             (-(cutout_size[0]/2-pad_spacing/2-length_right-gap_right-short_right-spacing_right), 0),
                             (0, track_right+2*self.overdev),
                             ((cutout_size[0]/2-pad_spacing/2-length_right-gap_right-short_right-spacing_right), 0)]
    right_track_points = self.append_points(right_track_raw_points)
    right_track = self.polyline_2D(right_track_points, name=name+"_track1", layer="TRACK") 
    
    left_track_raw_points = [(-cutout_size[0]/2+self.overdev,-track_left/2-self.overdev),
                             ((cutout_size[0]/2-pad_spacing/2-length_left-gap_left-short_left-spacing_left), 0),
                             (0, track_left+2*self.overdev),
                             (-(cutout_size[0]/2-pad_spacing/2-length_left-gap_left-short_left-spacing_left), 0)]
    left_track_points = self.append_points(left_track_raw_points)
    left_track = self.polyline_2D(left_track_points, name=name+"_track2", layer="TRACK") 
    
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
        right_short = self.polyline_2D(points, name=name+"_short1", layer="TRACK") 
        
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
        left_short = self.polyline_2D(points, name=name+"_short2", layer="TRACK") 
        
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
            self.set_current_coor([cutout_size[0]/2-self.overdev, track_right/2+gap_right+short_right+self.overdev], [1,0])
            right_quarter_up1 = self.draw_quarter_circle(name+'right_quarter_up1', 'TRACK', fillet_right1)
            self.set_current_coor([cutout_size[0]/2-self.overdev, -(track_right/2+gap_right+short_right)-self.overdev], [0,-1])
            right_quarter_down1 = self.draw_quarter_circle(name+'right_quarter_down1', 'TRACK', fillet_right1)
            right_short = self.unite([right_short, right_quarter_up1[1], right_quarter_down1[1]])
        else:
            fillet_right2 = (pad_size_right[1]/2-spacing_right-gap_right-track_right/2)/4+self.overdev
            self.set_current_coor( [cutout_size[0]/2-self.overdev, track_right/2+gap_right-self.overdev], [1,0],)
            right_quarter_up2 = self.draw_quarter_circle(name+'right_quarter_up2', 'TRACK', fillet_right2)
            self.set_current_coor([cutout_size[0]/2-self.overdev, -(track_right/2+gap_right)+self.overdev], [0,-1])
            right_quarter_down2 = self.draw_quarter_circle(name+'right_quarter_down2', 'TRACK' , fillet_right2)
            cutout = self.unite([cutout, right_quarter_up2[1], right_quarter_down2[1]], name+'nom1')
            
        if short_left!=0:
            self._fillet(track_left/2+gap_left+short_left-eps+self.overdev,[6,5], left_short)
            self._fillet(track_left/2+gap_left-eps-self.overdev,[2,1], left_short)
            
            fillet_left1= (pad_size_left[1]/2-spacing_left-short_left-gap_left-track_left/2)/4-self.overdev
            self.set_current_coor( [-cutout_size[0]/2+self.overdev,track_left/2+gap_left+short_left+self.overdev], [0,1])
            left_quarter_up1 = self.draw_quarter_circle(name+'left_quarter_up1', 'TRACK', fillet_left1)
            self.set_current_coor( [-cutout_size[0]/2+self.overdev,-(track_left/2+gap_left+short_left)-self.overdev], [-1,0])
            left_quarter_down1 = self.draw_quarter_circle(name+'left_quarter_down1', 'TRACK', fillet_left1)
            left_short = self.unite([left_short, left_quarter_up1[1], left_quarter_down1[1]])
        else:
            fillet_left2 = (pad_size_left[1]/2-spacing_left-gap_left-track_left/2)/4+self.overdev
            self.set_current_coor([-cutout_size[0]/2+self.overdev,track_left/2+gap_left-self.overdev], [0,1])
            left_quarter_up2 = self.draw_quarter_circle(name+'left_quarter_up2', 'TRACK', fillet_left2)
            self.set_current_coor([-cutout_size[0]/2+self.overdev,-(track_left/2+gap_left)+self.overdev], [-1,0])
            left_quarter_down2 = self.draw_quarter_circle(name+'left_quarter_down2', 'TRACK', fillet_left2)
            cutout = self.unite([cutout, left_quarter_up2[1], left_quarter_down2[1]], name+'nom2')
            
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
        added_gap_left = self.rect_center_2D([-cutout_size[0]/2, -(track_left+2*gap_left)/2+self.overdev], [self.overdev, (track_left+2*gap_left)-2*self.overdev], name=name+'_added_gap_left', layer='GAP')
        added_gap_right = self.rect_center_2D([cutout_size[0]/2, -(track_right+2*gap_right)/2+self.overdev], [-self.overdev, (track_right+2*gap_right)-2*self.overdev], name=name+'_added_gap_right', layer='GAP' )

        added_track_left = self.rect_center_2D([-cutout_size[0]/2, -(track_left)/2-self.overdev], [self.overdev, (track_left)+2*self.overdev], name=name+'_added_track_left', layer="TRACK")
        added_track_right = self.rect_center_2D([cutout_size[0]/2, -(track_right)/2-self.overdev], [-self.overdev, (track_right)+2*self.overdev], name=name+'_added_track_right', layer="TRACK")
        to_unite = to_unite+[added_track_left, added_track_right]
        
        cutout = self.unite([cutout, added_gap_left, added_gap_right], 'nom3')
        
    pads = self.unite(to_unite, name=name+'_pads')
    self.trackObjects.append(pads)
    
    self.gapObjects.append(cutout)
    
    if self.is_mask:
        self.maskObjects.append(mask)
    
    portOut1 = self.port(name+'portOut1', [cutout_size[0]/2,0], [1,0], track_right+2*self.overdev, gap_right-2*self.overdev)
    portOut2 = self.port(name+'portOut2', [-cutout_size[0]/2,0], [-1,0], track_left+2*self.overdev, gap_left-2*self.overdev)
    
    return [portOut1, portOut2, in_junction, out_junction]
    

print(__methods__)