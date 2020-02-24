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
layer_TRACK = 1
layer_GAP = 0
layer_RLC = 2
layer_MESH = 3
layer_MASK = 4
layer_Default = 10

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
    gap = self.rect_center_2D([0,0], [iLength, iTrack+2*iGap], name=name+'_gap', layer=layer_GAP)
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
    
    cutout = self.rect_center_2D([0,0], cutout_size-Vector([2*self.overdev, 2*self.overdev]), name=name+"_cutout", layer=layer_GAP)

    if self.is_mask:
        mask = self.rect_center_2D([0,0], cutout_size+Vector([track,track]),name=name+"_mask",layer=layer_MASK)
        self.maskObjects.append(mask)
    if not self.is_litho:
        mesh = self.rect_center_2D([0,0], cutout_size, name=name+"_mesh", layer=layer_MESH)
        self.mesh_zone(mesh, 2*track)

    right_pad = self.rect_center_2D(Vector(pad_spacing+pad_size[0],0)/2, pad_size+Vector([2*self.overdev, 2*self.overdev]), name=name+"_pad1", layer=layer_TRACK)
    left_pad = self.rect_center_2D(-Vector(pad_spacing+pad_size[0],0)/2, pad_size+Vector([2*self.overdev, 2*self.overdev]), name=name+"_pad2", layer=layer_TRACK)
    

    track_J=Jwidth*4.
    in_junction = self.port(name+'_in_jct', [-pad_spacing/2+self.overdev, 0], [1,0], track_J+2*self.overdev, 0)
    out_junction = self.port(name+'_out_jct', [pad_spacing/2-self.overdev, 0], [1,0], track_J+2*self.overdev, 0)
#        self.ports[name+'_in_jct'] = in_junction
#        self.ports[name+'_out_jct'] = out_junction
#        junction = self.connect_elt(name+'_junction', name+'_in_jct', name+'_out_jct')
    junction_pad = self._connect_JJ(name+'JJ', name+'_in_jct', name+'_out_jct', Jwidth+2*self.overdev, iInduct=Jinduc, fillet=None)
    pads = self.unite([right_pad, left_pad, junction_pad], name=name+'_pads')
    if fillet is not None:
        self._fillet(fillet+self.overdev,[19], pads)
        self._fillet(0.75*fillet-self.overdev,[18,17,14,13], pads)
        self._fillet(fillet+self.overdev,[12,11,10,9,8,7], pads)
        self._fillet(0.75*fillet-self.overdev,[6,5,2,1], pads)
        self._fillet(fillet+self.overdev,0, pads)
    if nport==1:
        self.port(name+'_1',cutout_size.px()/2, self.ori, track, gap)
        
    return pads, cutout, in_junction, out_junction

@register_method
@move
def draw_cylinder_transmon(self,
                   name, 
                   cutout_size,
                   pad_spacing,
                   pad_size,
                   Jwidth,
                   space,
                   spacing2,
                   tline_l,
                   tline_w,
                   thickness,
                   Jinduc,
                   fillet=None):
    cutout_size, pad_spacing, pad_size, Jwidth, space, thickness = parse_entry((cutout_size, pad_spacing, pad_size, Jwidth, space, thickness))
    insert_radius = 1.2/2*cutout_size[1]
    length = cutout_size[0]
    cylinder_side = self.cylinder([0,0,0], insert_radius, length, "X", name=name+"_cylinder", layer=layer_Default)
    
    sapphire_chip = self.box_corner([0,-cutout_size[1]/2-self.overdev/2,0], [cutout_size[0],cutout_size[1]+self.overdev, thickness], "3D", name=name+"_sapphire_chip")
    
    self.set_current_coor([0,0,thickness], [1,0])
    simple_transmon = self.draw_simple_transmon(name, 
                                               cutout_size,
                                               pad_spacing,
                                               pad_size,
                                               Jwidth,
                                               space,
                                               spacing2,
                                               tline_l,
                                               tline_w,
                                               Jinduc,
                                               fillet=None)
    return cylinder_side

@register_method
@move
def draw_simple_transmon(self,
                   name, 
                   cutout_size,
                   pad_spacing,
                   pad_size,
                   Jwidth,
                   space,
                   spacing2,
                   tline_l,
                   tline_w,
                   Jinduc,
                   fillet=None):

    cutout_size, pad_spacing, pad_size, Jwidth, space = parse_entry((cutout_size, pad_spacing, pad_size, Jwidth, space))
    cutout_size = Vector(cutout_size)
    pad_size = Vector(pad_size)
    capa_plasma=0
    iTrackJ=Jinduc

    
    cutout = self.rect_corner_2D([0,-cutout_size[1]/2+self.overdev], cutout_size-Vector([2*self.overdev, 2*self.overdev]), name=name+"_cutout", layer=layer_GAP)
    transmission_line = self.rect_center_2D(Vector(space+pad_size[0]+pad_spacing+spacing2+tline_l/2, 0), Vector([tline_l,tline_w]), name= name+"_transmission_line", layer=layer_GAP)
    self.assign_perfect_E(transmission_line)
#    if self.is_mask:
#        mask = self.rect_center_2D([0,0], cutout_size+Vector([track,track]),name=name+"_mask",layer=layer_MASK)
#        self.maskObjects.append(mask)
    right_pad = self.rect_corner_2D(Vector(space, -pad_size[1]/2), pad_size+Vector([2*self.overdev, 2*self.overdev]), name=name+"_pad1", layer=layer_TRACK)
    left_pad  = self.rect_corner_2D(Vector(space+pad_size[0]+pad_spacing, -pad_size[1]/2), pad_size+Vector([2*self.overdev, 2*self.overdev]), name=name+"_pad2", layer=layer_TRACK)
    JJ = self.rect_center_2D(Vector(space+pad_size[0]+pad_spacing/2,0), [pad_spacing, '5um'], name=name+'_lumped', layer=layer_RLC)
    plot_JJ = self.box_center([space+pad_size[0]+pad_spacing/2,0,0], [2*pad_spacing, '10um', pad_spacing], name=name+'_lumped', layer=layer_RLC, nonmodel = True)

    self.mesh_zone(JJ ,'0.01mm')

#    track_J=Jwidth*4.
#    in_junction = self.port(name+'_in_jct', [-pad_spacing/2+self.overdev, 0], [1,0], track_J+2*self.overdev, 0)
#    out_junction = self.port(name+'_out_jct', [pad_spacing/2-self.overdev, 0], [1,0], track_J+2*self.overdev, 0)
#        self.ports[name+'_in_jct'] = in_junction
#        self.ports[name+'_out_jct'] = out_junction
#        junction = self.connect_elt(name+'_junction', name+'_in_jct', name+'_out_jct')
#    flow_line = self.polyline_2D([[-pad_spacing/2+self.overdev, 0],[pad_spacing/2-self.overdev, 0]], closed = False, name=name+'_flowline', layer=layer_MESH)

#    self.assign_lumped_RLC(JJ,0, Jinduc, capa_plasma, flow_line)
    pads = self.unite([right_pad, left_pad], name=name+'_pads')
#    self.make_material(pads, "\"aluminum\"")
    self.assign_perfect_E(pads)
    
    if fillet is not None:
        self._fillet(fillet+self.overdev,[19], pads)
        self._fillet(0.75*fillet-self.overdev,[18,17,14,13], pads)
        self._fillet(fillet+self.overdev,[12,11,10,9,8,7], pads)
        self._fillet(0.75*fillet-self.overdev,[6,5,2,1], pads)
        self._fillet(fillet+self.overdev,0, pads)

    return [0,0], [pad_spacing, 10*iTrackJ], "x", 0, Jinduc, capa_plasma,

    

        
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
                    Jlength,
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
                          Jlength,
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
     Jlength,
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

    cutout = self.rect_center_2D([0,0], cutout_size-Vector([self.overdev*2, self.overdev*2]),name =name+"_cutout", layer=layer_Default)
    if not self.is_litho:
        mesh = self.rect_center_2D([0,0], cutout_size, name=name+"_mesh", layer=layer_MESH)
    if self.is_mask:
        mask = self.rect_center_2D([0,0], cutout_size+Vector([self.gap_mask,self.gap_mask])*2, name=name+"_mask", layer=layer_MASK)
    
    track_J=Jwidth*4.
    in_junction = self.port(name+'_in_jct', [-pad_spacing/2+self.overdev, 0], [1,0], track_J+2*self.overdev, 0)
    out_junction = self.port(name+'_out_jct', [+pad_spacing/2-self.overdev, 0], [-1,0], track_J+2*self.overdev, 0)
    
#        self.ports[name+'_in_jct'] = in_junction
#        self.ports[name+'_out_jct'] = out_junction
#        junction = self.connect_elt(name+'_junction', name+'_in_jct', name+'_out_jct')
    junction_pads = self._connect_JJ(name+'JJ', name+'_in_jct', name+'_out_jct', Jwidth+2*self.overdev, iLengthJ=Jlength, iInduct=Jinduc, fillet=None)
    
    raw_points_right = [(pad_spacing/2-self.overdev, -pad_size_right[1]/2-self.overdev),
                  (pad_size_right[0]+2*self.overdev, 0),
                  (0, pad_size_right[1]/2-spacing_right-short_right-gap_right-track_right/2+2*self.overdev),
                  (-pad_size_right[0]+length_right, 0),
                  (0, (spacing_right+short_right+gap_right+track_right/2)*2-2*self.overdev),
                  (pad_size_right[0]-length_right, 0),
                  (0, pad_size_right[1]/2-spacing_right-short_right-gap_right-track_right/2+2*self.overdev),
                  (-pad_size_right[0]-2*self.overdev, 0)]
    points_right = self.append_points(raw_points_right)
    right_pad = self.polyline_2D(points_right, name=name+"_pad1", layer=layer_TRACK)
    
    raw_points_left = [(-pad_spacing/2+self.overdev, -pad_size_left[1]/2-self.overdev),
                  (-pad_size_left[0]-2*self.overdev, 0),
                  (0, pad_size_left[1]/2-spacing_left-short_left-gap_left-track_left/2+2*self.overdev),
                  (pad_size_left[0]-length_left, 0),
                  (0, (spacing_left+short_left+gap_left+track_left/2)*2-2*self.overdev),
                  (-pad_size_left[0]+length_left, 0),
                  (0, pad_size_left[1]/2-spacing_left-short_left-gap_left-track_left/2+2*self.overdev),
                  (pad_size_left[0]+2*self.overdev, 0)]
    points_left = self.append_points(raw_points_left)
    left_pad = self.polyline_2D(points_left, name=name+"_pad2", layer=layer_TRACK) 
    
    right_track_raw_points = [(cutout_size[0]/2-self.overdev,-track_right/2-self.overdev),
                             (-(cutout_size[0]/2-pad_spacing/2-length_right-gap_right-short_right-spacing_right), 0),
                             (0, track_right+2*self.overdev),
                             ((cutout_size[0]/2-pad_spacing/2-length_right-gap_right-short_right-spacing_right), 0)]
    right_track_points = self.append_points(right_track_raw_points)
    right_track = self.polyline_2D(right_track_points, name=name+"_track1", layer=layer_TRACK) 
    
    left_track_raw_points = [(-cutout_size[0]/2+self.overdev,-track_left/2-self.overdev),
                             ((cutout_size[0]/2-pad_spacing/2-length_left-gap_left-short_left-spacing_left), 0),
                             (0, track_left+2*self.overdev),
                             (-(cutout_size[0]/2-pad_spacing/2-length_left-gap_left-short_left-spacing_left), 0)]
    left_track_points = self.append_points(left_track_raw_points)
    left_track = self.polyline_2D(left_track_points, name=name+"_track2", layer=layer_TRACK) 
    
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
        right_short = self.polyline_2D(points, name=name+"_short1", layer=layer_TRACK) 
        
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
        left_short = self.polyline_2D(points, name=name+"_short2", layer=layer_TRACK) 
        
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
            right_quarter_up1 = self.draw_quarter_circle(name+'right_quarter_up1', layer_TRACK, fillet_right1)
            self.set_current_coor([cutout_size[0]/2-self.overdev, -(track_right/2+gap_right+short_right)-self.overdev], [0,-1])
            right_quarter_down1 = self.draw_quarter_circle(name+'right_quarter_down1', layer_TRACK, fillet_right1)
            right_short = self.unite([right_short, right_quarter_up1, right_quarter_down1])
        else:
            fillet_right2 = (pad_size_right[1]/2-spacing_right-gap_right-track_right/2)/4+self.overdev
            self.set_current_coor( [cutout_size[0]/2-self.overdev, track_right/2+gap_right-self.overdev], [1,0],)
            right_quarter_up2 = self.draw_quarter_circle(name+'right_quarter_up2', layer_TRACK, fillet_right2)
            self.set_current_coor([cutout_size[0]/2-self.overdev, -(track_right/2+gap_right)+self.overdev], [0,-1])
            right_quarter_down2 = self.draw_quarter_circle(name+'right_quarter_down2', layer_TRACK , fillet_right2)
            cutout = self.unite([cutout, right_quarter_up2, right_quarter_down2], name+'nom1')
            
        if short_left!=0:
            self._fillet(track_left/2+gap_left+short_left-eps+self.overdev,[6,5], left_short)
            self._fillet(track_left/2+gap_left-eps-self.overdev,[2,1], left_short)
            
            fillet_left1= (pad_size_left[1]/2-spacing_left-short_left-gap_left-track_left/2)/4-self.overdev
            self.set_current_coor( [-cutout_size[0]/2+self.overdev,track_left/2+gap_left+short_left+self.overdev], [0,1])
            left_quarter_up1 = self.draw_quarter_circle(name+'left_quarter_up1', layer_TRACK, fillet_left1)
            self.set_current_coor( [-cutout_size[0]/2+self.overdev,-(track_left/2+gap_left+short_left)-self.overdev], [-1,0])
            left_quarter_down1 = self.draw_quarter_circle(name+'left_quarter_down1', layer_TRACK, fillet_left1)
            left_short = self.unite([left_short, left_quarter_up1, left_quarter_down1])
        else:
            fillet_left2 = (pad_size_left[1]/2-spacing_left-gap_left-track_left/2)/4+self.overdev
            self.set_current_coor([-cutout_size[0]/2+self.overdev,track_left/2+gap_left-self.overdev], [0,1])
            left_quarter_up2 = self.draw_quarter_circle(name+'left_quarter_up2', layer_TRACK, fillet_left2)
            self.set_current_coor([-cutout_size[0]/2+self.overdev,-(track_left/2+gap_left)+self.overdev], [-1,0])
            left_quarter_down2 = self.draw_quarter_circle(name+'left_quarter_down2', layer_TRACK, fillet_left2)
            cutout = self.unite([cutout, left_quarter_up2, left_quarter_down2], name+'nom2')
            
#        self._fillet(pad_size_right[0]/4+self.overdev,[7,0], right_pad)
#        self._fillet((pad_size_right[1]/2-spacing_right-short_right-gap_right-track_right/2)/4+self.overdev,[1,2,5,6], right_pad)
#        self._fillet(track_right/2+gap_right+short_right+spacing_right-eps-self.overdev,[5,6], right_pad)
#
#        self._fillet(pad_size_left[0]/4+self.overdev,[7,0], left_pad)
#        self._fillet((pad_size_left[1]/2-spacing_left-short_left-gap_left-track_left/2)/4+self.overdev,[1,2,5,6], left_pad)
#        self._fillet(track_left/2+gap_left+short_left+spacing_left-eps-self.overdev,[5,6], left_pad)
    
    if not self.is_litho:
        self.mesh_zone(mesh, 4*Jwidth)
    
    to_unite = [left_pad, right_pad, right_track, left_track, junction_pads]
    if short_right!=0:
        to_unite.append(right_short)
    if short_left!=0:
        to_unite.append(left_short)
    
    if self.is_overdev:
        added_gap_left = self.rect_center_2D([-cutout_size[0]/2, -(track_left+2*gap_left)/2+self.overdev], [self.overdev, (track_left+2*gap_left)-2*self.overdev], name=name+'_added_gap_left', layer=layer_GAP)
        added_gap_right = self.rect_center_2D([cutout_size[0]/2, -(track_right+2*gap_right)/2+self.overdev], [-self.overdev, (track_right+2*gap_right)-2*self.overdev], name=name+'_added_gap_right', layer=layer_GAP )

        added_track_left = self.rect_center_2D([-cutout_size[0]/2, -(track_left)/2-self.overdev], [self.overdev, (track_left)+2*self.overdev], name=name+'_added_track_left', layer=layer_TRACK)
        added_track_right = self.rect_center_2D([cutout_size[0]/2, -(track_right)/2-self.overdev], [-self.overdev, (track_right)+2*self.overdev], name=name+'_added_track_right', layer=layer_TRACK)
        to_unite = to_unite+[added_track_left, added_track_right]
        
#        cutout = self.unite([cutout, added_gap_left, added_gap_right], 'nom3')
        
    pads = self.unite(to_unite, name=name+'_pads')
    
        
    portOut2 = self.port(name+'_portOut2', [cutout_size[0]/2,0], [1,0], track_right+2*self.overdev, gap_right-2*self.overdev)
    portOut1 = self.port(name+'_portOut1', [-cutout_size[0]/2,0], [-1,0], track_left+2*self.overdev, gap_left-2*self.overdev)
    
    return [portOut1, portOut2, in_junction, out_junction]
    
@register_method
@move
def draw_capa_inline(self, name, iTrack, iGap, capa_length, pad_spacing, n_pad=1, iTrack_capa=None, iGap_capa=None, premesh=True, tight=False): #iGap_capa is added gap
    
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
        drawn_pads.append(self.rect_corner_2D([-capa_length/2,-iTrack/2-self.overdev], [capa_length/2-pad_spacing/2+self.overdev,iTrack+2*self.overdev], name=name+'_pad_left', layer=layer_TRACK))
        drawn_pads.append(self.rect_corner_2D([capa_length/2,-iTrack/2-self.overdev], [-(capa_length/2-pad_spacing/2)-self.overdev,iTrack+2*self.overdev], name=name+'_pad_right', layer=layer_TRACK))
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
            drawn_pads.append(self.rect_corner_2D([-capa_length/2+offset, curr_height-self.overdev], pad_size+Vector([self.overdev, 2*self.overdev]), name=name+"_pad"+str(ii), layer=layer_TRACK))
            drawn_pads.append(self.rect_corner_2D([-capa_length/2+offset+pad_spacing-self.overdev, curr_height+pad_width+pad_spacing-self.overdev], pad_size+Vector([self.overdev, 2*self.overdev]), name = name+"_pad"+str(ii)+'b', layer=layer_TRACK))
            curr_height = curr_height+2*(pad_width+pad_spacing)
        if n_pad%2!=0:
            drawn_pads.append(self.rect_corner_2D([-capa_length/2+offset, curr_height-self.overdev], pad_size+Vector([self.overdev, 2*self.overdev]), name=name+"_pad"+str(ii+1), layer=layer_TRACK))
        else:
            curr_height = curr_height-2*(pad_width+pad_spacing)
        if not usual:
#                drawn_pads.append(self.draw_rect(name+'_connect_left', self.coor([-capa_length/2-self.overdev, -iTrack_capa/2-self.overdev]), self.coor_vec([iTrack+2*self.overdev, curr_height+iTrack_capa/2+pad_width+2*self.overdev])))
            drawn_pads.append(self.rect_corner_2D([-capa_length/2-self.overdev, -iTrack_capa/2-self.overdev], [offset+2*self.overdev, curr_height+iTrack_capa/2+pad_width+2*self.overdev], name=name+'_connect_left', layer=layer_TRACK))
#                    drawn_pads.append(self.draw_rect(name+'_connect_right', self.coor([capa_length/2+self.overdev, -iTrack_capa/2+pad_width+pad_spacing-self.overdev]), self.coor_vec([-iTrack-2*self.overdev, curr_height+(iTrack_capa/2-2*pad_width-2*pad_spacing)+pad_width+2*self.overdev])))
            drawn_pads.append(self.rect_corner_2D([capa_length/2+self.overdev, -iTrack_capa/2+pad_width+pad_spacing-self.overdev], [-offset-2*self.overdev, curr_height+(iTrack_capa/2-2*pad_width-2*pad_spacing)+pad_width+2*self.overdev], name=name+'_connect_right', layer=layer_TRACK))
        else:
            pass
#                    drawn_pads.append(self.draw_rect(name+'_connect_right', self.coor([capa_length/2+self.overdev, -iTrack_capa/2+pad_width+pad_spacing-self.overdev]), self.coor_vec([-iTrack-2*self.overdev, curr_height+iTrack_capa/2+pad_width+2*self.overdev])))
            #TODO
            #drawn_pads.append(self.rect_corner_2D([capa_length/2+self.overdev, -iTrack_capa/2+pad_width+pad_spacing-self.overdev], [-iTrack-2*self.overdev, curr_height+iTrack_capa/2+pad_width+2*self.overdev], name=name+'_connect_right', layer=layer_TRACK ))
#                    
#                    drawn_pads.append(self.draw_rect(name+'_connect_right_b', self.coor([capa_length/2+self.overdev, -iTrack/2-self.overdev]), self.coor_vec([-iTrack-2*self.overdev, (iTrack/2+iTrack_capa/2)+2*self.overdev])))
#                    drawn_pads.append(self.draw_rect(name+'_connect_left_b', self.coor([-capa_length/2-self.overdev, iTrack/2+self.overdev]), self.coor_vec([iTrack+2*self.overdev, -(iTrack/2+iTrack_capa/2)-2*self.overdev])))
            #drawn_pads.append(self.rect_corner_2D([capa_length/2+self.overdev, -iTrack/2-self.overdev], [-iTrack-2*self.overdev, (iTrack/2+iTrack_capa/2)+2*self.overdev], name=name+'_connect_right_b', layer=layer_TRACK))
            #drawn_pads.append(self.rect_corner_2D([-capa_length/2-self.overdev, iTrack/2+self.overdev],[iTrack+2*self.overdev, -(iTrack/2+iTrack_capa/2)-2*self.overdev], name=name+'_connect_left_b', layer=layer_TRACK))

    if tight:
        portOut1 = self.port(name+'_outPort1',[capa_length/2,0], [1,0], iTrack_capa-2*pad_width-2*pad_spacing+2*self.overdev, iGap-(iTrack_capa-2*pad_width-2*pad_spacing-iTrack)/2-2*self.overdev)
    else:
        portOut1 = self.port(name+'_outPort1',[capa_length/2,0], [1,0], iTrack+2*self.overdev, iGap-2*self.overdev)

#            portOut2 = [self.pos-self.ori*capa_length/2, -self.ori, iTrack_capa+2*self.overdev, iGap-(iTrack_capa-iTrack)/2-2*self.overdev]
#        else:
    portOut2 = self.port(name+'_outPort2',[-capa_length/2, 0], [-1,0], iTrack+2*self.overdev, iGap-2*self.overdev)

    if self.is_overdev:
        drawn_pads.append(self.rect_corner_2D([capa_length/2, -iTrack/2], [-self.overdev, iTrack], name=name+'_added_right', layer=layer_TRACK))
        drawn_pads.append(self.rect_corner_2D([-capa_length/2, -iTrack/2], [self.overdev, iTrack], name=name+'_added_left', layer=layer_TRACK))

#    pads = self.unite(drawn_pads, name=name+'_pads')
    
#    self.trackObjects.append(pads)
    
    self.gapObjects.append(self.rect_center_2D([0,0], [capa_length, iTrack + 2*iGap+2*iGap_capa-2*self.overdev], name=name+"_gap", layer=layer_GAP))
    if self.is_mask:
        self.maskObjects.append(self.rect_center_2D([0,0], [capa_length, iTrack + 2*iGap+2*iGap_capa +2*self.gap_mask], name=name+"_mask", layer=layer_MASK))
    
    if premesh:
        if not self.is_litho:
            mesh = self.rect_center_2D([0,0], [capa_length, iTrack_capa+2*self.overdev], name=name+"_mesh", layer=layer_MESH)
            self.mesh_zone(mesh, pad_spacing)
#
#    self.ports[name+'_1'] = portOut1
#    self.ports[name+'_2'] = portOut2
            
    return portOut1, portOut2

@register_method
@move    
def draw_capa_interdigitated(self, name, iTrack, iGap, teeth_size, gap_size, N_period, fillet):
    '''
    '''
    iTrack, iGap = parse_entry((iTrack, iGap))
    fillet = parse_entry(fillet)
    teeth_size=parse_entry(teeth_size)
    gap_size=parse_entry(gap_size)
    teeth_size = Vector(teeth_size)
    portOut1 = self.port(name+'_outPort1',[0,teeth_size[0]+iTrack+iGap], [0,1], iTrack+2*self.overdev, iGap-2*self.overdev)
    portOut2 = self.port(name+'_outPort1',[0,-teeth_size[0]+iTrack+iGap], [0,-1], iTrack+2*self.overdev, iGap-2*self.overdev)


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
    #TODO Find append_absolute_point and append_points difference
    points = self.append_points(raw_points)
    connection = self.polyline_2D(points, closed=False, name=name+"_capagap", layer=layer_TRACK)
    
    
    self._fillet(connection, fillet)
    raw_points=[(-gap_size-teeth_size[0]+self.overdev,N_teeth*teeth_size[1]+self.overdev),(gap_size-teeth_size[0]-self.overdev,N_teeth*teeth_size[1]+self.overdev)]
    points=self.append_absolute_points(raw_points)
    capagap_starter = self.polyline2D(points, closed=False, name=name+'_width', layer=layer_GAP)
    
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
    points = self.append_points(raw_points)
    pads = self.polyline_2D(points, name+"_pads")
    #####Filets on edges of the capa
    self.fillet(fillet+self.overdev,[11,6,5,0], pads)
    
    pads_sub = self.subtract(pads, [capagap])
    
    #####Filets on edge
    self.fillet(fillet-self.overdev,[73,70], pads)

    #####Filets on last teeth
    self.fillet(0.5*fillet+self.overdev,[67,38,10])

    #####Filets on edge
    self.fillet(fillet-self.overdev,[7,4], pads)

    #####Filets on last teeth
    self.fillet(0.5*fillet+self.overdev,[1], pads)





    if not self.is_overdev:
        self.gapObjects.append(self.rect_center_2D([0,0], [2*teeth_size[0]+2*iTrack+2*iGap, 2*N_teeth*teeth_size[1]+2*iGap], name=name+"_gap", layer=layer_GAP))
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
        points = self.append_points(raw_points)
        self.gapObjects.append(self.polyline_2D(points, name+"_gap"))
#    if self.is_mask:

        
    if not self.is_litho:
        mesh=self.polyline_2D(points, name=name+"_mesh")
        self.modeler.assign_mesh_length(mesh, iTrack)

    self.trackObjects.append(pads_sub)

#
#
#    self.ports[name+'_1'] = portOut1
#    self.ports[name+'_2'] = portOut2
    return portOut1, portOut2