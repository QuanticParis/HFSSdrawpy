# -*- coding: utf-8 -*-
"""
Created on Mon Oct 28 16:18:36 2019

@author: Zaki
"""

import numpy as np

from ..utils import parse_entry, Vector
from ..parameters import TRACK, GAP, RLC, MESH, MASK, DEFAULT, eps

def create_port(self, widths=None, subnames=None, layers=None, offsets=0, name='port_0'):
    """
    Creates a port and draws a small triangle for each element of the port.
    This function does exactly the same thing as the Body method 'port' except
    that if 2 widths are provided and subnames, layers, and offsets are not,
    assumes a CPW port with first track then gap.

    Parameters
    ----------
    widths : float, 'VariableString' or list, optional
        Width of the different elements of the port. If None, assumes the
        creation of a constraint_port. The default is None.
    subnames : str or list, optional
        The cable's parts will be name cablename_subname. If None simply
        numbering the cable parts. The default is None.
    layers : int or list, optional
        Each layer is described by an int that is a python constant that one
        should import. If None, layer is the DEFAULT The default is
        None.
    offsets : float, 'VariableString' or list, optional
        Describes the offset of the cable part wrt the center of the cable.
        The default is 0.
    name : str, optional
        Name of the port.

    Returns
    -------
    'Port'
        Returns a Port object

    """

    if widths is not None and isinstance(widths, list) and len(widths)==2:
        if subnames is None and layers is None :
            subnames = ['track', 'gap']
            layers = [TRACK, GAP]
            offsets = [0, 0]
    port, = self.port(widths=widths, subnames=subnames, layers=layers,
                     offsets=offsets, name=name)
    return [port]

def draw_connector(self, track, gap, bond_length, pcb_track, pcb_gap,
                   slope=1, tr_line=True, name='connector_0'):
    '''
    Draws a CPW connector for inputs and outputs.

    Inputs:
    -------
    name : (str) should be different from other connector's name
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
    returns created entities with formalism [Port], [Entiy]
    '''

    track, gap, pcb_gap, pcb_track = parse_entry(track, gap, pcb_gap, pcb_track)
    bond_length, slope = parse_entry(bond_length, slope)

    adaptDist = (pcb_track/2-track/2)/slope

    points = [(pcb_gap-self.overdev, pcb_track/2+self.overdev),
              (pcb_gap+bond_length, pcb_track/2+self.overdev),
              (pcb_gap+bond_length+adaptDist, self.overdev+track/2),
              (pcb_gap+bond_length+adaptDist, -track/2-self.overdev),
              (pcb_gap+bond_length, -pcb_track/2-self.overdev),
              (pcb_gap-self.overdev, -pcb_track/2-self.overdev)]

    track_entity = self.polyline(points, layer=TRACK,
                                 name=name+'_track')

    points = [(pcb_gap/2+self.overdev, pcb_gap+pcb_track/2-self.overdev),
             (pcb_gap+bond_length, pcb_gap+pcb_track/2-self.overdev),
             (pcb_gap+bond_length+adaptDist, gap+track/2-self.overdev),
             (pcb_gap+bond_length+adaptDist, -gap-track/2+self.overdev),
             (pcb_gap+bond_length, -pcb_gap-pcb_track/2+self.overdev),
             (pcb_gap/2+self.overdev, -pcb_gap-pcb_track/2+self.overdev)]

    gap_entity = self.polyline(points, layer=GAP, name=name+'_gap')

    if self.is_mask:
        points =[(pcb_gap/2-self.gap_mask, pcb_gap+pcb_track/2+self.gap_mask),
                  (pcb_gap+bond_length, pcb_gap+pcb_track/2+self.gap_mask),
                  (pcb_gap+bond_length+adaptDist, gap+track/2+self.gap_mask),
                  (pcb_gap+bond_length+adaptDist, -gap-self.gap_mask),
                  (pcb_gap+bond_length, (pcb_gap)+(track-pcb_track)*0.5-self.gap_mask),
                  (pcb_gap/2-self.gap_mask, (pcb_gap)+(track-pcb_track)*0.5-self.gap_mask)]

        mask_entity = self.polyline(points, layer=MASK,
                                    name=name+"_mask")

    with self([adaptDist+pcb_gap+bond_length,0], [1,0]):
        portOut, = create_port(self, widths=[track+2*self.overdev,
                                             2*gap+track-2*self.overdev],
                               name=name)

    if tr_line:
        ohm = self.rect([pcb_gap/2+self.overdev, pcb_track/2+self.overdev],
                        [pcb_gap/2-2*self.overdev, -pcb_track-2*self.overdev],
                        layer=RLC, name=name+'_ohm')
        points = [(pcb_gap/2+self.overdev, 0), (pcb_gap-self.overdev, 0)]
        ohm.assign_lumped_RLC(points, ('50ohm', 0, 0))
        self.polyline(points, name=name+'_line', closed=False, layer=DEFAULT)

        ohm.assign_mesh_length(pcb_track/10)

    return [portOut]
