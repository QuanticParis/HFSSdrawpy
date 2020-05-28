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

def draw_connector(self, pcb_track, pcb_gap, bond_length,
                   tr_line=True, name='connector_0'):
    '''
    Draws a CPW connector for inputs and outputs.

    Inputs:
    -------
    name : (str) should be different from other connector's name
    iBondLength: (float) corresponds to dimension a in the drawing
    iLineTest (Bool): unclear, keep False

        ground plane
        +------+
        |      |
        |      |
        |   +--+
    iIn |   |    iOut
        |   +--+
        |      |
        |      |
        +------+

    Outputs:
    --------
    returns created entities with formalism [Port], [Entiy]
    '''

    pcb_gap, pcb_track = parse_entry(pcb_gap, pcb_track)
    bond_length = parse_entry(bond_length)

    # track
    self.rect([pcb_gap, pcb_track/2], [bond_length, -pcb_track],
                             layer=TRACK, name=name+'_track')

    # gap
    self.rect([pcb_gap/2, pcb_gap+pcb_track/2],
              [pcb_gap/2 + bond_length, -(2*pcb_gap+pcb_track)],
              layer=GAP, name=name+'_gap')

    if self.is_mask:
        self.rect([pcb_gap/2-self.is_mask, pcb_gap+pcb_track/2+self.is_mask],
                  [pcb_gap/2 + bond_length+2*self.is_mask, -(2*pcb_gap+pcb_track+2*self.is_mask)],
                  layer=MASK, name=name+'_mask')

    with self([pcb_gap+bond_length,0], [1,0]):
        portOut, = create_port(self, widths=[pcb_track,2*pcb_gap+pcb_track],
                               name=name)

    if tr_line:
        ohm = self.rect([pcb_gap/2, pcb_track/2],
                        [pcb_gap/2, -pcb_track],
                        layer=RLC, name=name+'_ohm')
        points = [(pcb_gap/2+self.overdev, 0), (pcb_gap-self.overdev, 0)]
        ohm.assign_lumped_RLC(points, ('50ohm', 0, 0))
        self.polyline(points, name=name+'_line', closed=False, layer=DEFAULT)

        ohm.assign_mesh_length(pcb_track/10)

    return [portOut]
