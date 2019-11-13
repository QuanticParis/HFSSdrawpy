# -*- coding: utf-8 -*-
"""
Created on Wed Nov 13 11:54:00 2019

@author: Zaki
"""
from hfss import parse_entry
from PythonModeler import PythonModeler
from CustomElement import Port
import matplotlib as plt
from matplotlib.patches import Circle, Wedge, Polygon
from matplotlib.collections import PatchCollection
import matplotlib.pyplot as plt

class ConnectElt2(PythonModeler):
    def __init__(self, chip):
        self.chip = chip
        self.interface = chip.interface
        self.coor_sys = chip.coor_sys
    def decorator(func):
        def decorated(*args, **kwargs):
            # parse the instructions of the user
            self1  = args[0]
            name = args[1]
            iIn = args[2]
            iOut = args[3]
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
    
            
                
                
    def append_points(self, coor_list):
        points = [coor_list[0]]

        for coor in coor_list[1:]:
            points.append((points[-1][0] + coor[0],points[-1][1] + coor[1]))
        return points
    
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
#        fig, ax = plt.subplots()
#        patches = []
#        patches.append(Polygon(points3))
#        p = PatchCollection(patches, alpha=0.4)
#        ax.add_collection(p)
        
        
#        TODO !!
#        if not self.is_litho:
#            self.draw(self.name+"_mesh", points)
#            self.modeler.assign_mesh_length(self.name+"_mesh",1/2*iLength)

        self.iIn = retIn
        self.iOut = retOut
#        return [retIn, retOut]