# -*- coding: utf-8 -*-
"""
Created on Mon Oct 28 16:21:38 2019

@author: Zaki
"""
from designer import Vector, Circuit, way, equalfloat, eps
from KeyElement import KeyElt
import numpy as np
from .hfss import parse_entry
from .hfss import VariableString

TOP = [0, 1]
DOWN = [0, -1]
RIGHT = [1, 0]
LEFT = [-1, 0]

POS = 0
ORI = 1
TRACK = 2
GAP = 3
    
class ConnectElt(KeyElt, Circuit):
    
    def __init__(self, name='connect_elt', iIn='iInt', iOut=None, layer=None):
        print(name)

        self.name = name
        self.layer = layer

        if isinstance(iIn, str):
            if iIn in self.ports:
                iInPort = parse_entry(self.ports[iIn])
                self.iIn = iIn # name of the in port
                self.pos = Vector(iInPort[POS]) # CONFLICT WITH SELF.POS OF KEY_ELT
                self.ori = Vector(iInPort[ORI])
                _, _, self.inTrack, self.inGap = iInPort
            elif iIn in self.ports_dc:
                iInPort = parse_entry(self.ports_dc[iIn])
                self.iIn = iIn # name of the in port
                self.pos = Vector(iInPort[POS])
                self.ori = Vector(iInPort[ORI])
                _, _, self.rel_posIn, self.widIn, self.cutIn, _ = iInPort
                self.multIn = int(iInPort[-1])
#                self.layerIn = layer
            else:
                raise ValueError('inPort %s does not exist' % iIn)
        #is the following really useful ?
        #it is when defining intermediate port for draw_adaptor
        elif isinstance(iIn, list):
            if iOut in self.ports: # if iIn is defined on the fly then iOut is well defined
                iInPort = parse_entry(iIn)
                self.iIn = 'iIn' # dummy name test
                self.pos = Vector(iInPort[POS])
                self.ori = Vector(iInPort[ORI])
                _, _, self.inTrack, self.inGap = iInPort
            elif iOut in self.ports_dc:
                iInPort = parse_entry(iIn)
                self.iIn = 'iIn' # dummy name test
                self.pos = Vector(iInPort[POS])
                self.ori = Vector(iInPort[ORI])
                _, _, self.rel_posIn, self.widIn, self.cutIn, _ = iInPort
                self.multIn = int(iInPort[-1])
        else:
            raise ValueError('iIn should be given a port name, a list or nothing')

        if isinstance(iOut, str):
            if iOut in self.ports:
                iOutPort = parse_entry(self.ports[iOut])
                self.isOut = True
                self.iOut = iOut
                self.posOut = Vector(iOutPort[POS])
                self.oriOut = Vector(iOutPort[ORI])
                _, _, self.outTrack, self.outGap = iOutPort
            elif iOut in self.ports_dc:
                iOutPort = parse_entry(self.ports_dc[iOut])
                self.isOut = True
                self.iOut = iOut
                self.posOut = Vector(iOutPort[POS])
                self.oriOut = Vector(iOutPort[ORI])
                _, _, self.rel_posOut, self.widOut, self.cutOut, _ = iOutPort
                self.multOut = int(iOutPort[-1])
#                self.layerOut = layer
            else:
                raise ValueError('outPort %s does not exist' % iOut)
        elif isinstance(iOut, list):
            if iIn in self.ports:
                iOutPort = parse_entry(iOut)
                self.outTrack, self.outGap = iOutPort[-2], iOutPort[-1]
                if len(iOut)>3:
                    self.posOut = Vector(iOutPort[POS])
                    self.oriOut = Vector(iOutPort[ORI])
                self.isOut = False
            elif iIn in self.ports_dc:
                iOutPort = parse_entry(iOut)
                self.rel_posOut, self.widOut, self.cutOut, self.multOut = iOutPort[-4], iOutPort[-3], iOutPort[-2], int(iOutPort[-1])
                if len(iOut)>4:
                    self.posOut = Vector(iOutPort[POS])
                    self.oriOut = Vector(iOutPort[ORI])
                self.isOut = False
        elif iOut is None:
            pass
        else:
            raise ValueError('iOut should be given a port name, a list or nothing')

    def rot(self, x, y):
        return Vector(x, y).rot(self.ori)

    def append_points(self, coor_list):
        points = [self.pos + self.rot(*coor_list[0])]
        for coor in coor_list[1:]:
            points.append(points[-1] + self.rot(*coor))
        return points

    def draw_capa(self, iLength, iWidth, iSize):
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
        retIn = [self.pos, -self.ori, self.inTrack, self.inGap]
        retOut = [self.pos+(self.inGap+self.outGap+iSize+2*iWidth)*self.ori,
                  self.ori, self.outTrack, self.outGap]

        points = self.append_points([(self.inGap+iWidth, 0),
                                     (0, -iLength/2),
                                     (-iWidth, 0),
                                     (0, iLength/2-self.inTrack/2),
                                     (-self.inGap, 0),
                                     (0, self.inTrack),
                                     (self.inGap, 0),
                                     (0, iLength/2-self.inTrack/2),
                                   (iWidth, 0)])
        self.trackObjects.append(self.draw(self.name+"_track1", points))

        points = self.append_points([(self.inGap+iWidth+iSize, 0),
                                     (0, -iLength/2),
                                     (+iWidth, 0),
                                     (0, iLength/2-self.outTrack/2),
                                     (+self.outGap, 0),
                                     (0, self.outTrack),
                                     (-self.outGap, 0),
                                     (0, iLength/2-self.outTrack/2),
                                     (-iWidth, 0)])
        self.trackObjects.append(self.draw(self.name+"_track2", points))

        points = self.append_points([(0, 0),
                                     (0, iLength/2+self.inGap),
                                     (self.inGap+iWidth+iSize/2, 0),
                                     (0, self.outGap-self.inGap),
                                     (self.outGap+iWidth+iSize/2, 0),
                                     (0, -iLength-2*self.outGap),
                                     (-(self.outGap+iWidth+iSize/2)),
                                     (0, self.outGap-self.inGap),
                                     (self.inGap+iWidth+iSize/2)])
        self.gapObjects.append(self.draw(self.name+"_gap", points))
        
        if not self.is_litho:
            self.draw(self.name+"_mesh", points)
            self.modeler.assign_mesh_length(self.name+"_mesh",1/2*iLength)

        self.iIn = retIn
        self.iOut = retOut
#        return [retIn, retOut]


#        CreateBondwire(name+"_bondwire", iIn)
    def find_slanted_path(self):
        
        iIn_pos = self.pos
        iIn_ori = self.ori
        iOut_pos = self.posOut
        iOut_ori = self.oriOut
        
        if iIn_ori.dot(iOut_ori)!=-1:
            raise ValueError('Cannot find slanted path: ports are not oriented correctly')
        else:
            dist = (iOut_pos-iIn_pos).dot(iIn_ori)
        
        pointA = iIn_pos+iIn_ori*dist/3
        pointB = iOut_pos+iOut_ori*dist/3
        
        self.to_bond.append([iIn_pos, pointA])
        self.to_bond.append([pointB, iOut_pos])
        return [iIn_pos, pointA, pointB, iOut_pos], dist/3

    def find_path(self, fillet, is_meander, to_meander, meander_length, meander_offset):
        
        iIn_pos = self.pos
        iIn_ori = self.ori
        iOut_pos = self.posOut
        iOut_ori = self.oriOut
        
        room_bonding = 0*100e-6 #SMPD MANU BOND SPACE
        print(str(1.1*fillet))
        
        pointA = iIn_pos+iIn_ori*room_bonding
        pointB = iOut_pos+iOut_ori*room_bonding

        point1 = iIn_pos+iIn_ori*(1.1*fillet+room_bonding)
        point2 = iOut_pos+iOut_ori*(1.1*fillet+room_bonding)

        def next_point(point1, point2, vec):
            choice1 = point1+vec*(point2-point1).dot(vec)
            choice2 = point1+vec.orth()*(point2-point1).dot(vec.orth())
            return [[point1, choice1, point2], [point1, choice2, point2]]


        points_choices = []
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
            prev_point = points[0]
            _points = [points[0]]
            vecs = []
            for point in points[1:]:
                if not equal_float(self.val(point)[0], self.val(prev_point)[0]) or not equal_float(self.val(point)[1], self.val(prev_point)[1]):
                    vec = self.val(point-prev_point)
                    length += self.val(vec).norm()
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

#        return_bonds(new_float_points, fillet, length_fillet, n_fillet)
#
#
#
#
#
#

#        return


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
        

    def draw_cable(self, fillet="0.3mm", is_bond=False, is_meander=False, to_meanders = [1,0,1,0,1,0,1,0,1,0], meander_length=0, meander_offset=0, is_mesh=False, constrains=[], reverse_adaptor=False, layer=None):
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
        
        
    def draw_slanted_cable(self, fillet=None, is_bond=False, is_mesh=False, constrains=[], reverse_adaptor=False, layer=None):
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

        '''
        if fillet is not None:
            fillet =parse_entry(fillet)
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
  
        
        
        tracks = []
        gaps = []
        masks = []
        to_add = ''
        
        points, calc_fillet = self.find_slanted_path()
        if fillet is None:
            fillet = calc_fillet
        connection = self.draw(self.name+'_track', points, closed=False)
# TODO Implement cable_lenght calculation for slanted cable
#        cable_length = self.length(points, 0, len(points)-1, fillet)+self.val(adaptor_length)#        print('{0}_length = {1:.3f} mm'.format(self.name, cable_length*1000))
        connection.fillets(fillet-eps)

        connection_gap = connection.copy(self.name+"_gap")

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
            print(tracks)
            names = [self.name+'_track', self.name+'_gap', self.name+'_mask']
            if track_adaptor is not None:
                names = [self.name+'_track_1', self.name+'_gap_1', self.name+'_mask_1']
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


    def draw_adaptor(self, iSlope=0.33):
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
        if not self.isOut:
            # calculate the output
            # do not forget to add the new port to dict
            adaptDist = abs(self.outTrack/2-self.inTrack/2)/iSlope
            outPort = [self.pos+self.ori*adaptDist, self.ori, self.outTrack, self.outGap]
            self.ports[self.iIn+'_bis'] = outPort
            self.__init__(self.name, self.iIn, self.iIn+'_bis')
        else:
            adaptDist = (self.pos-self.posOut).norm()



        points = self.append_points([(0, self.inTrack/2),
                                     (adaptDist, self.outTrack/2-self.inTrack/2),
                                     (0, -self.outTrack),
                                     (-adaptDist, self.outTrack/2-self.inTrack/2)])
        track = self.draw(self.name+"_track", points)
        self.trackObjects.append(track)

        points = self.append_points([(0, self.inGap+self.inTrack/2),
                                     (adaptDist, (self.outGap-self.inGap)+(self.outTrack-self.inTrack)/2),
                                     (0, -2*self.outGap-self.outTrack),
                                     (-adaptDist, (self.outGap-self.inGap)+(self.outTrack-self.inTrack)/2)])
        gap = self.draw(self.name+"_gap", points)
        self.gapObjects.append(gap)
        
        mask = None
        if self.is_mask:
            points = self.append_points([(0, self.gap_mask+self.inGap+self.inTrack/2),
                             (adaptDist, (self.outGap-self.inGap)+(self.outTrack-self.inTrack)/2),
                             (0, -2*self.gap_mask-2*self.outGap-self.outTrack),
                             (-adaptDist, (self.outGap-self.inGap)+(self.outTrack-self.inTrack)/2)])
            mask = self.draw(self.name+"_mask", points)
            self.maskObjects.append(mask)

        return self.iIn+'_bis', adaptDist, track, gap, mask
    
    def draw_dc_adaptor(self, iSlope=0.15):
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
        if not self.isOut:
            # calculate the output
            # do not forget to add the new port to dict
            # instead of super(), self would do the trick but we want to show that the variable is in the superclass
            # indeed the next part can be omitted if c.get_varible_value is placed when function is called
#            print(super(ConnectElt,self).get_variable_names())
            maybeString = [self.rel_posIn, self.rel_posOut, self.widIn, self.widOut]
            isnoString = []
            for kk in range(4):
                isnoS = []
                for ii in range(self.multOut):
                    if maybeString[kk][ii] in super().get_variable_names():
                        isnoS.append(parse_entry(super(ConnectElt,self).get_variable_value(maybeString[kk][ii])))
                    else:
                        isnoS.append(maybeString[kk][ii])
                isnoString.append(isnoS)
            rel_posIn, rel_posOut, widIn, widOut = isnoString  
            adaptDist = max([max([abs(widOut[ii]/2-widIn[ii]/2)/iSlope,\
                                  abs(self.rel_posOut[ii]/2-self.rel_posIn[ii]/2)/iSlope/2])\
                             for ii in range(self.multOut)])
            
            
            outPort = [self.pos+self.ori*adaptDist, self.ori, self.rel_posOut, self.widOut, self.cutOut, self.multOut] #DC type
            self.ports_dc[self.iIn+'_bis'] = outPort
            self.__init__(self.name, self.iIn, self.iIn+'_bis', layer=self.layer)
        else:
            adaptDist = (self.pos-self.posOut).norm()
        
        tracks = []
        for ii in range(len(self.rel_posOut)):
#            DEBUG
#            vec_center_in = Vector(self.pos)
#            vec_center_out = Vector(self.posOut)
#            vec_rel_in = Vector([0, self.rel_posIn[ii]])
#            vec_rel_out = Vector([0, self.rel_posOut[ii]])
#            print(self.ori)
#            vec_rel_in = vec_rel_in.rot(self.ori)
#            vec_rel_out = vec_rel_out.rot(-self.oriOut)
#            self.draw_rect_center(self.name+'subportin_'+str(ii), vec_center_in+vec_rel_in, Vector(['10um','10um']))
#            self.draw_rect_center(self.name+'subportout_'+str(ii), vec_center_out+vec_rel_out, Vector(['10um','10um']))
            
            points = self.append_absolute_points([(0, self.rel_posIn[ii]-self.widIn[ii]/2),
                                                  (0, self.rel_posIn[ii]+self.widIn[ii]/2),
                                                  (adaptDist, -self.rel_posOut[ii]+self.widOut[ii]/2),
                                                  (adaptDist, -self.rel_posOut[ii]-self.widOut[ii]/2)])
            tracks.append(self.draw(self.layer+'_'+self.name+"_track_"+str(ii), points))
        points = self.append_absolute_points([(0, -self.cutIn/2),
                                              (0, self.cutIn/2),
                                              (adaptDist, +self.cutOut/2), 
                                              (adaptDist, -self.cutOut/2)])
        cutout = self.draw(self.layer+'_'+self.name+'_cutout', points)
        
        if len(tracks)>1:
            tracks = self.unite(tracks, name=self.layer+'_'+self.name+'_track')
        else:
            tracks = self.rename(tracks[0], self.layer+'_'+self.name+'_track')
        self.layers[self.layer]['trackObjects'].append(tracks)
        self.layers[self.layer]['gapObjects'].append(cutout)

        return self.iIn+'_bis', adaptDist, tracks
    
    
    def draw_dc_cable(self, layer_name, fillet="0.5mm", constrains=[], reverse=True, iSlope=0.15):
        
        fillet = parse_entry(fillet)
        
        self.to_bond=[]
        adaptor_length=0
        track_adaptor = None
        
        if not self.multIn == self.multOut:
            raise ValueError('input and output ports do not have same multiplicity')
            
        testwid = [equal_float(self.val(self.widIn[ii]), self.val(self.widOut[ii])) for ii in range(self.multIn)]
        testpos = [equal_float(self.val(self.rel_posIn[ii]), self.val(self.rel_posOut[ii])) for ii in range(self.multIn)]
        if any(x==False for x in testwid) or any(x==False for x in testpos):
            if not reverse:
                adaptor = ConnectElt(self.name+'_adaptor', self.iOut, [self.rel_posIn, self.widIn, self.cutIn, self.multIn], layer=layer_name) #connect_elt dc type
                iOut, adaptor_length, track_adaptor = adaptor.draw_dc_adaptor() #, gap_adaptor
                self.__init__(self.name, self.iIn, iOut)
            else:
                adaptor = ConnectElt(self.name+'_adaptor', self.iIn, [self.rel_posOut, self.widOut, self.cutOut, self.multOut], layer=layer_name) #connect_elt dc type
                iIn, adaptor_length, track_adaptor = adaptor.draw_dc_adaptor(iSlope=iSlope) #, gap_adaptor
                self.ports_dc[iIn][2] = list(-np.array(self.ports_dc[iIn][2]))
                self.__init__(self.name, iIn, self.iOut)
        all_constrains = []
#        for constrain in constrains:
#            all_constrains.append([self.ports_dc[constrain][POS], -self.ports_dc[constrain][ORI]]+[self.ports_dc[constrain][ii] for ii in range(4)])
#            all_constrains.append([self.ports_dc[constrain][POS], self.ports_dc[constrain][ORI]+[self.ports_dc[constrain][ii] for ii in range(4)])
              
        cables = []
        port_names = [self.iIn]+all_constrains+[self.iOut]

        for ii in range(len(constrains)+1):
            if len(constrains)!=0:
                to_add = '_'+str(ii)
            else:
                to_add = ''
            self.__init__(self.name, *port_names[2*ii:2*ii+2])
            
            points = self.find_path(fillet, is_meander=False, to_meander=[], meander_length=0, meander_offset=0)
            connection = self.draw(self.name+'_track'+to_add, points, closed=False)
            connection.fillets(fillet-eps)
            
            connection_track = [connection.copy(layer_name+'_'+self.name+"_track"+to_add+str(ii)) for ii in range(self.multIn)]
            track_starter = [self.cable_starter('dc_track', index=ii) for ii in range(self.multIn)]
            connection_cutout = connection.copy(layer_name+'_'+self.name+'_cutout')
            cutout_starter = self.cable_starter('dc_cutout')
            for ii in range(len(self.rel_posIn)):
                cables.append(connection_track[ii].sweep_along_path(track_starter[ii]))
            cutout = connection_cutout.sweep_along_path(cutout_starter)
            
#            ''' what the hell the .copy routine ??? and .fillet ? '''
#            connection = []
#            for ii in range(self.multIn):
#                con = self.draw(layer_name+'_'+self.name+'_track'+to_add+'_'+str(ii), points, closed=False)
#                con.fillets(fillet-eps)
#                connection.append(con)
#            con = self.draw(layer_name+'_'+self.name+'_cutout', points, closed=False)
#            con.fillets(fillet-eps)
#            connection.append(con)   
#            track_starter = [self.cable_starter('dc_track', index=ii) for ii in range(self.multIn)]
#            cutout_starter = self.cable_starter('dc_cutout')
#            for ii in range(self.multIn):
#                cables.append(self.sweep_along_path(track_starter[ii], connection[ii]))
#            cutout = self.sweep_along_path(cutout_starter, connection[-1])

        if len(cables)>1:
            cable = self.unite(cables, layer_name+'_'+self.name+"_track")
        else:
            cable = self.rename(cables[0], layer_name+'_'+self.name+"_track")
        self.layers[layer_name]['trackObjects'].append(cable)
        self.layers[layer_name]['gapObjects'].append(cutout)    
            
    
    def _connect_JJ(self, iTrackJ, iInduct='1nH', fillet=None):
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
        iLength = (self.posOut-self.pos).norm()
        
        if 0:
            induc_H = self.val(iInduct)*1e-9
            print('induc'+str(induc_H))
            w_plasma = 2*np.pi*24*1e9
            capa_fF = 1/(induc_H*w_plasma**2)*10**15
            capa_plasma = self.set_variable(self.name+'_capa_plasma', str(capa_fF)+'fF')
            print(capa_fF)
            print(capa_plasma)
        else:
            capa_plasma = 0
        
        
        # No parsing needed, should not be called from outside
        self.pos = (self.pos+self.posOut)/2
        iTrack = self.inTrack # assume both track are identical
        
        adaptDist = iTrack/2-iTrackJ/2

        if self.val(adaptDist)>self.val(iLength/2-iTrackJ/2):
            raise ValueError('Increase iTrackJ %s' % self.name)


        raw_points = [(iTrackJ/2, iTrackJ/2),
                      ((iLength/2-iTrackJ/2-adaptDist), 0),
                      (adaptDist, (iTrack-iTrackJ)/2),
                      (0, -iTrack),
                      (-adaptDist, (iTrack-iTrackJ)/2),
                      (-(iLength/2-iTrackJ/2-adaptDist), 0)]
        points = self.append_points(raw_points)
        right_pad = self.draw(self.name+"_pad1", points)
        
            

        points = self.append_points(self.refy_points(raw_points))
        left_pad = self.draw(self.name+"_pad2", points)

        pads = self.unite([right_pad, left_pad], name=self.name+'_pads')
        
        if not self.is_litho:
            mesh = self.draw_rect_center(self.name+'_mesh', self.coor([0,0]), self.coor_vec([iLength, iTrack]))
            self.modeler.assign_mesh_length(mesh, iTrackJ/2)
    
            points = self.append_points([(iTrackJ/2,0),(-iTrackJ,0)])
            self.draw(self.name+'_line', points, closed=False)
    
            JJ = self.draw_rect_center(self.name, self.coor([0,0]), self.coor_vec([iTrackJ, iTrackJ]))
            self.assign_lumped_RLC(JJ, self.ori, (0, iInduct, capa_plasma))

        return pads
    
    def _connect_snails(self, squid_size, width_top, n_top, width_bot, n_bot, N, width_bridge):
        
        width_track = self.inTrack # assume both are equal
#        print(width_track)
        spacing = (self.posOut-self.pos).norm()
#        print(spacing)
        
        self.pos = (self.pos+self.posOut)/2
        
        width_snail = squid_size[0]+2*width_track
        
        tot_width = width_snail*N
        if np.abs(self.val(tot_width))>np.abs(self.val(spacing)):
            raise ValueError("cannot put all snail in given space")
        snails = []
        snails.append(self.draw_rect(self.name+'_pad_left', self.coor([-tot_width/2, -width_track/2]), self.coor_vec([-(spacing-tot_width)/2-5e-6,width_track])))
        snails.append(self.draw_rect(self.name+'_pad_right', self.coor([tot_width/2, -width_track/2]), self.coor_vec([(spacing-tot_width)/2+5e-6,width_track])))
        if N%2==1:
            x_pos=-(N//2)*width_snail
        else:
            x_pos=-(N//2-1/2)*width_snail

        for jj in range(int(N)):
            snail=[]
            snail.append(self.draw_rect(self.name+'_left', self.coor([x_pos-squid_size[0]/2-width_track/2, (-squid_size[1]/2-width_bot)-((-squid_size[1]/2-width_bot)+width_track/2)]), self.coor_vec([width_track/2, squid_size[1]+width_top+width_bot])))
            snail.append(self.draw_rect(self.name+'_right', self.coor([x_pos+squid_size[0]/2, -squid_size[1]/2-width_bot-((-squid_size[1]/2-width_bot)+width_track/2)]), self.coor_vec([width_track/2, squid_size[1]+width_top+width_bot])))
            for width, n, way in [[width_top, n_top, 1], [width_bot, n_bot, -1]]:
                if n==1:
                    snail.append(self.draw_rect(self.name+'_islandtop_left', self.coor([x_pos-squid_size[0]/2, way*squid_size[1]/2-((-squid_size[1]/2-width_bot)+width_track/2)]), self.coor_vec([(squid_size[0]-width_bridge)/2, way*width])))
                    snail.append(self.draw_rect(self.name+'_islandtop_right', self.coor([x_pos+squid_size[0]/2, way*squid_size[1]/2-((-squid_size[1]/2-width_bot)+width_track/2)]), self.coor_vec([-(squid_size[0]-width_bridge)/2, way*width])))
                elif n>1:
                    length_island = (squid_size[0]-n*width_bridge)/(n-1)
                    for ii in range(n_top-1):
                        snail.append(self.draw_rect(self.name+'_islandtop_'+str(ii), self.coor([x_pos-squid_size[0]/2+width_bridge+ii*(length_island+width_bridge), way*squid_size[1]/2-((-squid_size[1]/2-width_bot)+width_track/2)]), self.coor_vec([length_island, way*width])) )
            snail.append(self.draw_rect(self.name+'_connect_left', self.coor([x_pos-squid_size[0]/2-width_track/2, -width_track/2]), self.coor_vec([-width_track/2, width_track])))
            snail.append(self.draw_rect(self.name+'_connect_right', self.coor([x_pos+squid_size[0]/2+width_track/2, -width_track/2]), self.coor_vec([width_track/2, width_track])))
            
            self.unite(snail, name=self.name+'_snail_'+str(jj))
            x_pos = x_pos+width_snail
            
    def _connect_snails2(self, squid_size, width_top, n_top, width_bot, n_bot, N, width_bridge, spacing_bridge, litho='opt'):
        
        width_track = self.inTrack # assume both are equal
#        print(width_track)
        spacing = (self.posOut-self.pos).norm()
#        print(spacing)
        
        self.pos = (self.pos+self.posOut)/2
        
        width_snail = squid_size[0]+6*width_track #ZL
        
        tot_width = width_snail*N
        if np.abs(self.val(tot_width))>np.abs(self.val(spacing)):
            raise ValueError("cannot put all snail in given space")
        offset = 5.0e-6-0.630e-6
        overlap=0
        if litho=='elec':
            overlap = 5e-6
        snails = []
        snails.append(self.draw_rect(self.name+'_pad_left', self.coor([-tot_width/2, -width_track/2]), self.coor_vec([-(spacing-tot_width)/2-overlap,width_track])))
        snails.append(self.draw_rect(self.name+'_pad_right', self.coor([tot_width/2, -width_track/2]), self.coor_vec([(spacing-tot_width)/2+overlap,width_track])))
        if N%2==1:
            x_pos=-(N//2)*width_snail
        else:
            x_pos=-(N//2-1/2)*width_snail

        for jj in range(int(N)):
            snail=[]
            snail.append(self.draw_rect(self.name+'_left', self.coor([x_pos-squid_size[0]/2, -0.1e-6+(-squid_size[1]/2-width_bot+0.1e-6)-((-squid_size[1]/2-width_bot+0.1e-6)+width_track/2)-offset]), self.coor_vec([-3*width_track/2, squid_size[1]+width_top+width_bot+2*0.1e-6])))
            snail.append(self.draw_rect(self.name+'_right', self.coor([x_pos+squid_size[0]/2, -squid_size[1]/2-((-squid_size[1]/2)+width_track/2)-offset]), self.coor_vec([3*width_track/2, squid_size[1]+width_top+width_bot+0.1e-6])))
            for width, n, way in [[width_top, n_top, 1], [width_bot, n_bot, -1]]:
                if n==1:
                    snail.append(self.draw_rect(self.name+'_islandtop_left', self.coor([x_pos-squid_size[0]/2, way*squid_size[1]/2-((-squid_size[1]/2-width_bot-0.1e-6)+width_track/2)-offset]), self.coor_vec([(squid_size[0]-width_bridge)/2, way*width-2*0.1e-6])))
                    snail.append(self.draw_rect(self.name+'_islandtop_right', self.coor([x_pos+squid_size[0]/2, way*squid_size[1]/2-((-squid_size[1]/2-width_bot)+width_track/2)-offset]), self.coor_vec([-(squid_size[0]-width_bridge)/2, way*width])))
                if n==3: #TODO
                    snail.append(self.draw_rect(self.name+'_islandtop_left', self.coor([x_pos-squid_size[0]/2, way*squid_size[1]/2-((-squid_size[1]/2-width_bot+0.1e-6)+width_track/2)-offset]), self.coor_vec([(squid_size[0]-2*width_bridge-spacing_bridge)/2, way*width+2*0.1e-6])))
                    snail.append(self.draw_rect(self.name+'_islandtop_right', self.coor([x_pos+squid_size[0]/2, way*squid_size[1]/2-((-squid_size[1]/2-width_bot+0.1e-6)+width_track/2)-offset]), self.coor_vec([-(squid_size[0]-2*width_bridge-spacing_bridge)/2, way*width+2*0.1e-6])))
                    snail.append(self.draw_rect(self.name+'_island_middle_', self.coor([x_pos-spacing_bridge/2, way*squid_size[1]/2-((-squid_size[1]/2-width_bot)+width_track/2)-offset]), self.coor_vec([spacing_bridge, way*width])) )
                if n==4:
                    snail.append(self.draw_rect(self.name+'_islandtop_left', self.coor([x_pos-squid_size[0]/2, way*squid_size[1]/2-((-squid_size[1]/2-width_bot+0.1e-6)+width_track/2)-offset]), self.coor_vec([(squid_size[0]-3*width_bridge-2*spacing_bridge)/2, way*width+2*0.1e-6])))
                    snail.append(self.draw_rect(self.name+'_islandtop_right', self.coor([x_pos+squid_size[0]/2, way*squid_size[1]/2-((-squid_size[1]/2-width_bot)+width_track/2)-offset]), self.coor_vec([-(squid_size[0]-3*width_bridge-2*spacing_bridge)/2, way*width])))
                    snail.append(self.draw_rect(self.name+'_island_middle_left', self.coor([x_pos-spacing_bridge-width_bridge/2, way*squid_size[1]/2-((-squid_size[1]/2-width_bot)+width_track/2)-offset]), self.coor_vec([spacing_bridge, way*width])) )
                    snail.append(self.draw_rect(self.name+'_island_middle_right', self.coor([x_pos+spacing_bridge+width_bridge/2, way*squid_size[1]/2-((-squid_size[1]/2-width_bot+0.1e-6)+width_track/2)-offset]), self.coor_vec([-spacing_bridge, way*width+2*0.1e-6])) )
            snail.append(self.draw_rect(self.name+'_connect_left', self.coor([x_pos-squid_size[0]/2-3*width_track/2, -width_track/2]), self.coor_vec([-1.5*width_track, width_track]))) #ZL
            snail.append(self.draw_rect(self.name+'_connect_right', self.coor([x_pos+squid_size[0]/2+3*width_track/2, -width_track/2]), self.coor_vec([1.5*width_track, width_track])))
            
            self.unite(snail, name=self.name+'_snail_'+str(jj))
            x_pos = x_pos+width_snail
            
    def _connect_snails_Zaki(self, squid_size, width_top, n_top, width_bot, n_bot, N, width_bridge, spacing_bridge, litho='opt'):
        
        width_track = self.inTrack # assume both are equal
        print(width_track)
        spacing = (self.posOut-self.pos).norm()
        print(spacing)
        
        self.pos = (self.pos+self.posOut)/2
        
        width_snail = squid_size[0]+1*width_track #ZL
        
        tot_width = width_snail*N
        if np.abs(self.val(tot_width))>np.abs(self.val(spacing)):
            raise ValueError("cannot put all snail in given space")
        
        overlap=0
        if litho=='elec':
            overlap = 5e-6
        snails = []
        snails.append(self.draw_rect(self.name+'_pad_left', self.coor([-tot_width/2-width_track/2, -width_track/2]), self.coor_vec([-(spacing-tot_width)/2-overlap+width_track/2,width_track])))
        snails.append(self.draw_rect(self.name+'_pad_right', self.coor([tot_width/2+width_track/2, -width_track/2]), self.coor_vec([(spacing-tot_width)/2+overlap-width_track/2,width_track])))
        if N%2==1:
            x_pos=-(N//2)*width_snail
        else:
            x_pos=-(N//2-1/2)*width_snail

        for jj in range(int(N)):
            snail=[]
            snail.append(self.draw_rect(self.name+'_left',  self.coor([x_pos-squid_size[0]/2, -width_track/2]), self.coor_vec([-2*width_track/2, squid_size[1]+width_top+width_bot+2*0.1e-6]))) #ZL
            for width, n, way in [[width_top, n_top, 1], [width_bot, n_bot, -1]]:
                if n==1:
                    snail.append(self.draw_rect(self.name+'_islandtop_left', self.coor([x_pos-squid_size[0]/2, -width_track/2]), self.coor_vec([(squid_size[0]-width_bridge)/2, width+2*0.1e-6])))
                    snail.append(self.draw_rect(self.name+'_islandtop_right', self.coor([x_pos+squid_size[0]/2, -width_track/2+0.1e-6]), self.coor_vec([-(squid_size[0]-width_bridge)/2, width])))
                    #self.draw_rect(self.name+'_Lj_'+str(jj))
                if n==3: #TODO
                    snail.append(self.draw_rect(self.name+'_islandtop_left', self.coor([x_pos-squid_size[0]/2, squid_size[1]+width_top+width_bot-width_track/2+2*0.1e-6]), self.coor_vec([(squid_size[0]-2*width_bridge-spacing_bridge)/2, -(squid_size[1]+width_top)])))
                    snail.append(self.draw_rect(self.name+'_islandtop_right', self.coor([x_pos+squid_size[0]/2, squid_size[1]+width_top+width_bot-width_track/2+2*0.1e-6]), self.coor_vec([-(squid_size[0]-2*width_bridge-spacing_bridge)/2, -(squid_size[1]+width_top)-1*0.1e-6])))
                    snail.append(self.draw_rect(self.name+'_island_middle_', self.coor([x_pos-spacing_bridge/2, squid_size[1]+width_top+width_bot-width_track/2+0.1e-6]), self.coor_vec([spacing_bridge, -width])))
            #snail.append(self.draw_rect(self.name+'_connect_left', self.coor([x_pos-squid_size[0]/2-1*width_track/2, -width_track/2]), self.coor_vec([-1.5*width_track, width_track]))) #ZL
            #snail.append(self.draw_rect(self.name+'_connect_right', self.coor([x_pos+squid_size[0]/2+1*width_track/2, -width_track/2]), self.coor_vec([1.5*width_track, width_track])))

            #self.unite(snail, name=self.name+'_snail_'+str(jj))
            x_pos = x_pos+width_snail
        x_pos = x_pos-width_snail
        snail.append(self.draw_rect(self.name+'_right', self.coor([x_pos+squid_size[0]/2, -width_track/2]), self.coor_vec([2*width_track/2, squid_size[1]+width_top+width_bot+2*0.1e-6]))) #ZL

    def _connect_jct(self, width_bridge, n=1, spacing_bridge=0, assymetry=0.1e-6, overlap=5e-6, width_jct=None): #opt assymetry=0.25e-6
        limit_dose = 8e-6
        width = self.inTrack # assume both are equal
        spacing = (self.posOut-self.pos).norm()
        self.pos = (self.pos+self.posOut)/2
        n = int(n)
        
        tot_width = n*width_bridge+(n-1)*spacing_bridge
        
        self.draw_rect(self.name+'_left', self.coor([-tot_width/2-limit_dose,-width/2-assymetry]), self.coor_vec([-(spacing-tot_width)/2-overlap+limit_dose, width+2*assymetry]))
        self.draw_rect(self.name+'_left2', self.coor([-tot_width/2,-width/2-assymetry]), self.coor_vec([-limit_dose, width+2*assymetry]))

        print(n)
        if n%2==0:
            _width_right = width+2*assymetry
        else:
            _width_right = width
        if width_jct is not None:
            margin = 1e-6
            self.draw_rect(self.name+'_right2', self.coor([tot_width/2+margin,-_width_right/2]), self.coor_vec([limit_dose-margin, _width_right]))
            self.draw_rect(self.name+'_right3', self.coor([tot_width/2,-width_jct/2]), self.coor_vec([margin, width_jct]))
        else:
            self.draw_rect(self.name+'_right2', self.coor([tot_width/2,-_width_right/2]), self.coor_vec([limit_dose, _width_right]))

        self.draw_rect(self.name+'_right', self.coor([tot_width/2+limit_dose,-_width_right/2]), self.coor_vec([(spacing-tot_width)/2+overlap-limit_dose, _width_right]))

        x_pos = -(tot_width)/2+width_bridge
        for ii in range(n-1):
            if ii%2==1:
                _width_loc = width+2*assymetry
            else:
                _width_loc = width
            self.draw_rect(self.name+'_middle', self.coor([x_pos,-_width_loc/2]), self.coor_vec([spacing_bridge, _width_loc]))
            x_pos = x_
            pos+spacing_bridge+width_bridge
