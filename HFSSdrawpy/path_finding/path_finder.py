#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 14 18:18:44 2020

@author: Raphael
"""

from ..utils import VariableString, \
                             extract_value_unit, \
                             extract_value_dim, \
                             parse_entry, \
                             _val, val, \
                             way, \
                             Vector
import numpy as np


# useful function to find cable path
def next_point(point1, point2, vec):
    choice1 = point1+vec*(point2-point1).dot(vec)
    choice2 = point1+vec.orth()*(point2-point1).dot(vec.orth())
    return [[point1, choice1, point2], [point1, choice2, point2]]

def cost_f(x):
    if x==1:
        return 0
    elif x==0:
        return 1
    else:
        return 100

def right_left(points):
    # tells if given point is turning left or right
    vecs = []
    A = points[0]
    for B in points[1:]:
        vecs.append(way(val(B-A)))
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

def displace(points, rl, min_dist, displacement=0, offset=0, n_meander=-1):
    # points is simply a segment
    # rl indicates what type of corner we have on each end
    if np.abs(val(displacement))<val(min_dist)*1.1:
        print('here')
        displacement = min_dist*1.1

    elts = add_points(points, rl, min_dist, n_meander=n_meander)
    points, indices_corners, dist, ignore, n_add = elts

    new_points = [points[0]]
    parity = 1
    if indices_corners is not None:
        print("indices_corners is not None")
    #     for ii, B in enumerate(points[1:-1]):
    #         A = points[ii]
    #         AB= B-A
    #         vec = way(val(AB))
    #         if ii==0:
    #             parity = (-2*((indices_corners[0]-(ii+1))%2)+1)*(-rl[0])
    #         else:
    #             parity = -parity

    #         if ii+1 not in indices_corners:
    #             #regular point
    #             new_points[ii+1] = points[ii+1]+vec.orth()*parity*min_dist
    #         else:
    #             new_points[ii+1] = points[ii+1]+(vec.orth()*parity+vec).unit()*min_dist
    # else:

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
        vec=way(val(AB))
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

    return new_points, n_add

def add_points(points, rl, min_dist, n_meander):
    min_dist = min_dist*1.1
    n_points = len(points)
    if n_points!=2:
        print('n_points!=2')
    A = points[0]
    new_points =[A]
    B=points[-1]
    vec = way(val(B-A))
    AB = (B-A).norm()
    n_add = int(val(AB/min_dist))
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
    # new_points+=[A+vec*dist/2]
    new_points+=[A+vec*dist*(jj+1/2) for jj in range(n_add)]
    rl=None
    indices_corners=None

    print('n_add', n_add)

    if ignore:
        new_points.append(addedB)
    new_points.append(B)
    return new_points, indices_corners, dist, ignore, n_add



class Path(object):

    def __init__(self, name, port_in, port_out, fillet, points=[]):
        self.name = name
        self.port_in = port_in
        self.port_out = port_out
        self.fillet = fillet
        self.points = points

        if points==[]:
            in_pos = Vector(port_in.pos)
            in_ori = Vector(port_in.ori)
            out_pos = Vector(port_out.pos)
            out_ori = Vector(port_out.ori)
            room_bonding = 0*100e-6 #SMPD MANU BOND SPACE

            dist_y = (out_pos-in_pos).dot(in_ori.orth())
            self.slanted=False
            if in_ori.dot(out_ori)==-1 and abs(val(dist_y))<val(2*fillet):
                self.slanted=True  # need to do a slanted_path

            pointA = in_pos+in_ori*room_bonding  # first and last point
            pointB = out_pos+out_ori*room_bonding
            point1 = in_pos+in_ori*(1.1*fillet+room_bonding)  # after in
            point2 = out_pos+out_ori*(1.1*fillet+room_bonding)  # before out

            points_choices = []
            if in_ori.dot(out_ori)==-1:
                middle_point = (point1 + point2)/2

                choice_in = next_point(point1, middle_point, in_ori) #bon sens
                choice_out = next_point(point2, middle_point, out_ori) #à inverser
                for c_in in choice_in:
                    for c_out in choice_out:
                        points_choices.append([pointA, *c_in, *c_out[:-1][::-1], pointB])
            else:
                choice_in = next_point(point1, point2, in_ori)
                for c_in in choice_in:
                    points_choices.append([pointA, *c_in, pointB])

            final_choice= None
            cost=np.inf
            for ii, choice in enumerate(points_choices):
                new_cost, new_choice = self.clean(points=choice)
                if new_cost<cost:
                    final_choice = new_choice
                    cost = new_cost
            self.points = final_choice

    def __add__(self, other):
        assert(isinstance(other, Path))
        if self.points[-1] == other.points[0]:
            points = self.points + other.points[1:]
            return Path(self.name, self.port_in, other.port_out, self.fillet,
                        points=points)
        elif self.points[-1] == other.points[-1]:
            points = self.points + other.points[:-1][::-1]
            return Path(self.name, self.port_in, other.port_in, self.fillet,
                        points=points)
        elif self.points[0] == other.points[-1]:
            points = other.points + self.points[1:]
            return Path(self.name, other.port_in, self.port_out, self.fillet,
                        points=points)
        elif self.points[-1] == other.points[-1]:
            points = other.points + self.points[:-1][::-1]
            return Path(self.name, other.port_in, self.port_in, self.fillet,
                        points=points)
        else:
            raise ValueError('Added path do not coincide on one point')

    def clean(self, points=None):
        # if nothing is given, cleans the path and raise an exception if the path
        # is invalid
        # if points are given, cleans the points and returns the cost of the path
        # points : series of points that represent the cable path
        # returns : the cost (number of turns) and make sure each point is at a
        # turn ie that there are no consecutive segments in the same direction.

        if points is None:
            cleaning_only = True
            points = self.points.copy()
        else:
            cleaning_only = False

        prev_point = points[0]
        _points = [points[0]]
        vecs = []
        for ii, point in enumerate(points):
            if ii==0:
                pass
            else:
                if not point == prev_point:
                    _points.append(point)
                    vecs.append(way(val(point-prev_point)))
                    prev_point = point

        cost = 0
        points = _points.copy()
        new_points = [points[0]]
        prev_vec = vecs[0]
        for ii, curr_vec in enumerate(vecs[1:]):
            if curr_vec.dot(prev_vec)==0:
                new_points.append(points[ii+1])
            elif curr_vec.dot(prev_vec)==-1 and cleaning_only:
                msg = 'Provided path is invalid: the path goes back on \
                       itself, 180° corner'
                raise ValueError(msg)
            cost += cost_f(prev_vec.dot(curr_vec))
            prev_vec = curr_vec
        new_points.append(points[-1])

        if cleaning_only:
            self.points = new_points
        else:
            return cost, new_points

    # handling slanted parts
    def slanted(self):
        points = self.points.copy()
        if slanted:
            h = abs(dist_y/2)
            x = (fillet-h)/fillet
            d = h*x/(1-x**2)**0.5
            dist_x = (out_pos-in_pos).dot(in_ori)
            pointA = in_pos+in_ori*(dist_x/2-d)
            pointB = out_pos+out_ori*(dist_x/2-d)

            to_bond.append([in_pos, pointA])
            to_bond.append([pointB, out_pos])

            return [in_pos, pointA, pointB, out_pos]


    def to_bond(self):
        points = self.points
        fillet = self.fillet
        # TODO slanted path should not have bond
        bonding_segments = [[points[0]]]
        for ii, point in enumerate(points[1:-1]):
            index = ii+1
            p_vec = points[index-1]-point
            n_vec = points[index+1]-point
            bonding_segments[-1].append(point+way(val(p_vec))*fillet)
            bonding_segments.append([point+way(val(n_vec))*fillet])
        bonding_segments[-1].append(points[-1])
        return bonding_segments

    def meander(self, to_meander, meander_length, meander_offset): # to_meander is list of segments to be meander
        min_dist = 2*self.fillet
        points = self.points.copy()
        n_points = len(points)
        n_to_meander = len(to_meander)
        if n_points-1>n_to_meander:
            to_meander = to_meander+[0 for ii in range(n_points-1-n_to_meander)]
        else:
            to_meander = to_meander[:n_points-1]
        working_p, left_p, index_insertion = self.working_points(points,
                                                                 min_dist,
                                                                 to_meander)
        # # create adjustable variable for meander_length if it is not the case
        # if isinstance(meander_length, VariableString):
        #     meander_length_name = str(meander_length)
        # else:
        #     meander_length_name = self.name+'_meander_length'
        #     meander_length = VariableString(meander_length_name, value=1.1*min_dist)
        # VariableString.variables[meander_length_name] = 1.1*min_dist

        tot_add = 0 # number of added meanders
        if len(working_p) != 0:
            rl = right_left(working_p)
            working_ps = []

            for ii, isit in enumerate(to_meander):
                if isit!=0:
                    new_working_p, n_add = displace(working_p[ii:ii+2], rl[ii:ii+2], min_dist, displacement=meander_length, offset=meander_offset, n_meander=isit) # n_meander=-1 -> auto
                else:
                    new_working_p = working_p[ii:ii+2]
                    n_add = 0
                tot_add += n_add
                working_ps += new_working_p
            left_p[index_insertion:index_insertion] = working_ps

        # VariableString.variables[meander_length_name] = 1.1*self.fillet

        self.points = [self.points[0]] + left_p + [self.points[-1]]

    def working_points(self, points, min_dist, to_meander):
        min_dist = min_dist*1.1

        left_p_start=[points[0]]
        left_p_end=[points[-1]]

        index_start = 0
        for ii, point in enumerate(points[1:]):
            A = left_p_start[-1]
            B = point
            AB = B-A
            vec = way(val(B-A))
            if val(AB).norm() > val(min_dist):
                working_p_start = A+vec*min_dist/2
                index_start = ii+1
                break
            else:
                left_p_start.append(B)
                to_meander.pop(0)
        else:
            print('Warning: Could not find points to elongate cable %s' %self.name)
            left_p = left_p_start+left_p_end[::-1]
            return [], left_p, 0

        index_end = 0
        for ii, point in enumerate(points[::-1][1:]):
            A = left_p_end[-1]
            B = point
            AB = B-A
            vec = way(val(B-A))
            if val(AB).norm() > val(min_dist):
                working_p_end = A+vec*min_dist/2
                index_end = ii+1
                break
            else:
                left_p_end.append(B)
                to_meander.pop(-1)
        else:
            print('Warning: Could not find points to elongate cable %s' %self.name)
            left_p = left_p_start+left_p_end[::-1]
            return [], left_p, 0

        working_p = [working_p_start]+points[index_start:-index_end]+[working_p_end]
        index_insertion = len(left_p_start)
        left_p = left_p_start+left_p_end[::-1]

        return working_p, left_p, index_insertion

    def length(self):
        self.clean() # make sure each point is at a corner
        points = self.points
        value = 0
        corner = val(self.fillet*(2-np.pi/2))
        for ii in range(len(points)-1):
            value += val((points[ii+1]-points[ii]).norm())- corner
        value += corner
        return value

