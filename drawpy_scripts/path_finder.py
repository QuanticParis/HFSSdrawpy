#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 14 18:18:44 2020

@author: Raphael
"""

from .variable_string import VariableString, \
                             extract_value_unit, \
                             extract_value_dim, \
                             parse_entry, \
                             _val, val, \
                             var, \
                             Vector
class Path(object):

    def __init__(self, port_in, port_out, fillet):
        if value is not None:
            VariableString.store_variable(self, value)


    # useful function to find cable path
    def next_point(point1, point2, vec):
        choice1 = point1+vec*var((point2-point1).dot(vec))
        choice2 = point1+vec.orth()*var((point2-point1).dot(vec.orth()))
        return [[point1, choice1, point2], [point1, choice2, point2]]

    def check(points):
        # points : series of points that represent the cable path
        # returns : the cost (number of turns) and make sure each point is at a
        # turn ie that there are no consecutive segments in the same direction.
        length = 0
        prev_point = Vector(points[0])
        _points = [Vector(points[0])]
        vecs = []
        for point in points[1:]:
            if not ( equal_float(point[0], prev_point[0]) and equal_float(point[1], prev_point[1])):
               # TODO vec should already be a Vector
               vec = point-prev_point
               length += vec.norm()
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

    def cost_f(x):
        if x==1:
            return 0
        elif x==0:
            return 1
        else:
            return 100

    @move_port
    def find_simple_path(self, name, iInPort, iOutPort, fillet):

        iIn_pos = Vector(iInPort.pos)
        iIn_ori = Vector(iInPort.ori)
        iOut_pos = Vector(iOutPort.pos)
        iOut_ori = Vector(iOutPort.ori)
        room_bonding = 0*100e-6 #SMPD MANU BOND SPACE

        dist_y = (iOut_pos-iIn_pos).dot(iIn_ori.orth())
        slanted=False
        if iIn_ori.dot(iOut_ori)==-1 and abs(val(dist_y))<val(2*fillet):
            slanted=True  # need to do a slanted_path

        pointA = val(iIn_pos+iIn_ori*room_bonding)
        pointB = val(iOut_pos+iOut_ori*room_bonding)


        point1 = val(iIn_pos+iIn_ori*(1.1*fillet+room_bonding))
        point2 = val(iOut_pos+iOut_ori*(1.1*fillet+room_bonding))

        points_choices = Vector([])
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

        final_choice= None
        cost=np.inf
        for ii, choice in enumerate(points_choices):
            new_cost, new_choice, new_length = check(choice)

            if new_cost<cost:
                final_choice = new_choice
                cost = new_cost
                length = new_length

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
                vec = way(val(B-A))
                if val(AB).norm() > val(min_dist):
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
                    vec = way(val(B-A))
                    if val(AB).norm() > val(min_dist):
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

        def add_points(points, rl, min_dist, n_meander):
            min_dist = min_dist*1.1
            n_points = len(points)

            A = points[0]
            new_points =[]
            if n_points==2:
                new_points.append(A)
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
                new_points+=[A+vec*dist/2]
                new_points+=[A+vec*dist*(jj+1+1/2) for jj in range(n_add-1)]
                rl=None
                indices_corners=None

            if ignore:
                new_points.append(addedB)
            new_points.append(points[-1])
            return new_points, indices_corners, dist, ignore

        def displace(points, rl, min_dist, displacement=0, offset=0, n_meander=-1):
            if np.abs(val(displacement))<val(min_dist)*1.1:
                displacement = min_dist*1.1
            points, indices_corners, dist, ignore = add_points(points, rl, min_dist, n_meander=n_meander)
            new_points = [points[0]]
            parity = 1
            if indices_corners is not None:
                for ii, B in enumerate(points[1:-1]):
                    A = points[ii]
                    AB= B-A
                    vec = way(val(AB))
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

        if any(to_meander):
            min_dist = 2*fillet
            final_choice = meander(final_choice, min_dist, to_meander, meander_length, meander_offset)

        final_choice =  [val(iIn_pos)] + final_choice + [val(iOut_pos)]

# Needed to draw Manu bond
        def add_fillet_points(points, fillet):
            new_points = [points[0]]
            for ii, point in enumerate(points[1:-1]):
                index = ii+1
                p_vec = points[index-1]-point
                n_vec = points[index+1]-point
                new_points.append(point+way(val(p_vec))*fillet)
                new_points.append(point+way(val(n_vec))*fillet)
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

        _, final_choice, _ = check(final_choice)

        if slanted:
            h = abs(dist_y/2)
            x = (fillet-h)/fillet
            d = h*x/(1-x**2)**0.5
            dist_x = (iOut_pos-iIn_pos).dot(iIn_ori)
            pointA = iIn_pos+iIn_ori*(dist_x/2-d)
            pointB = iOut_pos+iOut_ori*(dist_x/2-d)

            to_bond.append([iIn_pos, pointA])
            to_bond.append([pointB, iOut_pos])

            return [iIn_pos, pointA, pointB, iOut_pos]
        else:
            to_bond_points = add_fillet_points(final_choice, fillet)
            for ii, point in enumerate(to_bond_points[::2]):
                points = [to_bond_points[2*ii], to_bond_points[2*ii+1]]
    #            print("points", points)
                # self.polyline_2D(points , closed=False, name='bef_test'+str(ii), layer=layer_Default)
                to_bond.append(points)
            return final_choice

    def find_meander_path(points, fillet, to_bond, to_meander, meander_length, meander_offset))

    def length(self, points, A, B, fillet): # A and B are integer point indices
#        for point in points:
#            print(val(point[0]), val(point[1]))
        if A<0 or A>=len(points):
            raise ValueError('First index should be within the point list')
        if B<0 or B>=len(points):
            raise ValueError('Second index should be within the point list')
        if A==B:
            return 0
        if A<B:
            value = 0
            for ii in range(B-A):
                value+=val((points[A+ii+1]-points[A+ii]).norm())
            return value-(B-A-1)*val(fillet*(2-np.pi/2))
        else:
            return self.length(points, B, A, fillet)
