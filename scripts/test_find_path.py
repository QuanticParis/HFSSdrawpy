# -*- coding: utf-8 -*-
"""
Created on Thu Apr 12 23:30:42 2018

@author: checkhov
"""

from designer import Vector
import matplotlib.pyplot as plt
import numpy as np
plt.close('all')
iIn_pos = Vector([0,0])
iIn_ori = Vector([-1,0])
iOut_pos = Vector([1, 1])
iOut_ori = Vector([0, -1])
start = [iIn_pos+ iIn_ori.orth()*0.05, iIn_pos- iIn_ori.orth()*0.05]
end = [iOut_pos+ iOut_ori.orth()*0.05, iOut_pos- iOut_ori.orth()*0.05]

points_fromIn = [iIn_pos]
points_fromOut = [iOut_pos]
point1 = iIn_pos+iIn_ori*0.1
point2 = iOut_pos+iOut_ori*0.1


def next_point(point1, point2, vec):
    choice1 = point1+(point2-point1).dot(vec)*vec
    choice2 = point1+(point2-point1).dot(vec.orth())*vec.orth()
    return [[point1, choice1, point2], [point1, choice2, point2]]


points_choices = []
if iIn_ori.dot(iOut_ori)==-1:
    middle_point = (point1 + point2)/2

    choice_in = next_point(point1, middle_point, iIn_ori) #bon sens
    choice_out = next_point(point2, middle_point, iOut_ori) #à inverser
#    print(choice_in)
#    print(choice_out)
    for c_in in choice_in:
        for c_out in choice_out:
            points_choices.append([iIn_pos, *c_in, *c_out[:-1][::-1], iOut_pos])
else:
    choice_in = next_point(point1, point2, iIn_ori)
    for c_in in choice_in:
        points_choices.append([iIn_pos, *c_in, iOut_pos])


def cost_f(x):
    if abs(x-1)<1e-2:
        return 0
    elif abs(x)<1e-2:
        return 1
    else:
        return 100

def check(points):
    points_array = np.array(points)
    vecs = points_array[:-1]-points_array[1:]

    cost = 0
    prev_vec = Vector(vecs[0]).unit()
    new_points=[]
    new_points.append(points[0])
    for ii, vec in enumerate(vecs[1:]):
        curr_vec = Vector(vec).unit()
        if curr_vec.dot(prev_vec)==0:
            new_points.append(points[ii+1])
        added_cost = cost_f(prev_vec.dot(curr_vec))
        cost += added_cost
        prev_vec = curr_vec
    new_points.append(points[-1])
    return cost, new_points

final_choice= None
cost=np.inf
for choice in points_choices:
    new_cost, new_points = check(choice)
    if new_cost<cost:
        final_choice = new_points
        cost = new_cost



###########
#def relative_append(A, vector_list):
#    new_vector_list = [A+vector_list[0]]
#    for vector in vector_list[1:]:
#        new_vector_list.append(new_vector_list[-1]+vector)
#    return new_vector_list
#
#def build_meander_list(fillet, vec, n_vec, n):
#    if n==0:
#        return [vec*fillet].extend(build_meander_list).extend([vec*fillet])
#    if n>0:
#        if n%2=0:
#            return build_meander_list(fillet, vec, n_vec, n-1)
#        else:
#            return build_meander_list(fillet, vec, -n_vec, n-1)
#
#fillet=0.1
#min_width = 4*fillet
#min_offset = 4*fillet*(np.pi/2-1) #off
#def meander(points, index, number, kind): #index = index segment on which you would like to make meanders
#    #kind will be start, middle, end
#    n_points = len(points)
#    A = points[index]
#    B = points[index+1]
#    AB = B-A
#    vec = way(AB)
#    if index==0:
#        prev_vec=None
#        next_vec = way(points[index+2]-B) # we assume we have at least one corner
#        points_to_append = relative_append(A, [vec*fillet, -next_vec*2*fillet, vec*2*fillet, next_vec*2*fillet, vec*fillet])
#        points[index+1:index+1] = points_to_append
#    elif index==n_points-1:
#        prev_vec = way(A-points[index-1])
#        next_vec=None
#    else:
#        prev_vec = way(A-points[index-1])
#        next_vec = way(points[index+2]-B)
###########

#def working_points(points, min_dist, to_meander):
#    min_dist = min_dist*1.1
#    working_p_start = []
#    working_p_end = []
#    left_p_start=[points[0]]
#    left_p_end=[points[-1]]
#    success=False
#    index_start = 0
#    for ii, point in enumerate(points[1:]):
#        A = left_p_start[-1]
#        B = point
#        AB = B-A
#        vec = way(B-A)
#        if AB.norm() > min_dist:
#            working_p_start.append(A+vec*min_dist/2)
#            success = True
#            index_start = ii+1
#            break
#        else:
#            left_p_start.append(B)
#            to_meander.pop(0)
#
#
#    if not success:
#        raise ValueError('Could not find points to elongate cable')
#
#    success=False
#    index_end = 0
#    for ii, point in enumerate(points[::-1][1:]):
#        A = left_p_end[-1]
#        B = point
#        AB = B-A
#        vec = way(B-A)
#        if AB.norm() > min_dist:
#            working_p_end.append(A+vec*min_dist/2)
#            success = True
#            index_end = ii+1
#            break
#        else:
#            left_p_end.append(B)
#            to_meander.pop(-1)
#
#    if not success:
#        raise ValueError('Could not find points to elongate cable')
#
#    working_p = working_p_start+points[index_start:-index_end]+working_p_end
#    index_insertion = len(left_p_start)
#    left_p = left_p_start+left_p_end[::-1]
#
#    return working_p, left_p, index_insertion
#
#def right_left(points):
#    vecs = []
#    A = points[0]
#    for B in points[1:]:
#        vecs.append(way(B-A))
#        A=B
##    print(points)
##    print(vecs)
#    vecA = vecs[0]
#    r_l = [0]
#    for vecB in vecs[1:]:
#        r_l.append(vecA.cross(vecB))
#        vecA=vecB
#    r_l.append(0)
#    return r_l
#
#def add_points(points, rl, min_dist):
#    min_dist = min_dist*1.1
#    n_points = len(points)
#    print(n_points)
#    A = points[0]
#    new_points =[]
#    if n_points==2:
#        new_points.append(A)
#        B=points[-1]
#        vec = way(B-A)
#        AB = (B-A).norm()
#        n_add = int(AB/min_dist)
#        if rl[0]*rl[1]==1 and n_add%2==0:
#            n_add=-1
#        if rl[0]*rl[1]==-1 and n_add%2==1:
#            n_add-=1
#        dist = AB/n_add
#        new_points+=[A+vec*dist/2]
#        new_points+=[A+vec*dist*(jj+1+1/2) for jj in range(n_add-1)]
#        rl=None
#        indices_corners=None
#    else:
#        indices_corners=[]
#        rl = right_left(points)
##        print(rl)
#        for ii, B in enumerate(points[1:]):
#            new_points.append(A)
#            vec = way(B-A)
#            AB = (B-A).norm()
#            if ii==0 or ii==n_points-2:
#                factor = 0.5
#            else:
#                factor = 1
#            n_add = int(AB/min_dist-factor)
#            if not(ii==0 or ii==n_points-2):
#                if rl[ii-1]*rl[ii]==-1 and n_add%2==1:
#                    n_add-=1
#                if rl[ii-1]*rl[ii]==1 and n_add%2==0:
#                    n_add-=1
#            print(n_add)
#            dist = AB/(n_add+factor)
#            if n_add>=1:
#                if ii==0:
#                    new_points+=[A+vec*dist/2]
#                    new_points+=[A+vec*dist*(jj+1+1/2) for jj in range(n_add-1)]
#                elif ii!=n_points-2:
#                    new_points+=[A+vec*dist*(jj+1) for jj in range(n_add)]
#                else:
#                    new_points+=[A+vec*dist*(jj+1) for jj in range(n_add)]
#            indices_corners.append(len(new_points))
#            A=B
#        indices_corners= indices_corners[:-1]
#    new_points.append(points[-1])
#    return new_points, indices_corners, dist
#
#def displace(points, rl, min_dist, displacement=0):
#    if displacement<min_dist:
#        displacement = min_dist*1.1
#    points, indices_corners, dist = add_points(points, rl, min_dist)
#    new_points = [points[0]]
#    parity = 1
#    if indices_corners is not None:
#        for ii, B in enumerate(points[1:-1]):
#            A = points[ii]
#            AB= B-A
#            vec = way(val(AB))
#            if ii==0:
#                parity = (-2*((indices_corners[0]-(ii+1))%2)+1)*(-rl[0])
#            else:
#                parity = -parity
#            print(parity)
#            if ii+1 not in indices_corners:
#                #regular point
#                new_points[ii+1] = points[ii+1]+vec.orth()*parity*min_dist
#            else:
#                new_points[ii+1] = points[ii+1]+(vec.orth()*parity+vec).unit()*min_dist
#    else:
#        if rl[0]!=0:
#            parity = -rl[0]
#            print(f'parity{parity}')
#        else:
#            parity = (2*(len(points)%2)-1) * (-rl[1])
#            print(f'parity{parity}')
#        for ii, B in enumerate(points[1:-1]):
#            A=points[ii]
#            AB=B-A
#            vec=way(val(AB))
#            new_points.append(points[ii+1]+vec.orth()*parity*displacement-vec*dist/2)
#            new_points.append(points[ii+1]+vec.orth()*parity*displacement+vec*dist/2)
#            parity = -parity
#        new_points.append(points[-1])
#
#
#    return new_points
#
#def meander(points, min_dist, to_meander): # to_meander is list of segments to be meander
#    n_points = len(points)
#    n_to_meander = len(to_meander)
#    if n_points-1>n_to_meander:
#        to_meander = to_meander+[0 for ii in range(n_points-1-n_to_meander)]
#    else:
#        to_meander = to_meander[:n_points-1]
#
#    working_p, left_p, index_insertion = working_points(points, min_dist, to_meander)
##    print(working_p, left_p, index_insertion)
#    rl = right_left(working_p)
#    print(to_meander)
#    working_ps = []
#    for ii, isit in enumerate(to_meander):
#        print(ii, isit)
#        if isit==1:
#            new_working_p = displace(working_p[ii:ii+2], rl[ii:ii+2], min_dist, displacement = min_dist)
#        else:
#            new_working_p = working_p[ii:ii+2]
#        working_ps += new_working_p
#
#    left_p[index_insertion:index_insertion] = working_ps
#    return left_p#, working_p
#
#
#fillet=0.07
#
#min_dist = 2*fillet
#
##print(final_choice)
#meander = meander(final_choice, min_dist, [0,1,0,5,2,3])
##print(meanders)
#
#
#fig, ax = plt.subplots(figsize = (6,12))
##for ii, meander in enumerate(meanders):
#ax.plot(np.array(start).T[0],np.array(start).T[1], color='red')
#ax.plot(np.array(end).T[0],np.array(end).T[1], color='red')
#ax.plot(np.array(meander).T[0],np.array(meander).T[1], 'o-')
#move_figure(fig, -1000, 10)