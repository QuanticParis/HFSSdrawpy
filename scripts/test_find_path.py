# -*- coding: utf-8 -*-
"""
Created on Thu Apr 12 23:30:42 2018

@author: checkhov
"""

from scripts.designer import Vector
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
plt.close('all')


def move_figure(f, x, y):
    """Move figure's upper left corner to pixel (x, y)"""
    backend = matplotlib.get_backend()
    if backend == 'TkAgg':
        f.canvas.manager.window.wm_geometry("+%d+%d" % (x, y))
    elif backend == 'WXAgg':
        f.canvas.manager.window.SetPosition((x, y))
    else:
        # This works for QT and GTK
        # You can also use window.setGeometry
        f.canvas.manager.window.move(x, y)
    plt.show()


def way(vec):
    if vec[1] != 0:
        if abs(vec[0]/vec[1])<1e-2:
            if vec[1]>0:
                return Vector(0,1)
            elif vec[1]<0:
                return Vector(0,-1)
    if vec[0] != 0 :
        if abs(vec[1]/vec[0])<1e-2:
            if vec[0]>0:
                return Vector(1,0)
            elif vec[0]<0:
                return Vector(-1,0)

def equal_float(float1, float2):
    if float1!=0:
        rel_diff = abs((float1-float2)/float1)
        if rel_diff<1e-5:
            return True
        else:
            return False
    elif float2!=0:
        rel_diff = abs((float1-float2)/float2)
        if rel_diff<1e-5:
            return True
        else:
            return False
    else:
        return True


iIn_pos = Vector([0,0])
iIn_ori = Vector([0, 1])
iOut_pos = Vector([0, 1])
iOut_ori = Vector([0, -1])
start = [iIn_pos+ iIn_ori.orth()*0.05, iIn_pos- iIn_ori.orth()*0.05]
end = [iOut_pos+ iOut_ori.orth()*0.05, iOut_pos- iOut_ori.orth()*0.05]


class Test(object):

    def val(self, nb):
        return nb

    def find_path(self, fillet, is_meander, to_meander, meander_length):

    #        iIn_pos = self.pos
    #        iIn_ori = self.ori
    #        iOut_pos = self.posOut
    #        iOut_ori = self.oriOut

    #        print(val(iIn_pos))
    #        print(val(iIn_ori))
    #        print(val(iOut_pos))
    #        print(val(iOut_ori))


            point1 = iIn_pos+iIn_ori*1.1*fillet
            point2 = iOut_pos+iOut_ori*1.1*fillet


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
                        points_choices.append([iIn_pos, *c_in, *c_out[:-1][::-1], iOut_pos])
            else:
                choice_in = next_point(point1, point2, iIn_ori)
                for c_in in choice_in:
                    points_choices.append([iIn_pos, *c_in, iOut_pos])


            def cost_f(x):
                if x==1:
                    return 0
                elif x==0:
                    return 1
                else:
                    return 100

            def check(points):
                print(points)
                length = 0
                prev_point = points[0]
                _points = [points[0]]
                vecs = []
                for point in points[1:]:
                    if not equal_float(self.val(point)[0], self.val(prev_point)[0]) or not equal_float(self.val(point)[1], self.val(prev_point)[1]):
                        vec = self.val(point-prev_point)
                        print(point, prev_point)
                        print(vec)
                        length += self.val(vec).norm()
                        vecs.append(way(vec))
                        prev_point = point
                        _points.append(point)
                cost = 0

                points = _points.copy()
                new_points = [points[0]]
                print(vecs)
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

            def displace(points, rl, min_dist, displacement=0, n_meander=-1):
                if self.val(displacement)<self.val(min_dist)*1.1:
                    displacement = min_dist*1.1
                print(n_meander)
                points, indices_corners, dist, ignore = add_points(points, rl, min_dist, n_meander=n_meander)
                print(points)
                print(dist)
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
                        print('rl '+str(-rl[1]*(rl[1]+1)+1))
                    print(ignore)
                    if ignore:
                        n_ignore=2
                        new_points.append(points[1])
                    else:
                        n_ignore=1
                    for ii, B in enumerate(points[n_ignore:-n_ignore]):
                        A=points[ii]
                        AB=B-A
                        vec=way(self.val(AB))
                        new_points.append(points[ii+n_ignore]+vec.orth()*parity*displacement-vec*dist/2)
                        new_points.append(points[ii+n_ignore]+vec.orth()*parity*displacement+vec*dist/2)
                        parity = -parity
                        print(points[ii+n_ignore])
                    if ignore:
                        new_points.append(points[-2])
                    new_points.append(points[-1])


                return new_points

            def meander(points, min_dist, to_meander, meander_length): # to_meander is list of segments to be meander
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
                        print(isit)
                        if isit!=0:
                            print('calling displace')
                            new_working_p = displace(working_p[ii:ii+2], rl[ii:ii+2], min_dist, displacement = meander_length, n_meander=isit) # n_meander=-1 -> auto
                        else:
                            new_working_p = working_p[ii:ii+2]
                        working_ps += new_working_p
#                        print(working_ps)

                    left_p[index_insertion:index_insertion] = working_ps
                return  left_p#left_p#,

            if is_meander:
                min_dist = 2*fillet
                final_choice = meander(final_choice, min_dist, to_meander, meander_length)



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


            _, final_choice, _ = check(final_choice)

            to_bond_points = add_fillet_points(final_choice, fillet)
            for ii, point in enumerate(to_bond_points[::2]):
    #            self.draw('bef_test', [to_bond_points[2*ii], to_bond_points[2*ii+1]], closed=False)
#                self.to_bond.append([to_bond_points[2*ii], to_bond_points[2*ii+1]])
                pass
#            to_bond_points= final_choice
            return final_choice, to_bond_points


dummy = Test()
meander, to_bond_points = dummy.find_path(fillet=0.1, is_meander=True, to_meander=[1,0,0,1,1,0,1], meander_length=0.3)


fig, ax = plt.subplots(figsize = (6,12))
#for ii, meander in enumerate(meanders):
ax.plot(np.array(start).T[0],np.array(start).T[1], color='red')
ax.plot(np.array(end).T[0],np.array(end).T[1], color='red')
ax.plot(np.array(meander).T[0],np.array(meander).T[1], 'o-')

for ii, point in enumerate(to_bond_points[:-2][::2]):
    ax.plot([to_bond_points[2*ii+1][0], to_bond_points[2*ii+2][0]],[to_bond_points[2*ii+1][1], to_bond_points[2*ii+2][1]], color='g')
#ax.axis('equal')
ax.set_xlim((-0.5,1.5))
ax.set_ylim((-0.5,1.5))
move_figure(fig, -1000, 10)




#
#points_fromIn = [iIn_pos]
#points_fromOut = [iOut_pos]
#point1 = iIn_pos+iIn_ori*0.1
#point2 = iOut_pos+iOut_ori*0.1
#
#
#
#
#
#
#def next_point(point1, point2, vec):
#    choice1 = point1+(point2-point1).dot(vec)*vec
#    choice2 = point1+(point2-point1).dot(vec.orth())*vec.orth()
#    return [[point1, choice1, point2], [point1, choice2, point2]]
#
#
#points_choices = []
#if iIn_ori.dot(iOut_ori)==-1:
#    middle_point = (point1 + point2)/2
#
#    choice_in = next_point(point1, middle_point, iIn_ori) #bon sens
#    choice_out = next_point(point2, middle_point, iOut_ori) #à inverser
##    print(choice_in)
##    print(choice_out)
#    for c_in in choice_in:
#        for c_out in choice_out:
#            points_choices.append([iIn_pos, *c_in, *c_out[:-1][::-1], iOut_pos])
#else:
#    choice_in = next_point(point1, point2, iIn_ori)
#    for c_in in choice_in:
#        points_choices.append([iIn_pos, *c_in, iOut_pos])
#
#
#def cost_f(x):
#    if abs(x-1)<1e-2:
#        return 0
#    elif abs(x)<1e-2:
#        return 1
#    else:
#        return 100
#
#def check(points):
#    points_array = np.array(points)
#    vecs = points_array[:-1]-points_array[1:]
#
#    cost = 0
#    prev_vec = Vector(vecs[0]).unit()
#    new_points=[]
#    new_points.append(points[0])
#    for ii, vec in enumerate(vecs[1:]):
#        curr_vec = Vector(vec).unit()
#        if curr_vec.dot(prev_vec)==0:
#            new_points.append(points[ii+1])
#        added_cost = cost_f(prev_vec.dot(curr_vec))
#        cost += added_cost
#        prev_vec = curr_vec
#    new_points.append(points[-1])
#    return cost, new_points
#
#final_choice= None
#cost=np.inf
#for choice in points_choices:
#    new_cost, new_points = check(choice)
#    if new_cost<cost:
#        final_choice = new_points
#        cost = new_cost
#
#
#
############
##def relative_append(A, vector_list):
##    new_vector_list = [A+vector_list[0]]
##    for vector in vector_list[1:]:
##        new_vector_list.append(new_vector_list[-1]+vector)
##    return new_vector_list
##
##def build_meander_list(fillet, vec, n_vec, n):
##    if n==0:
##        return [vec*fillet].extend(build_meander_list).extend([vec*fillet])
##    if n>0:
##        if n%2=0:
##            return build_meander_list(fillet, vec, n_vec, n-1)
##        else:
##            return build_meander_list(fillet, vec, -n_vec, n-1)
##
##fillet=0.1
##min_width = 4*fillet
##min_offset = 4*fillet*(np.pi/2-1) #off
##def meander(points, index, number, kind): #index = index segment on which you would like to make meanders
##    #kind will be start, middle, end
##    n_points = len(points)
##    A = points[index]
##    B = points[index+1]
##    AB = B-A
##    vec = way(AB)
##    if index==0:
##        prev_vec=None
##        next_vec = way(points[index+2]-B) # we assume we have at least one corner
##        points_to_append = relative_append(A, [vec*fillet, -next_vec*2*fillet, vec*2*fillet, next_vec*2*fillet, vec*fillet])
##        points[index+1:index+1] = points_to_append
##    elif index==n_points-1:
##        prev_vec = way(A-points[index-1])
##        next_vec=None
##    else:
##        prev_vec = way(A-points[index-1])
##        next_vec = way(points[index+2]-B)
############
#
##def working_points(points, min_dist, to_meander):
##    min_dist = min_dist*1.1
##    working_p_start = []
##    working_p_end = []
##    left_p_start=[points[0]]
##    left_p_end=[points[-1]]
##    success=False
##    index_start = 0
##    for ii, point in enumerate(points[1:]):
##        A = left_p_start[-1]
##        B = point
##        AB = B-A
##        vec = way(B-A)
##        if AB.norm() > min_dist:
##            working_p_start.append(A+vec*min_dist/2)
##            success = True
##            index_start = ii+1
##            break
##        else:
##            left_p_start.append(B)
##            to_meander.pop(0)
##
##
##    if not success:
##        raise ValueError('Could not find points to elongate cable')
##
##    success=False
##    index_end = 0
##    for ii, point in enumerate(points[::-1][1:]):
##        A = left_p_end[-1]
##        B = point
##        AB = B-A
##        vec = way(B-A)
##        if AB.norm() > min_dist:
##            working_p_end.append(A+vec*min_dist/2)
##            success = True
##            index_end = ii+1
##            break
##        else:
##            left_p_end.append(B)
##            to_meander.pop(-1)
##
##    if not success:
##        raise ValueError('Could not find points to elongate cable')
##
##    working_p = working_p_start+points[index_start:-index_end]+working_p_end
##    index_insertion = len(left_p_start)
##    left_p = left_p_start+left_p_end[::-1]
##
##    return working_p, left_p, index_insertion
##
##def right_left(points):
##    vecs = []
##    A = points[0]
##    for B in points[1:]:
##        vecs.append(way(B-A))
##        A=B
###    print(points)
###    print(vecs)
##    vecA = vecs[0]
##    r_l = [0]
##    for vecB in vecs[1:]:
##        r_l.append(vecA.cross(vecB))
##        vecA=vecB
##    r_l.append(0)
##    return r_l
##
##def add_points(points, rl, min_dist):
##    min_dist = min_dist*1.1
##    n_points = len(points)
##    print(n_points)
##    A = points[0]
##    new_points =[]
##    if n_points==2:
##        new_points.append(A)
##        B=points[-1]
##        vec = way(B-A)
##        AB = (B-A).norm()
##        n_add = int(AB/min_dist)
##        if rl[0]*rl[1]==1 and n_add%2==0:
##            n_add=-1
##        if rl[0]*rl[1]==-1 and n_add%2==1:
##            n_add-=1
##        dist = AB/n_add
##        new_points+=[A+vec*dist/2]
##        new_points+=[A+vec*dist*(jj+1+1/2) for jj in range(n_add-1)]
##        rl=None
##        indices_corners=None
##    else:
##        indices_corners=[]
##        rl = right_left(points)
###        print(rl)
##        for ii, B in enumerate(points[1:]):
##            new_points.append(A)
##            vec = way(B-A)
##            AB = (B-A).norm()
##            if ii==0 or ii==n_points-2:
##                factor = 0.5
##            else:
##                factor = 1
##            n_add = int(AB/min_dist-factor)
##            if not(ii==0 or ii==n_points-2):
##                if rl[ii-1]*rl[ii]==-1 and n_add%2==1:
##                    n_add-=1
##                if rl[ii-1]*rl[ii]==1 and n_add%2==0:
##                    n_add-=1
##            print(n_add)
##            dist = AB/(n_add+factor)
##            if n_add>=1:
##                if ii==0:
##                    new_points+=[A+vec*dist/2]
##                    new_points+=[A+vec*dist*(jj+1+1/2) for jj in range(n_add-1)]
##                elif ii!=n_points-2:
##                    new_points+=[A+vec*dist*(jj+1) for jj in range(n_add)]
##                else:
##                    new_points+=[A+vec*dist*(jj+1) for jj in range(n_add)]
##            indices_corners.append(len(new_points))
##            A=B
##        indices_corners= indices_corners[:-1]
##    new_points.append(points[-1])
##    return new_points, indices_corners, dist
##
##def displace(points, rl, min_dist, displacement=0):
##    if displacement<min_dist:
##        displacement = min_dist*1.1
##    points, indices_corners, dist = add_points(points, rl, min_dist)
##    new_points = [points[0]]
##    parity = 1
##    if indices_corners is not None:
##        for ii, B in enumerate(points[1:-1]):
##            A = points[ii]
##            AB= B-A
##            vec = way(val(AB))
##            if ii==0:
##                parity = (-2*((indices_corners[0]-(ii+1))%2)+1)*(-rl[0])
##            else:
##                parity = -parity
##            print(parity)
##            if ii+1 not in indices_corners:
##                #regular point
##                new_points[ii+1] = points[ii+1]+vec.orth()*parity*min_dist
##            else:
##                new_points[ii+1] = points[ii+1]+(vec.orth()*parity+vec).unit()*min_dist
##    else:
##        if rl[0]!=0:
##            parity = -rl[0]
##            print(f'parity{parity}')
##        else:
##            parity = (2*(len(points)%2)-1) * (-rl[1])
##            print(f'parity{parity}')
##        for ii, B in enumerate(points[1:-1]):
##            A=points[ii]
##            AB=B-A
##            vec=way(val(AB))
##            new_points.append(points[ii+1]+vec.orth()*parity*displacement-vec*dist/2)
##            new_points.append(points[ii+1]+vec.orth()*parity*displacement+vec*dist/2)
##            parity = -parity
##        new_points.append(points[-1])
##
##
##    return new_points
##
##def meander(points, min_dist, to_meander): # to_meander is list of segments to be meander
##    n_points = len(points)
##    n_to_meander = len(to_meander)
##    if n_points-1>n_to_meander:
##        to_meander = to_meander+[0 for ii in range(n_points-1-n_to_meander)]
##    else:
##        to_meander = to_meander[:n_points-1]
##
##    working_p, left_p, index_insertion = working_points(points, min_dist, to_meander)
###    print(working_p, left_p, index_insertion)
##    rl = right_left(working_p)
##    print(to_meander)
##    working_ps = []
##    for ii, isit in enumerate(to_meander):
##        print(ii, isit)
##        if isit==1:
##            new_working_p = displace(working_p[ii:ii+2], rl[ii:ii+2], min_dist, displacement = min_dist)
##        else:
##            new_working_p = working_p[ii:ii+2]
##        working_ps += new_working_p
##
##    left_p[index_insertion:index_insertion] = working_ps
##    return left_p#, working_p
##
##
##fillet=0.07
##
##min_dist = 2*fillet
##
###print(final_choice)
##meander = meander(final_choice, min_dist, [0,1,0,5,2,3])
###print(meanders)
#
#
