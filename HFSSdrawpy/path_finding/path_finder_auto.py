import numpy as np

from ..utils import Vector, val, way, way_approx, equal_float, parse_entry
from .test_finding import find_path

# useful function to find cable path
def next_point(point1, point2, vec):
    choice1 = point1 + vec * (point2 - point1).dot(vec)
    choice2 = point1 + vec.orth() * (point2 - point1).dot(vec.orth())
    return [[point1, choice1, point2], [point1, choice2, point2]]


def cost_f(x):
    if x == 1:
        return 0
    elif x == 0:
        return 1
    else:
        return 100


def right_left(points):
    # tells if given point is turning left or right
    vecs = []
    A = points[0]
    for B in points[1:]:
        vecs.append(way(val(B - A)))
        A = B
    #    print(points)
    #    print(vecs)
    vecA = vecs[0]
    r_l = [0]
    for vecB in vecs[1:]:
        r_l.append(vecA.scalar_cross(vecB))
        vecA = vecB
    r_l.append(0)
    return r_l


def displace(points, rl, min_dist, displacement=0, offset=0, n_meander=-1):
    # points is simply a segment
    # rl indicates what type of corner we have on each end
    if np.abs(val(displacement)) < val(min_dist) * 1.1:
        displacement = min_dist * 1.1

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

    if rl[0] != 0:
        parity = -rl[0]
    else:
        parity = (2 * (len(points) % 2) - 1) * (-rl[1] * (rl[1] + 1) + 1)
    if ignore:
        n_ignore = 2
        new_points.append(points[1])
    else:
        n_ignore = 1
    for ii, B in enumerate(points[n_ignore:-n_ignore]):
        A = points[ii]
        AB = B - A
        vec = way(val(AB))
        if parity == 1:
            new_points.append(
                points[ii + n_ignore]
                + vec.orth() * parity * (displacement + offset)
                - vec * dist / 2
            )
            new_points.append(
                points[ii + n_ignore]
                + vec.orth() * parity * (displacement + offset)
                + vec * dist / 2
            )
        else:
            new_points.append(
                points[ii + n_ignore]
                + vec.orth() * parity * (displacement - offset)
                - vec * dist / 2
            )
            new_points.append(
                points[ii + n_ignore]
                + vec.orth() * parity * (displacement - offset)
                + vec * dist / 2
            )
        parity = -parity
    if ignore:
        new_points.append(points[-2])
    new_points.append(points[-1])

    return new_points, n_add


def add_points(points, rl, min_dist, n_meander):
    min_dist = min_dist * 1.1
    n_points = len(points)
    if n_points != 2:
        print("n_points!=2")
    A = points[0]
    new_points = [A]
    B = points[-1]
    vec = way(val(B - A))
    AB = (B - A).norm()
    n_add = int(val(AB / min_dist))
    if rl[0] * rl[1] == 1 and n_add % 2 == 0:
        n_add -= 1
    if rl[0] * rl[1] == -1 and n_add % 2 == 1:
        n_add -= 1
    if n_meander == -1 or n_meander >= n_add:
        dist = AB / n_add
        ignore = False
    elif n_meander < n_add:
        n_add = n_meander
        centerAB = (A + B) / 2
        addedA = centerAB - vec * n_add / 2 * min_dist
        addedB = centerAB + vec * n_add / 2 * min_dist
        dist = min_dist
        A = addedA
        new_points.append(addedA)
        ignore = True
    # new_points+=[A+vec*dist/2]
    new_points += [A + vec * dist * (jj + 1 / 2) for jj in range(n_add)]
    rl = None
    indices_corners = None

    # print('n_add', n_add)

    if ignore:
        new_points.append(addedB)
    new_points.append(B)
    return new_points, indices_corners, dist, ignore, n_add


def ori_from_path(path):
    pt0 = path[0]
    path_ori = []
    for pt in path[1:]:
        path_ori.append(way(val(pt-pt0)))
        pt0 = pt
    return path_ori

def left_right(path_ori):
    ori0 = path_ori[0]
    lr = [0] # starting point is not a 90° corner
    for ori in path_ori[1:]:
        lr.append(ori[0]*ori0[1]-ori[1]*ori0[0])
        ori0 = ori
    lr.append(0) # ending point is not a 90° corner
    return lr

def wonder_stretch(link, fillet, n, added):
    if link<(n*2+added)*fillet: # do not fit
        if link>(n*2+added-2)*fillet: # does it fit once stretched
            return True, n
        else:
            while link<(n*2+added-2)*fillet:
                n = n-1
            return True, n
    else:
        return False, n

def meander_index(meander):
    return getattr(meander, 'index')

class Meander(object):
    def __init__(self, index, length=None, n=1, offset='0mm', stretch=True, pos=0):
        length, offset = parse_entry(length, offset)
        
        self.index = index # which cable segment should host the meander
        self.n = n # number of meander on a cable section
        self.length = length # meander_length
        self.offset = offset # meander_offset
        self.pos = pos # should be 0, 0.5, 1 depending on where the meander should
                        # be along the segment
        self.stretch = stretch # should the meander run along all the segment
        
class PathAuto(object):
    def __init__(self, name, port_in, port_out, fillet, path=None,
                 mid=0.5):
        self.name = name
        self.port_in = port_in
        self.port_out = port_out
        self.fillet = fillet
        self.thetas = None
        self.lr = None
        self.path_ori = None
        self.length = None
        
        if path is None:
            self.path = find_path(port_in.pos, port_in.ori, 
                                     port_out.pos, port_out.ori,
                                     fillet=fillet, 
                                     mid=mid)
        else:
            self.path = path
                
        # nob
        # ensemble of points, for each the direction for shifting 
        # and the ratio between shifting and elongation
        # the type of shift could be provided (shift or meander for instance)
    
    def analyse_path(self):
        self.thetas = None
        self.path_ori = ori_from_path(self.path)
        self.lr = left_right(self.path_ori)
        self.length = self.compute_length()
        
        # print(self.path_ori)
        # print(self.lr)
        # print(self.thetas)
        
        # compute where can be added meanders
        # compute where we can offset
        
        return self.length
    
    
    def meander(self, meanders, target_length=None):
        fillet = 1.1*self.fillet
        fillet_val = val(fillet)
        
        if target_length is not None:
            if val(self.length) > val(target_length):
                print('Target length is smaller than bare cable length')
        
        path_copy = self.path.copy()
        meanders.sort(key=meander_index, reverse=True)
        n_adj = 0 # number of nobs to tune the length
        for meander in meanders:
            # TODO order meanders in wrt to their index
            if meander.length is None:
                meander_length = 2*self.fillet
            else:
                meander_length = meander.length
                
            ii = meander.index

            ori = self.path_ori[ii]
            link = self.links[ii]
            if ori == Vector(0, 0):
                ori = way_approx(path_copy[ii+1]-path_copy[ii])
            
            # # check corner integration
            # if self.lr[ii]==0 and meander.pos==0:
            #     print('Cannot integrate meander in a non 90° corner')
            #     print('Meander is positionned in the center')
            #     meander.pos = 0.5
            # if self.lr[ii+1]==0 and meander.pos==1:
            #     print('Cannot integrate meander in a non 90° corner')
            #     print('Meander is positionned in the center')
            #     meander.pos = 0.5
            
            # check room
            if not meander.stretch:
                if abs(self.lr[ii]*self.lr[ii+1])==1:
                    if (meander.pos==0) or (meander.pos==1):
                        meander.stretch , meander.n = wonder_stretch(link, fillet_val, meander.n, 2)
                        if not meander.stretch:
                            if meander.pos==0:
                                adapt0 = []
                                adapt1 = [path_copy[ii]+2*meander.n*fillet*ori, path_copy[ii+1]]
                                pt0 = path_copy[ii]
                                pt1 = adapt1[0]
                                prev_path_ori = self.path_ori[ii-1]
                            else:
                                adapt0 = [path_copy[ii], path_copy[ii+1]-2*meander.n*fillet*ori]
                                adapt1 = []
                                pt0 = adapt0[-1]
                                pt1 = path_copy[ii+1]
                                prev_path_ori = self.path_ori[ii+1]*(-1)**meander.n
                    else:
                        meander.stretch , meander.n = wonder_stretch(link, fillet_val, meander.n, 4)
                        if not meander.stretch:
                            adapt0 = [path_copy[ii], (path_copy[ii+1]+path_copy[ii])/2-meander.n*fillet*ori]
                            adapt1 = [(path_copy[ii+1]+path_copy[ii])/2+meander.n*fillet*ori, path_copy[ii+1]]
                            pt0 = adapt0[-1]
                            pt1 = adapt1[0]
                            prev_path_ori = ori.orth()
                
                if (self.lr[ii]==0) ^ (self.lr[ii+1]==0):
                    if self.lr[ii]==0:
                        if (meander.pos==0) or (meander.pos==0.5):
                            meander.stretch , meander.n = wonder_stretch(link, fillet_val, meander.n, 3)
                            if not meander.stretch:
                                if meander.pos==0:
                                    adapt0 = [path_copy[ii], path_copy[ii]+fillet*ori]
                                    adapt1 = [path_copy[ii]+(2*meander.n+1)*fillet*ori, path_copy[ii+1]]
                                    pt0 = adapt0[-1]
                                    pt1 = adapt1[0]
                                    prev_path_ori = ori.orth()
                                else:
                                    adapt0 = [path_copy[ii], (path_copy[ii+1]+path_copy[ii])/2-(2*meander.n+1)/2*fillet*ori]
                                    adapt1 = [(path_copy[ii+1]+path_copy[ii])/2+(2*meander.n-1)/2*fillet*ori, path_copy[ii+1]]
                                    pt0 = adapt0[-1]
                                    pt1 = adapt1[0]
                                    prev_path_ori = ori.orth()
                        else:
                            meander.stretch , meander.n = wonder_stretch(link, fillet_val, meander.n, 1)
                            if not meander.stretch:
                                adapt0 = [path_copy[ii], path_copy[ii+1]-2*meander.n*fillet*ori]
                                adapt1 = []
                                pt0 = adapt0[-1]
                                pt1 = path_copy[ii+1]
                                prev_path_ori = self.path_ori[ii+1]*(-1)**meander.n
                    else:
                        if (meander.pos==1) or (meander.pos==0.5):
                            meander.stretch , meander.n = wonder_stretch(link, fillet_val, meander.n, 3)
                            if not meander.stretch:
                                if meander.pos==1:
                                    adapt0 = [path_copy[ii], path_copy[ii+1]-(2*meander.n+1)*fillet*ori]
                                    adapt1 = [path_copy[ii+1]-fillet*ori, path_copy[ii+1]]
                                    pt0 = adapt0[-1]
                                    pt1 = adapt1[0]
                                    prev_path_ori = ori.orth()
                                else:
                                    adapt0 = [path_copy[ii], (path_copy[ii+1]+path_copy[ii])/2-(2*meander.n-1)/2*fillet*ori]
                                    adapt1 = [(path_copy[ii+1]+path_copy[ii])/2+(2*meander.n+1)/2*fillet*ori, path_copy[ii+1]]
                                    pt0 = adapt0[-1]
                                    pt1 = adapt1[0]
                                    prev_path_ori = ori.orth()
                        else:
                            meander.stretch , meander.n = wonder_stretch(link, fillet_val, meander.n, 1)
                            if not meander.stretch:
                                adapt0 = []
                                adapt1 = [path_copy[ii]+2*meander.n*fillet*ori, path_copy[ii+1]]
                                pt0 = path_copy[ii]
                                pt1 = path_copy[ii]+2*meander.n*fillet*ori
                                prev_path_ori = self.path_ori[ii-1]
                                
                if (self.lr[ii]==0) and (self.lr[ii+1]==0):
                    meander.stretch , meander.n = wonder_stretch(link, fillet_val, meander.n, 2)
                    if not meander.stretch:
                        if meander.pos==0:
                            adapt0 = [path_copy[ii], path_copy[ii]+fillet*ori]
                            adapt1 = [path_copy[ii]+(2*meander.n+1)*fillet*ori, path_copy[ii+1]]
                            pt0 = adapt0[-1]
                            pt1 = adapt1[0]
                            prev_path_ori = ori.orth()
                        elif meander.pos==0.5:
                            adapt0 = [path_copy[ii], (path_copy[ii+1]+path_copy[ii])/2-meander.n*fillet*ori]
                            adapt1 = [(path_copy[ii+1]+path_copy[ii])/2+meander.n*fillet*ori, path_copy[ii+1]]
                            pt0 = adapt0[-1]
                            pt1 = adapt1[0]
                            prev_path_ori = ori.orth()
                        else:
                            adapt0 = [path_copy[ii], path_copy[ii+1]-(2*meander.n+1)*fillet*ori]
                            adapt1 = [path_copy[ii+1]-fillet*ori, path_copy[ii+1]]
                            pt0 = adapt0[-1]
                            pt1 = adapt1[0]
                            prev_path_ori = ori.orth()     
                            
                            
            if meander.stretch:
                if abs(self.lr[ii]*self.lr[ii+1])==1:
                    # ensure stretch will work
                    if (self.lr[ii]*self.lr[ii+1]==-1) and (meander.n%2==1):
                        meander.n += 1
                    if (self.lr[ii]*self.lr[ii+1]==1) and (meander.n%2==0):
                        meander.n += 1
                    while link<(meander.n*2)*fillet_val: # only case where we want to adjust the meander number 2 by 2
                        meander.n -= 2
                    adapt0 = []
                    adapt1 = []
                    pt0 = path_copy[ii]
                    pt1 = path_copy[ii+1]
                    prev_path_ori = self.path_ori[ii-1]
                    if meander.length is None:
                        if meander.n <= 1:
                            meander_length = 0
                        else:
                            meander_length = self.fillet
                if (self.lr[ii]==0) ^ (self.lr[ii+1]==0):
                    while link<(meander.n*2 + 1)*fillet_val:
                        meander.n -= 1      
                    if self.lr[ii]==0:
                        adapt0 = [path_copy[ii], path_copy[ii] + ori*fillet]
                        adapt1 = []
                        pt0 = adapt0[-1]
                        pt1 = path_copy[ii+1]
                        prev_path_ori = self.path_ori[ii+1]*(-1)**meander.n
                    else:
                        adapt0 = []
                        adapt1 = [path_copy[ii+1]-ori*fillet, path_copy[ii+1]]
                        pt0 = path_copy[ii]
                        pt1 = adapt1[0]
                        prev_path_ori = self.path_ori[ii-1]
                if (self.lr[ii]==0) and (self.lr[ii+1]==0):
                    # TODO, in the slated case, one should compute projection instead of link
                    while link<(meander.n*2 + 2)*fillet_val:
                        meander.n -= 1   
                    adapt0 = [path_copy[ii], path_copy[ii] + ori*fillet]
                    adapt1 = [path_copy[ii+1]-ori*fillet, path_copy[ii+1]]
                    pt0 = adapt0[-1]
                    pt1 = adapt1[0]
                    prev_path_ori = ori.orth()
                    
            if meander.n>0:
                _points = [(pt1-pt0)*(jj/meander.n)+pt0 for jj in range(meander.n+1) for _ in (0, 1)][1:-1]
                points = [point + prev_path_ori*(meander_length*(-1)**(jj//2)+meander.offset) for jj, point in enumerate(_points)]
                if meander.length is None:
                    adjustable_points = [prev_path_ori*(-1)**(jj//2) for jj, point in enumerate(_points)]
                    adjustable_path = [None]*len(path_copy)
                    adj_adapt0 = [None]*len(adapt0)
                    adj_adapt1 = [None]*len(adapt1)
                    n_adj += len(adjustable_points)
                path_copy = path_copy[:ii] + adapt0 + points + adapt1 + path_copy[ii+2:]
                if meander.length is None:
                    adjustable_path = adjustable_path[:ii] + adj_adapt0 + adjustable_points + adj_adapt1 + adjustable_path[ii+2:]
            else:
                print('Could not place a meander')
        self.path = path_copy

        self.analyse_path() # make the adjustment parametrisable i.e. compute length in a symbolic manner

        
        if target_length is not None:
            if val(self.length) > val(target_length):
                print('Target length is smaller than bare cable length')
            else:
                if n_adj==0:
                    print('No nob was given to tune the cable length')
                else:
                    shift = (target_length-self.length)/n_adj
    
                    for jj, adj in enumerate(adjustable_path):
                        if adj is not None:
                            path_copy[jj] = path_copy[jj] + adj*shift
                            
                    self.path = path_copy
            
                                                   
            
            
            
            
                
        # check if enough room to put meanders
        # this room depends on lr and pos
        # create the points
        # reanalyse_path or update self.stuff more smartly
        
        
        pass

    def check_path(self):

        __path = np.diff(path, axis=0)
        __path = np.sum(__path, axis=1)
        if abs(__path[0])<fillet-1e-5:
            print('start', ii, jj)
            print(__path)
        if abs(__path[-1])<fillet-1e-5:
            print('end', ii, jj)
            print(__path)
        if np.any(np.abs(__path[1:-1])<2*fillet-1e-5):
            print('mid', ii, jj)
            print(__path)

    def adjust_path():
        
        adjust_links = []
        if len(lr)>1:
            index = 0
            if len(shifts)<len(lr):
                shifts += [0]*(len(lr)-len(shifts))
            for kk in range(len(lr)-1):
                if lr[kk]==lr[kk+1]:
                    if shifts[index] is None:
                        shift = 0
                        adjust_links.append(kk)
                    else:
                        shift = shifts[index]
                    path[kk+1] = path[kk+1]+shift*path_ori[kk]
                    path[kk+2] = path[kk+2]+shift*path_ori[kk]
                    index += 1
        true_fillet = val(fillet)/1.1
        length_path, thetas = compute_length(path, true_fillet)
        print("%.3f"%(length_path*1000))
        if target_length is not None:
            if val(target_length)<val(length_path):
                raise Exception('Target_length (%.3f mm) is smaller than minimal length (%.3f mm)'%(1000*val(target_length), 1000*val(length_path)))
            elif sum(x is None for x in shifts)==0:
                raise Exception('No adjustable parameter was given, please use None inputs')
            elif len(adjust_links)==0:
                raise Exception('No adjustable link in cable')
            else:
                shift = (target_length-length_path)/2/len(adjust_links)
                for kk in adjust_links:
                    path[kk+1] = path[kk+1]+shift*path_ori[kk]
                    path[kk+2] = path[kk+2]+shift*path_ori[kk]
                
                length_path, thetas = compute_length(path, true_fillet)
                print("%.3f"%(length_path*1000))
                
    def __add__(self, other):
        assert isinstance(other, PathAuto)
        assert self.path[-1] == other.path[0]
        path = self.path[:-1] + other.path[1:]
        return PathAuto(self.name, 
                        self.port_in, 
                        other.port_out,
                        self.fillet,
                        path=path)
    
    def to_bond(self):
        bonding_segments=[]
        for ii, ori in enumerate(self.path_ori):
            if ori[0]!=0 or ori[1]!=0:
                bonding_segments.append([self.path[ii]+ori*self.fillet, self.path[ii+1]-ori*self.fillet])
        return bonding_segments

    # def meander(
    #     self, to_meander, meander_length, meander_offset
    # ):  # to_meander is list of segments to be meander
    #     min_dist = 2 * self.fillet
    #     points = self.path.copy()
    #     n_points = len(points)
    #     n_to_meander = len(to_meander)
    #     if n_points - 1 > n_to_meander:
    #         to_meander = to_meander + [0 for ii in range(n_points - 1 - n_to_meander)]
    #     else:
    #         to_meander = to_meander[: n_points - 1]
    #     working_p, left_p, index_insertion = self.working_points(points, min_dist, to_meander)
    #     # # create adjustable variable for meander_length if it is not the case
    #     # if isinstance(meander_length, VariableString):
    #     #     meander_length_name = str(meander_length)
    #     # else:
    #     #     meander_length_name = self.name+'_meander_length'
    #     #     meander_length = VariableString(meander_length_name, value=1.1*min_dist)
    #     # VariableString.variables[meander_length_name] = 1.1*min_dist

    #     tot_add = 0  # number of added meanders
    #     if len(working_p) != 0:
    #         rl = right_left(working_p)
    #         working_ps = []

    #         for ii, isit in enumerate(to_meander):
    #             if isit != 0:
    #                 new_working_p, n_add = displace(
    #                     working_p[ii : ii + 2],
    #                     rl[ii : ii + 2],
    #                     min_dist,
    #                     displacement=meander_length,
    #                     offset=meander_offset,
    #                     n_meander=isit,
    #                 )  # n_meander=-1 -> auto
    #             else:
    #                 new_working_p = working_p[ii : ii + 2]
    #                 n_add = 0
    #             tot_add += n_add
    #             working_ps += new_working_p
    #         left_p[index_insertion:index_insertion] = working_ps

    #     # VariableString.variables[meander_length_name] = 1.1*self.fillet

    #     self.path = [self.path[0]] + left_p + [self.path[-1]]

    def working_points(self, points, min_dist, to_meander):
        min_dist = min_dist * 1.1

        left_p_start = [points[0]]
        left_p_end = [points[-1]]

        index_start = 0
        for ii, point in enumerate(points[1:]):
            A = left_p_start[-1]
            B = point
            AB = B - A
            vec = way(val(B - A))
            if val(AB).norm() > val(min_dist):
                working_p_start = A + vec * min_dist / 2
                index_start = ii + 1
                break
            else:
                left_p_start.append(B)
                to_meander.pop(0)
        else:
            #print("Warning: Could not find points to elongate cable %s" % self.name)
            left_p = left_p_start + left_p_end[::-1]
            return [], left_p, 0

        index_end = 0
        for ii, point in enumerate(points[::-1][1:]):
            A = left_p_end[-1]
            B = point
            AB = B - A
            vec = way(val(B - A))
            if val(AB).norm() > val(min_dist):
                working_p_end = A + vec * min_dist / 2
                index_end = ii + 1
                break
            else:
                left_p_end.append(B)
                to_meander.pop(-1)
        else:
            print("Warning: Could not find points to elongate cable %s" % self.name)
            left_p = left_p_start + left_p_end[::-1]
            return [], left_p, 0

        working_p = [working_p_start] + points[index_start:-index_end] + [working_p_end]
        index_insertion = len(left_p_start)
        left_p = left_p_start + left_p_end[::-1]

        return working_p, left_p, index_insertion

    def compute_length(self):
        if self.thetas is None:
            # angles have not already been computed
            self.thetas = []
            flag = True
        else:
            # should compute the angles
            flag = False
    
        fillet_val = val(self.fillet)
        # corner = val(self.fillet * (2 - np.pi / 2))
        link0 = 0
        value = 0
        links = []
        for ii in range(len(self.path) - 1):
            # TODO prefer better norm so that no sqrt and square since
            # horizontal/vertical segments can be computed by abs
            link = val((self.path[ii + 1] - self.path[ii]).norm())
            value += link
            links.append(link)
            if ii>0:
                if flag:
                    if self.lr[ii-1]==0:
                        dot_prod = (val(self.path[ii + 1] - self.path[ii])).dot(val(self.path[ii] - self.path[ii-1]))
                        theta = np.arccos(dot_prod/link0/link)/(np.pi/2) # in units of pi/2
                        value -= (2*np.tan(np.pi/2*theta/2)-np.pi/2*theta)*fillet_val
                    elif abs(self.lr[ii-1])==1:
                        theta = 1
                        value -= (2-np.pi/2)*fillet_val
                    else:
                        raise Exception()
                    self.thetas.append(theta) # in units of pi/2
                else:
                    theta = self.thetas[ii-1]
                    value -= (2*np.tan(np.pi/2*theta/2)-np.pi/2*theta)*fillet_val
            link0 = link
        self.length = value
        self.links = links
        return value
