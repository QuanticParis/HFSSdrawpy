import numpy as np

from ..utils import Vector, val, way, way_approx, equal_float, parse_entry
from .test_finding import find_path


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
        
        # TODO check that length > offset + 2*fillet, however I don't know how 
        # to retrieve fillet value here. 
        
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
        # print(self.length*1e3)
        # print(val(self.compute_length_analytic())*1e3)
        
        # print(self.path_ori)
        # print(self.lr)
        # print(self.thetas)
        
        # compute where can be added meanders
        # compute where we can offset
        
        return self.length
    
    
    def meander(self, meanders, target_length=None, editable_in_hfss=False):
        fillet = 1.1*self.fillet
        fillet_val = val(fillet)
                
        if target_length is not None:
            if val(self.length) > val(target_length):
                print('Target length is smaller than bare cable length')
        
        path_copy = self.path.copy()
        adjustable_path = [[None]]*len(path_copy)
        adjustable_offset = [[None]]*len(path_copy)
        meanders.sort(key=meander_index, reverse=True)
        n_adj = 0 # number of nobs to tune the length
        
        for meander in meanders:
            # TODO order meanders in wrt to their index
            if meander.length is None:
                meander_length = 2*self.fillet
            else:
                meander_length = meander.length
                
            ii = meander.index
            if ii < len(self.path_ori):
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
                    if meander.length is not None:
                        points = [point + prev_path_ori*(meander_length*(-1)**(jj//2)+meander.offset) for jj, point in enumerate(_points)]
                    else:
                        points = [point + prev_path_ori*(meander_length*(-1)**(jj//2)+0*meander.offset) for jj, point in enumerate(_points)]
                    
                    if meander.length is None:
                        adjustable_points = [[prev_path_ori*(-1)**(jj//2)] for jj, point in enumerate(_points)]
                        offset_points = [[prev_path_ori*meander.offset] for jj, point in enumerate(_points)]
                        adj_adapt0 = [[None]]*len(adapt0)
                        adj_adapt1 = [[None]]*len(adapt1)
                        n_adj += len(adjustable_points)
                    path_copy = path_copy[:ii] + adapt0 + points + adapt1 + path_copy[ii+2:]
                    if meander.length is None:
                        adjusted =  adj_adapt0 + adjustable_points + adj_adapt1
                        adjusted_offset = adj_adapt0 + offset_points + adj_adapt1
                        adjustable_path = adjustable_path[:ii] + [adjustable_path[ii]+adjusted[0]]+ adjusted[1:-1] + [adjustable_path[ii+1]+adjusted[-1]] + adjustable_path[ii+2:]
                        adjustable_offset = adjustable_offset[:ii] + [adjustable_offset[ii]+adjusted_offset[0]]+ adjusted_offset[1:-1] + [adjustable_offset[ii+1]+adjusted_offset[-1]] + adjustable_offset[ii+2:]
                else:
                    print('Could not place a meander in segment %d.'%ii)
            else:
                print('Ignored meander %d (only %d segments)'%(ii, len(self.path_ori)))
        self.path = path_copy
        
        self.analyse_path()
        
        if target_length is not None:
            if editable_in_hfss: # only work for 90° angles
                if not (sum([abs(lr) for lr in self.lr])) == (len(self.lr)-2):
                    raise ValueError('Slanted cable cannot be adjusted a posteriori in HFSS')
                actual_length = self.compute_length_analytic()
            else:
                actual_length = self.length
            
            if val(actual_length) > val(target_length):
                print('Target length is smaller than bare cable length.')
            else:
                if n_adj==0:
                    print('No nob was given to tune the cable length.')
                else:
                    shift = (target_length-actual_length)/n_adj
                    for jj, (_adj, _off) in enumerate(zip(adjustable_path, adjustable_offset)):
                        for adj, off in zip(_adj, _off):
                            if adj is not None:
                                path_copy[jj] = path_copy[jj] + adj*shift
                            if off is not None:
                                path_copy[jj] = path_copy[jj] + off
                    self.path = path_copy  
                
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
    
    def compute_length_analytic(self):
        value = 0
        for ii in range(len(self.path) - 1):
            link = (self.path[ii + 1] - self.path[ii]).dot(self.path_ori[ii])
            value += link
            if ii>0:
                value -= (2-np.pi/2)*self.fillet
        return value
