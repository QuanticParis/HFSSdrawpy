# -*- coding: utf-8 -*-
"""
Created on Sun Apr  3 20:45:34 2022

@author: dione
"""

import matplotlib.pyplot as plt
import numpy as np

from HFSSdrawpy.utils import Vector, val, way, equal_float

plt.close('all')

def sign(x):
    if x>=0:
        return 1
    else:
        return -1

def outside_margin(val, fillet_val):
    return (val<=-2*fillet_val or val>=2*fillet_val)

def plot(pos, ori, color='C0'):
    ax.plot([pos[0]], [pos[1]], 'o', color=color)
    point = pos + fillet*ori
    ax.plot([pos[0], point[0]], [pos[1], point[1]], color=color)
    
def plot_path(path, color='C0'):
    path = np.array(path).T
    ax.plot(path[0], path[1], color=color)
    
def solve(in_pos, in_ori, out_pos, out_ori):
    ax, ay = in_pos[0], in_pos[1]
    ux, uy = in_ori[0], in_ori[1]
    bx, by = out_pos[0], out_pos[1]
    vx, vy = out_ori[0], out_ori[1]
    
    # ax + l ux = bx + m vx
    # ay + l uy = by + m vy
    # ax*vy - bx*vy = m vx vy - l ux vy
    # -ay*vx + by*vx = -m vy vx + l uy vx
    if not equal_float(val(uy*vx-ux*vy), 0):
        l = (ax*vy - bx*vy -ay*vx + by*vx)/(uy*vx-ux*vy)
        m = (ax*uy - bx*uy -ay*ux + by*ux)/(vx*uy-vy*ux)
        return l, m, val(l), val(m)
    else:
        return None, None, None, None
    
def find_path(in_pos, in_ori, out_pos, out_ori, fillet=0.1, 
              mid=0.5):
    
    fillet = 1.1*fillet 
    fillet_val = val(fillet)

    if mid==0:
        mid = 1e-5
    if mid==1:
        mid = 1-1e-5
    
    in_pos = Vector(in_pos)
    in_ori = Vector(in_ori)
    out_pos = Vector(out_pos)
    out_ori = Vector(out_ori)
    
    path_in = [in_pos]
    path_in_ori = [in_ori]
    path_out = [out_pos]
    path_out_ori = [out_ori]
    
    # slanted case:
    align = val(in_ori.dot(out_ori))
    along = (out_pos-in_pos).dot(in_ori)
    offset = val((out_pos-in_pos).dot(in_ori.orth()))
    if equal_float(align, 1.0) and equal_float(abs(offset), 0) and val(along)>0:
        success=True
        path = path_in + path_out[::-1]
    elif equal_float(align, 1.0) and abs(offset)<=2*fillet_val and val(along)>0:
        print('SLANTED')
        path_in.append(in_pos+along/3*in_ori)
        path_out.append(out_pos-along/3*out_ori)
        success = True
        path = path_in + path_out[::-1]
    else:
        # solve natural intersect
        l, m, l_val, m_val = solve(path_in[-1], path_in_ori[-1], 
              path_out[-1], path_out_ori[-1])    
        if l is not None:
            if l_val>=fillet_val and m_val<=-fillet_val:
                # valid path
                success = True
                path_in.append(path_in[-1]+l*path_in_ori[-1])
            else:
                success = False
        else:
            success = False
        
        if not success:
            # test in ext
            proj = (path_out[-1]-path_in[-1]).dot(path_in_ori[-1])
            if val(proj)>=2*fillet_val:
                ext_in = mid*(proj - 2*fillet) + fillet
            else:
                ext_in = fillet
            
            path_in.append(path_in[-1]+ext_in*path_in_ori[-1])
            path_in_ori.append(path_in_ori[-1].orth())
            # plot(path_in[-1], path_in_ori[-1], color='C2')
            
            l, m, l_val, m_val = solve(path_in[-1], path_in_ori[-1], 
                         path_out[-1], path_out_ori[-1])  
            
            if l is not None:
                if m_val<=-fillet_val and outside_margin(l_val, fillet_val): # U turn = 2 fillets
                    # valid path
                    success = True
                    path_in.append(path_in[-1]+l*path_in_ori[-1])
                else:
                    success = False
            else:
                success = False
                
            if not success:
                # test out ext
                # here
                path_in.pop()
                path_in_ori.pop()
                # here
                
                proj = (path_out[-1]-path_in[-1]).dot(path_out_ori[-1])
                if val(proj)>3*fillet_val:
                    ext_out = (1-mid)*(proj - 3*fillet) + fillet
                else:
                    ext_out = fillet
                    
                path_out.append(path_out[-1]-ext_out*path_out_ori[-1])
                path_out_ori.append(path_out_ori[-1].orth())
                # plot(path_out[-1], path_out_ori[-1], color='C3')
                
                l, m, l_val, m_val = solve(path_in[-1], path_in_ori[-1], 
                             path_out[-1], path_out_ori[-1])  
                
                if l is not None: 
                    if l_val>=fillet_val and outside_margin(m_val, fillet_val):
                        success = True
                        path_out.append(path_out[-1]+m*path_out_ori[-1])
                    else:
                        success = False
                else:
                    success = False
                    
                if not success:
                    # test in out ext
                    
                    # change middle since now turns on out end
                    proj = (path_out[-1]-path_in[-1]).dot(path_in_ori[-1])
                    if val(proj)>=3*fillet_val:
                        ext_in = mid*(proj - 3*fillet) + fillet
                    else:
                        ext_in = fillet
                        
                    path_in.append(path_in[-1]+ext_in*path_in_ori[-1])
                    path_in_ori.append(path_in_ori[-1].orth())
                    # already plotted
                    
                    # out already present
                    
                    l, m, l_val, m_val = solve(path_in[-1], path_in_ori[-1], 
                                 path_out[-1], path_out_ori[-1]) 
                    if l is not None:
                        outside_l = outside_margin(l_val, fillet_val) 
                        outside_m = outside_margin(m_val, fillet_val)
                        if outside_l and outside_m:
                            success = True
                            path_in.append(path_in[-1]+l*path_in_ori[-1])
                        elif outside_l and not outside_m:
                            # modify in point 
                            path_in.pop()
                            path_out.append(path_out[-1]+2*fillet*path_in_ori[-2])
                            path_in.append(path_out[-1]-l*path_in_ori[-1])
                            success = True
                            # success = False
                        elif not outside_l and outside_m:
                            # modify out point
                            path_out.pop()
                            path_in.append(path_in[-1]-2*fillet*path_out_ori[-2])
                            path_out.append(path_in[-1]-m*path_out_ori[-1])
                            success = True
                            # sucess = False
                        else:
                            success = False
                    else:
                        success = False
                        
                    if not success:
                        # test in in out ext
                        # test in in probably equivalent to in out
                        
                        # make sure we go towards out
                        ori = path_in_ori.pop()
                        proj = ori.dot(path_out[-2]-path_in[-2])
                        ori = ori*sign(proj)
                        path_in_ori.append(ori)
                        
                        proj = (path_out[-1]-path_in[-1]).dot(path_in_ori[-1])
                        if val(proj)>4*fillet_val:
                            ext_in = mid*(proj - 4*fillet) + 2*fillet
                        else:
                            ext_in = fillet

                        path_in.append(path_in[-1]+ext_in*path_in_ori[-1])
                        path_in_ori.append(path_in_ori[-1].orth())
                        # plot(path_in[-1], path_in_ori[-1], color='C4')
                        
                        l, m, l_val, m_val = solve(path_in[-1], path_in_ori[-1], 
                                     path_out[-1], path_out_ori[-1])  
                        
                        if l is not None:
                            outside_l = outside_margin(l_val, fillet_val) 
                            outside_m = outside_margin(m_val, fillet_val)
                            if outside_l and outside_m:
                                # valid path
                                success = True
                                path_in.append(path_in[-1]+l*path_in_ori[-1])
                            elif outside_l and not outside_m:
                                # modify in point 
                                path_in.pop()
                                path_out.append(path_out[-1]+2*fillet*path_in_ori[-2])
                                path_in.append(path_out[-1]-l*path_in_ori[-1])
                                success = True
                                # success = False
                            elif not outside_l and outside_m:
                                # modify out point
                                # print(ii, jj)
                                path_out.pop()
                                path_in.append(path_in[-1]-2*fillet*path_out_ori[-2])
                                path_out.append(path_in[-1]-m*path_out_ori[-1])
                                success = True
                                # sucess = False
                            else:
                                success = False
                        else:
                            success = False
                            
                        if not success:
                            print('Did not manage to find path')
                            # raise Exception('Did not manage to find path')

        path = path_in + path_out[::-1]
            
    return path



if __name__ == '__main__':
    in_pos = np.array([0, 0])
    in_ori = np.array([1, 0])
    
    fillet = 0.075
    # fillet = 0.15
    fillet = 0.075
    
    fig_tot, ax_tot = plt.subplots(4, 4, figsize=(18, 18))
    
    out_poss = [[1, 1], [-1, 1], [-1, -1], [1, -1]]
    out_poss = [[0.1, 1], [-0.1, 1], [-0.1, -1], [0.1, -1]]
    out_poss = [[1, 0.1], [-1, 0.1], [-1, -0.1], [1, -0.1]]

    
    out_oris = [[1, 0], [0, 1], [-1, 0], [0, -1]]
    for ii, out_pos in enumerate(out_poss):
        for jj, out_ori in enumerate(out_oris):
            out_pos = np.array(out_pos)
            out_ori = np.array(out_ori)
            ax = ax_tot[ii, jj]
            
            plot(in_pos, in_ori, color='C0')
            plot(out_pos, out_ori, color='C1')
    
            path = find_path(in_pos, in_ori, out_pos, out_ori, 
                             fillet=fillet, mid=0.5)
            plot_path(path)
            
            ax.set_xlim(-1.5, 1.5)
            ax.set_ylim(-1.5, 1.5)
            ax.set_aspect('equal')



