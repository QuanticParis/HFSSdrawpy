'''
TODO

- debug comments
- automate Bonds
- delete "mm", instead, select the units in HFSS
- define in classes
class CircuitElement(object):
    def __init__(self, gates, params):

    def get_gates(self):
        returns list of gates (inputs and outpts)

class Capacitor(CircuitElement):
    def draw(self):
        draws the element

- LitExp: check if variable exists, updates value if exists, create a new one if not.
- Create drawing script which can start from blank file
- Premesh
- add finiding the limits of the drawing

'''

'''
Assumption

- only assign lumpRLC to rectangles
- assume to have list and not tuples for polylines
- TODO, do not do Lj+'1nH' for now, but can do bx+'1mm'
'''
eps = 0.00000001

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

'''
Throughout the code, a single frame is used: the default one
+----------+
|          |
|          |
^          |
y          |
+x>--------+
'''

TOP = [0, 1]
DOWN = [0, -1]
RIGHT = [1, 0]
LEFT = [-1, 0]

POS = 0
ORI = 1
TRACK = 2
GAP = 3

from .hfss import parse_entry
from .hfss import VariableString
import numpy as np


class Vector(list):

    def __init__(self, vec, vec_y=None):
        if vec_y is not None:
            vec = [vec, vec_y]
        super().__init__(parse_entry(vec))

    def check(self, elt):
        return isinstance(elt, list) or isinstance(elt, tuple)

    def check_nb(self, nb):
        return isinstance(nb, float) or isinstance(nb, int) or isinstance(nb, VariableString)

    def __add__(self, other):
        if self.check(other):
            return Vector([self[0]+other[0], self[1]+other[1]])
        else:
            raise TypeError('Could not perform add operation')

    def __radd__(self, other):
        return self + other

    def __sub__(self, other):
        if self.check(other):
            return Vector([self[0]-other[0], self[1]-other[1]])
        else:
            raise TypeError('Could not perform sub operation')

    def __neg__(self):
        return Vector([-self[0], -self[1]])

    def __rsub__(self, other):
        return -self + other

    def __mul__(self, other):
        if self.check(other):
            return Vector([self[0]*other[0], self[1]*other[1]])
        elif self.check_nb(other):
            return Vector([self[0]*other, self[1]*other])
        else:
            raise TypeError('Could not perform mul operation')

    def __rmul__(self, other):
        return self * other

    def __truediv__(self, other):
        if self.check(other):
            return Vector([self[0]/other[0], self[1]/other[1]])
        elif self.check_nb(other):
            return Vector([self[0]/other, self[1]/other])
        else:
            raise TypeError('Could not perform div operation')

    def __rtruediv__(self, other):
        if self.check(other):
            return Vector([other[0]/self[0], other[1]/self[1]])
        elif self.check_nb(other):
            return Vector([other/self[0], other/self[1]])
        else:
            raise TypeError('Could not perform rdiv operation')

    def dot(self, other):
        if self.check(other):
            return self[0]*other[0]+self[1]*other[1]
        else:
            raise TypeError('Could not perform dot operation')

    def cross(self, other):
        if self.check(other):
            return self[0]*other[1]-self[1]*other[0]
        else:
            raise TypeError('Could not perform dot operation')

    def norm(self):
        return (self[0]**2+self[1]**2)**0.5

    def abs(self):
        return Vector([abs(self[0]), abs(self[1])])

    def unit(self):
        norm = self.norm()
        return Vector([self[0]/norm, self[1]/norm])

    def orth(self):
        return Vector([-self[1], self[0]])

    def rot(self, other):
        '''
        Inputs:
        -------
        other: vector

        Returns
        -------
        vector: rotation around z of self by an angle given by other w.r.t to x
        '''
        if self.check(other):
            unitOther = other.unit()
            return Vector([self.dot(unitOther.refx()), self.dot(unitOther.orth().refy())])
        else:
            raise TypeError('Could not perform rdiv operation')

    def px(self):
        return Vector([self[0], 0])

    def py(self):
        return Vector([0, self[1]])

    def refx(self, offset=0):
        return Vector([self[0], -self[1]+2*offset])

    def refy(self, offset=0):
        return Vector([-self[0]+2*offset, self[1]])


class Circuit(object):

    modeler = None
    design = None

    trackObjects, gapObjects, bondwireObjects, maskObjects = [], [], [], []
    layers = {}
    ports = {}
    ports_dc = {}
    all_points = []
    all_points_val = []

    def __init__(self, design=None, modeler=None):
        if design is not None and modeler is not None:
            Circuit.modeler = modeler
            Circuit.design = design

            self.modeler.set_units('mm')
        '''
        trackObjects (list of strings) e.g. ['readout_track', 'transmon_pad']
        gapObjects (list of strings) e.g. ['readout_box', 'transmon_box']
        One starts from a ground plane
        1. Substract the gapObject
        2. Add the trackObject

        throught the code, one populates these two lists with objects, and at the very end
        apply the above points 1 and 2
        '''
    
    def new_layer(self, name):
        '''
        name is a string
        '''
        self.layers[name] = {}
        self.layers[name]['trackObjects'] = []
        self.layers[name]['gapObjects'] = []
        
        
    def val(self, name): # to use if you want to compare two litt expressions
        if isinstance(name, list):
            name_list = []
            for elt in name:
                name_list.append(self.design.eval_var_str(elt))
            if isinstance(name, Vector):
                return Vector(name_list)
            else:
                return name_list
        else:
            return self.design.eval_var_str(name)

    def set_variable(self, name, iVal):
        """
        Inputs:
        -------
        name (str): name of the variable in HFSS e.g. 'chip_length'
        iVal (str, VarStr, float): value of the variable if str will try to analyse the unit
        Returns:
        --------
        takes create a design variable with name name and value iVal
        """
        variable = self.design.set_variable(name, iVal)
        self.__dict__[name] = variable
        return variable
    
    def get_variable_value(self, name):
        val = self.design.get_variable_value(name)
        return val
    
    def get_variable_names(self):
        names = self.design.get_variable_names()
        return names

    def available_ports(self):
        return self.ports.keys()

    def assign_lumped_RLC(self, iObj, ori, iVal, iSuff="LumpRLC"):
        """
        Inputs:
        -------
        iObj (Rect): name of the object e.g. 'transmon_junction'
        ori (Vector): (x0, y0) direction of the current line
        iVal (str or VarStr,)*3: resistance, inductance and capacitance value
                                        e.g. (R1, Lj, '50fF')
        iSuff (str): suffix to add to iObj name
        Returns:
        --------
        takes existing name object and assigns boundary condition lumped RLC, L=iVal nH, R=0, C=0, along current_line
        """
        r, l, c = iVal
        if self.val(ori[1])==0:
            axis='X'
        elif self.val(ori[0])==0:
            axis='Y'
        else:
            raise ValueError('Expect good orientation to assign_RLC to rectangle')

        iObj.make_rlc_boundary(axis, r=r, l=l, c=c, name=iSuff)


    def assign_lumped_RLC_ArbitratryLine(self, iObj, iVal, start, end, iSuff="LumpRLC"):
        """
        Inputs:
        -------
        iObj (Rect): name of the object e.g. 'transmon_junction'
        iVal (str or VarStr,)*3: resistance, inductance and capacitance value
                                        e.g. (R1, Lj, '50fF')
        start (str or VarStr,)*3: 3D coordinate of the start of the current line
        end (str or VarStr,)*3: 3D coordinate of the end of the current line
        iSuff (str): suffix to add to iObj name
        Returns:
        --------
        takes existing name object and assigns boundary condition lumped RLC, L=iVal nH, R=0, C=0, along current_line
        """
        r, l, c = iVal
        
        start = [self.modeler.eval_var_str(s, unit="meter") for s in start]
        end = [self.modeler.eval_var_str(s, unit="meter") for s in end]
        
        name = str(iObj)+'_'+iSuff
        
        self.modeler._make_lumped_rlc(r, l, c, start, end, ["Objects:=", [iObj]], name=name)


    def assign_perfE(self, iObj, iSuff='PerfE'):
        '''
        Assigns boundary condition PerfE of name name to object iObject

        Inputs:
        -------
        iObject: name of the object e.g. transmon_pad
        name: name of the perfect E e.g. PerfE1

        '''
        self.modeler.assign_perfect_E(iObj, name=iSuff)

    def assign_perfect_E_faces(self, iObj):
        # very peculiar to Si cavity
        return self.modeler.assign_perfect_E_faces(iObj)
        
        
    def make_material(self, iObj, material='vacuum'):
        '''
        Changes the material type of a 3D object

        Inputs:
        -------
        iObject: name of the object
        material: 'pec', 'vacuum'... mind the syntax 

        '''
        if material is 'vacuum':
            self.modeler.make_material(iObj, material="\"vacuum\"")
        if material is 'pec':
            self.modeler.make_material(iObj, material="\"pec\"")
        if material is 'Rogers':
            self.modeler.make_material(iObj, material="\"Rogers RO4003 (tm)\"")
       
    def delete(self, iObj):
        self.modeler.delete(iObj)



    def unite(self, iObjs, name=None):
        '''
        Performs the unions of the elements of iObjects

        Input:
        iObjects (list of strings): list of object names e.g. ['transmon_pad','transmon_pad_extension']

        Returns:
        string: name of the merged object: iObjects[0]
        '''
        if len(iObjs) > 1:
            iObj = self.modeler.unite(iObjs, name=name)
        else:
            print('Watch out: trying to unit a single element list --> renaming can be an issue')
            iObj = iObjs[0]
        return iObj
    
    
    def intersect(self, iObjs, keep_originals=False):
        '''
        Performs the intersection of the elements of iObjects

        Input:
        iObjects (list of strings): list of object names e.g. ['transmon_pad','transmon_pad_extension']

        Returns:
        string: name of the merged object: iObjects[0]
        '''
        iObj = self.modeler.intersect(iObjs, keep_originals)
        return iObj


    def subtract(self, iObjBlank, iObjTools, keep_originals=False):
        '''
        suBTracts iObjectTools from iObjectBlanc

        Inputs:
        -------
        iObjectBlanc (string) : HFSS object name e.g. 'ground_plane'
        iObjectTools (list) : HFSS object name e.g. ['readout_box', 'res_gap']

        '''
        iObjBlank = self.modeler.subtract(iObjBlank, iObjTools, keep_originals=keep_originals)
        return iObjBlank
    
    def duplicate_along_line(self, iObject, iVector, n=2):
        '''
        copies iObjand moves the copy by iVector
        Inputs:
        -------
        iObject (string) : HFSS object name e.g. 'ground_plane'
        iVector (list) : list of iVector's coordinates
        '''
        while len(iVector) < 3:
            iVector.append(0)
            self.modeler.duplicate_along_line(iObject, iVector, n=n)

    def rename(self, iObj, name):
        '''
        rename iObj by "name"

        Inputs:
        -------
        name (string) : the new name
        iObj (string) : HFSS object name e.g. ['readout_box', 'res_gap']

        '''
        name = self.modeler.rename_obj(iObj, name)
        return name
    
    
    def copy(self, iObject, name=None, bruteforce=False):
        '''
        copies iObj

        Inputs:
        -------
        iObject (string) : HFSS object name e.g. 'ground_plane'
        name (string) : the new name

        '''
        if not bruteforce:
            new_Obj = self.modeler.copy(iObject)
            if name is not None:
                name = self.rename(new_Obj, name)
    #            self.__dict__[name] = name
        else:
            vertices = self.get_vertex_ids(iObject)
            print(vertices)
            new_Obj = self.modeler.draw_polyline(vertices)
        return new_Obj
    
    
    def duplicate_along_line(self, iObject, iVector):
        '''
        copies iObjand moves the copy by iVector

        Inputs:
        -------
        iObject (string) : HFSS object name e.g. 'ground_plane'
        iVector (list) : list of iVector's coordinates

        '''
        while len(iVector) < 3:
            iVector.append(0)
        self.modeler.duplicate_along_line(iObject, iVector)

    def translate(self, iObject, iVector):
        '''
        moves iObj by iVector

        Inputs:
        -------
        iObject (string) : HFSS object name e.g. 'ground_plane'
        iVector (list) : list of iVector's STRING coordinates

        '''
        while len(iVector) < 3:
            print('size of translation vector artificially increased')
            iVector.append('0mm')
        iObject = self.modeler.translate(iObject, iVector)
        return iObject
        
    
    def sweep_along_path(self, iObject, path_length):
        '''
        creates 3D structure from 2D one

        Inputs:
        -------
        iObject (string) : HFSS object name e.g. 'ground_plane'
        path_lenght (list) : height of 3D structure

        '''
        iObject = self.modeler._sweep_along_path(path_length, iObject)
        return iObject
    
    def thicken_sheet(self, iSheet, iThickness): 
        '''
        creates 3D structure from 2D one (new and more specific version of the sweep_along_path function)

        Inputs:
        -------
        iSheet (string) : HFSS object name of a 2D structure e.g. 'ground_plane'
        iThickness (list) : height of 3D structure

        '''
        iThickness = parse_entry((iThickness,))
        self.modeler.separate_bodies(iSheet)
        iSheet_matched = self.modeler.get_matched_object_name(iSheet)
#        print(iSheet_matched)
        self.modeler.sweep_along_vector(iSheet_matched, ['0','0',iThickness])
        self.unite(iSheet_matched)
          
        
    def mirrorZ(self, iObject):
        '''
        transforms iObject through a mirror operation along Z
        this function is design to correct bugs when usin the thicken_sheet function: some sheets are drawn
            in such a way that they thicken in the wrong direction, hence need for flipping with mirrorZ
        
        '''
        iObject = self.modeler.mirrorZ(iObject)
        

    def homothetic(self, iObjects, incOdec, ratio='0.4um', N=6):
        '''
        Inputs:
        -------
        iObjects (list)  : list of HFSS object names
        incOdec          : 'increase' or 'decrease' the surface
        ratio            : ratio of the transformation (string with HFSS-readable unit)
        N                : resolution of the process
        
        Outputs:
        -------
        
        Remarks:
        -------
        not stable for N>6 but 6 is enough if ratio/width_track is not too big
        ''' 
        
        theta_list = (2*np.pi/N)*np.linspace(0,N-1,N)
        for ii, obj in enumerate(iObjects):
            subObjects = [obj] 
            for t, theta in enumerate(theta_list):
                u_x = parse_entry(ratio)*np.cos(theta)
                u_y = parse_entry(ratio)*np.sin(theta)
                self.duplicate_along_line(obj, [u_x, u_y])
                subObjects.append(obj+'_'+str(t+1)) #careful: if iObjects has already been copied, _1 might already exist
            if incOdec is 'increase':
                iObjects[ii] = self.unite(subObjects)
            if incOdec is 'decrease':
                iObjects[ii] = self.intersect(subObjects)
            
            
            
    # main function
    def draw(self, iObj, iPoints, closed=True): #assume everything is in the plane z=0
        '''
        Inputs
        ------
        name : object name in HFSS e.g. transmon_pad
        iPoints : list of tuples [(x0, y0), (x1, y1), ..., (xn, yn)]

        '''
        points = []
        points_to_append = []
        points_to_append_val = []
        
        prev_val_x=0
        prev_val_y=0
        for iPoint in iPoints:
            val_x = self.val(iPoint[0])
            val_y = self.val(iPoint[1])
            x = iPoint[0]
            y = iPoint[1]
#            print('{} == {}'.format(val_x, prev_val_x))
#            print('{} == {}\n'.format(val_y, prev_val_y))
            if points == []:
                points.append([x, y, 0])
                points_to_append.append([x, y])
                points_to_append_val.append([val_x, val_y])
                prev_val_x = val_x
                prev_val_y = val_y
            elif not equal_float(prev_val_x, val_x) or not equal_float(prev_val_y, val_y):
                points.append([x, y, 0])
                points_to_append.append([x, y])
                points_to_append_val.append([val_x, val_y])
                prev_val_x = val_x
                prev_val_y = val_y
            else:
                pass
#                print('Warning: Found two overlapping points while creating the polyline, supressed one')
        iObj = self.modeler.draw_polyline(points, closed=closed, name=iObj)
        self.all_points += points_to_append
        self.all_points_val += points_to_append_val
        return iObj

    def draw_wirebond_port(self, iObj, iIn): # draw wire bond from port name
        iIn = self.ports[iIn]
        pos, ori = Vector(iIn[POS]), Vector(iIn[ORI])
        track, gap = iIn[TRACK], iIn[GAP]

        width = (track+2*gap)*1.5
        iObj = self.modeler.draw_wirebond(pos, ori, width, name = iObj, material='perfect conductor', solve_inside=False)
        self.bondwireObjects.append(iObj)
        return iObj

    def draw_wirebond(self, iObj, pos, ori, width): # draw wire bond from dimensions
        iObj = self.modeler.draw_wirebond(pos, ori, width, name = iObj, material='perfect conductor', solve_inside=False)
        self.bondwireObjects.append(iObj)
        return iObj

    def draw_box(self, name, pos, iSize, iMaterial='vaccum'):
        box = self.modeler.draw_box_corner(pos, iSize, material=iMaterial, name=name)
        self.__dict__[box] = box
        return box

    def draw_cylinder(self, name, pos, iSize, axis, iMaterial='vaccum'):
        box = self.modeler.draw_cylinder(pos, iSize[0], iSize[1], axis, material=iMaterial, name=name)
        self.__dict__[box] = box
        return box
    
    def draw_trapeze(self, name, pos, z, size, height, angle=54.74*np.pi/180):
        # pos, is central vector position
        # z is position of basis
        # size is the diagonal vector of the basis
        #   assume size is positive
        # height is the height of the pyramid, can be positive or negative
        pos, z, size, height = parse_entry((pos, z, size, height))
        pos = Vector(pos)
        size = Vector(size)
        rect_base = self.draw_rect_center(name, pos, size, z=z)
        size_top = size - Vector([height/np.sign(self.val(height))*2/np.tan(angle)*np.sign(self.val(size[0])), 
                                  height/np.sign(self.val(height))*2/np.tan(angle)*np.sign(self.val(size[1]))])
        rect_top = self.draw_rect_center('rect2', pos, size_top, z=z+height)
        
        pyramid = self.modeler.connect_faces(rect_base, rect_top)
        return pyramid
    
    def draw_disk(self, name, pos, iSize, axis):
        disk = self.modeler.draw_disk(pos, iSize, axis, name=name)
        self.__dict__[disk] = disk
        return disk

    def draw_rect_center(self, name, pos, iSize, z=0):
        pos = [pos[0], pos[1], z]
        size = [iSize[0], iSize[1], 0]
        rect = self.modeler.draw_rect_center(pos, size, name=name)
        self.__dict__[rect] = rect
        corner1 = pos+iSize/2
        corner2 = pos-iSize/2
        self.all_points += [[corner1[0], corner1[1]], 
                            [corner2[0], corner2[1]]]
        self.all_points_val += [[self.val(corner1[0]), self.val(corner1[1])], 
                                [self.val(corner2[0]), self.val(corner2[1])]]
        return rect

    def draw_rect(self, name, pos, iSize, z=0):
        pos = [pos[0], pos[1], z]
        size = [iSize[0], iSize[1], 0]
        rect = self.modeler.draw_rect_corner(pos, size, name=name)
        self.__dict__[rect] = rect
        corner1 = pos
        corner2 = pos+iSize
        self.all_points += [[corner1[0], corner1[1]], 
                            [corner2[0], corner2[1]]]
        self.all_points_val += [[self.val(corner1[0]), self.val(corner1[1])], 
                                [self.val(corner2[0]), self.val(corner2[1])]]
        return rect
    
    def create_object_from_face(self, name):
        face = self.modeler.create_object_from_face(name)
        return face
    
    def get_extent(self, margin=0):
        margin = parse_entry(margin)
        
        points = self.all_points
        points_val = np.array(self.all_points_val)
        
        bottom = points[np.argmin(points_val[:,1])][1]-margin
        top = points[np.argmax(points_val[:,1])][1]+margin
        left = points[np.argmin(points_val[:,0])][0]-margin
        right = points[np.argmax(points_val[:,0])][0]+margin
        
        return bottom, top, left, right

    def key_elt(self, name='key_elt_0', pos=[0,0], ori=[1,0]):
        obj = KeyElt(name, pos=pos, ori=ori)
        self.__dict__[name] = obj
        return obj

    def connect_elt(self, name='connect_elt_0', iIn='iInt', iOut=None, layer=None):
        obj = ConnectElt(name=name, iIn=iIn, iOut=iOut, layer=layer)
        self.__dict__[name] = obj
        return obj    
    
    def split_dc_ports(self, port, subgroups, subnames, iGap=None):
        ''' 
        Given a n-multiport 'port', it returns a set of len(subgroups) multiports
        of size subgroups[ii] and names subnames[ii]
        '''
        if len(subgroups) is not len(subnames):
            raise ValueError('number of sub-groups and names do not match')

        oldport = self.ports_dc[port]
        for ii in range(len(subgroups)):
            if subgroups[ii] is 1:
                if iGap is not None:
                    vec_center = Vector([oldport[0][0], oldport[0][1]])
                    vec_rel = Vector([0, oldport[2][ii]])
                    vec_rel = vec_rel.rot(oldport[1])
                    newport = [vec_center+vec_rel, oldport[1], oldport[3][ii], iGap]
                    self.ports[port+'_'+subnames[ii]] = newport
                else:
                    print(port+'_'+subnames[ii]+' defined as a single dc port')
                    newport = [oldport[0], \
                               oldport[1], \
                               oldport[2][sum(subgroups[:ii]):sum(subgroups[:ii+1])], \
                               oldport[3][sum(subgroups[:ii]):sum(subgroups[:ii+1])], \
                               oldport[4][sum(subgroups[:ii]):sum(subgroups[:ii+1])], \
                               subgroups[ii]]
                    self.ports_dc[port+'_'+subnames[ii]] = newport
            else:
                newport = [oldport[0], \
                           oldport[1], \
                           oldport[2][sum(subgroups[:ii]):sum(subgroups[:ii+1])], \
                           oldport[3][sum(subgroups[:ii]):sum(subgroups[:ii+1])], \
                           oldport[4][sum(subgroups[:ii]):sum(subgroups[:ii+1])], \
                           subgroups[ii]]
                self.ports_dc[port+'_'+subnames[ii]] = newport

    def draw_dc_port(self, name):
        '''
        test function to locate and see the orientation of a multiple dc port
        '''
        port = self.ports_dc[name]
        for ii in range(port[-1]):
            vec_center = Vector([port[0][0], port[0][1]])
            vec_rel = Vector([0, port[3][ii]])
            vec_rel = vec_rel.rot(port[1])
            self.draw_rect_center(name+'subport_'+str(ii), vec_center+vec_rel, Vector(['100um','100um']))
    
    def create_dc_layout(self, chip_width, chip_length, pos_neg, homothetic=None):
        
        z_shift = [str(0.5*ii)+'mm' for ii in range(len(self.layers.keys()))]
        for ii, layer in enumerate(self.layers.keys()):
            if len(self.layers[layer]['gapObjects']) > 1:
                gapObject = self.unite(self.layers[layer]['gapObjects'])#, name=layer+'_gapObject')
            else:
                if len(self.layers[layer]['gapObjects']) is not 0:
                    gapObject = self.rename(self.layers[layer]['gapObjects'][0], layer+'_gapObject')
                else:
                    print('no gapObjects defined for the layer '+layer)
                    gapObject = None
            if len(self.layers[layer]['trackObjects']) > 1:
                trackObject = self.unite(self.layers[layer]['trackObjects'])#, name=layer+'_trackObject')
            else:
                if len(self.layers[layer]['trackObjects']) is not 0:
                    trackObject = self.rename(self.layers[layer]['trackObjects'][0], layer+'_trackObject')
                else:
                    print('no trackObjects defined for the layer '+layer)
                    trackObject = None
                    
            if pos_neg is 'neg':
                ground_plane = self.draw_rect(layer+'_groud_plane', [0,0], [chip_width, chip_length])
                negatif = self.draw_rect(layer+'_negatif', [0, 0], [chip_width, chip_length])
                if gapObject and trackObject is not None:
                    self.subtract(ground_plane, [gapObject])
                    self.subtract(negatif, [ground_plane]+[trackObject])
                elif gapObject is None: #should do the trick for the pits defined as trackO when creating etch drawing
                    self.subtract(negatif, [trackObject])
                if homothetic is 'increase':
                    layout = self.homothetic([negatif], 'decrease', self.overetch, N=6)
                if homothetic is 'decrease':
                    layout = self.homothetic([negatif], 'increase', self.overetch, N=6)
            if pos_neg is 'pos':
                if gapObject is not None:
                    self.delete(gapObject)
#                HOMOTHETIC HAS NO OUTPUT
                if trackObject is not None:
                    if homothetic is 'increase':
                        layout = self.homothetic(trackObject, 'increase', self.overetch, N=6)
                    if homothetic is 'decrease':
                        layout = self.homothetic(trackObject, 'decrease', self.overetch, N=6)
                    else:
                        layout = trackObject
            if layout:
                self.translate(layout, ['0mm', '0mm', z_shift[ii]])


class KeyElt(Circuit):

    pcb_track = parse_entry('300um')
    pcb_gap = parse_entry('200um')
    is_mask = False
    gap_mask = parse_entry('20um')
    overdev = parse_entry('0um')
    is_overdev = False
    is_litho = False

#    @property
#    def pcb_track(self):
#        return self._pcb_track
#
#    @pcb_track.setter
#    def pcb_track(self, new_pcb_track):
#        self._pcb_track = parse_entry(new_pcb_track)
#
#    @property
#    def pcb_gap(self):
#        return self._pcb_gap
#
#    @pcb_gap.setter
#    def pcb_gap(self, new_pcb_gap):
#        self._pcb_gap = parse_entry(new_pcb_gap)



    def __init__(self, name, pos=[0,0], ori=[1,0]):
        pos, ori = parse_entry((pos, ori))
        self.name = name
        self.pos = Vector(pos)
        self.ori = Vector(ori)
        print(self.name)

    def rot(self, x, y=0):
#        if isinstance(x, Vector):
#            return(Vector(x))
#        else:
        return Vector(x, y).rot(self.ori)

    def append_points(self, coor_list):
        points = [self.pos + self.rot(*coor_list[0])]
        for coor in coor_list[1:]:
            points.append(points[-1] + self.rot(*coor))
        return points
    
    def append_absolute_points(self, coor_list):
        points=[]
        for coor in coor_list:
            points.append(self.pos + self.rot(*coor))
        return points

    def refx_points(self, coor_list, offset=0, absolute=False):
        points=[]
        for ii, coor in enumerate(coor_list):
            if ii>0 and not absolute:
                offset=0
            points.append(Vector(*coor).refx(offset))
        return points

    def refy_points(self, coor_list, offset=0, absolute=False):
        points=[]
        for ii, coor in enumerate(coor_list):
            if ii>0 and not absolute:
                offset=0
            points.append(Vector(*coor).refy(offset))
        return points
    
    def move_points(self, coor_list, move, absolute=False):
        points=[]
        for ii, coor in enumerate(coor_list):
            if ii>0 and not absolute:
                move=[0,0]
            points.append(Vector(*coor)+Vector(*move))
        return points


    def coor(self, vec): # Change of coordinate for a point
        return self.rot(*vec)+self.pos

    def coor_vec(self, vec): # Change of coordinate for a vector
        return self.rot(*vec)
    
    def create_port(self, iTrack=0, iGap=0):
        iTrack, iGap = parse_entry((iTrack, iGap))
        portOut = [self.coor([0,0]), self.coor_vec([1,0]), iTrack+2*self.overdev, iGap-2*self.overdev]
        self.ports[self.name] = portOut
    
    def create_dc_port(self, layer, cut, rel_pos, wid):
        portOut = [self.coor([0,0]), self.coor_vec([1,0]), cut, rel_pos, wid, len(rel_pos)]
        self.ports_dc[layer+'_'+self.name] = portOut
     
        
    def draw_fluxtrappingholes(self, chip_length,chip_width,margin,hole_spacing,hole_size):

        chip_length,chip_width,margin,hole_spacing=parse_entry((chip_length,chip_width,margin,hole_spacing))
        i=0
        Nx=int(self.val(chip_length-2*margin)/self.val(hole_spacing)/2)
        Ny=int(self.val(chip_width-2*margin)/self.val(hole_spacing)/2)
        holes=[]
        for p in range(-Nx,Nx+1):
            for q in range(-Ny,Ny+1):
                y=chip_length/2+p*hole_spacing
                x=chip_width/2+q*hole_spacing
                holes.append(self.draw_rect_center(self.name+'_hole'+str(i), [x,y], [hole_size, hole_size]))
                i+=1
        self.unite(holes)
        
    def draw_connector(self, iTrack, iGap, iBondLength, iSlope=1, pcbTrack=None, pcbGap=None, tr_line=True):
        '''
        Draws a CPW connector for inputs and outputs.

        Inputs:
        -------

        iBondLength: (float) corresponds to dimension a in the drawing
        iSlope (float): 1 is 45 degrees to adpat lengths between Bonding pad and iout
        iLineTest (Bool): unclear, keep False

            ground plane
            +------+
            |       \
            |        \
            |   +-a-+\|
        iIn |   |.....| iOut
            |   +---+/|
            |        /
            |       /
            +------+

        Outputs:
        --------
        returns iIn and recalculated iOut
        '''
        iTrack, iGap = parse_entry((iTrack, iGap))
        iBondLength, iSlope = parse_entry((iBondLength, iSlope))
        
        if pcbGap is not None:
            pcbGap = parse_entry(pcbGap)
            self.pcb_gap = pcbGap
            
        if pcbTrack is not None:
            pcbTrack = parse_entry(pcbTrack)
            self.pcb_track = pcbTrack
            
        adaptDist = (self.pcb_track/2-iTrack/2)/iSlope

        portOut = [self.coor([adaptDist+self.pcb_gap+iBondLength,0]), self.coor_vec([1,0]), iTrack+2*self.overdev, iGap-2*self.overdev]
#        print(self.pos, self.ori)
#        print(adaptDist)
#        print(self.pos+self.ori*(adaptDist+iGap+iBondLength), self.ori)
        points = self.append_points([(self.pcb_gap-self.overdev, self.pcb_track/2+self.overdev),
                                     (iBondLength+self.overdev, 0),
                                     (adaptDist, iTrack/2-self.pcb_track/2),
                                     (0, -iTrack-2*self.overdev),
                                     (-adaptDist, iTrack/2-self.pcb_track/2),
                                     (-iBondLength-self.overdev, 0)])
        self.trackObjects.append(self.draw(self.name+"_track", points))

        points = self.append_points([(self.pcb_gap/2+self.overdev, self.pcb_gap+self.pcb_track/2-self.overdev),
                             (self.pcb_gap/2+iBondLength-self.overdev, 0),
                             (adaptDist, (iGap-self.pcb_gap)+(iTrack-self.pcb_track)*0.5),
                             (0, -2*iGap-iTrack+2*self.overdev),
                             (-adaptDist, (iGap-self.pcb_gap)+(iTrack-self.pcb_track)*0.5),
                             (-(self.pcb_gap/2+iBondLength)+self.overdev, 0)])
        self.gapObjects.append(self.draw(self.name+"_gap", points))

        if self.is_mask:
            points = self.append_points([(self.pcb_gap/2-self.gap_mask, self.pcb_gap+self.pcb_track/2+self.gap_mask),
                             (self.pcb_gap/2+iBondLength+self.gap_mask, 0),
                             (adaptDist, (iGap-self.pcb_gap)+(iTrack-self.pcb_track)*0.5),
                             (0, -2*iGap-iTrack-2*self.gap_mask),
                             (-adaptDist, (iGap-self.pcb_gap)+(iTrack-self.pcb_track)*0.5),
                             (-(self.pcb_gap/2+iBondLength)-self.gap_mask, 0)])
            self.maskObjects.append(self.draw(self.name+"_mask", points))


        if not self.is_litho and tr_line:
            points = self.append_points([(self.pcb_gap/2+self.overdev, self.pcb_track/2+self.overdev),
                                         (self.pcb_gap/2-2*self.overdev, 0),
                                         (0, -self.pcb_track-2*self.overdev),
                                         (-self.pcb_gap/2+2*self.overdev, 0)])
            ohm = self.draw(self.name+"_ohm", points)
            self.assign_lumped_RLC(ohm, self.ori, ('50ohm',0,0))
            self.modeler.assign_mesh_length(ohm, self.pcb_track/10)
#            self.trackObjects.append(ohm)
            points = self.append_points([(self.pcb_gap/2+self.overdev,0),(self.pcb_gap/2-2*self.overdev,0)])
            self.draw(self.name+'_line', points, closed=False)

        self.ports[self.name] = portOut


    def draw_JJ(self, iTrack, iGap, iTrackJ, iLength, iInduct='1nH', fillet=None):
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
        iTrack, iGap, iTrackJ, iLength = parse_entry((iTrack, iGap, iTrackJ, iLength))

        portOut1 = [self.coor([iLength/2,0]), self.coor_vec([1,0]), iTrack, iGap]
        portOut2 = [self.coor([-iLength/2,0]), self.coor_vec([-1,0]), iTrack, iGap]
        self.ports[self.name+'_1'] = portOut1
        self.ports[self.name+'_2'] = portOut2

        junction = self.connect_elt(self.name, self.name+'_1', self.name+'_2')
        pads = junction._connect_JJ(iTrackJ, iInduct=iInduct, fillet=None)
        self.trackObjects.append(pads)
        
        self.gapObjects.append(self.draw_rect_center(self.name+"_cutout", self.coor([0,0]), self.coor_vec([iLength, iTrack+2*iGap])))

    def draw_IBM_tansmon(self,
                       cutout_size,
                       pad_spacing,
                       pad_size,
                       Jwidth,
                       track,
                       gap,
                       Jinduc,
                       nport=1,
                       fillet=None):

        cutout_size, pad_spacing, pad_size, Jwidth, track, gap = parse_entry((cutout_size, pad_spacing, pad_size, Jwidth, track, gap))
        cutout_size = Vector(cutout_size)
        pad_size = Vector(pad_size)

        cutout = self.draw_rect_center(self.name+"_cutout", self.coor([0,0]), self.coor_vec(cutout_size-Vector([2*self.overdev, 2*self.overdev])))

        if self.is_mask:
            self.maskObjects.append(self.draw_rect_center(self.name+"_mask", self.coor([0,0]), self.coor_vec(cutout_size+Vector([track,track]))))
        
        if not self.is_litho:
            mesh = self.draw_rect_center(self.name+"_mesh", self.coor([0,0]), self.coor_vec(cutout_size))
            self.modeler.assign_mesh_length(mesh, 2*track, suff='')

        track_J=Jwidth*4.
        in_junction = [self.coor([-pad_spacing/2+self.overdev, 0]), self.coor_vec([1,0]), track_J+2*self.overdev, 0]
        out_junction = [self.coor([pad_spacing/2-self.overdev, 0]), self.coor_vec([-1,0]), track_J+2*self.overdev, 0]
        junction = self.connect_elt(self.name+'_junction', in_junction, out_junction)
        junction_pads = junction._connect_JJ(Jwidth+2*self.overdev, iInduct=Jinduc, fillet=None)

        right_pad = self.draw_rect_center(self.name+"_pad1", self.coor(Vector(pad_spacing+pad_size[0],0)/2), self.coor_vec(pad_size+Vector([2*self.overdev, 2*self.overdev])))
        left_pad = self.draw_rect_center(self.name+"_pad2", self.coor(-Vector(pad_spacing+pad_size[0],0)/2), self.coor_vec(pad_size+Vector([2*self.overdev, 2*self.overdev])))
        

        pads = self.unite([right_pad, left_pad, junction_pads], name=self.name+'_pads')
        
        if fillet is not None:
            pads.fillet(fillet+self.overdev,[19])
            pads.fillet(0.75*fillet-self.overdev,[18,17,14,13])
            pads.fillet(fillet+self.overdev,[12,11,10,9,8,7])
            pads.fillet(0.75*fillet-self.overdev,[6,5,2,1])
            pads.fillet(fillet+self.overdev,0)
      
        self.trackObjects.append(pads)

        if nport==1:
            self.ports[self.name+'_1'] = [self.coor(cutout_size.px()/2), self.ori, track, gap]
        elif nport==2:
            self.ports[self.name+'_1'] = [self.coor(cutout_size.px()/2), self.ori, track, gap]
            self.ports[self.name+'_2'] = [self.coor(-cutout_size.px()/2), -self.ori, track, gap]
        elif nport==3:
            self.ports[self.name+'_1a'] = [self.coor(cutout_size.px()/2+pad_size.py()/2), self.ori, track, gap]
            self.ports[self.name+'_1b'] = [self.coor(cutout_size.px()/2-pad_size.py()/2), self.ori, track, gap]
            self.ports[self.name+'_2'] = [self.coor(-cutout_size.px()/2), -self.ori, track, gap]

        elif nport==4:
            self.ports[self.name+'_1a'] = [self.coor(cutout_size.px()/2+pad_size.py()/2), self.ori, track, gap]
            self.ports[self.name+'_1b'] = [self.coor(cutout_size.px()/2-pad_size.py()/2), self.ori, track, gap]
            self.ports[self.name+'_2a'] = [self.coor(-cutout_size.px()/2+pad_size.py()/2), -self.ori, track, gap]
            self.ports[self.name+'_2b'] = [self.coor(-cutout_size.px()/2-pad_size.py()/2), -self.ori, track, gap]
        elif nport==5:
            self.ports[self.name+'_1'] = [self.coor(cutout_size.px()/2), self.ori, track+2*self.overdev, gap-2*self.overdev]
            self.ports[self.name+'_2'] = [self.coor(-cutout_size.px()/2), -self.ori, track+2*self.overdev, gap-2*self.overdev]
            self.ports[self.name+'_3a'] = [self.coor(Vector(pad_spacing/2+pad_size[0],-cutout_size[1]/2)), -self.ori.orth() ,track/2+2*self.overdev,gap/2-2*self.overdev]
            if self.is_overdev:
                sub_1 = self.draw_rect(self.name + '_sub_1', self.coor([cutout_size[0]/2, -track/2-gap+self.overdev]), self.coor_vec([-self.overdev, track+2*gap-2*self.overdev]))
                sub_2 = self.draw_rect(self.name + '_sub_2', self.coor([-cutout_size[0]/2, -track/2-gap+self.overdev]), self.coor_vec([self.overdev, track+2*gap-2*self.overdev]))
                sub_3a = self.draw_rect(self.name + '_sub_3a', self.coor([-track/2-gap+self.overdev+pad_spacing/2+pad_size[0],-cutout_size[1]/2]), self.coor_vec([track+2*gap-2*self.overdev, self.overdev]))

        if self.is_overdev:
            cutout = self.unite([cutout, sub_1, sub_2, sub_3a])
            
        self.gapObjects.append(cutout)
        
#        self.draw_rect_center(self.name+"check1", self.coor(cutout_size.px()/2)+self.ori*pad_size[0]/20, self.rot(*pad_size)/10)
#        self.draw_rect_center(self.name+"check2", self.coor(-cutout_size.px()/2)-self.ori*pad_size[0]/20, self.rot(*pad_size)/10)
#        self.draw_rect_center(self.name+"check3", self.coor(Vector(pad_spacing/2+pad_size[0],cutout_size[1]/2))+self.ori.orth()*pad_size[1]/20, self.rot(*pad_size)/10)


    def draw_ZR_transmon(self,
                        cutout_size,
                        pad_spacing,
                        pad_size_right,
                        track_right,
                        gap_right,
                        length_right,
                        spacing_right,
                        short_right,
                        Jwidth,
                        Jinduc,
                        pad_size_left=None,
                        track_left=None,
                        gap_left=None,
                        length_left=None,
                        spacing_left=None,
                        short_left=None,
                        fillet=None):
        
        # Short should be 0 for no short

        parsed = parse_entry((cutout_size,
                              pad_spacing,
                              pad_size_right,
                              track_right,
                              gap_right,
                              length_right,
                              spacing_right,
                              short_right,
                              Jwidth,
                              Jinduc,
                              pad_size_left,
                              track_left,
                              gap_left,
                              length_left, 
                              spacing_left, 
                              short_left))
                        
        (cutout_size,
         pad_spacing,
         pad_size_right,
         track_right,
         gap_right,
         length_right,
         spacing_right,
         short_right,
         Jwidth,
         Jinduc,
         pad_size_left,
         track_left,
         gap_left,
         length_left,
         spacing_left,
         short_left) = parsed
         
        if pad_size_left is None:
            pad_size_left=pad_size_right
        if track_left is None:
            track_left=track_right
        if gap_left is None:
            gap_left=gap_right
        if length_left is None:
            length_left=length_right
        if spacing_left is None:
            spacing_left=spacing_right
        if short_left is None:
            short_left=short_right
        
        cutout_size = Vector(cutout_size)
        pad_size_left = Vector(pad_size_left)
        pad_size_right = Vector(pad_size_right)

        cutout = self.draw_rect_center(self.name+"_cutout", self.coor([0,0]), self.coor_vec(cutout_size-Vector([self.overdev*2, self.overdev*2])))
        if not self.is_litho:
            mesh = self.draw_rect_center(self.name+"_mesh", self.coor([0,0]), self.coor_vec(cutout_size))
        if self.is_mask:
            mask = self.draw_rect_center(self.name+"_mask", self.coor([0,0]), self.coor_vec(cutout_size+Vector([self.gap_mask,self.gap_mask])*2))
        
        track_J=Jwidth*4.
        in_junction = [self.coor([-pad_spacing/2+self.overdev, 0]), self.coor_vec([1,0]), track_J+2*self.overdev, 0]
        out_junction = [self.coor([pad_spacing/2-self.overdev, 0]), self.coor_vec([-1,0]), track_J+2*self.overdev, 0]
        junction = self.connect_elt(self.name+'_junction', in_junction, out_junction)
        junction_pads = junction._connect_JJ(Jwidth+2*self.overdev, iInduct=Jinduc, fillet=None)

        raw_points = [(pad_spacing/2-self.overdev, -pad_size_right[1]/2-self.overdev),
                      (pad_size_right[0]+2*self.overdev, 0),
                      (0, pad_size_right[1]/2-spacing_right-short_right-gap_right-track_right/2+2*self.overdev),
                      (-pad_size_right[0]+length_right, 0),
                      (0, (spacing_right+short_right+gap_right+track_right/2)*2-2*self.overdev),
                      (pad_size_right[0]-length_right, 0),
                      (0, pad_size_right[1]/2-spacing_right-short_right-gap_right-track_right/2+2*self.overdev),
                      (-pad_size_right[0]-2*self.overdev, 0)]
        points = self.append_points(raw_points)
        right_pad = self.draw(self.name+"_pad1", points) 
            
        raw_points = [(-pad_spacing/2+self.overdev, -pad_size_left[1]/2-self.overdev),
                      (-pad_size_left[0]-2*self.overdev, 0),
                      (0, pad_size_left[1]/2-spacing_left-short_left-gap_left-track_left/2+2*self.overdev),
                      (pad_size_left[0]-length_left, 0),
                      (0, (spacing_left+short_left+gap_left+track_left/2)*2-2*self.overdev),
                      (-pad_size_left[0]+length_left, 0),
                      (0, pad_size_left[1]/2-spacing_left-short_left-gap_left-track_left/2+2*self.overdev),
                      (pad_size_left[0]+2*self.overdev, 0)]
        points = self.append_points(raw_points)
        left_pad = self.draw(self.name+"_pad2", points)
        
        right_track_raw_points = [(cutout_size[0]/2-self.overdev,-track_right/2-self.overdev),
                                 (-(cutout_size[0]/2-pad_spacing/2-length_right-gap_right-short_right-spacing_right), 0),
                                 (0, track_right+2*self.overdev),
                                 ((cutout_size[0]/2-pad_spacing/2-length_right-gap_right-short_right-spacing_right), 0)]
        right_track = self.draw(self.name+"_track1", self.append_points(right_track_raw_points))
        
        left_track_raw_points = [(-cutout_size[0]/2+self.overdev,-track_left/2-self.overdev),
                                 ((cutout_size[0]/2-pad_spacing/2-length_left-gap_left-short_left-spacing_left), 0),
                                 (0, track_left+2*self.overdev),
                                 (-(cutout_size[0]/2-pad_spacing/2-length_left-gap_left-short_left-spacing_left), 0)]
        left_track = self.draw(self.name+"_track2", self.append_points(left_track_raw_points))
        
        if short_right!=0:
            raw_points = [(cutout_size[0]/2-self.overdev,-track_right/2-gap_right+self.overdev),
                          (-(cutout_size[0]/2-pad_spacing/2-length_right-short_right-spacing_right)+2*self.overdev, 0),
                          (0, 2*gap_right+track_right-2*self.overdev),
                          ((cutout_size[0]/2-pad_spacing/2-length_right-short_right-spacing_right)-2*self.overdev, 0),
                          (0, short_right+2*self.overdev),
                          (-(cutout_size[0]/2-pad_spacing/2-length_right-spacing_right), 0),
                          (0, -(2*gap_right+track_right+2*short_right)-2*self.overdev),
                          ((cutout_size[0]/2-pad_spacing/2-length_right-spacing_right), 0)]
            points = self.append_points(raw_points)
            right_short = self.draw(self.name+"_short1", points) 
            
        if short_left!=0:
            raw_points = [(-cutout_size[0]/2+self.overdev,-track_left/2-gap_left+self.overdev),
                          ((cutout_size[0]/2-pad_spacing/2-length_left-short_left-spacing_left)-2*self.overdev, 0),
                          (0, 2*gap_left+track_left-2*self.overdev),
                          (-(cutout_size[0]/2-pad_spacing/2-length_left-short_left-spacing_left)+2*self.overdev, 0),
                          (0, short_left+2*self.overdev),
                          ((cutout_size[0]/2-pad_spacing/2-length_left-spacing_left), 0),
                          (0, -(2*gap_left+track_left+2*short_left)-2*self.overdev),
                          (-(cutout_size[0]/2-pad_spacing/2-length_left-spacing_left), 0)]
            points = self.append_points(raw_points)
            left_short = self.draw(self.name+"_short2", points) 
            
        if fillet is not None:
            right_track.fillet(track_right/2-eps+self.overdev,[2,1])
            left_track.fillet(track_left/2-eps+self.overdev,[2,1])
            cutout.fillet(cutout_size[1]/6-self.overdev,[0,1,2,3])
            if not self.is_litho:
                mesh.fillet(cutout_size[1]/6,[0,1,2,3])
            
            if self.is_mask:
                mask.fillet(cutout_size[1]/6+self.gap_mask,[0,1,2,3])
            
            if short_right!=0:
                right_short.fillet(track_right/2+gap_right+short_right-eps+self.overdev,[6,5])
                right_short.fillet(track_right/2+gap_right-eps-self.overdev,[2,1])
                
                right_quarter_up = self.draw_quarter_circle('quarter_up1', (pad_size_right[1]/2-spacing_right-short_right-gap_right-track_right/2)/4-self.overdev, [cutout_size[0]/2-self.overdev, track_right/2+gap_right+short_right+self.overdev], ori=[-1,1])
                right_quarter_down = self.draw_quarter_circle('quarter_down1', (pad_size_right[1]/2-spacing_right-short_right-gap_right-track_right/2)/4-self.overdev, [cutout_size[0]/2-self.overdev, -(track_right/2+gap_right+short_right)-self.overdev], ori=[-1,-1])
                right_short = self.unite([right_short, right_quarter_up, right_quarter_down])
            else:
                right_quarter_up = self.draw_quarter_circle('quarter_up1', (pad_size_right[1]/2-spacing_right-gap_right-track_right/2)/4+self.overdev, [cutout_size[0]/2-self.overdev, track_right/2+gap_right-self.overdev], ori=[1,1])
                right_quarter_down = self.draw_quarter_circle('quarter_down1', (pad_size_right[1]/2-spacing_right-gap_right-track_right/2)/4+self.overdev, [cutout_size[0]/2-self.overdev, -(track_right/2+gap_right)+self.overdev], ori=[1,-1])
                cutout = self.unite([cutout, right_quarter_up, right_quarter_down])
                
            if short_left!=0:
                left_short.fillet(track_left/2+gap_left+short_left-eps+self.overdev,[6,5])
                left_short.fillet(track_left/2+gap_left-eps-self.overdev,[2,1])
                
                left_quarter_up = self.draw_quarter_circle('quarter_up2', (pad_size_left[1]/2-spacing_left-short_left-gap_left-track_left/2)/4-self.overdev, [-cutout_size[0]/2+self.overdev,track_left/2+gap_left+short_left+self.overdev], ori=[1,1])
                left_quarter_down = self.draw_quarter_circle('quarter_down2', (pad_size_left[1]/2-spacing_left-short_left-gap_left-track_left/2)/4-self.overdev, [-cutout_size[0]/2+self.overdev,-(track_left/2+gap_left+short_left)-self.overdev], ori=[1,-1])
                left_short = self.unite([left_short, left_quarter_up, left_quarter_down])
            else:
                left_quarter_up = self.draw_quarter_circle('quarter_up1', (pad_size_left[1]/2-spacing_left-gap_left-track_left/2)/4+self.overdev, [-cutout_size[0]/2+self.overdev,track_left/2+gap_left-self.overdev], ori=[-1,1])
                left_quarter_down = self.draw_quarter_circle('quarter_down1', (pad_size_left[1]/2-spacing_left-gap_left-track_left/2)/4+self.overdev, [-cutout_size[0]/2+self.overdev,-(track_left/2+gap_left)+self.overdev], ori=[-1,-1])
                cutout = self.unite([cutout, left_quarter_up, left_quarter_down])
                
            right_pad.fillet(pad_size_right[0]/4+self.overdev,[7,0])
            right_pad.fillet((pad_size_right[1]/2-spacing_right-short_right-gap_right-track_right/2)/4+self.overdev,[1,2,5,6])
            right_pad.fillet(track_right/2+gap_right+short_right+spacing_right-eps-self.overdev,[5,6])

            left_pad.fillet(pad_size_left[0]/4+self.overdev,[7,0])
            left_pad.fillet((pad_size_left[1]/2-spacing_left-short_left-gap_left-track_left/2)/4+self.overdev,[1,2,5,6])
            left_pad.fillet(track_left/2+gap_left+short_left+spacing_left-eps-self.overdev,[5,6])
        
        if not self.is_litho:
            self.modeler.assign_mesh_length(mesh, 4*Jwidth)
        
        to_unite = [right_pad, left_pad, right_track, left_track, junction_pads]
        if short_right!=0:
            to_unite.append(right_short)
        if short_left!=0:
            to_unite.append(left_short)
        
        if self.is_overdev:
            added_gap_left = self.draw_rect(self.name+'_added_gap_left', self.coor([-cutout_size[0]/2, -(track_left+2*gap_left)/2+self.overdev]), self.coor_vec([self.overdev, (track_left+2*gap_left)-2*self.overdev]))
            added_gap_right = self.draw_rect(self.name+'_added_gap_right', self.coor([cutout_size[0]/2, -(track_right+2*gap_right)/2+self.overdev]), self.coor_vec([-self.overdev, (track_right+2*gap_right)-2*self.overdev]))

            added_track_left = self.draw_rect(self.name+'_added_track_left', self.coor([-cutout_size[0]/2, -(track_left)/2-self.overdev]), self.coor_vec([self.overdev, (track_left)+2*self.overdev]))
            added_track_right = self.draw_rect(self.name+'_added_track_right', self.coor([cutout_size[0]/2, -(track_right)/2-self.overdev]), self.coor_vec([-self.overdev, (track_right)+2*self.overdev]))
            to_unite = to_unite+[added_track_left, added_track_right]
            
            cutout = self.unite([cutout, added_gap_left, added_gap_right])
            
        pads = self.unite(to_unite, name=self.name+'_pads')
        self.trackObjects.append(pads)
        
        self.gapObjects.append(cutout)
        
        if self.is_mask:
            self.maskObjects.append(mask)
        
        portOut1 = [self.coor([cutout_size[0]/2,0]), self.coor_vec([1,0]), track_right+2*self.overdev, gap_right-2*self.overdev]
        self.ports[self.name+'_1'] = portOut1
        portOut2 = [self.coor([-cutout_size[0]/2,0]), self.coor_vec([-1,0]), track_left+2*self.overdev, gap_left-2*self.overdev]
        self.ports[self.name+'_2'] = portOut2

    def draw_quarter_circle(self, name, fillet, coor, ori=Vector([1,1])):
        ori=Vector(ori)
        temp = self.draw_rect(self.name+"_"+name, self.coor(coor), self.coor_vec(ori*2*fillet))
        temp_fillet = self.draw_rect(self.name+"_"+name+'f', self.coor(coor), self.coor_vec(ori*2*fillet))
        temp_fillet.fillet(fillet, 0)
        
        quarter = self.subtract(temp, [temp_fillet])
        return quarter
    
    def mesh_zone(self, zone_size, mesh_length):
        zone_size, mesh_length = parse_entry((zone_size, mesh_length))
        if not self.is_litho:
            zone = self.draw_rect_center(self.name, self.coor([0, 0]), self.coor_vec(zone_size))
            self.modeler.assign_mesh_length(zone, mesh_length)
            
    def cutout(self, zone_size):
        zone_size = parse_entry(zone_size)
        self.maskObjects.append(self.draw_rect_center(self.name, self.coor([0, 0]), self.coor_vec(zone_size)))
        
        
    def draw_capa(self, iTrack, iGap, pad_spacing, pad_size, half=False):
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


              +--+  +--+
              |  |  |  |
            +-+  |  |  +-+
        iIn |    |  |    | iOut
            +-+  |  |  +-+
              |  |  |  |
              +--+  +--+
        '''
        _pos = self.pos
        if half:
            self.pos = self.pos - self.ori*(pad_spacing/2+pad_size[0]+iGap)
            self.pos2 = self.pos
        iTrack, iGap, pad_spacing, pad_size = parse_entry((iTrack, iGap, pad_spacing, pad_size))
        pad_size = Vector(pad_size)

        portOut1 = [self.pos+self.ori*(pad_spacing/2+pad_size[0]+iGap), self.ori, iTrack+2*self.overdev, iGap-2*self.overdev]
        if not half:
            portOut2 = [self.pos-self.ori*(pad_spacing/2+pad_size[0]+iGap), -self.ori, iTrack+2*self.overdev, iGap-2*self.overdev]

        raw_points = [(pad_spacing/2-self.overdev, pad_size[1]/2+self.overdev),
                      (pad_size[0]+2*self.overdev, 0),
                      (0, -(pad_size[1]-iTrack)/2),
                      (iGap-self.overdev, 0),
                      (0, -iTrack-2*self.overdev),
                      (-iGap+self.overdev, 0),
                      (0, -(pad_size[1]-iTrack)/2),
                      (-pad_size[0]-2*self.overdev, 0)]
        points = self.append_points(raw_points)
        right_pad = self.draw(self.name+"_pad1", points)

        if not half:
            points = self.append_points(self.refy_points(raw_points))
            left_pad = self.draw(self.name+"_pad2", points)

        padlist = [right_pad]
        if not half:
            padlist.append(left_pad)
        pads = self.unite(padlist, name=self.name+'_pads')
        self.trackObjects.append(pads)

        jj = 2 if not half else 1
        pos_cutout = self.pos if not half else self.pos+self.ori*(pad_spacing/4+pad_size[0]/2+iGap/2)
        cutout = self.draw_rect_center(self.name+"_cutout", pos_cutout, self.coor_vec([pad_spacing*jj/2 + jj*pad_size[0]+jj*iGap-jj*self.overdev, pad_size[1] + 2*iGap-2*self.overdev]))
        
        if self.is_overdev:
            sub_1 = self.draw_rect(self.name + '_sub_1', self.coor([pad_spacing/2+pad_size[0]+iGap, -iTrack/2-iGap+self.overdev]), self.coor_vec([-self.overdev, iTrack+2*iGap-2*self.overdev]))
            if not half:
                sub_2 = self.draw_rect(self.name + '_sub_1', self.coor([-pad_spacing/2-pad_size[0]-iGap, -iTrack/2-iGap+self.overdev]), self.coor_vec([self.overdev, iTrack+2*iGap-2*self.overdev]))
            cutout_list = [cutout, sub_1]
            if not half:
                cutout_list.append(sub_2)
            cutout = self.unite(cutout_list)
        
        self.gapObjects.append(cutout)
        if self.is_mask:
            self.maskObjects.append(self.draw_rect_center(self.name+"_mask", self.coor([0,0]), self.coor_vec([pad_spacing + 2*pad_size[0]+4*iGap, pad_size[1] + 4*iGap])))
        if not self.is_litho:
            self.draw(self.name+"_mesh", points)
            self.modeler.assign_mesh_length(self.name+"_mesh",iTrack)

        self.ports[self.name+'_1'] = portOut1
        if not half:
            self.ports[self.name+'_2'] = portOut2

        self.pos = _pos
        
    def draw_capa_inline(self, iTrack, iGap, capa_length, pad_spacing, n_pad=1, iTrack_capa=None, iGap_capa=None, premesh=True, tight=False): #iGap_capa is added gap
        
        iTrack, iGap, capa_length, pad_spacing, n_pad = parse_entry((iTrack, iGap, capa_length, pad_spacing, n_pad))
        if iTrack_capa is not None:
            iTrack_capa = parse_entry((iTrack_capa))
            usual=False
        else:
            iTrack_capa = iTrack
            usual=True
            
        if iGap_capa is not None:
            iGap_capa = parse_entry((iGap_capa))
        else:
            iGap_capa=0

        drawn_pads = []
        if n_pad==1:
            drawn_pads.append(self.draw_rect(self.name+'_pad_left', self.coor([-capa_length/2,-iTrack/2-self.overdev]), self.coor_vec([capa_length/2-pad_spacing/2+self.overdev,iTrack+2*self.overdev])))
            drawn_pads.append(self.draw_rect(self.name+'_pad_right', self.coor([capa_length/2,-iTrack/2-self.overdev]), self.coor_vec([-(capa_length/2-pad_spacing/2)-self.overdev,iTrack+2*self.overdev])))
        else:
            pad_width = (iTrack_capa-(n_pad-1)*pad_spacing)/n_pad
            if not usual:
#                pad_length = capa_length-pad_spacing-2*iTrack
                pad_length = capa_length-pad_spacing-2*pad_spacing
#                offset = iTrack
                offset = pad_spacing
            else:
                pad_length = capa_length-pad_spacing
                offset = 0
            curr_height = -iTrack_capa/2
            pad_size = Vector([pad_length, pad_width])
            for ii in range(int(n_pad/2)):
                drawn_pads.append(self.draw_rect(self.name+"_pad"+str(ii), self.coor([-capa_length/2+offset, curr_height-self.overdev]), self.coor_vec(pad_size+Vector([self.overdev, 2*self.overdev]))))
                drawn_pads.append(self.draw_rect(self.name+"_pad"+str(ii)+'b', self.coor([-capa_length/2+offset+pad_spacing-self.overdev, curr_height+pad_width+pad_spacing-self.overdev]), self.coor_vec(pad_size+Vector([self.overdev, 2*self.overdev]))))
                curr_height = curr_height+2*(pad_width+pad_spacing)
            if n_pad%2!=0:
                drawn_pads.append(self.draw_rect(self.name+"_pad"+str(ii+1), self.coor([-capa_length/2+offset, curr_height-self.overdev]), self.coor_vec(pad_size+Vector([self.overdev, 2*self.overdev]))))
            else:
                curr_height = curr_height-2*(pad_width+pad_spacing)
            if not usual:
#                drawn_pads.append(self.draw_rect(self.name+'_connect_left', self.coor([-capa_length/2-self.overdev, -iTrack_capa/2-self.overdev]), self.coor_vec([iTrack+2*self.overdev, curr_height+iTrack_capa/2+pad_width+2*self.overdev])))
                drawn_pads.append(self.draw_rect(self.name+'_connect_left', self.coor([-capa_length/2-self.overdev, -iTrack_capa/2-self.overdev]), self.coor_vec([offset+2*self.overdev, curr_height+iTrack_capa/2+pad_width+2*self.overdev])))

                if n_pad%2!=0:
#                    drawn_pads.append(self.draw_rect(self.name+'_connect_right', self.coor([capa_length/2+self.overdev, -iTrack_capa/2+pad_width+pad_spacing-self.overdev]), self.coor_vec([-iTrack-2*self.overdev, curr_height+(iTrack_capa/2-2*pad_width-2*pad_spacing)+pad_width+2*self.overdev])))
                    drawn_pads.append(self.draw_rect(self.name+'_connect_right', self.coor([capa_length/2+self.overdev, -iTrack_capa/2+pad_width+pad_spacing-self.overdev]), self.coor_vec([-offset-2*self.overdev, curr_height+(iTrack_capa/2-2*pad_width-2*pad_spacing)+pad_width+2*self.overdev])))
                else:
#                    drawn_pads.append(self.draw_rect(self.name+'_connect_right', self.coor([capa_length/2+self.overdev, -iTrack_capa/2+pad_width+pad_spacing-self.overdev]), self.coor_vec([-iTrack-2*self.overdev, curr_height+iTrack_capa/2+pad_width+2*self.overdev])))
                    drawn_pads.append(self.draw_rect(self.name+'_connect_right', self.coor([capa_length/2+self.overdev, -iTrack_capa/2+pad_width+pad_spacing-self.overdev]), self.coor_vec([-offset-2*self.overdev, curr_height+iTrack_capa/2+pad_width+2*self.overdev])))
#                    
#                    drawn_pads.append(self.draw_rect(self.name+'_connect_right_b', self.coor([capa_length/2+self.overdev, -iTrack/2-self.overdev]), self.coor_vec([-iTrack-2*self.overdev, (iTrack/2+iTrack_capa/2)+2*self.overdev])))
#                    drawn_pads.append(self.draw_rect(self.name+'_connect_left_b', self.coor([-capa_length/2-self.overdev, iTrack/2+self.overdev]), self.coor_vec([iTrack+2*self.overdev, -(iTrack/2+iTrack_capa/2)-2*self.overdev])))
                    drawn_pads.append(self.draw_rect(self.name+'_connect_right_b', self.coor([capa_length/2+self.overdev, -iTrack/2-self.overdev]), self.coor_vec([-offset-2*self.overdev, (iTrack/2+iTrack_capa/2)+2*self.overdev])))
                    drawn_pads.append(self.draw_rect(self.name+'_connect_left_b', self.coor([-capa_length/2-self.overdev, iTrack/2+self.overdev]), self.coor_vec([offset+2*self.overdev, -(iTrack/2+iTrack_capa/2)-2*self.overdev])))
    
        if tight:
            portOut1 = [self.pos+self.ori*capa_length/2, self.ori, iTrack_capa-2*pad_width-2*pad_spacing+2*self.overdev, iGap-(iTrack_capa-2*pad_width-2*pad_spacing-iTrack)/2-2*self.overdev]
        else:
            portOut1 = [self.pos+self.ori*capa_length/2, self.ori, iTrack+2*self.overdev, iGap-2*self.overdev]

#            portOut2 = [self.pos-self.ori*capa_length/2, -self.ori, iTrack_capa+2*self.overdev, iGap-(iTrack_capa-iTrack)/2-2*self.overdev]
#        else:
        portOut2 = [self.pos-self.ori*capa_length/2, -self.ori, iTrack+2*self.overdev, iGap-2*self.overdev]

        if self.is_overdev:
            drawn_pads.append(self.draw_rect(self.name+'_added_right', self.coor([capa_length/2, -iTrack/2]), self.coor_vec([-self.overdev, iTrack])))
            drawn_pads.append(self.draw_rect(self.name+'_added_left', self.coor([-capa_length/2, -iTrack/2]), self.coor_vec([self.overdev, iTrack])))

        pads = self.unite(drawn_pads, name=self.name+'_pads')
        
        self.trackObjects.append(pads)
        
        self.gapObjects.append(self.draw_rect_center(self.name+"_cutout", self.coor([0,0]), self.coor_vec([capa_length, iTrack + 2*iGap+2*iGap_capa-2*self.overdev])))
        if self.is_mask:
            self.maskObjects.append(self.draw_rect_center(self.name+"_mask", self.coor([0,0]), self.coor_vec([capa_length, iTrack + 2*iGap+2*iGap_capa +2*self.gap_mask])))
        
        if premesh:
            if not self.is_litho:
                self.draw_rect_center(self.name+"_mesh", self.coor([0,0]), self.coor_vec([capa_length, iTrack_capa+2*self.overdev]))
                self.modeler.assign_mesh_length(self.name+"_mesh",pad_spacing)

        self.ports[self.name+'_1'] = portOut1
        self.ports[self.name+'_2'] = portOut2
        
    def draw_capa_interdigitated(self, iTrack, iGap, teeth_size,gap_size, N_period, fillet):
        '''
        '''
        iTrack, iGap = parse_entry((iTrack, iGap))
        fillet = parse_entry(fillet)
        teeth_size=parse_entry(teeth_size)
        gap_size=parse_entry(gap_size)
        teeth_size = Vector(teeth_size)
        portOut1 = [self.pos+self.ori*(teeth_size[0]+iTrack+iGap), self.ori, iTrack+2*self.overdev, iGap-2*self.overdev]
        portOut2 = [self.pos-self.ori*(teeth_size[0]+iTrack+iGap), -self.ori, iTrack+2*self.overdev, iGap-2*self.overdev]


        N_teeth=2*N_period+1    
        raw_points = [(teeth_size[0], -N_teeth*teeth_size[1]-self.overdev)]
        raw_points.append((teeth_size[0], (-N_teeth+1)*teeth_size[1]))
        for i in range(-N_teeth+1,N_teeth-1,4):
            raw_points.append((-teeth_size[0], i*teeth_size[1]))
            raw_points.append((-teeth_size[0], (i+2)*teeth_size[1]))
            raw_points.append((teeth_size[0], (i+2)*teeth_size[1]))
            raw_points.append((teeth_size[0], (i+4)*teeth_size[1]))
        raw_points.append((-teeth_size[0], (N_teeth-1)*teeth_size[1]))
        raw_points.append((-teeth_size[0], N_teeth*teeth_size[1]+self.overdev))

        points = self.append_absolute_points(raw_points)
        connection = self.draw(self.name+"_capagap", points, closed=False)
        
        
        connection.fillets(fillet)
        raw_points=[(-gap_size-teeth_size[0]+self.overdev,N_teeth*teeth_size[1]+self.overdev),(gap_size-teeth_size[0]-self.overdev,N_teeth*teeth_size[1]+self.overdev)]
        points=self.append_absolute_points(raw_points)
        capagap_starter = self.draw(self.name+'_width', points, closed=False)
        
        capagap = connection.sweep_along_path(capagap_starter)
        
   
        
        raw_points = [(-teeth_size[0]-iTrack-self.overdev, -N_teeth*teeth_size[1]-self.overdev),
                      (-teeth_size[0]-iTrack-self.overdev,-iTrack/2-self.overdev),
                      (-teeth_size[0]-iTrack-iGap,-iTrack/2-self.overdev),
                      (-teeth_size[0]-iTrack-iGap, iTrack/2+self.overdev),
                      (-teeth_size[0]-iTrack-self.overdev, iTrack/2+self.overdev),
                      (-teeth_size[0]-iTrack-self.overdev, N_teeth*teeth_size[1]+self.overdev),
                      (teeth_size[0]+iTrack+self.overdev,  N_teeth*teeth_size[1]+self.overdev),
                      (teeth_size[0]+iTrack+self.overdev,iTrack/2+self.overdev),
                      (teeth_size[0]+iTrack+iGap,iTrack/2+self.overdev),
                      (teeth_size[0]+iTrack+iGap,-iTrack/2-self.overdev),
                      (teeth_size[0]+iTrack+self.overdev, -iTrack/2-self.overdev),
                      (teeth_size[0]+iTrack+self.overdev, -N_teeth*teeth_size[1]-self.overdev)]
        points = self.append_absolute_points(raw_points)
        pads = self.draw(self.name+"_pads", points)
        #####Filets on edges of the capa
        pads.fillet(fillet+self.overdev,11)
        pads.fillet(fillet+self.overdev,6)
        pads.fillet(fillet+self.overdev,5)
        pads.fillet(fillet+self.overdev,0)
        
        pads_sub = self.subtract(pads, [capagap])
        #print(pads_sub.vertices())
        
        #####Filets on edge
        pads.fillet(fillet-self.overdev,73)
        pads.fillet(fillet-self.overdev,70)

        #####Filets on last teeth
        pads.fillet(0.5*fillet+self.overdev,67)
        pads.fillet(0.5*fillet+self.overdev,38)
        pads.fillet(0.5*fillet+self.overdev,10)

        #####Filets on edge
        pads.fillet(fillet-self.overdev,7)
        pads.fillet(fillet-self.overdev,4)

        #####Filets on last teeth
        pads.fillet(0.5*fillet+self.overdev,1)





        if not self.is_overdev:
            self.gapObjects.append(self.draw_rect_center(self.name+"_gap", self.coor([0,0]), self.coor_vec([2*teeth_size[0]+2*iTrack+2*iGap, 2*N_teeth*teeth_size[1]+2*iGap])))
        else:
            raw_points = [(-teeth_size[0]-iTrack-iGap, iTrack/2+iGap-self.overdev),
                          (-teeth_size[0]-iTrack-iGap+self.overdev, iTrack/2+iGap-self.overdev),
                          (-teeth_size[0]-iTrack-iGap+self.overdev, N_teeth*teeth_size[1]+iGap-self.overdev),
                          (teeth_size[0]+iTrack+iGap-self.overdev, N_teeth*teeth_size[1]+iGap-self.overdev),
                          (teeth_size[0]+iTrack+iGap-self.overdev, iTrack/2+iGap-self.overdev),
                          (teeth_size[0]+iTrack+iGap, iTrack/2+iGap-self.overdev),
                          (teeth_size[0]+iTrack+iGap, -iTrack/2-iGap+self.overdev),
                          (teeth_size[0]+iTrack+iGap-self.overdev, -iTrack/2-iGap+self.overdev),
                          (teeth_size[0]+iTrack+iGap-self.overdev, -N_teeth*teeth_size[1]-iGap+self.overdev),
                          (-teeth_size[0]-iTrack-iGap+self.overdev, -N_teeth*teeth_size[1]-iGap+self.overdev),
                          (-teeth_size[0]-iTrack-iGap+self.overdev, -iTrack/2-iGap+self.overdev),
                          (-teeth_size[0]-iTrack-iGap, -iTrack/2-iGap+self.overdev)]
            points = self.append_absolute_points(raw_points)
            self.gapObjects.append(self.draw(self.name+"_gap", points))
        if self.is_mask:
            self.maskObjects.append(self.draw_rect_center(self.name+"_mask", self.coor([0,0]), self.coor_vec([2*teeth_size[0]+2*iTrack+4*iGap, 2*N_teeth*teeth_size[1]+4*iGap])))

            
        if not self.is_litho:
            self.draw(self.name+"_mesh", points)
            self.modeler.assign_mesh_length(self.name+"_mesh",iTrack)

        self.trackObjects.append(pads_sub)



        self.ports[self.name+'_1'] = portOut1
        self.ports[self.name+'_2'] = portOut2
#    

    def draw_squid(self, iTrack, iGap, squid_size, iTrackPump, iGapPump, iTrackSquid=None, iTrackJ=None, Lj_down='1nH', Lj_up=None,  typePump='down', doublePump=False, iSlope=1, iSlopePump=0.5, fillet=None): #for now assume left and right tracks are the same width
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
        if iTrackSquid is None:
            iTrackSquid = iTrack/4
        if iTrackJ is None:
            iTrackJ = iTrackSquid/2
        if Lj_up is None:
            Lj_up = Lj_down
        iTrack, iGap, squid_size, iTrackPump, iGapPump, iTrackSquid, iTrackJ = parse_entry((iTrack, iGap, squid_size, iTrackPump, iGapPump, iTrackSquid, iTrackJ))
        squid_size = Vector(squid_size)

        adapt_dist = squid_size[1]/2 #if slope ==1 !!!!

        #adapta
        raw_points_a = [(squid_size[0]/2+adapt_dist,0),
                        (0, iTrack/2),
                        (-adapt_dist, squid_size[1]/2+iTrackSquid-iTrack/2),
                        (-squid_size[0]/2+iTrackSquid/2, 0),
                        (0, -iTrackSquid),
                        (squid_size[0]/2-iTrackSquid/2, 0)]
        points_a = self.append_points(raw_points_a)
        track_a = self.draw(self.name+"_track_a", points_a)

        raw_points_b = self.refx_points(raw_points_a)
        points_b = self.append_points(raw_points_b)
        track_b = self.draw(self.name+"_track_b", points_b)

        raw_points_c = self.refy_points(raw_points_a)
        points_c = self.append_points(raw_points_c)
        track_c = self.draw(self.name+"_track_c", points_c)

        raw_points_d = self.refy_points(raw_points_b)
        points_d = self.append_points(raw_points_d)
        track_d = self.draw(self.name+"_track_d", points_d)

        #junction up
        print(self.ori)
        in_junction_up = [self.coor([-iTrackSquid/2,squid_size[1]/2+iTrackSquid/2]), self.coor_vec([1,0]), iTrackSquid, 0]
        out_junction_up = [self.coor([iTrackSquid/2,squid_size[1]/2+iTrackSquid/2]), self.coor_vec([-1,0]), iTrackSquid, 0]
        junction = self.connect_elt(self.name+'_junction_up', in_junction_up, out_junction_up)
        junction_pads_up = junction._connect_JJ(iTrackJ, iInduct=Lj_up, fillet=fillet)
        #junction down
        in_junction_down = [self.coor([-iTrackSquid/2,-squid_size[1]/2-iTrackSquid/2]), self.coor_vec([1,0]), iTrackSquid, 0]
        out_junction_down = [self.coor([iTrackSquid/2,-squid_size[1]/2-iTrackSquid/2]), self.coor_vec([-1,0]), iTrackSquid, 0]
        junction = self.connect_elt(self.name+'_junction_down', in_junction_down, out_junction_down)
        junction_pads_down = junction._connect_JJ(iTrackJ, iInduct=Lj_up, fillet=fillet)

        right_track = self.draw_rect_center(self.name+"_added_track1", self.coor([2*(squid_size[0]/2+adapt_dist),0]), self.coor_vec([squid_size[0]+2*adapt_dist, iTrack]))
        left_track = self.draw_rect_center(self.name+"_added_track2", self.coor([-2*(squid_size[0]/2+adapt_dist),0]), self.coor_vec([squid_size[0]+2*adapt_dist, iTrack]))

        squid = self.unite([right_track, left_track, track_a, track_b, track_c, track_d, junction_pads_down, junction_pads_up], name=self.name)
        self.trackObjects.append(squid)


        if fillet is not None:
            fillet=parse_entry(fillet)
            squid.fillet(fillet,[32,31,26,25,24,19,18,16,11,10,9,4,3,0])

        adapt_dist_pump = 4*iTrackPump#(4*iTrackPump - 2*iTrackSquid)/2/iSlopePump


        self.gapObjects.append(self.draw_rect_center(self.name+"_added_gap", self.coor([0,0]), self.coor_vec([3*(squid_size[0]+2*adapt_dist), iTrack+2*iGap])))
        if self.is_mask:
            self.maskObjects.append(self.draw_rect_center(self.name+"_added_gap", self.coor([0,0]), self.coor_vec([3*(squid_size[0]+2*adapt_dist), iTrack+3*iGap])))

        raw_points_adapt_pump_a = [(3/2*(squid_size[0]+2*adapt_dist),-iTrack/2-iGap-iTrackSquid),
                                         (-(squid_size[0]+2*adapt_dist)+iTrackSquid,0),
                                         (iTrackPump/2-iTrackSquid, -adapt_dist_pump),
                                         (iGapPump, 0)]
        if self.is_mask:
            raw_points_adapt_pump_mask = [(3/2*(squid_size[0]+2*adapt_dist)+iGapPump,-iTrack/2-iGap-iTrackSquid-iGapPump),
                                          (0,iTrackSquid+iGapPump),
                                             (-(squid_size[0]+2*adapt_dist)+iTrackSquid-iGapPump-iTrackPump/2,0),
                                             (iTrackPump/2-iTrackSquid, -adapt_dist_pump-iTrackSquid),
                                             (3*iGapPump+iTrackPump/2, 0)]
            self.maskObjects.append(self.draw(self.name+"_cutout_pump_a_mask", self.append_points(raw_points_adapt_pump_mask)))
            raw_points_adapt_pump_mask_b = self.refy_points(raw_points_adapt_pump_mask, offset = (squid_size[0]+2*adapt_dist)/2)
            self.maskObjects.append(self.draw(self.name+"_cutout_pump_b_mask", self.append_points(raw_points_adapt_pump_mask_b)))

        if typePump == 'up' or typePump == 'Up':
            raw_points_adapt_pump_a = self.refx_points(raw_points_adapt_pump_a)
            ori_pump = [0,1]
            pos_pump = [(squid_size[0]+2*adapt_dist)/2, iTrack/2+iGap+iTrackSquid+adapt_dist_pump]
        elif typePump =='down' or typePump == 'Down':
            ori_pump = [0,-1]
            pos_pump = [(squid_size[0]+2*adapt_dist)/2, -iTrack/2-iGap-iTrackSquid-adapt_dist_pump]
        else:
            raise ValueError("typePump should be 'up' or 'down', given %s" % typePump)
        points_adapt_pump_a = self.append_points(raw_points_adapt_pump_a)
        cutout_pump_a=self.draw(self.name+"_cutout_pump_a", points_adapt_pump_a)
        self.gapObjects.append(cutout_pump_a)
        
        raw_points_adapt_pump_b = self.refy_points(raw_points_adapt_pump_a, offset = (squid_size[0]+2*adapt_dist)/2)
        points_adapt_pump_b = self.append_points(raw_points_adapt_pump_b)
        cutout_pump_b=self.draw(self.name+"_cutout_pump_b", points_adapt_pump_b)
        self.gapObjects.append(cutout_pump_b)
        
        if fillet is not None:
            cutout_pump_a.fillet(fillet/2,1)
            cutout_pump_a.fillet(fillet/4,0)
            cutout_pump_b.fillet(fillet/2,1)
            cutout_pump_b.fillet(fillet/4,0)
        if doublePump:
            raw_points_adapt_pump_c = self.refx_points(raw_points_adapt_pump_a)
            points_adapt_pump_c = self.append_points(raw_points_adapt_pump_c)
            self.gapObjects.append(self.draw(self.name+"_cutout_pump_c", points_adapt_pump_c))

            raw_points_adapt_pump_d = self.refx_points(raw_points_adapt_pump_b)
            points_adapt_pump_d = self.append_points(raw_points_adapt_pump_d)
            self.gapObjects.append(self.draw(self.name+"_cutout_pump_d", points_adapt_pump_d))

        portOut1 = [self.pos+self.ori*(squid_size[0]/2+adapt_dist)*3, self.ori, iTrack, iGap]
        portOut2 = [self.pos-self.ori*(squid_size[0]/2+adapt_dist)*3, -self.ori, iTrack, iGap]
        self.ports[self.name+'_1'] = portOut1
        self.ports[self.name+'_2'] = portOut2

        if doublePump:
            portOutpump1 = [self.coor(pos_pump), self.coor_vec(ori_pump), iTrackPump, iGapPump]
            self.ports[self.name+'_pump1'] = portOutpump1
            portOutpump2 = [self.coor(Vector(pos_pump).refx()), -self.coor_vec(ori_pump), iTrackPump, iGapPump]
            self.ports[self.name+'_pump2'] = portOutpump2
        else:
            portOutpump1 = [self.coor(pos_pump), self.coor_vec(ori_pump), iTrackPump, iGapPump]
            self.ports[self.name+'_pump'] = portOutpump1


    def draw_snails(self, iTrack, iGap, array_room, array_offset, iTrackPump, iGapPump, mode='litho', iTrackMinPump=None, iTrackSnail=None, fillet=None, N_snails=1, snail_dict={'loop_width':20e-6, 'loop_length':20e-6, 'length_big_junction':10e-6, 'length_small_junction':2e-6, 'bridge':1e-6, 'bridge_spacing':1e-6}, L_eq = '1nH'): #for now assume left and right tracks are the same width
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
        if iTrackSnail is None:
            iTrackSnail = iTrack/10
        if iTrackMinPump is None:
            iTrackMinPump=iTrackSnail

        iTrack, iGap, array_room, array_offset, iTrackPump, iGapPump, iTrackMinPump, iTrackSnail = parse_entry((iTrack, iGap, array_room, array_offset, iTrackPump, iGapPump, iTrackMinPump, iTrackSnail))

        adapt_dist = iTrack/2 #if slope ==1 !!!!
        
        if self.val(array_offset)>=0:
            raw_points_a = [(array_room/2+adapt_dist-iTrack,-iTrack/2),
                            (iTrack,0),
                            (0, iTrack),
                            (-adapt_dist, -iTrack/2+iTrackSnail/2+array_offset),
                            (0, -iTrackSnail),
                            (-iTrack, 0)]
        else: 
            array_offset = -array_offset
            raw_points_a = [(array_room/2+adapt_dist-iTrack,-iTrack/2),
                (iTrack,0),
                (0, iTrack),
                (-adapt_dist, -iTrack/2+iTrackSnail/2+array_offset),
                (0, -iTrackSnail),
                (-iTrack, 0)]
            raw_points_a = self.refx_points(raw_points_a)
            array_offset = -array_offset
            
        
        points_a = self.append_points(raw_points_a)
        track_a = self.draw(self.name+"_track_a", points_a)
        
        raw_points_c = self.refy_points(raw_points_a)
        points_c = self.append_points(raw_points_c)
        track_c = self.draw(self.name+"_track_c", points_c)

        #snail array
        if mode=='litho':
            in_array = [self.coor([array_room/2, array_offset]), -self.ori, iTrackSnail, 0]
            out_array = [self.coor([-array_room/2, array_offset]), self.ori, iTrackSnail, 0]      
            snail_array = self.connect_elt(self.name+'_array', out_array, in_array)
            snail_track = snail_array._connect_snails2([snail_dict['loop_width'], snail_dict['loop_length']], snail_dict['length_big_junction'], 4, snail_dict['length_small_junction'], 1, N_snails, snail_dict['bridge'], snail_dict['bridge_spacing'])
        
        if 0:
            in_array = [self.coor([array_room/2, array_offset]), -self.ori, iTrackSnail, 0]
            out_array = [self.coor([-array_room/2, array_offset]), self.ori, iTrackSnail, 0]      
            snail_array = self.connect_elt(self.name+'_array', in_array, out_array)
            snail_track = snail_array._connect_JJ(iTrack, iInduct=L_eq)
        
        if mode=='equivalent':
            connect_left = self.draw_rect(self.name+"_left", self.coor([array_room/2,array_offset-iTrackSnail/2]), self.coor_vec([-iTrackSnail, iTrackSnail]))
            connect_right = self.draw_rect(self.name+"_right", self.coor([-array_room/2,array_offset-iTrackSnail/2]), self.coor_vec([iTrackSnail, iTrackSnail]))
            array_eq = self.draw_rect_center(self.name+"_array_eq", self.coor([0,array_offset]), self.coor_vec([array_room-2*iTrackSnail, iTrackSnail]))
            self.assign_lumped_RLC(array_eq, self.ori, (0, L_eq, 0))
            points = self.append_points([(-array_room/2+iTrackSnail,array_offset),(array_room-2*iTrackSnail,0)])
            self.draw(self.name+'_array_eq_line', points, closed=False)
            
#        #junction up
#        print(self.ori)
#        in_junction_up = [self.coor([-iTrackSnail/2,squid_size[1]/2+iTrackSnail/2]), self.coor_vec([1,0]), iTrackSnail, 0]
#        out_junction_up = [self.coor([iTrackSnail/2,squid_size[1]/2+iTrackSnail/2]), self.coor_vec([-1,0]), iTrackSnail, 0]
#        junction = self.connect_elt(self.name+'_junction_up', in_junction_up, out_junction_up)
#        junction_pads_up = junction._connect_JJ(iTrackJ, iInduct=Lj_up, fillet=fillet)
#        
#        #junction down
#        in_junction_down = [self.coor([-iTrackSnail/2,-squid_size[1]/2-iTrackSnail/2]), self.coor_vec([1,0]), iTrackSnail, 0]
#        out_junction_down = [self.coor([iTrackSnail/2,-squid_size[1]/2-iTrackSnail/2]), self.coor_vec([-1,0]), iTrackSnail, 0]
#        junction = self.connect_elt(self.name+'_junction_down', in_junction_down, out_junction_down)
#        junction_pads_down = junction._connect_JJ(iTrackJ, iInduct=Lj_up, fillet=fillet)

        right_track = self.draw_rect_center(self.name+"_added_track1", self.coor([2*(array_room/2+adapt_dist),0]), self.coor_vec([array_room+2*adapt_dist, iTrack+2*self.overdev]))
        left_track = self.draw_rect_center(self.name+"_added_track2", self.coor([-2*(array_room/2+adapt_dist),0]), self.coor_vec([array_room+2*adapt_dist, iTrack+2*self.overdev]))
        if mode=='equivalent':
            squid = self.unite([right_track, left_track, track_a, track_c, connect_left, connect_right], name=self.name+'_temp')
        else:
            squid = self.unite([right_track, left_track, track_a, track_c], name=self.name+'_temp')
        if fillet is not None:
            squid.fillet(iTrack/2,[0, 3, 8, 12])
            
        if 0:
            squid = self.unite([squid, snail_track], name=self.name)
        self.trackObjects.append(squid)

        adapt_dist_pump = 4*iTrackPump
        
        gaps=[]
        masks=[]

        gaps.append(self.draw_rect_center(self.name+"_added_gap", self.coor([0,0]), self.coor_vec([3*(array_room+2*adapt_dist), iTrack+2*iGap-2*self.overdev])))
        if self.is_mask:
            masks.append(self.draw_rect_center(self.name+"_added_gap", self.coor([0,0]), self.coor_vec([3*(array_room+2*adapt_dist), iTrack+2*iGap+2*self.gap_mask])))

        raw_points_adapt_pump_a = [(3/2*(array_room+2*adapt_dist)-self.overdev,-iTrack/2-iGap-iTrackMinPump-self.overdev),
                                   (-((array_room+2*adapt_dist))+iTrackMinPump+2*self.overdev,0),
                                   (iTrackPump/2-iTrackMinPump, -adapt_dist_pump+self.overdev),
                                   (iGapPump-2*self.overdev, 0)]
        if self.is_mask:
            raw_points_adapt_pump_a_mask = [(3/2*(array_room+2*adapt_dist)+self.gap_mask,-iTrack/2-iGap-iTrackMinPump+self.gap_mask),
                                            (-((array_room+2*adapt_dist))+iTrackMinPump-2*self.gap_mask,0),
                                            (iTrackPump/2-iTrackMinPump-(iTrackPump/2-self.gap_mask), -adapt_dist_pump-self.gap_mask),
                                            (iGapPump+2*self.gap_mask+(iTrackPump/2-self.gap_mask), 0)]
        typePump='down'
        doublePump=False

        if typePump == 'up' or typePump == 'Up':
            raw_points_adapt_pump_a = self.refx_points(raw_points_adapt_pump_a)
            if self.is_mask:
                raw_points_adapt_pump_a_mask = self.refx_points(raw_points_adapt_pump_a_mask)
            ori_pump = [0,1]
            pos_pump = [((array_room+2*adapt_dist))/2, iTrack/2+iGap+iTrackMinPump+adapt_dist_pump]
        elif typePump =='down' or typePump == 'Down':
            ori_pump = [0,-1]
            pos_pump = [((array_room+2*adapt_dist))/2, -iTrack/2-iGap-iTrackMinPump-adapt_dist_pump]
        else:
            raise ValueError("typePump should be 'up' or 'down', given %s" % typePump)
        points_adapt_pump_a = self.append_points(raw_points_adapt_pump_a)
        cutout_pump_a=self.draw(self.name+"_cutout_pump_a", points_adapt_pump_a)
        gaps.append(cutout_pump_a)
        if self.is_mask:
            points_adapt_pump_a_mask = self.append_points(raw_points_adapt_pump_a_mask)
            masks.append(self.draw(self.name+"_cutout_pump_a_mask", points_adapt_pump_a_mask))
            
        raw_points_adapt_pump_b = self.refy_points(raw_points_adapt_pump_a, offset = ((array_room+2*adapt_dist))/2)
        points_adapt_pump_b = self.append_points(raw_points_adapt_pump_b)
        cutout_pump_b=self.draw(self.name+"_cutout_pump_b", points_adapt_pump_b)
        gaps.append(cutout_pump_b)
        if self.is_mask:
            raw_points_adapt_pump_b_mask = self.refy_points(raw_points_adapt_pump_a_mask, offset = ((array_room+2*adapt_dist))/2)
            points_adapt_pump_b_mask = self.append_points(raw_points_adapt_pump_b_mask)
            masks.append(self.draw(self.name+"_cutout_pump_b_mask", points_adapt_pump_b_mask))
        

        if doublePump:
            raw_points_adapt_pump_c = self.refx_points(raw_points_adapt_pump_a)
            points_adapt_pump_c = self.append_points(raw_points_adapt_pump_c)
            gaps.append(self.draw(self.name+"_cutout_pump_c", points_adapt_pump_c))
            if self.is_mask:
                raw_points_adapt_pump_c_mask = self.refx_points(raw_points_adapt_pump_a_mask)
                points_adapt_pump_c_mask = self.append_points(raw_points_adapt_pump_c_mask)
                masks.append(self.draw(self.name+"_cutout_pump_c_mask", points_adapt_pump_c_mask))

            raw_points_adapt_pump_d = self.refx_points(raw_points_adapt_pump_b)
            points_adapt_pump_d = self.append_points(raw_points_adapt_pump_d)
            gaps.append(self.draw(self.name+"_cutout_pump_d", points_adapt_pump_d))
            if self.is_mask:
                raw_points_adapt_pump_d_mask = self.refx_points(raw_points_adapt_pump_b_mask)
                points_adapt_pump_d_mask = self.append_points(raw_points_adapt_pump_d_mask)
                masks.append(self.draw(self.name+"_cutout_pump_d_mask", points_adapt_pump_d_mask))

        gaps = self.unite(gaps, self.name+'_cutout')
        self.gapObjects.append(gaps)
        if self.is_mask:
            masks = self.unite(masks, self.name+'_mask')
            self.maskObjects.append(masks)
            
        portOut1 = [self.coor([(3*(array_room+2*adapt_dist))/2,0]), self.ori, iTrack+2*self.overdev, iGap-2*self.overdev]
        portOut2 = [self.coor([-(3*(array_room+2*adapt_dist))/2,0]), -self.ori, iTrack+2*self.overdev, iGap-2*self.overdev]
        self.ports[self.name+'_1'] = portOut1
        self.ports[self.name+'_2'] = portOut2

        if doublePump:
            portOutpump1 = [self.coor(pos_pump), self.coor_vec(ori_pump), iTrackPump+2*self.overdev, iGapPump-2*self.overdev]
            self.ports[self.name+'_pump1'] = portOutpump1
            portOutpump2 = [self.coor(Vector(pos_pump).refx()), -self.coor_vec(ori_pump), iTrackPump+2*self.overdev, iGapPump-2*self.overdev]
            self.ports[self.name+'_pump2'] = portOutpump2
        else:
            portOutpump1 = [self.coor(pos_pump), self.coor_vec(ori_pump), iTrackPump+2*self.overdev, iGapPump-2*self.overdev]
            self.ports[self.name+'_pump'] = portOutpump1


    def draw_snails_Zaki(self, iTrack, iGap, array_room, array_offset, iTrackPump,
                    iGapPump, iTrackSnail=None, fillet=None, N_snails=1,
                    snail_dict={'loop_width':20e-6, 'loop_length':20e-6,
                                'length_big_junction':10e-6,
                                'length_small_junction':2e-6, 'bridge':1e-6,
                                'bridge_spacing':1e-6}, L_eq = '1nH'): #for now assume left and right tracks are the same width
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
        if iTrackSnail is None:
            iTrackSnail = iTrack/10

        iTrack, iGap, array_room, array_offset, iTrackPump, iGapPump, iTrackSnail = parse_entry((iTrack, iGap, array_room, array_offset, iTrackPump, iGapPump, iTrackSnail))

        adapt_dist = iTrack/2 #if slope ==1 !!!!

        #adapta
        raw_points_a = [(array_room/2+adapt_dist,-iTrack/2),
                        (0, iTrack),
                        (-adapt_dist, -iTrack/2+iTrackSnail/2+array_offset),
                        (0, -iTrackSnail)]
        
        points_a = self.append_points(raw_points_a)
        track_a = self.draw(self.name+"_track_a", points_a)
        
        raw_points_c = self.refy_points(raw_points_a)
        points_c = self.append_points(raw_points_c)
        track_c = self.draw(self.name+"_track_c", points_c)

        #snail array
        if L_eq is None:
            in_array = [self.coor([array_room/2, array_offset]), -self.ori, iTrackSnail, 0]
            out_array = [self.coor([-array_room/2, array_offset]), self.ori, iTrackSnail, 0]      
            snail_array = self.connect_elt(self.name+'_array', in_array, out_array)
            snail_track = snail_array._connect_snails_Zaki([snail_dict['loop_width'], snail_dict['loop_length']], snail_dict['length_big_junction'], 3, snail_dict['length_small_junction'], 1, N_snails, snail_dict['bridge'], snail_dict['bridge_spacing'])
        
        if L_eq is not None:
            array_eq = self.draw_rect_center(self.name+"_array_eq", self.coor([0,array_offset]), self.coor_vec([array_room, iTrackSnail]))
            self.assign_lumped_RLC(array_eq, self.ori, (0, L_eq, 0))
            points = self.append_points([(-array_room/2,array_offset),(array_room,0)])
            self.draw(self.name+'_array_eq_line', points, closed=False)
            
#        #junction up
#        print(self.ori)
#        in_junction_up = [self.coor([-iTrackSnail/2,squid_size[1]/2+iTrackSnail/2]), self.coor_vec([1,0]), iTrackSnail, 0]
#        out_junction_up = [self.coor([iTrackSnail/2,squid_size[1]/2+iTrackSnail/2]), self.coor_vec([-1,0]), iTrackSnail, 0]
#        junction = self.connect_elt(self.name+'_junction_up', in_junction_up, out_junction_up)
#        junction_pads_up = junction._connect_JJ(iTrackJ, iInduct=Lj_up, fillet=fillet)
#        
#        #junction down
#        in_junction_down = [self.coor([-iTrackSnail/2,-squid_size[1]/2-iTrackSnail/2]), self.coor_vec([1,0]), iTrackSnail, 0]
#        out_junction_down = [self.coor([iTrackSnail/2,-squid_size[1]/2-iTrackSnail/2]), self.coor_vec([-1,0]), iTrackSnail, 0]
#        junction = self.connect_elt(self.name+'_junction_down', in_junction_down, out_junction_down)
#        junction_pads_down = junction._connect_JJ(iTrackJ, iInduct=Lj_up, fillet=fillet)

        right_track = self.draw_rect_center(self.name+"_added_track1", self.coor([2*(array_room/2+adapt_dist),0]), self.coor_vec([array_room+2*adapt_dist, iTrack]))
        left_track = self.draw_rect_center(self.name+"_added_track2", self.coor([-2*(array_room/2+adapt_dist),0]), self.coor_vec([array_room+2*adapt_dist, iTrack]))

        squid = self.unite([right_track, left_track, track_a, track_c], name=self.name)
        self.trackObjects.append(squid)


        if fillet is not None:
            squid.fillet(iTrack/2,[0, 3, 7, 10])

        adapt_dist_pump = 4*iTrackPump#(4*iTrackPump - 2*iTrackSnail)/2/iSlopePump


        self.gapObjects.append(self.draw_rect_center(self.name+"_added_gap", self.coor([0,0]), self.coor_vec([3*(array_room+2*adapt_dist), iTrack+2*iGap])))
        if self.is_mask:
            self.maskObjects.append(self.draw_rect_center(self.name+"_added_gap", self.coor([0,0]), self.coor_vec([3*(array_room+2*adapt_dist), iTrack+3*iGap])))

        raw_points_adapt_pump_a = [(3/2*(array_room+2*adapt_dist),-iTrack/2-iGap-iTrackSnail),
                                         (-(array_room+2*adapt_dist)+iTrackSnail,0),
                                         (iTrackPump/2-iTrackSnail, -adapt_dist_pump),
                                         (iGapPump, 0)]
        if self.is_mask:
            raw_points_adapt_pump_mask = [(3/2*(array_room+2*adapt_dist)+iGapPump,-iTrack/2-iGap-iTrackSnail-iGapPump),
                                          (0,iTrackSnail+iGapPump),
                                             (-(array_room+2*adapt_dist)+iTrackSnail-iGapPump-iTrackPump/2,0),
                                             (iTrackPump/2-iTrackSnail, -adapt_dist_pump-iTrackSnail),
                                             (3*iGapPump+iTrackPump/2, 0)]
            self.maskObjects.append(self.draw(self.name+"_cutout_pump_a_mask", self.append_points(raw_points_adapt_pump_mask)))
            raw_points_adapt_pump_mask_b = self.refy_points(raw_points_adapt_pump_mask, offset = (array_room+2*adapt_dist)/2)
            self.maskObjects.append(self.draw(self.name+"_cutout_pump_b_mask", self.append_points(raw_points_adapt_pump_mask_b)))

        typePump='up'
        doublePump=False

        if typePump == 'up' or typePump == 'Up':
            raw_points_adapt_pump_a = self.refx_points(raw_points_adapt_pump_a)
            ori_pump = [0,1]
            pos_pump = [(array_room+2*adapt_dist)/2, iTrack/2+iGap+iTrackSnail+adapt_dist_pump]
        elif typePump =='down' or typePump == 'Down':
            ori_pump = [0,-1]
            pos_pump = [(array_room+2*adapt_dist)/2, -iTrack/2-iGap-iTrackSnail-adapt_dist_pump]
        else:
            raise ValueError("typePump should be 'up' or 'down', given %s" % typePump)
        points_adapt_pump_a = self.append_points(raw_points_adapt_pump_a)
        cutout_pump_a=self.draw(self.name+"_cutout_pump_a", points_adapt_pump_a)
        self.gapObjects.append(cutout_pump_a)
        
        raw_points_adapt_pump_b = self.refy_points(raw_points_adapt_pump_a, offset = (array_room+2*adapt_dist)/2)
        points_adapt_pump_b = self.append_points(raw_points_adapt_pump_b)
        cutout_pump_b=self.draw(self.name+"_cutout_pump_b", points_adapt_pump_b)
        self.gapObjects.append(cutout_pump_b)
        
        if fillet is not None and False:    
            cutout_pump_a.fillet(fillet/2,1)
            cutout_pump_a.fillet(fillet/4,0)
            cutout_pump_b.fillet(fillet/2,1)
            cutout_pump_b.fillet(fillet/4,0)
        if doublePump:
            raw_points_adapt_pump_c = self.refx_points(raw_points_adapt_pump_a)
            points_adapt_pump_c = self.append_points(raw_points_adapt_pump_c)
            self.gapObjects.append(self.draw(self.name+"_cutout_pump_c", points_adapt_pump_c))

            raw_points_adapt_pump_d = self.refx_points(raw_points_adapt_pump_b)
            points_adapt_pump_d = self.append_points(raw_points_adapt_pump_d)
            self.gapObjects.append(self.draw(self.name+"_cutout_pump_d", points_adapt_pump_d))

        portOut1 = [self.pos+self.ori*(array_room/2+adapt_dist)*3, self.ori, iTrack, iGap]
        portOut2 = [self.pos-self.ori*(array_room/2+adapt_dist)*3, -self.ori, iTrack, iGap]
        self.ports[self.name+'_1'] = portOut1
        self.ports[self.name+'_2'] = portOut2

        if doublePump:
            portOutpump1 = [self.coor(pos_pump), self.coor_vec(ori_pump), iTrackPump, iGapPump]
            self.ports[self.name+'_pump1'] = portOutpump1
            portOutpump2 = [self.coor(Vector(pos_pump).refx()), -self.coor_vec(ori_pump), iTrackPump, iGapPump]
            self.ports[self.name+'_pump2'] = portOutpump2
        else:
            portOutpump1 = [self.coor(pos_pump), self.coor_vec(ori_pump), iTrackPump, iGapPump]
            self.ports[self.name+'_pump'] = portOutpump1
			

    def draw_snails_sym(self, iTrack, iGap, array_room, iTrackPump, iGapPump, approach='20um', mode='litho', iTrackMinPump=None, iTrackSnail=None, fillet=None, N_snails=1, snail_dict={'loop_width':20e-6, 'loop_length':20e-6, 'length_big_junction':10e-6, 'length_small_junction':2e-6, 'bridge':1e-6, 'bridge_spacing':1e-6}, L_eq = '1nH'): #for now assume left and right tracks are the same width
        '''
        --------

        '''
        if iTrackSnail is None:
            iTrackSnail = iTrack/10
        if iTrackMinPump is None:
            iTrackMinPump=iTrackSnail

        iTrack, iGap, array_room, iTrackPump, iGapPump, iTrackMinPump, iTrackSnail, approach = parse_entry((iTrack, iGap, array_room, iTrackPump, iGapPump, iTrackMinPump, iTrackSnail, approach))

        adapt_dist = iTrack/2 #if slope ==1 !!!!
        
        track_a = self.draw_rect(self.name+"_track_a", self.coor([array_room/2,-iTrackSnail/2]), self.coor_vec([adapt_dist, iTrackSnail]))
        track_c = self.draw_rect(self.name+"_track_c", self.coor([-array_room/2,-iTrackSnail/2]), self.coor_vec([-adapt_dist, iTrackSnail]))
        

        #snail array
        if mode=='litho':
            in_array = [self.coor([array_room/2, 0]), -self.ori, iTrackSnail, 0]
            out_array = [self.coor([-array_room/2, 0]), self.ori, iTrackSnail, 0]      
            snail_array = self.connect_elt(self.name+'_array', out_array, in_array)
            snail_track = snail_array._connect_snails2([snail_dict['loop_width'], snail_dict['loop_length']], snail_dict['length_big_junction'], 4, snail_dict['length_small_junction'], 1, N_snails, snail_dict['bridge'], snail_dict['bridge_spacing'])
        
        if 0:
            in_array = [self.coor([array_room/2, 0]), -self.ori, iTrackSnail, 0]
            out_array = [self.coor([-array_room/2, 0]), self.ori, iTrackSnail, 0]      
            snail_array = self.connect_elt(self.name+'_array', in_array, out_array)
            snail_track = snail_array._connect_JJ(iTrack, iInduct=L_eq)
        
        if mode=='equivalent':
            connect_left = self.draw_rect(self.name+"_left", self.coor([array_room/2,-iTrackSnail/2]), self.coor_vec([-iTrackSnail, iTrackSnail]))
            connect_right = self.draw_rect(self.name+"_right", self.coor([-array_room/2,-iTrackSnail/2]), self.coor_vec([iTrackSnail, iTrackSnail]))
            if not self.is_litho:
                array_eq = self.draw_rect_center(self.name+"_array_eq", self.coor([0,0]), self.coor_vec([array_room-2*iTrackSnail, iTrackSnail]))
                self.assign_lumped_RLC(array_eq, self.ori, (0, L_eq, 0))
                points = self.append_points([(-array_room/2+iTrackSnail,0),(array_room-2*iTrackSnail,0)])
                self.draw(self.name+'_array_eq_line', points, closed=False)
            
#        #junction up
#        print(self.ori)
#        in_junction_up = [self.coor([-iTrackSnail/2,squid_size[1]/2+iTrackSnail/2]), self.coor_vec([1,0]), iTrackSnail, 0]
#        out_junction_up = [self.coor([iTrackSnail/2,squid_size[1]/2+iTrackSnail/2]), self.coor_vec([-1,0]), iTrackSnail, 0]
#        junction = self.connect_elt(self.name+'_junction_up', in_junction_up, out_junction_up)
#        junction_pads_up = junction._connect_JJ(iTrackJ, iInduct=Lj_up, fillet=fillet)
#        
#        #junction down
#        in_junction_down = [self.coor([-iTrackSnail/2,-squid_size[1]/2-iTrackSnail/2]), self.coor_vec([1,0]), iTrackSnail, 0]
#        out_junction_down = [self.coor([iTrackSnail/2,-squid_size[1]/2-iTrackSnail/2]), self.coor_vec([-1,0]), iTrackSnail, 0]
#        junction = self.connect_elt(self.name+'_junction_down', in_junction_down, out_junction_down)
#        junction_pads_down = junction._connect_JJ(iTrackJ, iInduct=Lj_up, fillet=fillet)

        right_track = self.draw_rect_center(self.name+"_added_track1", self.coor([2*(array_room/2+adapt_dist)+iTrackMinPump/2,0]), self.coor_vec([array_room+2*adapt_dist+iTrackMinPump, iTrack+2*self.overdev]))
        left_track = self.draw_rect_center(self.name+"_added_track2", self.coor([-2*(array_room/2+adapt_dist)-iTrackMinPump/2,0]), self.coor_vec([array_room+2*adapt_dist+iTrackMinPump, iTrack+2*self.overdev]))
        if mode=='equivalent':
            squid = self.unite([right_track, left_track, track_a, track_c, connect_left, connect_right], name=self.name+'_temp')
        else:
            squid = self.unite([right_track, left_track, track_a, track_c], name=self.name+'_temp')
            
        self.trackObjects.append(squid)

        adapt_dist_pump = 4*iTrackPump
        
        gaps=[]
        masks=[]
        
        raw_points_gap = [(-(3*(array_room+2*adapt_dist)+2*iTrackMinPump)/2,iTrack/2+iGap-self.overdev),
                          (iTrackMinPump-self.overdev, 0),
                          (0, -approach),
                          (array_room+2*adapt_dist+2*adapt_dist+array_room+2*self.overdev,0),
                          (0,approach),
                          (array_room+2*adapt_dist+iTrackMinPump-self.overdev,0),
                          (0, -(iTrack/2+iGap-self.overdev)),
                          (-(3*(array_room+2*adapt_dist)+2*iTrackMinPump),0)]
        points_gap = self.append_points(raw_points_gap)
        gap1 = self.draw(self.name+"_gap_1", points_gap)
        gap1.fillet(2*iTrackMinPump,[2,3])
        gaps.append(gap1)
        points_gap = self.append_points(self.refy_points(self.refx_points(raw_points_gap)))
        gap2 = self.draw(self.name+"_gap_2", points_gap)
        gap2.fillet(2*iTrackMinPump,[2,3])
        gaps.append(gap2)
        
#        gaps.append(self.draw_rect_center(self.name+"_added_gap", self.coor([0,0]), self.coor_vec([3*(array_room+2*adapt_dist)+2*iTrackMinPump, iTrack+2*iGap-2*self.overdev])))
        if self.is_mask:
            masks.append(self.draw_rect_center(self.name+"_added_gap", self.coor([0,0]), self.coor_vec([3*(array_room+2*adapt_dist), iTrack+2*iGap+2*self.gap_mask])))

        raw_points_adapt_pump_a = [(3/2*(array_room+2*adapt_dist)-self.overdev-iTrackMinPump,-iTrack/2-iGap-iTrackMinPump-self.overdev+approach),
                                   (-((array_room+2*adapt_dist))+2*iTrackMinPump+2*self.overdev,0),
                                   (iTrackPump/2-iTrackMinPump, -adapt_dist_pump+self.overdev),
                                   (iGapPump-2*self.overdev, 0)]
        if self.is_mask:
            raw_points_adapt_pump_a_mask = [(3/2*(array_room+2*adapt_dist)+self.gap_mask,-iTrack/2-iGap-iTrackMinPump+self.gap_mask),
                                            (-((array_room+2*adapt_dist))+iTrackMinPump-2*self.gap_mask,0),
                                            (iTrackPump/2-iTrackMinPump-(iTrackPump/2-self.gap_mask), -adapt_dist_pump-self.gap_mask),
                                            (iGapPump+2*self.gap_mask+(iTrackPump/2-self.gap_mask), 0)]



        

        points_adapt_pump_a = self.append_points(raw_points_adapt_pump_a)
        cutout_pump_a=self.draw(self.name+"_cutout_pump_a", points_adapt_pump_a)
        cutout_pump_a.fillet(iTrackMinPump,[0,1])
        gaps.append(cutout_pump_a)
        if self.is_mask:
            points_adapt_pump_a_mask = self.append_points(raw_points_adapt_pump_a_mask)
            masks.append(self.draw(self.name+"_cutout_pump_a_mask", points_adapt_pump_a_mask))
            
        raw_points_adapt_pump_b = self.refy_points(raw_points_adapt_pump_a, offset = ((array_room+2*adapt_dist))/2)
        points_adapt_pump_b = self.append_points(raw_points_adapt_pump_b)
        cutout_pump_b=self.draw(self.name+"_cutout_pump_b", points_adapt_pump_b)
        cutout_pump_b.fillet(iTrackMinPump,[0,1])
        gaps.append(cutout_pump_b)
        if self.is_mask:
            raw_points_adapt_pump_b_mask = self.refy_points(raw_points_adapt_pump_a_mask, offset = ((array_room+2*adapt_dist))/2)
            points_adapt_pump_b_mask = self.append_points(raw_points_adapt_pump_b_mask)
            masks.append(self.draw(self.name+"_cutout_pump_b_mask", points_adapt_pump_b_mask))
        
        ori_pump = [0,-1]
        pos_pump = [((array_room+2*adapt_dist))/2, -iTrack/2-iGap-iTrackMinPump-adapt_dist_pump+approach]
        
        raw_points_adapt_pump_c = self.refy_points(self.refx_points(raw_points_adapt_pump_a))
        points_adapt_pump_c = self.append_points(raw_points_adapt_pump_c)
        cutout_pump_c = self.draw(self.name+"_cutout_pump_c", points_adapt_pump_c)
        cutout_pump_c.fillet(iTrackMinPump,[0,1])
        gaps.append(cutout_pump_c)
        if self.is_mask:
            raw_points_adapt_pump_c_mask = self.refy_points(self.refx_points(raw_points_adapt_pump_a_mask))
            points_adapt_pump_c_mask = self.append_points(raw_points_adapt_pump_c_mask)
            masks.append(self.draw(self.name+"_cutout_pump_c_mask", points_adapt_pump_c_mask))

        raw_points_adapt_pump_d = self.refy_points(self.refx_points(raw_points_adapt_pump_b))
        points_adapt_pump_d = self.append_points(raw_points_adapt_pump_d)
        cutout_pump_d = self.draw(self.name+"_cutout_pump_d", points_adapt_pump_d)
        cutout_pump_d.fillet(iTrackMinPump,[0,1])
        gaps.append(cutout_pump_d)
        if self.is_mask:
            raw_points_adapt_pump_d_mask = self.refy_points(self.refx_points(raw_points_adapt_pump_b_mask))
            points_adapt_pump_d_mask = self.append_points(raw_points_adapt_pump_d_mask)
            masks.append(self.draw(self.name+"_cutout_pump_d_mask", points_adapt_pump_d_mask))
        
#        if fillet is not None:
#            cutout_pump_a.fillet(fillet/2,1)
#            cutout_pump_a.fillet(fillet/4,0)
#            cutout_pump_b.fillet(fillet/2,1)
#            cutout_pump_b.fillet(fillet/4,0)

        gaps = self.unite(gaps, self.name+'_cutout')
        self.gapObjects.append(gaps)
        if self.is_mask:
            masks = self.unite(masks, self.name+'_mask')
            self.maskObjects.append(masks)
            
        portOut1 = [self.coor([(3*(array_room+2*adapt_dist))/2+iTrackMinPump,0]), self.ori, iTrack+2*self.overdev, iGap-2*self.overdev]
        portOut2 = [self.coor([-(3*(array_room+2*adapt_dist))/2-iTrackMinPump,0]), -self.ori, iTrack+2*self.overdev, iGap-2*self.overdev]
        self.ports[self.name+'_1'] = portOut1
        self.ports[self.name+'_2'] = portOut2


        portOutpump1 = [self.coor(pos_pump), self.coor_vec(ori_pump), iTrackPump+2*self.overdev, iGapPump-2*self.overdev]
        self.ports[self.name+'_pump_1'] = portOutpump1
        portOutpump2 = [self.coor(Vector(pos_pump).refx().refy()), -self.coor_vec(ori_pump), iTrackPump+2*self.overdev, iGapPump-2*self.overdev]
        self.ports[self.name+'_pump_2'] = portOutpump2
#        else:
#            portOutpump1 = [self.coor(pos_pump), self.coor_vec(ori_pump), iTrackPump+2*self.overdev, iGapPump-2*self.overdev]
#            self.ports[self.name+'_pump'] = portOutpump1
            
    def draw_squid_protect(self, iTrack, iGap, squid_size, iTrackPump, iGapPump, iTrackSquid=None, iTrackJ=None, Lj_down='1nH', Lj_up=None,  typePump='down', doublePump=False, iSlope=1, iSlopePump=0.5, fillet=None): #for now assume left and right tracks are the same width
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
            
        if iTrackSquid is None:
            iTrackSquid = iTrack/4
        if iTrackJ is None:
            iTrackJ = iTrackSquid/2
        if Lj_up is None:
            Lj_up = Lj_down
        iTrack, iGap, squid_size, iTrackPump, iGapPump, iTrackSquid, iTrackJ = parse_entry((iTrack, iGap, squid_size, iTrackPump, iGapPump, iTrackSquid, iTrackJ))
        squid_size = Vector(squid_size)

        adapt_dist = squid_size[1]/2 #if slope ==1 !!!!

        #adapta
        raw_points_out = [(squid_size[0]/2+iTrackSquid,squid_size[1]/2+iTrackSquid),
                        (-(squid_size[0]+2*iTrackSquid), 0),
                        (0, -(squid_size[1]/2+iTrackSquid)+2*iTrackSquid),
                        (-squid_size[0]*2/3, 0),
                        (0, -4*iTrackSquid),
                        (squid_size[0]*2/3, 0),
                        (0, -(squid_size[1]/2+iTrackSquid)+2*iTrackSquid),
                        ((squid_size[0]+2*iTrackSquid), 0),]
        points_out = self.append_points(raw_points_out)
        track_out = self.draw(self.name+"_track_out", points_out)
        
        raw_points_in = [(squid_size[0]/2,squid_size[1]/2),
                        (-(squid_size[0]), 0),
                        (0, -(squid_size[1]/2)+iTrackSquid),
                        (-squid_size[0]/3, 0),
                        (0, iTrackSquid),
                        (-iTrackSquid,0),
                        (0,-iTrackSquid),
                        (-squid_size[0]/3+iTrackSquid,0),
                        (0, -2*iTrackSquid),
                        (squid_size[0]/3-iTrackSquid, 0),
                        (0, -iTrackSquid),
                        (iTrackSquid,0),
                        (0,iTrackSquid),
                        (squid_size[0]/3, 0),
                        (0, -(squid_size[1]/2)+iTrackSquid),
                        ((squid_size[0]), 0),]
        points_in = self.append_points(raw_points_in)
        track_in = self.draw(self.name+"_track_in", points_in)
        
        track_out=self.subtract(track_out, [track_in])

        #junction up
        print(self.ori)
        in_junction_up = [self.coor([-squid_size[0]/2-squid_size[0]/3-iTrackSquid,3/2*iTrackSquid]), self.coor_vec([1,0]), iTrackSquid, 0]
        out_junction_up = [self.coor([-squid_size[0]/2-squid_size[0]/3,3/2*iTrackSquid]), self.coor_vec([-1,0]), iTrackSquid, 0]
        junction = self.connect_elt(self.name+'_junction_up', in_junction_up, out_junction_up)
        junction_pads_up = junction._connect_JJ(iTrackJ, iInduct=Lj_up, fillet=fillet)
        #junction down
        in_junction_down = [self.coor([-squid_size[0]/2-squid_size[0]/3-iTrackSquid,-3/2*iTrackSquid]), self.coor_vec([1,0]), iTrackSquid, 0]
        out_junction_down = [self.coor([-squid_size[0]/2-squid_size[0]/3,-3/2*iTrackSquid]), self.coor_vec([-1,0]), iTrackSquid, 0]
        junction = self.connect_elt(self.name+'_junction_down', in_junction_down, out_junction_down)
        junction_pads_down = junction._connect_JJ(iTrackJ, iInduct=Lj_up, fillet=fillet)

        right_track = self.draw_rect(self.name+"_added_track1", self.coor([squid_size[0]/2+iTrackSquid,-iTrack/2-self.overdev]), self.coor_vec([squid_size[0]+iTrackSquid, iTrack+2*self.overdev]))
        left_track = self.draw_rect(self.name+"_added_track2", self.coor([-squid_size[0]/2-iTrackSquid-squid_size[0]*2/3,-iTrack/2-self.overdev]), self.coor_vec([-squid_size[0]-iTrackSquid+squid_size[0]*2/3, iTrack+2*self.overdev]))

        squid = self.unite([right_track, left_track, track_out, junction_pads_down, junction_pads_up], name=self.name)
        self.trackObjects.append(squid)


        adapt_dist_pump = 4*iTrackPump#(4*iTrackPump - 2*iTrackSquid)/2/iSlopePump


        gaps=[]
        masks=[]
        
        gaps.append(self.draw_rect_center(self.name+"_added_gap", self.coor([0,0]), self.coor_vec([3*squid_size[0]+4*iTrackSquid, iTrack+2*iGap-2*self.overdev])))
        if self.is_mask:
            masks.append(self.draw_rect_center(self.name+"_mask", self.coor([0,0]), self.coor_vec([3*squid_size[0]+4*iTrackSquid, iTrack+2*iGap+2*self.gap_mask])))

        raw_points_adapt_pump_a = [(3/2*squid_size[0]+2*iTrackSquid-self.overdev,-iTrack/2-iGap-iTrackSquid-self.overdev),
                                   (-(squid_size[0])+2*self.overdev,0),
                                   (iTrackPump/2-iTrackSquid, -adapt_dist_pump+self.overdev),
                                   (iGapPump-2*self.overdev, 0)]
        if self.is_mask:
            raw_points_adapt_pump_a_mask = [(3/2*squid_size[0]+2*iTrackSquid+self.gap_mask,-iTrack/2-iGap-iTrackSquid+self.gap_mask),
                                            (-(squid_size[0])-2*self.gap_mask,0),
                                            (iTrackPump/2-iTrackSquid-(iTrackPump/2-self.gap_mask), -adapt_dist_pump-self.gap_mask),
                                            (iGapPump+2*self.gap_mask+(iTrackPump/2-self.gap_mask), 0)]

        if typePump == 'up' or typePump == 'Up':
            raw_points_adapt_pump_a = self.refx_points(raw_points_adapt_pump_a)
            if self.is_mask:
                raw_points_adapt_pump_a_mask = self.refx_points(raw_points_adapt_pump_a_mask)
            ori_pump = [0,1]
            pos_pump = [(squid_size[0])/2+iTrackSquid, iTrack/2+iGap+iTrackSquid+adapt_dist_pump]
        elif typePump =='down' or typePump == 'Down':
            ori_pump = [0,-1]
            pos_pump = [(squid_size[0])/2+iTrackSquid, -iTrack/2-iGap-iTrackSquid-adapt_dist_pump]
        else:
            raise ValueError("typePump should be 'up' or 'down', given %s" % typePump)
        points_adapt_pump_a = self.append_points(raw_points_adapt_pump_a)
        cutout_pump_a=self.draw(self.name+"_cutout_pump_a", points_adapt_pump_a)
        gaps.append(cutout_pump_a)
        if self.is_mask:
            points_adapt_pump_a_mask = self.append_points(raw_points_adapt_pump_a_mask)
            masks.append(self.draw(self.name+"_cutout_pump_a_mask", points_adapt_pump_a_mask))
            
        raw_points_adapt_pump_b = self.refy_points(raw_points_adapt_pump_a, offset = (squid_size[0])/2+iTrackSquid)
        points_adapt_pump_b = self.append_points(raw_points_adapt_pump_b)
        cutout_pump_b=self.draw(self.name+"_cutout_pump_b", points_adapt_pump_b)
        gaps.append(cutout_pump_b)
        if self.is_mask:
            raw_points_adapt_pump_b_mask = self.refy_points(raw_points_adapt_pump_a_mask, offset = (squid_size[0])/2+iTrackSquid)
            points_adapt_pump_b_mask = self.append_points(raw_points_adapt_pump_b_mask)
            masks.append(self.draw(self.name+"_cutout_pump_b_mask", points_adapt_pump_b_mask))
        
        if fillet is not None:
            cutout_pump_a.fillet(fillet/2,1)
            cutout_pump_a.fillet(fillet/4,0)
            cutout_pump_b.fillet(fillet/2,1)
            cutout_pump_b.fillet(fillet/4,0)
        if doublePump:
            raw_points_adapt_pump_c = self.refx_points(raw_points_adapt_pump_a)
            points_adapt_pump_c = self.append_points(raw_points_adapt_pump_c)
            gaps.append(self.draw(self.name+"_cutout_pump_c", points_adapt_pump_c))
            if self.is_mask:
                raw_points_adapt_pump_c_mask = self.refx_points(raw_points_adapt_pump_a_mask)
                points_adapt_pump_c_mask = self.append_points(raw_points_adapt_pump_c_mask)
                masks.append(self.draw(self.name+"_cutout_pump_c_mask", points_adapt_pump_c_mask))

            raw_points_adapt_pump_d = self.refx_points(raw_points_adapt_pump_b)
            points_adapt_pump_d = self.append_points(raw_points_adapt_pump_d)
            gaps.append(self.draw(self.name+"_cutout_pump_d", points_adapt_pump_d))
            if self.is_mask:
                raw_points_adapt_pump_d_mask = self.refx_points(raw_points_adapt_pump_b_mask)
                points_adapt_pump_d_mask = self.append_points(raw_points_adapt_pump_d_mask)
                masks.append(self.draw(self.name+"_cutout_pump_d_mask", points_adapt_pump_d_mask))

        gaps = self.unite(gaps, self.name+'_cutout')
        self.gapObjects.append(gaps)
        if self.is_mask:
            masks = self.unite(masks, self.name+'_mask')
            self.maskObjects.append(masks)
            
        portOut1 = [self.coor([(3*squid_size[0]+4*iTrackSquid)/2,0]), self.ori, iTrack+2*self.overdev, iGap-2*self.overdev]
        portOut2 = [self.coor([-(3*squid_size[0]+4*iTrackSquid)/2,0]), -self.ori, iTrack+2*self.overdev, iGap-2*self.overdev]
        self.ports[self.name+'_1'] = portOut1
        self.ports[self.name+'_2'] = portOut2

        if doublePump:
            portOutpump1 = [self.coor(pos_pump), self.coor_vec(ori_pump), iTrackPump+2*self.overdev, iGapPump-2*self.overdev]
            self.ports[self.name+'_pump1'] = portOutpump1
            portOutpump2 = [self.coor(Vector(pos_pump).refx()), -self.coor_vec(ori_pump), iTrackPump+2*self.overdev, iGapPump-2*self.overdev]
            self.ports[self.name+'_pump2'] = portOutpump2
        else:
            portOutpump1 = [self.coor(pos_pump), self.coor_vec(ori_pump), iTrackPump, iGapPump]
            self.ports[self.name+'_pump'] = portOutpump1
            
    def draw_T(self, iTrack, iGap):
        
        if not self.is_overdev or self.val(self.overdev<0):
            cutout = self.draw_rect_center(self.name+'_cutout', self.coor([0,self.overdev/2]), self.coor_vec([2*iGap+iTrack, 2*iGap+iTrack-self.overdev]))
            self.gapObjects.append(cutout)
        else:
            points = self.append_points([(-(iGap+iTrack/2),-iTrack/2-iGap+self.overdev),
                             (0, 2*iGap+iTrack-2*self.overdev),
                             ((iGap+iTrack/2)*2, 0),
                             (0, -(2*iGap+iTrack)+2*self.overdev),
                             (-self.overdev, 0),
                             (0, -self.overdev), 
                             (-iTrack-2*iGap+2*self.overdev, 0),
                             (0, self.overdev)])
            cutout = self.draw(self.name+'_cutout', points)
            self.gapObjects.append(cutout)
        
        if self.is_mask:
            mask = self.draw_rect(self.name+'_mask', self.coor([-iGap-iTrack/2,-iGap-iTrack/2]), self.coor_vec([2*iGap+iTrack, 2*iGap+iTrack+self.gap_mask]))
            self.maskObjects.append(mask)
            
        points = self.append_points([(-(iGap+iTrack/2),-iTrack/2-self.overdev),
                                     (0, iTrack+2*self.overdev),
                                     ((iGap+iTrack/2)*2, 0),
                                     (0, -iTrack-2*self.overdev),
                                     (-iGap+self.overdev, 0),
                                     (0, -iGap+self.overdev), 
                                     (-iTrack-2*self.overdev, 0),
                                     (0, iGap-self.overdev)])
        track = self.draw(self.name+'_track', points)
        if self.val(iGap)<self.val(iTrack):
            fillet=iGap
        else:
            fillet=iTrack
        track.fillet(fillet-eps,[4,7])
        
        self.trackObjects.append(track)
        
        portOut1 = [self.coor([iTrack/2+iGap, 0]), self.coor_vec([1,0]), iTrack+2*self.overdev, iGap-2*self.overdev]
        self.ports[self.name+'_1'] = portOut1
        portOut2 = [self.coor([-(iTrack/2+iGap), 0]), self.coor_vec([-1,0]), iTrack+2*self.overdev, iGap-2*self.overdev]
        self.ports[self.name+'_2'] = portOut2
        portOut3 = [self.coor([0, -(iTrack/2+iGap)]), self.coor_vec([0,-1]), iTrack+2*self.overdev, iGap-2*self.overdev]
        self.ports[self.name+'_3'] = portOut3


    def draw_inductance(self, iTrack, iGap, all_length, ind_length, mode='litho', L_eq = '1nH'): #for now assume left and right tracks are the same width
        '''
        --------

        '''

        iTrack, iGap, all_length, ind_length = parse_entry((iTrack, iGap, all_length, ind_length))

        #snail array
#        if mode=='litho':
#            in_array = [self.coor([all_length/2, 0]), -self.ori, iTrackSnail, 0]
#            out_array = [self.coor([-all_length/2, 0]), self.ori, iTrackSnail, 0]      
#            snail_array = self.connect_elt(self.name+'_array', out_array, in_array)
#            snail_track = snail_array._connect_snails2([snail_dict['loop_width'], snail_dict['loop_length']], snail_dict['length_big_junction'], 4, snail_dict['length_small_junction'], 1, N_snails, snail_dict['bridge'], snail_dict['bridge_spacing'])
#        
#        if 0:
#            in_array = [self.coor([array_room/2, 0]), -self.ori, iTrackSnail, 0]
#            out_array = [self.coor([-array_room/2, 0]), self.ori, iTrackSnail, 0]      
#            snail_array = self.connect_elt(self.name+'_array', in_array, out_array)
#            snail_track = snail_array._connect_JJ(iTrack, iInduct=L_eq)
        
        if mode=='equivalent':
            connect_left = self.draw_rect(self.name+"_left", self.coor([ind_length/2,-iTrack/2]), self.coor_vec([all_length/2-ind_length/2, iTrack]))
            connect_right = self.draw_rect(self.name+"_right", self.coor([-ind_length/2,-iTrack/2]), self.coor_vec([-(all_length/2-ind_length/2), iTrack]))
            if not self.is_litho:
                array_eq = self.draw_rect_center(self.name+"_eq", self.coor([0,0]), self.coor_vec([ind_length, iTrack]))
                self.assign_lumped_RLC(array_eq, self.ori, (0, L_eq, 0))
                points = self.append_points([(-ind_length/2,0),(ind_length,0)])
                self.draw(self.name+'_eq_line', points, closed=False)
 

        connect = self.unite([connect_left, connect_right])
        self.trackObjects.append(connect)
        
        gap = self.draw_rect(self.name+"_cutout", self.coor([-all_length/2,-iTrack/2-iGap]), self.coor_vec([all_length, iTrack+2*iGap]))
        self.gapObjects.append(gap)
        if self.is_mask:
            mask = self.draw_rect(self.name+"_mask", self.coor([-all_length/2,-iTrack/2-iGap-self.gap_mask]), self.coor_vec([all_length, iTrack+2*iGap+2*self.gap_mask]))
            self.maskObjects.append(mask)
            
        portOut1 = [self.coor([all_length/2,0]), self.ori, iTrack+2*self.overdev, iGap-2*self.overdev]
        portOut2 = [self.coor([-all_length/2,0]), -self.ori, iTrack+2*self.overdev, iGap-2*self.overdev]
        self.ports[self.name+'_1'] = portOut1
        self.ports[self.name+'_2'] = portOut2
        
        
    def draw_halfJRM(self, iTrack, iGap, all_length, mode='litho', LS = '0.1nH',
                     Lind ='1nH', LJ='5nH'): #for now assume left and right tracks are the same width
        '''
        --------

        '''

        iTrack, iGap, all_length = parse_entry((iTrack, iGap, all_length))
        T = iTrack # Track
        A = all_length # legnth of halfJRM element top to bottom
        C = iGap # Gap. Width of halfJRM element is T+2*C
        L = C/3 # length of lumped inductor
        IJ = (C-L)/2 # Island of metal linking junctions to ground C = L+2*IJ
        CLR = (A-3*T)/4
        I = (A-2*CLR-3*T-2*L)/4 # Island of metal linking shunt junctions

        #snail array
#        if mode=='litho':
#            in_array = [self.coor([all_length/2, 0]), -self.ori, iTrackSnail, 0]
#            out_array = [self.coor([-all_length/2, 0]), self.ori, iTrackSnail, 0]      
#            snail_array = self.connect_elt(self.name+'_array', out_array, in_array)
#            snail_track = snail_array._connect_snails2([snail_dict['loop_width'], snail_dict['loop_length']], snail_dict['length_big_junction'], 4, snail_dict['length_small_junction'], 1, N_snails, snail_dict['bridge'], snail_dict['bridge_spacing'])
#        
#        if 0:
#            in_array = [self.coor([array_room/2, 0]), -self.ori, iTrackSnail, 0]
#            out_array = [self.coor([-array_room/2, 0]), self.ori, iTrackSnail, 0]      
#            snail_array = self.connect_elt(self.name+'_array', in_array, out_array)
#            snail_track = snail_array._connect_JJ(iTrack, iInduct=L_eq)

        raw_points_center = [(I+T/2, T/2),(-I-T/2, T/2),(-I-T/2, -T/2),(I+T/2, -T/2)]
        raw_points_right = [(3*I/2+L+T+I/2+T/2, T/2), (3*I/2+L+T-I/2-T/2, T/2),
                            (3*I/2+L+T-I/2-T/2, -T/2), (3*I/2+L+T+I/2+T/2, -T/2)]
        raw_points_left = self.refy_points(raw_points_right, absolute=True)
        raw_points_connector_top = [(T/2, T/2), (T/2, T/2+C), (-T/2, T/2+C),(-T/2, T/2)]
        raw_points_center1 =[(T/2, -T/2),(-T/2, -T/2),(-T/2, -T/2-IJ),(T/2, -T/2-IJ)]
        raw_points_center2 =self.move_points(raw_points_center1, [0, -IJ-L], absolute=True)
        raw_points_left1 =self.move_points(raw_points_center1, [-L-T-2*I, 0], absolute=True)
        raw_points_left2 =self.move_points(raw_points_center2, [-L-T-2*I, 0], absolute=True)
        raw_points_right1 =self.move_points(raw_points_center1, [L+T+2*I, 0], absolute=True)
        raw_points_right2 =self.move_points(raw_points_center2, [L+T+2*I, 0], absolute=True)

        x0 = -T/2-I-L-I-T
        raw_points_connector_left = [(x0, T/2),(x0-CLR, T/2),(x0-CLR, -T/2),(x0,-T/2)]
        raw_points_connector_right = self.refy_points(raw_points_connector_left, absolute=True)
        
        island_center = self.draw(self.name+"_center", self.append_absolute_points(raw_points_center))
        island_right = self.draw(self.name+"_right",self.append_absolute_points( raw_points_right))
        island_left = self.draw(self.name+"_left", self.append_absolute_points(raw_points_left))
        island_center1 = self.draw(self.name+"_center1", self.append_absolute_points(raw_points_center1))
        island_center2 = self.draw(self.name+"_center2", self.append_absolute_points(raw_points_center2))
        island_left1 = self.draw(self.name+"_left1", self.append_absolute_points(raw_points_left1))
        island_left2 = self.draw(self.name+"_left2", self.append_absolute_points(raw_points_left2))
        island_right1 = self.draw(self.name+"_right1", self.append_absolute_points(raw_points_right1))
        island_right2 = self.draw(self.name+"_right2", self.append_absolute_points(raw_points_right2))
        connector_top = self.draw(self.name+"_conn_top", self.append_absolute_points(raw_points_connector_top))
        connector_left = self.draw(self.name+"_conn_left", self.append_absolute_points(raw_points_connector_left))
        connector_right = self.draw(self.name+"_conn_right", self.append_absolute_points(raw_points_connector_right))
        
        if not self.is_litho:
            array_eqS1 = self.draw_rect(self.name+"_LS1", self.coor([T/2+I, -T/2]), self.coor_vec([L, T]))
            array_eqS2 = self.draw_rect(self.name+"_LS2", self.coor([-T/2-I-L, -T/2]), self.coor_vec([L, T]))
            array_eqL = self.draw_rect(self.name+"_L", self.coor([-T/2, -T/2-IJ]), self.coor_vec([T, -L]))
            array_eqJ1 = self.draw_rect(self.name+"_LJ1", self.coor([2*I+L+T-T/2, -T/2-IJ]), self.coor_vec([T, -L]))
            array_eqJ2 = self.draw_rect(self.name+"_LJ2", self.coor([-2*I-L-T-T/2, -T/2-IJ]), self.coor_vec([T, -L]))
            self.assign_lumped_RLC(array_eqS1, self.ori, (0, LS, 0))
            self.assign_lumped_RLC(array_eqS2, self.ori, (0, LS, 0))
            self.assign_lumped_RLC(array_eqL, self.ori, (0, Lind, 0))
            self.assign_lumped_RLC(array_eqJ1, self.ori, (0, LJ, 0))
            self.assign_lumped_RLC(array_eqJ2, self.ori, (0, LJ, 0))
            pointsS1 = [(I+T/2,0),(I+T/2+L,0)]
            self.draw(self.name+'_eq_lineS1', self.append_absolute_points(pointsS1), closed=False)
            pointsS2 = self.refy_points(pointsS1, absolute=True)
            self.draw(self.name+'_eq_lineS2', self.append_absolute_points(pointsS2), closed=False)
            pointsL = [(0,-T/2-IJ),(0,-T/2-IJ-L)]
            self.draw(self.name+'_eq_lineL', self.append_absolute_points(pointsL), closed=False)
            pointsJ1 = self.move_points(pointsL, [T+2*I+L, 0], absolute=True)
            self.draw(self.name+'_eq_lineJ1', self.append_absolute_points(pointsJ1), closed=False)
            pointsJ2 = self.move_points(pointsL, [-T-2*I-L, 0], absolute=True)
            self.draw(self.name+'_eq_lineJ2', self.append_absolute_points(pointsJ2), closed=False)
 

        connect = self.unite([island_center, island_right, island_left,
                              connector_top, island_center1, island_center2,
                              island_left1,island_left2,
                              island_right1,island_right2,
                              connector_left, connector_right])
        self.trackObjects.append(connect)
        

        gap = self.draw_rect(self.name+"_cutout", self.coor([-A/2,-T/2-C]), self.coor_vec([A, T+2*C]))
        self.gapObjects.append(gap)
#        if self.is_mask:
#            mask = self.draw_rect(self.name+"_mask", self.coor([-all_length/2,-iTrack/2-iGap-self.gap_mask]), self.coor_vec([all_length, iTrack+2*iGap+2*self.gap_mask]))
#            self.maskObjects.append(mask)
#            
        portT = [self.coor([0, T/2+C]), self.coor_vec([0,1]), iTrack+2*self.overdev, iGap-2*self.overdev]
        portL = [self.coor([-A/2,0]), self.coor_vec([-1,0]), iTrack+2*self.overdev, iGap-2*self.overdev]
        portR = [self.coor([A/2,0]), self.coor_vec([1,0]), iTrack+2*self.overdev, iGap-2*self.overdev]
        portB = [self.coor([0,-T/2-IJ-L-IJ]), self.coor_vec([0,-1]), iTrack+2*self.overdev, iGap-2*self.overdev]
        self.ports[self.name+'_T'] = portT
        self.ports[self.name+'_L'] = portL
        self.ports[self.name+'_R'] = portR
        self.ports[self.name+'_B'] = portB


    def draw_T_join(self, iTrack, iGap):
        
        iTrack, iGap = parse_entry((iTrack, iGap))
        
        if not self.is_overdev or self.val(self.overdev<0):
            cutout = self.draw_rect_center(self.name+'_cutout', self.coor([0,self.overdev/2]), self.coor_vec([4*iGap+2*iTrack, 2*iGap+iTrack-self.overdev]))
            self.gapObjects.append(cutout)
        else:
            points = self.append_points([(-(iGap*2+iTrack),-iTrack/2-iGap+self.overdev),
                             (0, 2*iGap+iTrack-2*self.overdev),
                             ((iGap*2+iTrack)*2, 0),
                             (0, -(2*iGap+iTrack)+2*self.overdev),
                             (-self.overdev, 0),
                             (0, -self.overdev), 
                             (-iTrack*2-4*iGap+2*self.overdev, 0),
                             (0, self.overdev)])
            cutout = self.draw(self.name+'_cutout', points)
            self.gapObjects.append(cutout)
        
        if self.is_mask:
            mask = self.draw_rect(self.name+'_mask', self.coor([-iGap*2-iTrack,-iGap-iTrack/2]), self.coor_vec([4*iGap+2*iTrack, 2*iGap+iTrack+self.gap_mask]))
            self.maskObjects.append(mask)
            
        points = self.append_points([(-(iGap+iTrack/2)*2,-iTrack/2-self.overdev),
                                     (0, iTrack+2*self.overdev),
                                     ((iGap+iTrack/2)*4, 0),
                                     (0, -iTrack-2*self.overdev),
                                     (-iGap*2+self.overdev, 0),
                                     (0, -iGap+self.overdev), 
                                     (-iTrack*2-2*self.overdev, 0),
                                     (0, iGap-self.overdev)])
        track = self.draw(self.name+'_track', points)
        if self.val(iGap)<self.val(iTrack):
            fillet=iGap
        else:
            fillet=iTrack
        track.fillet(fillet-eps,[4,7])
        
        self.trackObjects.append(track)
        
        portOut1 = [self.coor([iTrack+iGap*2, 0]), self.coor_vec([1,0]), iTrack+2*self.overdev, iGap-2*self.overdev]
        self.ports[self.name+'_1'] = portOut1
        portOut2 = [self.coor([-(iTrack+iGap*2), 0]), self.coor_vec([-1,0]), iTrack+2*self.overdev, iGap-2*self.overdev]
        self.ports[self.name+'_2'] = portOut2
        portOut3 = [self.coor([0, -(iTrack/2+iGap)]), self.coor_vec([0,-1]), iTrack*2+2*self.overdev, iGap*2-2*self.overdev]
        self.ports[self.name+'_3'] = portOut3

        
    def draw_fluxline(self, iTrack, iGap, length, track_flux, slope=0.5, sym='center', return_spacing=0, return_depth=0, opposite=False):
        is_fillet=True
        iTrack, iGap, length, track_flux, return_spacing, return_depth = parse_entry((iTrack, iGap, length, track_flux, return_spacing, return_depth))
        if sym == 'center':
            offset_gap_down = 0
            offset_track_down = 0
            offset_gap_up = 0
            offset_track_up = 0
            adapt_length = (iTrack/2-track_flux)/slope
        else:
            adapt_length = (iTrack/2-track_flux+length/2)/slope
            if sym=='down':
                offset_gap_down = iGap+2*track_flux
                offset_track_down = length/2+track_flux
                offset_gap_up = 0
                offset_track_up = 0
            elif sym=='up':
                offset_gap_down = 0
                offset_track_down = 0
                offset_gap_up = iGap+2*track_flux
                offset_track_up = length/2+track_flux
            else:
                raise ValueError("sym should be either True, 'down' or 'up'.")
        
        points = self.append_absolute_points([(adapt_length, iGap+iTrack/2),
                                              (adapt_length-track_flux, iGap+iTrack/2),
                                              (track_flux/2, length/2+offset_gap_up),
                                              (track_flux/2, -length/2-offset_gap_down),
                                              (adapt_length-track_flux, -(iGap+iTrack/2)),
                                              (adapt_length, -(iGap+iTrack/2))])
                

        if not sym == 'center':
            if return_spacing != 0:
                if not opposite:
                    test1 = self.draw_rect(self.name+'_test1', self.coor([return_depth, length/2+offset_gap_up]), self.coor_vec([-(return_spacing+track_flux/2+return_depth), return_spacing+track_flux]))
                else:
                    test1 = self.draw_rect(self.name+'_test1', self.coor([return_depth, -(length/2+offset_gap_up)]), self.coor_vec([-(return_spacing+track_flux/2+return_depth), -(return_spacing+track_flux)]))
                test1.fillet(return_spacing+track_flux+eps, 2)
                self.gapObjects.append(test1)                    
                    
            
            
        gap = self.draw(self.name+'_gap', points)
        
        
        track = self.draw(self.name+'_track_guide', points, closed=False)
        gapext = self.draw(self.name+'_gapext_guide', points, closed=False)
        fillet = eps
        fillet2 = track_flux+eps
        if is_fillet:
            track.fillet(fillet2, [1,4])
            gapext.fillet(fillet2, [1,4])
            gap.fillet(fillet2,[1,4])
            
            track.fillet(fillet, [3,4])
            gapext.fillet(fillet, [3,4])
            gap.fillet(fillet,[3,4])

        points_starter = self.append_absolute_points([(adapt_length, iGap+iTrack/2),
                                              (adapt_length, iGap+iTrack/2+track_flux)])
        track_starter = self.draw(self.name+'_track', points_starter, closed=False)
        gapext_starter = self.draw(self.name+'_gapext', points_starter, closed=False)
        
        
        track = track.sweep_along_path(track_starter)
        gapext = gapext.sweep_along_path(gapext_starter)
        
        gap = self.unite([gap, gapext])
        
        if self.is_mask:
            mask = self.draw(self.name+'_mask', points)
            maskext = self.draw(self.name+'_mask_ext_guide', points, closed=False)
            points_starter_mask = self.append_absolute_points([(adapt_length, iGap+iTrack/2),
                                              (adapt_length, iGap+iTrack/2+self.gap_mask)])
            maskext_starter = self.draw(self.name+'_maskext', points_starter_mask, closed=False)
            maskext = maskext.sweep_along_path(maskext_starter)
            
            mask = self.unite([mask, maskext])
            self.maskObjects.append(mask)

        points = self.append_absolute_points([(adapt_length, iTrack/2),
                                              (adapt_length-track_flux, iTrack/2),
                                              (track_flux/2, track_flux-offset_track_down+offset_track_up),
                                              (track_flux/2, -track_flux-offset_track_down+offset_track_up),
                                              (adapt_length-track_flux, -iTrack/2),
                                              (adapt_length, -iTrack/2)])
                  
        track_adapt = self.draw(self.name+'_track_adapt', points)
        if is_fillet:
            track_adapt.fillet(fillet2, [1,4])
        
        track = self.unite([track, track_adapt])
#            track.fillet(fillet, [15,20])
#            track.fillet(fillet2, [17,20])
        
        self.gapObjects.append(gap)
        self.trackObjects.append(track)

        portOut = [self.coor([adapt_length, 0]), self.coor_vec([1,0]), iTrack+2*self.overdev, iGap-2*self.overdev]
        self.ports[self.name] = portOut
        

    def draw_double_fluxline(self, iTrack, iGap, length, track_flux, slope=0.7, sym='center', return_spacing=0, return_depth=0, opposite=False):
        is_fillet = True
        iTrack, iGap, length, track_flux, return_spacing, return_depth = parse_entry((iTrack, iGap, length, track_flux, return_spacing, return_depth))

        offset_gap_down = 0
        offset_track_down = 0
        offset_gap_up = 0
        offset_track_up = 0
        adapt_length = (iTrack/2-track_flux/2)/slope
        small_gap = length/4-3*track_flux/2/2
        big_length = 3*iTrack + 4*iGap
        
        points = self.append_absolute_points([(adapt_length, big_length/2),
                                              (adapt_length-track_flux,big_length/2),
                                              (track_flux/2, length/2+offset_gap_up),
                                              (track_flux/2, -length/2-offset_gap_down),
                                              (adapt_length-track_flux, -big_length/2),
                                              (adapt_length, -big_length/2)])
                      
                    
        gap = self.draw(self.name+'_gap', points)
        
        track = self.draw(self.name+'_track_guide', points, closed=False)
        
        gapext = self.draw(self.name+'_gapext_guide', points, closed=False)


        points_starter = self.append_absolute_points([(adapt_length, big_length/2),
                                              (adapt_length, big_length/2+track_flux)])
        track_starter = self.draw(self.name+'_track', points_starter, closed=False)
        gapext_starter = self.draw(self.name+'_gapext', points_starter, closed=False)
        
        
        track = track.sweep_along_path(track_starter)
        gapext = gapext.sweep_along_path(gapext_starter)
        
        gap = self.unite([gap, gapext])
        
        if self.is_mask:
            mask = self.draw(self.name+'_mask', points)
            maskext = self.draw(self.name+'_mask_ext_guide', points, closed=False)
            points_starter_mask = self.append_absolute_points([(adapt_length, iGap+iTrack/2),
                                              (adapt_length, iGap+iTrack/2+self.gap_mask)])
            maskext_starter = self.draw(self.name+'_maskext', points_starter_mask, closed=False)
            maskext = maskext.sweep_along_path(maskext_starter)
            
            mask = self.unite([mask, maskext])
            self.maskObjects.append(mask)
            
        raw_points_top = [(adapt_length,big_length/2-iGap),
                          (adapt_length-track_flux, big_length/2-iGap),
                          (track_flux/2, small_gap+3*track_flux/2),
                          (track_flux/2, small_gap+3*track_flux/2-2/2*track_flux),
                          (adapt_length-track_flux, big_length/2-iGap-iTrack),
                          (adapt_length, big_length/2-iGap-iTrack)]
            
        raw_points_middle = [(adapt_length, iTrack/2),
                          (adapt_length-track_flux, iTrack/2),
                          (track_flux/2, track_flux/2-offset_track_down+offset_track_up),
                          (track_flux/2, -track_flux/2-offset_track_down+offset_track_up),
                          (adapt_length-track_flux, -iTrack/2),
                          (adapt_length, -iTrack/2)]
        
        raw_points_bottom = self.refx_points(raw_points_top)

        points_top = self.append_absolute_points(raw_points_top)
        points_middle = self.append_absolute_points(raw_points_middle)
        points_bottom = self.append_absolute_points(raw_points_bottom)
                  
        track_adapt_top = self.draw(self.name+'_track_adapt_top', points_top)
        track_adapt_middle = self.draw(self.name+'_track_adapt_middle',
                                       points_middle)
        track_adapt_bottom = self.draw(self.name+'_track_adapt_bottom',
                                       points_bottom)
        

        
        track = self.unite([track, track_adapt_top,
                            track_adapt_middle,
                            track_adapt_bottom])
#            track.fillet(fillet, [15,20])
#            track.fillet(fillet2, [17,20])

        fillet2 = 2*track_flux/2+eps
        if is_fillet:
            track.fillet(fillet2, [1, 2, 3, 4, 7,8,9,10,13,14,15,16,19,20,21,22,25,26,27,28])
            gap.fillet(fillet2, [2,3,4, 5])
        
        self.gapObjects.append(gap)
        self.trackObjects.append(track)

        portOut_top = [self.coor([adapt_length, iTrack/2+iGap+iTrack/2]), self.coor_vec([1,0]), iTrack+2*self.overdev, iGap-2*self.overdev]
        portOut_bottom = [self.coor([adapt_length, -iTrack/2-iGap-iTrack/2]), self.coor_vec([1,0]), iTrack+2*self.overdev, iGap-2*self.overdev]
        self.ports[self.name+'_top'] = portOut_top
        self.ports[self.name+'_bottom'] = portOut_bottom
 
        
    def draw_end_cable(self, iTrack, iGap, typeEnd = 'open', fillet=None):
        iTrack, iGap = parse_entry((iTrack, iGap))
        if typeEnd=='open' or typeEnd=='Open':
            cutout = self.draw_rect(self.name+'_cutout', self.coor([iGap,-(iTrack+2*iGap)/2+self.overdev]), self.coor_vec([-iGap+self.overdev, iTrack+2*iGap-2*self.overdev]))
            if fillet is not None:
                if abs(self.ori[0])==1:
                    cutout.fillet(iGap-self.overdev-eps,[2,1])
                else:
                    cutout.fillet(iGap-self.overdev-eps,[2,3])
            self.gapObjects.append(cutout)
            portOut = [self.coor([iGap, 0]), self.ori, iTrack+2*self.overdev, iGap-2*self.overdev]
            
            if self.is_overdev:
                track = self.draw_rect(self.name+'_track', self.coor([iGap,-iTrack/2-self.overdev]), self.coor_vec([-self.overdev, iTrack+2*self.overdev]))
                self.trackObjects.append(track)
            if self.is_mask:
                mask = self.draw_rect(self.name+'_mask', self.coor([-self.gap_mask,-(iTrack+2*iGap)/2-self.gap_mask]), self.coor_vec([iGap+self.gap_mask, iTrack+2*iGap+2*self.gap_mask]))
                if fillet is not None:
                    if abs(self.ori[0])==1:
                        mask.fillet(iGap+self.gap_mask-eps,[0,3])
                    else:
                        mask.fillet(iGap+self.gap_mask-eps,[0,1])
                self.maskObjects.append(mask)
                
        elif typeEnd=='short' or typeEnd=='Short':
            cutout1 = self.draw_rect(self.name+'_cutout1', self.coor([iGap/2,-(iTrack/2+iGap)+self.overdev]), self.coor_vec([-iGap/2+self.overdev, iGap-2*self.overdev]))
            cutout2 = self.draw_rect(self.name+'_cutout2', self.coor([iGap/2,(iTrack/2+iGap)-self.overdev]), self.coor_vec([-iGap/2+self.overdev, -iGap+2*self.overdev]))
            if fillet is not None:
                if abs(self.ori[0])==1:
                    cutout1.fillet(iGap/2-self.overdev-eps,[2,1])
                    cutout2.fillet(iGap/2-self.overdev-eps,[2,1])
                else:
                    cutout1.fillet(iGap/2-self.overdev-eps,[2,3])
                    cutout2.fillet(iGap/2-self.overdev-eps,[2,3])
            cutout = self.unite([cutout1, cutout2], self.name+'_cutout')
            self.gapObjects.append(cutout)
            portOut = [self.coor([iGap/2, 0]), self.ori, iTrack+2*self.overdev, iGap-2*self.overdev]
            
            if self.is_mask:
                mask = self.draw_rect(self.name+'_mask', self.coor([-self.gap_mask,-(iTrack+2*iGap)/2-self.gap_mask]), self.coor_vec([iGap/2+self.gap_mask, iTrack+2*iGap+2*self.gap_mask]))
                if fillet is not None:
                    if abs(self.ori[0])==1:
                        mask.fillet(iGap/2+self.gap_mask-eps,[0,3])
                    else:
                        mask.fillet(iGap/2+self.gap_mask-eps,[0,1])
                self.maskObjects.append(mask)
        else:
            raise ValueError("typeEnd should be 'open' or 'short', given %s" % typeEnd)
        self.ports[self.name] = portOut


    def draw_alignement_mark(self, iSize, iXdist, iYdist):
        iXdist, iYdist, iSize=parse_entry((iXdist, iYdist, iSize))
        raw_points = [(iXdist,iYdist),
                      (iSize/2,iSize/2),
                      (-iSize,0)]
        mark1=self.draw(self.name+"_mark_a", self.append_points(raw_points))
        
        raw_points = [(iXdist,iYdist),
                              (-iSize/2,-iSize/2),
                              (iSize,0)]
        mark2=self.draw(self.name+"_mark_b", self.append_points(raw_points))
        self.gapObjects.append(self.unite([mark1,mark2]))
        
    def draw_alignement_mark_r(self, size, disp, suff=''):
        size, disp = parse_entry((size, disp))
        raw_points = [(*disp,),
                      (size/2,size/2),
                      (0,-size)]        
        marks = []
        marks.append(self.draw(self.name+'_'+suff+"a", self.append_points(raw_points)))
        marks.append(self.draw(self.name+'_'+suff+"b", self.append_points(self.refy_points(raw_points, disp[0]))))
        self.gapObjects.append(self.unite([*marks]))
        if self.is_mask:
            raw_points_mask = [(disp[0]-self.gap_mask,disp[1]),
                               (size/2+2*self.gap_mask,size/2+self.gap_mask*2),
                               (0,-size-self.gap_mask*4)]
            marks_mask = []
            marks_mask.append(self.draw(self.name+suff+"a_mask", self.append_points(raw_points_mask)))
            marks_mask.append(self.draw(self.name+suff+"b_mask", self.append_points(self.refy_points(raw_points_mask, disp[0]))))
            self.maskObjects.append(self.unite(marks_mask))
            
    def draw_alignement_marks(self, size, disp, dict_except=None):
        size, disp = parse_entry((size, disp))
        disp = Vector(disp)
        moves = []
        directions_name = ['NE', 'NO', 'SE', 'SO']        
        directions_exp = [[1,1], [-1,1], [1,-1], [-1,-1]]
        if dict_except is not None:
            for direction_name in directions_name:
                try:
                    moves.append(parse_entry(dict_except[direction_name]))
                except Exception:
                    moves.append([0,0])
        else:
            moves = [[0,0] for ii in range(4)]
        
        for move, direction_exp, direction_name in zip(moves, directions_exp, directions_name):
            self.draw_alignement_mark_r(size, disp*direction_exp+move, suff=direction_name)      
        
    def draw_dose_test_Nb(self, pad_size, pad_spacing, array):
        pad_size, pad_spacing = parse_entry((pad_size, pad_spacing))
        pad_size = Vector(pad_size)
        
        pads = []
        pos_x = pad_spacing
        pos_y = pad_spacing
        for jj in range(array[1]):
            for ii in range(array[0]):
                pads.append(self.draw_rect(self.name+'_%d_%d'%(ii, jj), self.coor([pos_x, pos_y]), self.coor_vec(pad_size)))
                pos_x += pad_size[0]+pad_spacing
                pads.append(self.draw_rect(self.name+'_%d_%d'%(ii, jj), self.coor([pos_x, pos_y]), self.coor_vec(pad_size)))
                pos_x += pad_size[0]+2*pad_spacing
            pos_y += pad_size[1]+pad_spacing
            pos_x = pad_spacing
                
        cutout = self.draw_rect(self.name+'_cutout', self.coor([0, 0]), self.coor_vec([pad_size[0]*2*array[0]+pad_spacing*3*array[0],pad_size[1]*array[1]+pad_spacing*(array[1]+1)]))
        if self.is_mask:
            mask = self.draw_rect(self.name+'_cutout', self.coor([-self.gap_mask, -self.gap_mask]), self.coor_vec([pad_size[0]*2*array[0]+pad_spacing*3*array[0]+2*self.gap_mask,pad_size[1]*array[1]+pad_spacing*(array[1]+1)+2*self.gap_mask]))
            self.maskObjects.append(mask)
        cutout = self.subtract(cutout, pads)
        self.gapObjects.append(cutout)

    def draw_dose_test(self, pad_size, pad_spacing, iTrack, bridge, N, bridge_spacing, length_big_junction, length_small_junction):
        pad_size, pad_spacing, iTrack, bridge, bridge_spacing, length_big_junction, length_small_junction = parse_entry((pad_size, pad_spacing, iTrack, bridge, bridge_spacing, length_big_junction, length_small_junction))
        pad_size = Vector(pad_size)
        
        self.draw_rect(self.name+'_left', self.coor([-pad_spacing/2, -pad_size[1]/2]), self.coor_vec([-pad_size[0],pad_size[1]]))
        self.draw_rect(self.name+'_right', self.coor([pad_spacing/2, pad_size[1]/2]), self.coor_vec([pad_size[0],-pad_size[1]]))
        
        portOut1 = [self.coor([pad_spacing/2, 0]), -self.ori, iTrack, 0]
        portOut2 = [self.coor([-pad_spacing/2, 0]), self.ori, iTrack, 0]
        self.ports[self.name+'_1'] = portOut1
        self.ports[self.name+'_2'] = portOut2
        
        in_array = portOut2
        out_array = portOut1
        snail_array = self.connect_elt(self.name+'_junction', in_array, out_array)
        snail_track = snail_array._connect_snails2([20e-6,20e-6], length_big_junction, 3, length_small_junction, 1, N, bridge, bridge_spacing)#(squid_size, width_top, n_top, width_bot, n_bot, N, width_bridge)
    
    def draw_dose_test_snail_Zaki(self, pad_size, pad_spacing, iTrack, N_snails=1,
                                  snail_dict={'loop_width':20e-6, 'loop_length':20e-6,
                                                'length_big_junction':10e-6,
                                                'length_small_junction':2e-6, 'bridge':1e-6,
                                                'bridge_spacing':1e-6}):
        pad_size, pad_spacing, iTrack = parse_entry((pad_size, pad_spacing, iTrack))
        pad_size = Vector(pad_size)
        
        self.draw_rect(self.name+'_left', self.coor([-pad_spacing/2, -pad_size[1]/2]), self.coor_vec([-pad_size[0],pad_size[1]]))
        self.draw_rect(self.name+'_right', self.coor([pad_spacing/2, pad_size[1]/2]), self.coor_vec([pad_size[0],-pad_size[1]]))
        
        portOut1 = [self.coor([pad_spacing/2, 0]), -self.ori, iTrack, 0]
        portOut2 = [self.coor([-pad_spacing/2, 0]), self.ori, iTrack, 0]
        self.ports[self.name+'_1'] = portOut1
        self.ports[self.name+'_2'] = portOut2
        
        in_array = portOut2
        out_array = portOut1
        snail_array = self.connect_elt(self.name+'_junction', in_array, out_array)
        snail_track = snail_array._connect_snails_Zaki([snail_dict['loop_width'], snail_dict['loop_length']], snail_dict['length_big_junction'], 3, snail_dict['length_small_junction'], 1, N_snails, snail_dict['bridge'], snail_dict['bridge_spacing'])#(squid_size, width_top, n_top, width_bot, n_bot, N, width_bridge)
        
    def draw_dose_test_junction(self, pad_size, pad_spacing, width, width_bridge, n_bridge=1, spacing_bridge=0, alternate_width=True):
        pad_size, pad_spacing,width, spacing_bridge, width_bridge = parse_entry((pad_size, pad_spacing, width, spacing_bridge, width_bridge))
        pad_size = Vector(pad_size)
        if self.val(width)<1.5e-6 and n_bridge==1:
            width_jct = width
            width = 1.5e-6
        else: 
            width_jct=None
        
        self.draw_rect(self.name+'_left', self.coor([-pad_spacing/2, -pad_size[1]/2]), self.coor_vec([-pad_size[0],pad_size[1]]))
        self.draw_rect(self.name+'_right', self.coor([pad_spacing/2, pad_size[1]/2]), self.coor_vec([pad_size[0],-pad_size[1]]))
        
        portOut1 = [self.coor([pad_spacing/2, 0]), -self.ori, width, 0]
        portOut2 = [self.coor([-pad_spacing/2, 0]), self.ori, width, 0]
        self.ports[self.name+'_1'] = portOut1
        self.ports[self.name+'_2'] = portOut2
        
        in_array = portOut2
        out_array = portOut1
        jcts = self.connect_elt(self.name+'_junction', in_array, out_array)
        print(n_bridge)
        if alternate_width:
            jcts._connect_jct(width_bridge, n=n_bridge, spacing_bridge=spacing_bridge, width_jct=width_jct)
        else:
            jcts._connect_jct(width_bridge, n=n_bridge, spacing_bridge=spacing_bridge, assymetry=0, width_jct=width_jct)
    
    def draw_cav_Si(self, size, chip_height, tunnel_length, size_plot=None, height_plot=None, n=8, track='42um', gap='25um', ports_senw = ['500um',0,0,0], margin_litho='10um'):
        # can only make square cavity
        # size : intra size of the cavity : 4mm->15GHz
        # chip_heigth : silicon height
        # tunnel length : length of the silicon tunnel between 2 via
        # size_plot : base size of the central plot
        # height_plot : thickness of the capacitance at the top of the plot
        # n : nb of vias
        # track, gap : if we have an access port, track and gap of this port 
        #   currently all port have same gap and track
        # ports_senw = list giving the length of the cable coming from each port ; south, east, north, west
        #   0 means there should be no port
        # ports will be numeroted from 1 to n in the trigo way
        # margin_litho : width of the Si3N4 mask between the vias
        
        ports_senw_bool = [bool(ii) for ii in ports_senw]
        size, chip_height, tunnel_length, track, gap, margin_litho, ports_senw = parse_entry((size, chip_height, tunnel_length, track, gap, margin_litho, ports_senw))
#        self.draw_rect(self.name+'_rect', self.coor([-size/2, -size/2]), self.coor_vec([size,size]))
        access_width = track+2*gap+2*margin_litho
        if size_plot is None:
            size_plot=size/2
        else:
            size_plot = parse_entry((size_plot))
        
        if height_plot is None:
            height_plot=chip_height/2
        else:
            height_plot = parse_entry((height_plot))
            
        # bottom row
        width_trapeze = (size-(n+1)*margin_litho)/n
        track = '42um'
        gap = '25um'
        name_s, name_e, name_n, name_w = self.name+'_s', self.name+'_e', self.name+'_n', self.name+'_w'
        port_nb = 0
        for ii, name in enumerate([name_s, name_e, name_n, name_w]):
            if ports_senw_bool[ii]:
                port_nb += 1 
                _width_trapeze = width_trapeze-access_width/n
                x_pos = -size/2+margin_litho+_width_trapeze/2
                _n = n//2
                self.draw_trapeze(name, self.coor([x_pos, -size/2-tunnel_length/2]), 0, self.coor_vec([_width_trapeze, tunnel_length]), -chip_height/2)
                self.draw_trapeze(name+'_'+str('bis'), self.coor([x_pos, -size/2-tunnel_length/2]), -chip_height, self.coor_vec([_width_trapeze, tunnel_length]), chip_height/2)
                self.unite([name, name+'_'+str('bis')])

                self.duplicate_along_line(name, self.coor_vec([_width_trapeze+margin_litho, 0]), n=_n)
                unite_list = [name]+[name+'_'+str(ii+1) for ii in range(_n-1)]
                self.unite(unite_list)
                
                self.duplicate_along_line(name, self.coor_vec([_n*(_width_trapeze+margin_litho)+access_width, 0]), n=2)
                unite_list = [name, name+'_'+str(_n)]
                self.unite(unite_list)
                
                portOut = [self.coor([0, -size/2-tunnel_length]), self.coor_vec([0,-1]), track, gap]
                self.ports[self.name+'_port_'+str(port_nb)] = portOut
                
                end = self.key_elt(self.name+'_end_'+str(ii), self.coor([0, -size/2-tunnel_length+ports_senw[ii]]), self.coor_vec([0,-1]))
                end.draw_end_cable(track, gap, fillet=gap)
                cable = self.connect_elt(self.name+'_cable_'+str(ii), self.name+'_end_'+str(ii), self.name+'_port_'+str(port_nb))
                cable.draw_cable(is_bond=False)
            else:
                x_pos = -size/2+margin_litho+width_trapeze/2
        
                self.draw_trapeze(name, self.coor([x_pos, -size/2-tunnel_length/2]), 0, self.coor_vec([width_trapeze, tunnel_length]), -chip_height/2)
                self.draw_trapeze(name+'_'+str('bis'), self.coor([x_pos, -size/2-tunnel_length/2]), -chip_height, self.coor_vec([width_trapeze, tunnel_length]), chip_height/2)
                self.unite([name, name+'_'+str('bis')])
                
                self.duplicate_along_line(name, self.coor_vec([width_trapeze+margin_litho, 0]), n=n)
                unite_list = [name]+[name+'_'+str(ii+1) for ii in range(n-1)]
                self.unite(unite_list)
                
            self.ori = self.ori.rot(Vector([0,1]))

        name_c = self.name+'_c'
        self.draw_trapeze(name_c, self.coor([-size/2-tunnel_length/2, -size/2-tunnel_length/2]), 0, self.coor_vec([tunnel_length, tunnel_length]), -chip_height/2)
        self.draw_trapeze(name_c+'_bis', self.coor([-size/2-tunnel_length/2, -size/2-tunnel_length/2]), -chip_height, self.coor_vec([tunnel_length, tunnel_length]), chip_height/2)
        self.unite([name_c, name_c+'_bis'])
        
        self.duplicate_along_line(name_c, self.coor_vec([size+tunnel_length, 0]), n=2)
        unite_list = [name_c, name_c+'_1']
        self.unite(unite_list)
        
        self.duplicate_along_line(name_c, self.coor_vec([0, size+tunnel_length]), n=2)
        unite_list = [name_c, name_c+'_2']
        self.unite(unite_list)
        
        name_p = name+'_p'
        self.draw_trapeze(name_p, self.coor([0, 0]), -chip_height, self.coor_vec([size_plot, size_plot]), chip_height-height_plot)


        unite_list = [name_s, name_e, name_n, name_w, name_c, name_p]
        self.unite(unite_list)
        self.name = name_s
        
    def size_dc_gap(self, length, positions, widths, border):
        
        length, positions, widths, border = parse_entry((length, positions, widths, border))
        
        low = np.argmin(positions)
        high = np.argmax(positions)
        pos_cutout = self.coor([-length/2, positions[low]-widths[low]/2-border])
        width_cutout = positions[high]-positions[low]+widths[high]/2+widths[low]/2+2*border
        vec_cutout = self.coor_vec([length, width_cutout])
      
        return pos_cutout, width_cutout, vec_cutout
        
    
    def draw_dc_Nrect(self, layer_name, length, rel_pos, widths, border='10um'):
        '''
        REL_POS HAS TO BE A LIST
        '''

        if isinstance(widths, list):
            if len(widths) is not len(rel_pos):
                raise ValueError('widths and rel_pos do not have same length')
        elif isinstance(widths, str):
            widths = [widths for ii in range(len(rel_pos))]
        else:
            raise ValueError('widths should be a string or a list of strings')
        mult = len(rel_pos)
        
        rel_pos, length, widths, border = parse_entry((rel_pos, length, widths, border))
        
        rect = []
        for ii in range(mult):
            rect.append(self.draw_rect(layer_name+'_'+self.name+'_track_'+str(ii), \
                                        self.coor([-length/2, rel_pos[ii]-widths[ii]/2]), \
                                        self.coor_vec([length, widths[ii]])))
        if len(rect) > 1:
            rect = self.unite(rect, name=layer_name+'_'+self.name+'_track')
        else:
            rect = self.rename(rect[0], layer_name+'_'+self.name+'_track')
        self.layers[layer_name]['trackObjects'].append(rect)
        
        pos_cutout, width_cutout, vec_cutout = self.size_dc_gap(length, rel_pos, widths, border)
        cutout = self.draw_rect(layer_name+'_'+self.name+'_cutout', pos_cutout, vec_cutout)
        self.layers[layer_name]['gapObjects'].append(cutout)
        
        portOut1 = [self.coor([length/2, 0]), self.ori, rel_pos, widths, width_cutout, mult] #portOut 
        portOut2 = [self.coor([-length/2, 0]), -self.ori, rel_pos[::-1], widths, width_cutout, mult] #portIn
        self.ports_dc[self.name+'_1'] = portOut1
        self.ports_dc[self.name+'_2'] = portOut2
        
    def draw_dc_test_Nrect(self, layer_name, N, length, rel_pos, widths, border='10um'):
        '''
<<<<<<< HEAD
        Poorly coded: rel_pos is the center of one edge but we use draw_rect instead of draw_rect_center
=======
>>>>>>> simplified routines draw rect
        REL_POS HAS TO BE A LIST
        '''

        if isinstance(widths, list):
            if len(widths) is not len(rel_pos):
                raise ValueError('widths and rel_pos do not have same length')
        elif isinstance(widths, str):
            widths = [widths for ii in range(len(rel_pos))]
        else:
            raise ValueError('widths should be a string or a list of strings')
        mult = len(rel_pos)
        
        rel_pos, length, widths, border = parse_entry((rel_pos, length, widths, border))
        
        rect = []
        for kk in range(N):
            y0 = parse_entry((str(kk*100)+'um'))
            for ii in range(mult):
                rect.append(self.draw_rect(layer_name+'_'+self.name+'_track_'+str(ii), \
                                            self.coor([-length/2, y0+rel_pos[ii]-widths[ii]/2]), \
                                            self.coor_vec([length, widths[ii]])))
        if len(rect) > 1:
            rect = self.unite(rect, name=layer_name+'_'+self.name+'_track')
        else:
            rect = self.rename(rect[0], layer_name+'_'+self.name+'_track')
        self.layers[layer_name]['trackObjects'].append(rect)
        
#        pos_cutout, width_cutout, vec_cutout = self.size_dc_gap(length, rel_pos, widths, border)
#        cutout = self.draw_rect(layer_name+'_'+self.name+'_cutout', pos_cutout, vec_cutout)
#        self.layers[layer_name]['gapObjects'].append(cutout)
        
#        portOut1 = [self.coor([length/2, 0]), self.ori, rel_pos, widths, width_cutout, mult] #portOut 
#        portOut2 = [self.coor([-length/2, 0]), -self.ori, list(-np.array(rel_pos)), widths, width_cutout, mult] #portIn
#        self.ports_dc[self.name+'_1'] = portOut1
#        self.ports_dc[self.name+'_2'] = portOut2
        
        
    def draw_dc_pad(self, layer_name, iTrack, iGap, xlength='250um', ylength='250um'):
        
        iTrack, iGap, xlength, ylength = parse_entry((iTrack, iGap, xlength, ylength))
        
        pad = self.draw_rect('pad',\
                             self.coor([-xlength/2, -ylength/2]),\
                             self.coor_vec([xlength, ylength]))
        padout = self.draw_rect('padout',\
                                self.coor([xlength/2, -iTrack/2]), \
                                self.coor_vec([2*iGap, iTrack]))
        pad = self.unite([pad, padout], name=layer_name+'_'+self.name+'_track')
        pad_gap = self.draw_rect(layer_name+'_'+self.name+'_gap',\
                                 self.coor([-xlength/2-2*iGap, -ylength/2-2*iGap]), \
                                 self.coor_vec([xlength+4*iGap, ylength+4*iGap]))

        self.layers[layer_name]['trackObjects'].append(pad)
        self.layers[layer_name]['gapObjects'].append(pad_gap)

        portOut = [self.coor([xlength/2+2*iGap, 0]), self.ori, iTrack, iGap]
        self.ports[self.name] = portOut
        
        
    def draw_dc_pads(self, layer_name, rel_pos, widths, gaps, xlength='250um', ylength='250um'):

        if isinstance(widths, list):
            if len(widths) is not len(rel_pos):
                raise ValueError('width and rel_pos lists do not have same length')
        elif isinstance(widths, str):
            widths = [widths for ii in range(len(rel_pos))]
        else:
            raise ValueError('width should be a string or a list of strings')
            
        if isinstance(gaps, list):
            if len(gaps) is not len(rel_pos):
                raise ValueError('gap and rel_pos lists do not have same length')
        elif isinstance(gaps, str):
            gaps = [gaps for ii in range(len(rel_pos))]
        else:
            raise ValueError('gap should be a string or a list of strings')
            
        mult = len(rel_pos)
        
        rel_pos, widths, gaps, xlength, ylength = parse_entry((rel_pos, widths, gaps, xlength, ylength))
        
        pads = []
        pad_gaps = []
        for ii in range(mult):
            pad = self.draw_rect('pad'+str(ii), \
                                 self.coor([-xlength/2, rel_pos[ii]-ylength/2]), \
                                 self.coor_vec([xlength, ylength]))
            padout = self.draw_rect('pad_out'+str(ii), \
                                    self.coor([xlength/2, rel_pos[ii]-widths[ii]/2]), \
                                    self.coor_vec([4*gaps[ii], widths[ii]+2*gaps[ii]]))
            pads.append(self.unite([pad, padout]))
            pad_gap = self.draw_rect('gap'+str(ii), \
                                     self.coor([-xlength/2-3*gaps[ii], rel_pos[ii]-ylength/2-3*gaps[ii]]), \
                                     self.coor_vec([xlength+6*gaps[ii], ylength+6*gaps[ii]]))
            pad_gapout = self.draw_rect('gap_out'+str(ii), \
                                        self.coor([xlength/2+3*gaps[ii], rel_pos[ii]-widths[ii]/2-gaps[ii]]), \
                                        self.coor_vec([gaps[ii], widths[ii]+2*gaps[ii]]))
            pad_gaps.append(self.unite([pad_gap, pad_gapout]))
        if len(pads) > 1:
            self.unite(pads, name=layer_name+'_'+self.name+'_track')
            self.unite(pad_gaps, name=layer_name+'_'+self.name+'_gap')
        else:
            self.rename(pads[0], layer_name+'_'+self.name+'_track')
            self.rename(pad_gaps[0], layer_name+'_'+self.name+'_gap')
        self.layers[layer_name]['trackObjects'].append(pads[0])
        self.layers[layer_name]['gapObjects'].append(pad_gaps[0])
        
        cutout_list = [widths[ii]+2*gaps[ii] for ii in range(mult)]
        track_list = [widths[ii] for ii in range(mult)]

        portOut = [self.coor([xlength/2+4*gaps[ii], 0]), self.ori, cutout_list, rel_pos, track_list, mult]
        self.ports_dc[self.name] = portOut
      
    def draw_alignment_marks(self, layer_name, isLayer63=True):
        
        w_thin, l_thin, w_large, l_large, square = '0.1um', '4um', '1um', '3um', '12um'
        w_thin, l_thin, w_large, l_large, square = parse_entry((w_thin, l_thin, w_large, l_large, square))
        marks = []
        if isLayer63 is True:
            self.new_layer('layer63')
            squares = []
        count = -1
        for x0 in ['-42um', '42um']:
            count += 1
            for y0 in ['-42um', '42um']:
                count += 1
                mark = []
                x0, y0 = parse_entry((x0, y0))
                mark.append(self.draw_rect('thin_x', self.coor([x0-l_thin/2, y0-w_thin/2]), self.coor_vec([l_thin, w_thin])))
                mark.append(self.draw_rect('thin_y', self.coor([x0-w_thin/2, y0-l_thin/2]), self.coor_vec([w_thin, l_thin])))
                mark.append(self.draw_rect('large_right', self.coor([x0+l_thin/2, y0-w_large/2]), self.coor_vec([l_large, w_large])))
                mark.append(self.draw_rect('large_left', self.coor([x0-l_thin/2, y0-w_large/2]), self.coor_vec([-l_large, w_large])))
                mark.append(self.draw_rect('large_top', self.coor([x0-w_large/2, y0+l_thin/2]), self.coor_vec([w_large, l_large])))
                mark.append(self.draw_rect('large_bot', self.coor([x0-w_large/2, y0-l_thin/2]), self.coor_vec([w_large, -l_large])))
                mark = self.unite(mark, name=layer_name+'_'+self.name+'_alignement_mark_'+str(count))
                marks.append(mark)
                if isLayer63 is True:
                    squares.append(self.draw_rect_center('layer63', self.coor([x0, y0]), self.coor_vec([square, square])))
                    
        marks = self.unite(marks, name=layer_name+'_'+self.name+'_alignement_mark')
        self.layers[layer_name]['trackObjects'].append(marks)
        if isLayer63 is True:
            squares = self.unite(squares, name='layer63_'+self.name)
            
        
    def draw_pits(self, rel_pos, length, width):
        '''
        rel_pos gives the position of the fine structure
        length is the length along the CNT region
        width is in the transverse direction, spans the whole comb
        '''
        self.new_layer('pits')
        rel_pos, length, width = parse_entry((rel_pos, length, width))
        
        print(rel_pos)
        dist2CNT = max([abs(pos) for pos in rel_pos])
        dist2CNT += parse_entry('10um')
        space = parse_entry('40um')
        close_pit_left = self.draw_rect('close_pit_left', \
                                        self.coor([-length/2, dist2CNT]), \
                                        self.coor_vec([length, space]))
        close_pit_right = self.draw_rect('close_pit_right', \
                                         self.coor([-length/2, -dist2CNT]), \
                                         self.coor_vec([length, -space]))
        pit_left = self.draw_rect('pit_left', \
                                  self.coor([-length, dist2CNT+space]), \
                                  self.coor_vec([2*length, width]))
        pit_right = self.draw_rect('pit_right', \
                                    self.coor([-length, -dist2CNT-space]), \
                                    self.coor_vec([2*length, -width]))
        pit_left = self.unite([pit_left, close_pit_left], name='pits_'+self.name+'_left')
        pit_right = self.unite([pit_right, close_pit_right], name='pits_'+self.name+'_right')
        self.layers['pits']['trackObjects'].append(pit_left)
        self.layers['pits']['trackObjects'].append(pit_right)
        
        
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

    