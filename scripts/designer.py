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


from hfss import VariableString
from hfss import parse_entry

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
        print(self, other)
        print("breakpoint1")
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

class Design(object):
    def __init__(self, name):
        self.name = name
        
class Object(object):
    def __init__(self, model, material):
        self.model = model # string with "Sheet" or "Solids"
        self.material = material
        self.boundaries = []
        self.subObject = []
    
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
            
            
    
class Chip(Design):

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
        POS = 0
        ORI = 1
        TRACK = 2
        GAP = 3
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

   
    