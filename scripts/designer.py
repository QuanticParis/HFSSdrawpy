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

'''

'''
Assumption

- only assign lumpRLC to rectangles
- assume to have list and not tuples for polylines
- TODO, do not do Lj+'1nH' for now, but can do bx+'1mm'
'''

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
        super().__init__(vec)

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
    ports = {}

    def __init__(self, design, modeler):

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


    def unite(self, iObjs, name=None):
        '''
        Performs the unions of the elements of iObjects

        Input:
        iObjects (list of strings): list of object names e.g. ['transmon_pad','transmon_pad_extension']

        Returns:
        string: name of the merged object: iObjects[0]
        '''
        iObj = self.modeler.unite(iObjs, name=name)
        return iObj


    def assign_perfE(self, iObj, iSuff='PerfE'):
        '''
        Assigns boundary condition PerfE of name name to object iObject

        Inputs:
        -------
        iObject: name of the object e.g. transmon_pad
        name: name of the perfect E e.g. PerfE1

        '''
        self.modeler.assign_perfect_E(iObj, name=iSuff)


    def subtract(self, iObjBlank, iObjTools):
        '''
        suBTracts iObjectTools from iObjectBlanc

        Inputs:
        -------
        iObjectBlanc (string) : HFSS object name e.g. 'ground_plane'
        iObjectTools (list) : HFSS object name e.g. ['readout_box', 'res_gap']

        '''
        iObjBlank = self.modeler.subtract(iObjBlank, iObjTools, keep_originals=False)
        return iObjBlank


    # main function
    def draw(self, iObj, iPoints, closed=True): #assume everything is in the plane z=0
        '''
        Inputs
        ------
        name : object name in HFSS e.g. transmon_pad
        iPoints : list of tuples [(x0, y0), (x1, y1), ..., (xn, yn)]

        '''
        points = []
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
                prev_val_x = val_x
                prev_val_y = val_y
            elif not equal_float(prev_val_x, val_x) or not equal_float(prev_val_y, val_y):
                points.append([x, y, 0])
                prev_val_x = val_x
                prev_val_y = val_y
            else:
                pass
#                print('Warning: Found two overlapping points while creating the polyline, supressed one')
        iObj = self.modeler.draw_polyline(points, closed=closed, name=iObj)
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

    def draw_rect_center(self, name, pos, iSize):
        pos = [pos[0], pos[1], 0]
        size = [iSize[0], iSize[1], 0]
        rect = self.modeler.draw_rect_center(pos, size, name=name)
        self.__dict__[rect] = rect
        return rect

    def draw_rect(self, name, pos, iSize):
        pos = [pos[0], pos[1], 0]
        size = [iSize[0], iSize[1], 0]
        rect = self.modeler.draw_rect_corner(pos, size, name=name)
        self.__dict__[rect] = rect
        return rect

    def key_elt(self, name='key_elt_0', pos=[0,0], ori=[1,0]):
        obj = KeyElt(name, pos=pos, ori=ori)
        self.__dict__[name] = obj
        return obj

    def connect_elt(self, name='connect_elt_0', iIn='iInt', iOut=None):
        obj = ConnectElt(name=name, iIn=iIn, iOut=iOut)
        self.__dict__[name] = obj
        return obj

class KeyElt(Circuit):

    pcb_track = parse_entry('300um')
    pcb_gap = parse_entry('200um')
#
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

    def refx_points(self, coor_list, offset=0):
        points=[]
        for ii, coor in enumerate(coor_list):
            if ii>0:
                offset=0
            points.append(Vector(*coor).refx(offset))
        return points

    def refy_points(self, coor_list, offset=0):
        points=[]
        for ii, coor in enumerate(coor_list):
            if ii>0:
                offset=0
            points.append(Vector(*coor).refy(offset))
        return points

    def coor(self, vec): # Change of coordinate for a point
        return self.rot(*vec)+self.pos

    def coor_vec(self, vec): # Change of coordinate for a vector
        return self.rot(*vec)

    def draw_connector(self, iTrack, iGap, iBondLength, iSlope=1, iLineTest=False):
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

        adaptDist = (self.pcb_track/2-iTrack/2)/iSlope

        portOut = [self.pos+self.ori*(adaptDist+self.pcb_gap+iBondLength), self.ori, iTrack, iGap]
#        print(self.pos, self.ori)
#        print(adaptDist)
#        print(self.pos+self.ori*(adaptDist+iGap+iBondLength), self.ori)
        points = self.append_points([(self.pcb_gap, self.pcb_track/2),
                                     (iBondLength, 0),
                                     (adaptDist, iTrack/2-self.pcb_track/2),
                                     (0, -iTrack),
                                     (-adaptDist, iTrack/2-self.pcb_track/2),
                                     (-iBondLength, 0)])
        self.trackObjects.append(self.draw(self.name+"_track", points))

        points = self.append_points([(self.pcb_gap/2, self.pcb_gap+self.pcb_track/2),
                             (self.pcb_gap/2+iBondLength, 0),
                             (adaptDist, (iGap-self.pcb_gap)+(iTrack-self.pcb_track)*0.5),
                             (0, -2*iGap-iTrack),
                             (-adaptDist, (iGap-self.pcb_gap)+(iTrack-self.pcb_track)*0.5),
                             (-(self.pcb_gap/2+iBondLength), 0)])
        self.gapObjects.append(self.draw(self.name+"_gap", points))

        if(iLineTest==False):
            points = self.append_points([(self.pcb_gap/2, 0.5*self.pcb_track),
                                         (self.pcb_gap/2, 0),
                                         (0, -self.pcb_track),
                                         (-self.pcb_gap/2, 0)])
            ohm = self.draw(self.name+"_ohm", points)
            self.assign_lumped_RLC(ohm, self.ori, ('50ohm',0,0))

#        self.iIn = retIn
#        self.iOut = retOut
        self.ports[self.name] = portOut


    def draw_JJ(self, iTrack, iGap, iTrackJ, iLength, iInduct='1nH', add_port=True):
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

        adaptDist = iTrack/2-iTrackJ/2

        if self.val(adaptDist)>self.val(iLength/2-iTrackJ/2):
            raise ValueError('Increase iTrackJ %s' % self.name)

        portOut1 = [self.pos+self.ori*iLength/2, self.ori, iTrack, iGap]
        portOut2 = [self.pos-self.ori*iLength/2, -self.ori, iTrack, iGap]

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
        self.trackObjects.append(pads)

        mesh = self.draw_rect_center(self.name+'_mesh', self.coor([0,0]), self.coor_vec([iLength, iTrack]))
        self.modeler.assign_mesh_length(mesh, iTrackJ/2, suff='')

        points = self.append_points([(iTrackJ/2,0),(-iTrackJ,0)])
        line = self.draw(self.name+'_line', points, closed=False)

        if add_port:
            self.gapObjects.append(self.draw_rect_center(self.name+"_cutout", self.coor([0,0]), self.coor_vec([iLength, iTrack+2*iGap])))

        JJ = self.draw_rect_center(self.name, self.coor([0,0]), self.coor_vec([iTrackJ, iTrackJ]))
        self.assign_lumped_RLC(JJ, self.ori, (0, iInduct, 0))

        if add_port:
            self.ports[self.name+'_1'] = portOut1
            self.ports[self.name+'_2'] = portOut2

        return pads

    def draw_IBM_tansmon(self,
                       cutout_size,
                       pad_spacing,
                       pad_size,
                       Jwidth,
                       track,
                       gap,
                       Jinduc,
                       nport=1,
                       fillet=None,
                       maskbox=False):

        cutout_size, pad_spacing, pad_size, Jwidth, track, gap = parse_entry((cutout_size, pad_spacing, pad_size, Jwidth, track, gap))
        cutout_size = Vector(cutout_size)
        pad_size = Vector(pad_size)

        self.gapObjects.append(self.draw_rect_center(self.name+"_cutout", self.coor([0,0]), self.coor_vec(cutout_size)))
        
        if maskbox:
            self.gapObjects.append(self.draw_rect_center(self.name+"_cutout", self.coor([0,0]), self.coor_vec(cutout_size)))
        
        mesh = self.draw_rect_center(self.name+"_mesh", self.coor([0,0]), self.coor_vec(cutout_size))
        self.modeler.assign_mesh_length(mesh, 2*track, suff='')


        track_J=Jwidth*4.
        pos_junction = self.coor([0,0])
        junction = self.key_elt(self.name+'_junction', pos_junction, self.ori)
        junction_pad = junction.draw_JJ(track_J, (cutout_size[1]-track_J)/2, Jwidth, pad_spacing, iInduct=Jinduc, add_port=False)
        self.trackObjects.pop()

        right_pad = self.draw_rect_center(self.name+"_pad1", self.coor(Vector(pad_spacing+pad_size[0],0)/2), self.coor_vec(pad_size))
        left_pad = self.draw_rect_center(self.name+"_pad2", self.coor(-Vector(pad_spacing+pad_size[0],0)/2), self.coor_vec(pad_size))
        
        if fillet is not None:
            right_pad.fillet(fillet,3)
            right_pad.fillet(fillet,2)
            right_pad.fillet(fillet,1)
            right_pad.fillet(fillet,0)
            
            left_pad.fillet(fillet,3)
            left_pad.fillet(fillet,2)
            left_pad.fillet(fillet,1)
            left_pad.fillet(fillet,0)

        pads = self.unite([right_pad, left_pad, junction_pad], name=self.name+'_pads')
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
            self.ports[self.name+'_1'] = [self.coor(cutout_size.px()/2), self.ori, track, gap]
            self.ports[self.name+'_2'] = [self.coor(-cutout_size.px()/2), -self.ori, track, gap]
            self.ports[self.name+'_3a'] = [self.coor(Vector(pad_spacing/2+pad_size[0],-cutout_size[1]/2)), -self.ori.orth(),track/2,gap/2]

#        self.draw_rect_center(self.name+"check1", self.coor(cutout_size.px()/2)+self.ori*pad_size[0]/20, self.rot(*pad_size)/10)
#        self.draw_rect_center(self.name+"check2", self.coor(-cutout_size.px()/2)-self.ori*pad_size[0]/20, self.rot(*pad_size)/10)
#        self.draw_rect_center(self.name+"check3", self.coor(Vector(pad_spacing/2+pad_size[0],cutout_size[1]/2))+self.ori.orth()*pad_size[1]/20, self.rot(*pad_size)/10)

    def draw_capa(self, iTrack, iGap, pad_spacing, pad_size):
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
        iIn |    ||      | iOut
            +-+  |  |  +-+
              |  |  |  |
              +--+  +--+
        '''
        iTrack, iGap, pad_spacing, pad_size = parse_entry((iTrack, iGap, pad_spacing, pad_size))
        pad_size = Vector(pad_size)

        portOut1 = [self.pos+self.ori*(pad_spacing/2+pad_size[0]+iGap), self.ori, iTrack, iGap]
        portOut2 = [self.pos-self.ori*(pad_spacing/2+pad_size[0]+iGap), -self.ori, iTrack, iGap]

        raw_points = [(pad_spacing/2, pad_size[1]/2),
                      (pad_size[0], 0),
                      (0, -(pad_size[1]-iTrack)/2),
                      (iGap, 0),
                      (0, -iTrack),
                      (-iGap, 0),
                      (0, -(pad_size[1]-iTrack)/2),
                      (-pad_size[0], 0)]
        points = self.append_points(raw_points)
        right_pad = self.draw(self.name+"_pad1", points)

        points = self.append_points(self.refy_points(raw_points))
        left_pad = self.draw(self.name+"_pad2", points)

        pads = self.unite([right_pad, left_pad], name=self.name+'_pads')
        self.trackObjects.append(pads)

        self.gapObjects.append(self.draw_rect_center(self.name+"_cutout", self.coor([0,0]), self.coor_vec([pad_spacing + 2*pad_size[0]+2*iGap, pad_size[1] + 2*iGap])))

        self.draw(self.name+"_mesh", points)
        self.modeler.assign_mesh_length(self.name+"_mesh",iTrack)

        self.ports[self.name+'_1'] = portOut1
        self.ports[self.name+'_2'] = portOut2
        
        
    def draw_capa_interdigitated(self, iTrack, iGap, teeth_size,gap_size, N_period, fillet):
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
        iTrack, iGap = parse_entry((iTrack, iGap))
        fillet = parse_entry(fillet)
        teeth_size=parse_entry(teeth_size)
        gap_size=parse_entry(gap_size)
        teeth_size = Vector(teeth_size)
        portOut1 = [self.pos+self.ori*(teeth_size[0]+iTrack+iGap), self.ori, iTrack, iGap]
        portOut2 = [self.pos-self.ori*(teeth_size[0]+iTrack+iGap), -self.ori, iTrack, iGap]


        N_teeth=2*N_period+1
        raw_points = [(teeth_size[0], -N_teeth*teeth_size[1])]
        raw_points.append((teeth_size[0], (-N_teeth+1)*teeth_size[1]))
        for i in range(-N_teeth+1,N_teeth-1,4):
            raw_points.append((-teeth_size[0], i*teeth_size[1]))
            raw_points.append((-teeth_size[0], (i+2)*teeth_size[1]))
            raw_points.append((teeth_size[0], (i+2)*teeth_size[1]))
            raw_points.append((teeth_size[0], (i+4)*teeth_size[1]))
        raw_points.append((-teeth_size[0], (N_teeth-1)*teeth_size[1]))
        raw_points.append((-teeth_size[0], N_teeth*teeth_size[1]))


#        raw_points = [(teeth_size[0], -3*teeth_size[1]),
#                      (teeth_size[0], -2*teeth_size[1]),
#                      (-teeth_size[0], -2*teeth_size[1]),
#                      (-teeth_size[0], 0*teeth_size[1]),
#                      (teeth_size[0], 0*teeth_size[1]),
#                      (teeth_size[0], 2*teeth_size[1]),
#                      (-teeth_size[0], 2*teeth_size[1]),
#                      (-teeth_size[0], 3*teeth_size[1])]
        points = self.append_absolute_points(raw_points)
        connection = self.draw(self.name+"_capagap", points, closed=False)
        
        
        connection.fillets(fillet)
        raw_points=[(-gap_size-teeth_size[0],N_teeth*teeth_size[1]),(gap_size-teeth_size[0],N_teeth*teeth_size[1])]
        points=self.append_absolute_points(raw_points)
        capagap_starter = self.draw(self.name+'_width', points, closed=False)
        
        capagap = connection.sweep_along_path(capagap_starter)
        
   
        
        raw_points = [(-teeth_size[0]-iTrack, -N_teeth*teeth_size[1]),
                      (-teeth_size[0]-iTrack,-iTrack/2),
                      (-teeth_size[0]-iTrack-iGap,-iTrack/2),
                      (-teeth_size[0]-iTrack-iGap, iTrack/2),
                      (-teeth_size[0]-iTrack, iTrack/2),
                      (-teeth_size[0]-iTrack, N_teeth*teeth_size[1]),
                      (teeth_size[0]+iTrack,  N_teeth*teeth_size[1]),
                      (teeth_size[0]+iTrack,iTrack/2),
                      (teeth_size[0]+iTrack+iGap,iTrack/2),
                      (teeth_size[0]+iTrack+iGap,-iTrack/2),
                      (teeth_size[0]+iTrack, -iTrack/2),
                      (teeth_size[0]+iTrack, -N_teeth*teeth_size[1])]
        points = self.append_absolute_points(raw_points)
        pads = self.draw(self.name+"_pads", points)
        pads.fillet(fillet,11)
        pads.fillet(fillet,6)
        pads.fillet(fillet,5)
        pads.fillet(fillet,0)
        
        pads_sub = self.subtract(pads, [capagap])
        
        
        raw_points = [(-teeth_size[0]-iTrack-iGap, -N_teeth*teeth_size[1]-iGap),
                      (teeth_size[0]+iTrack+iGap, -N_teeth*teeth_size[1]-iGap),
                      (teeth_size[0]+iTrack+iGap, N_teeth*teeth_size[1]+iGap),
                      (-teeth_size[0]-iTrack-iGap, N_teeth*teeth_size[1]+iGap)]
        points = self.append_absolute_points(raw_points)
        gap = self.draw(self.name+"_gap", points)
        #gap.fillet(fillet,3)
        #gap.fillet(fillet,2)
        #gap.fillet(fillet,1)
        #gap.fillet(fillet,0)
        self.draw(self.name+"_mesh", points)
        self.modeler.assign_mesh_length(self.name+"_mesh",iTrack)


        self.trackObjects.append(pads_sub)

        self.gapObjects.append(gap)
        #self.gapObjects.append(capagap)
        #self.draw(self.name+"_mesh", points)
        #self.modeler.assign_mesh_length(self.name+"_mesh",1/2*pad_size[1])

        self.ports[self.name+'_1'] = portOut1
        self.ports[self.name+'_2'] = portOut2
    

    def draw_squid(self, iTrack, iGap, squid_size, iTrackPump, iGapPump, iTrackSquid=None, iTrackJ=None, Lj_down='1nH', Lj_up=None,  typePump='down', doublePump=False, iSlope=1, iSlopePump=0.5): #for now assume left and right tracks are the same width
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
        pos_junction_up = self.coor(squid_size.py()/2+Vector([0,1])*iTrackSquid/2)
        junction = self.key_elt(self.name+'_junction_up', pos_junction_up, self.ori)
        junction_pads_up = junction.draw_JJ(iTrackSquid, iTrackSquid, iTrackJ, iTrackSquid, iInduct=Lj_up, add_port=False)
        self.trackObjects.pop()
        #junction down
        pos_junction_down = self.coor((squid_size.py()/2+Vector([0,1])*iTrackSquid/2).refx())
        junction = self.key_elt(self.name+'_junction_down', pos_junction_down, self.ori)
        junction_pads_down = junction.draw_JJ(iTrackSquid, iTrackSquid, iTrackJ, iTrackSquid, iInduct=Lj_down, add_port=False)
        self.trackObjects.pop()

        right_track = self.draw_rect_center(self.name+"_added_track1", self.coor([2*(squid_size[0]/2+adapt_dist),0]), self.coor_vec([squid_size[0]+2*adapt_dist, iTrack]))
        left_track = self.draw_rect_center(self.name+"_added_track2", self.coor([-2*(squid_size[0]/2+adapt_dist),0]), self.coor_vec([squid_size[0]+2*adapt_dist, iTrack]))

        squid = self.unite([right_track, left_track, track_a, track_b, track_c, track_d, junction_pads_down, junction_pads_up], name=self.name)
        self.trackObjects.append(squid)


        self.gapObjects.append(self.draw_rect_center(self.name+"_added_gap", self.coor([0,0]), self.coor_vec([3*(squid_size[0]+2*adapt_dist), iTrack+2*iGap])))

        adapt_dist_pump = (2*iTrackPump - 2*iTrackSquid)/2/iSlopePump

        raw_points_adapt_pump_a = [(3/2*(squid_size[0]+2*adapt_dist),-iTrack/2-iGap-iTrackSquid),
                                 (-(squid_size[0]+2*adapt_dist)+iTrackSquid,0),
                                 (iTrackPump/2-iTrackSquid, -adapt_dist_pump),
                                 (iGapPump, 0)]
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
        self.gapObjects.append(self.draw(self.name+"_cutout_pump_a", points_adapt_pump_a))

        raw_points_adapt_pump_b = self.refy_points(raw_points_adapt_pump_a, offset = (squid_size[0]+2*adapt_dist)/2)
        points_adapt_pump_b = self.append_points(raw_points_adapt_pump_b)
        self.gapObjects.append(self.draw(self.name+"_cutout_pump_b", points_adapt_pump_b))

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
            portOutpump2 = [self.coor(pos_pump.refx()), -self.coor_vec(ori_pump), iTrackPump, iGapPump]
            self.ports[self.name+'_pump2'] = portOutpump2
        else:
            portOutpump1 = [self.coor(pos_pump), self.coor_vec(ori_pump), iTrackPump, iGapPump]
            self.ports[self.name+'_pump'] = portOutpump1

    def draw_end_cable(self, iTrack, iGap, typeEnd = 'open'):

        self.gapObjects.append(self.draw_rect_center(self.name+"_cutout", self.coor([iGap/2,0]), self.coor_vec([iGap, iTrack+2*iGap])))
        if typeEnd=='open' or typeEnd=='Open':
            portOut = [self.pos+self.ori*iGap, self.ori, iTrack, iGap]
        elif typeEnd=='short' or typeEnd=='Short':
            portOut = [self.pos, self.ori, iTrack, iGap]
        else:
            raise ValueError("typeEnd should be 'open' or 'short', given %s" % typeEnd)
        self.ports[self.name] = portOut



class ConnectElt(KeyElt, Circuit):

    def __init__(self, name='connect_elt', iIn='iInt', iOut=None):
        print(name)
#        print(iIn)
#        print(iOut)
        self.name = name
        self.iIn = iIn # name of the in port

        if iIn in self.ports:
            iInPort = parse_entry(self.ports[iIn])
            self.pos = Vector(iInPort[POS])
            self.ori = Vector(iInPort[ORI])
            _, _, self.inTrack, self.inGap = iInPort
        else:
            raise ValueError('inPort %s does not exist' % iIn)

        if isinstance(iOut, str):
            if iOut in self.ports:
                iOutPort = parse_entry(self.ports[iOut])
                self.isOut = True
                self.iOut = iOut
                self.posOut = Vector(iOutPort[POS])
                self.oriOut = Vector(iOutPort[ORI])
                _, _, self.outTrack, self.outGap = iOutPort
            else:
                raise ValueError('outPort %s does not exist' % iOut)
        elif isinstance(iOut, list):
            iOutPort = parse_entry(iOut)
            self.outTrack, self.outGap = iOutPort[-2], iOutPort[-1]
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

        self.draw(self.name+"_mesh", points)
        self.modeler.assign_mesh_length(self.name+"_mesh",1/2*iLength)

        self.iIn = retIn
        self.iOut = retOut
#        return [retIn, retOut]

    def draw_half_capa(self, iLength, iWidth, iGap, add_gap=False,fillet=None):
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

        Outputs:
        --------
        retIn: same as iIn, with flipped vector
        retOut: calculated output port to match all input dimensions

            igap iWidth
                 +--+
                 |  |
            +----+  | iLength
        iIn |       |
            +----+  |
                 |  |
                 +--+
        '''
        iLength, iWidth, iGap = parse_entry((iLength, iWidth, iGap))
        self.ori = -self.ori

        points = self.append_points([(0, self.inTrack/2),
                                     (iGap, 0),
                                     (0, (iLength-self.inTrack)/2),
                                     (iWidth, 0),
                                     (0, -iLength),
                                     (-iWidth, 0),
                                     (0, (iLength-self.inTrack)/2),
                                     (-iGap, 0)])
        halfcapa=self.draw(self.name+"_pad", points)
        if fillet is not None:
            halfcapa.fillet(fillet,6)
            halfcapa.fillet(fillet,5)
            halfcapa.fillet(fillet,4)
            halfcapa.fillet(fillet,3)
            halfcapa.fillet(fillet,2)
            halfcapa.fillet(fillet,1)
            
        self.trackObjects.append(halfcapa)

#        CreateBondwire(name+"_bondwire", iIn)

    def find_path(self, fillet, is_meander, to_meander, meander_length):

        iIn_pos = self.pos
        iIn_ori = self.ori
        iOut_pos = self.posOut
        iOut_ori = self.oriOut

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

        def displace(points, rl, min_dist, displacement=0, n_meander=-1):
            if self.val(displacement)<self.val(min_dist)*1.1:
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
                    new_points.append(points[ii+n_ignore]+vec.orth()*parity*displacement-vec*dist/2)
                    new_points.append(points[ii+n_ignore]+vec.orth()*parity*displacement+vec*dist/2)
                    parity = -parity
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
                    if isit!=0:
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
        if A<0 or A>=len(points):
            raise ValueError('First index should be within the point list')
        if B<0 or B>=len(points):
            raise ValueError('Second index should be within the point list')
        if A==B:
            return 0
        if abs(A-B)==1:
            if A<B:
#                if A%2==1:
#                    return self.val(fillet*np.pi/2)
#                else:
#                    return self.val((points[A]-points[B]).norm())

                if A==0 or B==len(points)-1:
                    return self.val((points[A]-points[B]).norm())-self.val(fillet*(2-np.pi/2))
                else:
                    return self.val((points[A]-points[B]).norm())-self.val(fillet*(1-np.pi/4))
            else:
                return self.length(points, B, A, fillet)
        if abs(A-B)>1:
            if A<B:
                return self.length(points, A, B-1, fillet) + self.length(points, B-1, B, fillet)
            else:
                return self.length(points, B, A, fillet)


#    def length_exp(self, points, A, B, fillet): # A and B are integer point indices
#        if A<0 or A>=len(points):
#            raise ValueError('First index should be within the point list')
#        if B<0 or B>=len(points):
#            raise ValueError('Second index should be within the point list')
#        if A==B:
#            return 0
#        if abs(A-B)==1:
#            if A<B:
#                if A==0 or B==len(points)-1:
#  
#                    return (points[A]-points[B]).norm()-(2-np.pi/2)*fillet
#                else:
#                    return (points[A]-points[B]).norm()-(1-np.pi/4)*fillet
#            else:
#                return self.length(points, B, A, fillet)
#        if abs(A-B)>1:
#            if A<B:
#                return self.length(points, A, B-1, fillet) + self.length(points, B-1, B, fillet)
#            else:
#                return self.length(points, B, A, fillet)

    def cable_starter(self, width = 'track'): # width can also be 'gap'
        if width=='track' or width=='Track':
            points = self.append_points([(0, self.inTrack/2),
                                         (0, -self.inTrack)])
        elif width=='gap' or width=='Gap':
            points = self.append_points([(0, self.inGap+self.inTrack/2),
                                         (0, -2*self.inGap-self.inTrack)])
        return self.draw(self.name+'_width'+width, points, closed=False)


    def draw_cable(self, fillet="0.3mm", is_bond=True, is_meander=False, to_meander = [1,0,1,0,1,0,1,0,1,0], meander_length=0, ismesh=False):
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
        fillet, meander_length=parse_entry((fillet, meander_length))
#        inPos, inOri = Vector(iIn[POS]), Vector(iIn[ORI])
#        outPos, outOri = Vector(iOut[POS]), Vector(iOut[ORI])
#        _, _, track, gap = iIn
        self.to_bond=[]
        adaptor_length=0
        track_adaptor = None
        if(self.val(self.inTrack) != self.val(self.outTrack) or self.val(self.inGap) != self.val(self.outGap)):
            if self.val(self.inTrack+self.inGap) > self.val(self.outTrack+self.outGap):
                adaptor = ConnectElt(self.name+'_adaptor', self.iOut, [self.inTrack, self.inGap])
                iOut, adaptor_length, track_adaptor = adaptor.draw_adaptor()
                self.__init__(self.name, self.iIn, iOut)
            else:
                adaptor = ConnectElt(self.name+'_adaptor', self.iIn, [self.outTrack, self.outGap])
                iIn, adaptor_length, track_adaptor = adaptor.draw_adaptor()
                self.__init__(self.name, iIn, self.iOut)

        points = self.find_path(fillet, is_meander, to_meander, meander_length)

        connection = self.draw(self.name+"_track", points, closed=False)



        cable_length = self.length(points, 0, len(points)-1, fillet)+self.val(adaptor_length)
        print('{0}_length = {1:.3f} mm'.format(self.name, cable_length*1000))
#        if self.name=='capa_cable':
#            print(points)
        

        connection.fillets(fillet)

        connection_gap = connection.copy(self.name+"_gap")

        track_starter = self.cable_starter('track')
        gap_starter = self.cable_starter('gap')

        track = connection.sweep_along_path(track_starter)
        
        if ismesh is True:
            self.modeler.assign_mesh_length(track,2*self.inTrack)

        if track_adaptor is not None:
            self.trackObjects.pop()
            to_unite = [track, track_adaptor]
            track = self.unite(to_unite)
        self.trackObjects.append(track)

        self.gapObjects.append(connection_gap.sweep_along_path(gap_starter))

        if is_bond:
            self.draw_bond((self.inTrack+self.inGap*2)*1.5)


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

        self.gapObjects.append(self.draw(self.name+"_gap", points))

        return self.iIn+'_bis', adaptDist, track
