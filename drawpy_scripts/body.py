import numpy as np
from functools import wraps

from .utils import Vector, \
                   parse_entry, \
                   check_name, \
                   find_last_list, \
                   find_corresponding_list, \
                   to_move, \
                   _val, \
                   val, \
                   way

from .model_entity import ModelEntity
from .python_modeler import PythonModeler
from .path_finder import Path

from .parameters import layer_TRACK, \
                        layer_GAP, \
                        layer_MASK, \
                        layer_Default, \
                        layer_RLC, \
                        layer_PORT, \
                        layer_MESH

class Port():
    instances_to_move = None
    dict_instances  = {}

    def __init__(self, name, pos, ori, widths, subnames, layers, offsets, constraint_port, key='name'):
        if not (isinstance(key, Port) or key is None):
            name = check_name(self.__class__, name)
        self.name = name
        self.pos = Vector(pos)
        self.ori = Vector(ori)
        self.constraint_port = constraint_port
        self.save = None
        if not constraint_port:
            self.widths = parse_entry(widths)
            self.subnames = subnames
            self.layers = layers
            self.offsets = parse_entry(offsets)
            self.N = len(widths)
        else:
            self.widths = widths
            self.subnames = subnames
            self.layers = layers
            self.offsets = offsets
            self.N = 0

        if Port.instances_to_move is not None:
            find_last_list(Port.instances_to_move).append(self)
        if key=='name':  # normal initialisation
            self.dict_instances[name] = self

            # create a reversed version of the port that can be called by either
            # port.r or 'port_name.r'
            reversed_ori = -self.ori
            reversed_offsets = None
            if self.offsets is not None:
                reversed_offsets = []
                for ii in range(self.N):
                    reversed_offsets.append(-self.offsets[ii])
            self.r = Port(self.name+'_r', self.pos, reversed_ori, self.widths,
                          self.subnames, self.layers, reversed_offsets,
                          self.constraint_port, key=self)

        elif isinstance(key, Port):  # reverse initialisation, key is the previous port
            self.dict_instances[name] = self
            self.r = key
        else:
            pass  # when the port is only a float eval do not add it in dict

    def __str__(self):
        return self.name

    def __repr__(self):
        return self.name

    @staticmethod
    def reset():
        Port.instances_to_move = []
        Port.dict_instances  = {}

    @classmethod
    def print_instances(cls):
        for instance_name in cls.dict_instances:
            print(instance_name)#, cls.dict_instances[instance_name])

    def compare(self, other, pm):
        points = []

        adapt_dist = pm.set_variable(1e-5, name=self.name+'_adapt')
        max_diff = 0
        for ii in range(self.N):
            if self.layers[ii]!=other.layers[ii]:
                raise ValueError('Tried to connect ports form different \
                                 layers: %s != %s'%(self.layers[ii],
                                                    other.layers[ii]))
            width1 = self.widths[ii]
            width2 = other.widths[ii]

            offset1 = self.offsets[ii]
            offset2 = -other.offsets[ii]

            if width1!=width2 or offset1!=offset2:
                # need adaptor
                points.append([Vector(0, offset1+width1/2).rot(self.ori)+self.pos,
                               Vector(adapt_dist, offset2+width2/2).rot(self.ori)+self.pos,
                               Vector(adapt_dist, offset2-width2/2).rot(self.ori)+self.pos,
                               Vector(0, offset1-width1/2).rot(self.ori)+self.pos])
            max_diff = max(max_diff, abs(_val(offset1+width1/2-(offset2+width2/2))),
                           abs(_val(offset2-width2/2-(offset1-width1/2))))

        adapt_dist = pm.set_variable(2*max_diff, name=self.name+'_adapt')

        if len(points) != 0:
            self.save = {'pos':self.pos, 'widths':self.widths,
                         'offsets':self.offsets}
            self.pos = self.pos + Vector(2*max_diff, 0).rot(self.ori)
            self.widths = other.widths
            self.offsets = [-offset for offset in other.offsets]

            self.r.save = {'pos':self.r.pos, 'widths':self.r.widths,
                         'offsets':self.r.offsets}
            self.r.pos = self.pos
            self.r.widths = self.widths
            self.r.offsets = other.offsets
        return points, 2*max_diff

    def val(self):
        _widths = []
        _offsets = []
        for ii in range(self.N):
            width = self.widths[ii]
            offset = self.offsets[ii]
            _widths.append(_val(width))
            _offsets.append(_val(offset))

        _pos = []
        for coor in self.pos:
            _pos.append(_val(coor))
        _pos = Vector(_pos)

        _ori = []
        for coor in self.ori:
            _ori.append(_val(coor))
        _ori = Vector(_ori)

        return Port(self.name, _pos, _ori, _widths, self.subnames,
                    self.layers, _offsets, self.constraint_port, key=None)

    def revert(self):
        if self.save is not None:
            self.pos = self.save['pos']
            self.widths = self.save['widths']
            self.offsets = self.save['offsets']

            self.r.pos = self.save['pos']
            self.r.widths = self.save['widths']
            reversed_offsets = []
            for ii in range(self.N):
                reversed_offsets.append(-self.offsets[ii])
            self.r.offsets = reversed_offsets

    def bond_params(self):
        y_max = -np.infty
        y_max_val = -np.infty
        y_min = np.infty
        y_min_val = np.infty
        for ii in range(self.N):
            # widths should not be negative
            _y_max_val = _val(self.offsets[ii]+self.widths[ii]/2)
            _y_min_val = _val(self.offsets[ii]-self.widths[ii]/2)
            if _y_max_val > y_max_val:
                y_max = self.offsets[ii]+self.widths[ii]/2
                y_max_val = _y_max_val
            if _y_min_val < y_min_val:
                y_min = self.offsets[ii]-self.widths[ii]/2
                y_min_val = _y_min_val
        return y_max, y_min

class Body(PythonModeler):

    dict_instances = {}
    def __init__(self, pm=None, coor_sys=None, rel_coor=None, ref_name='Global'): #network
        # Note: for now coordinate systems are not reactualized at each run
        if rel_coor is None:
            rel_coor = [[0, 0, 0],  # origin
                        [1, 0, 0],  # new_x
                        [0, 1, 0]]  # new_y
        else:
            rel_coor = parse_entry(rel_coor)

        pm.interface.create_coor_sys(coor_sys=coor_sys, rel_coor=rel_coor,
                                     ref_name=ref_name)

        self.pm = pm
        self.coor_sys = coor_sys # also called body_name
                                 # a body is equivalent to a coor_sys
        self.rel_coor = rel_coor
        self.ref_name = ref_name
        self.interface = pm.interface
        self.mode = pm.mode # 'hfss' or 'gds'
        self.dict_instances[coor_sys] = self  # coor_sys is the name
        self.list_entities = []
        self.list_ports = []
        self.cursors = [] # tuple to escape list parsing

    def __call__(self, pos, ori):
        pos, ori = parse_entry(pos, ori)
        if len(pos)==2:
            pos.append(0)
        self.cursors.append((pos, ori))
        return self

    def __enter__(self):
        #1 We need to keep track of the entities created during the execution of a function
        self.list_entities = to_move(ModelEntity) # save "indentation level"
        self.list_ports = to_move(Port)
        return self

    def __exit__(self, *exc):
        #4 We move the entity that were created by the last function
        list_entities_new = find_last_list(ModelEntity.instances_to_move)
        list_ports_new = find_last_list(Port.instances_to_move)
        pos, angle = self.cursors[-1]

        #5 We move the entities_to_move with the right operation
        if len(list_entities_new)>0:
            self.rotate(list_entities_new, angle=angle)
            self.translate(list_entities_new, vector=[pos[0], pos[1], pos[2]])

        if len(list_ports_new)>0:
            self.rotate_port(list_ports_new, angle)
            self.translate_port(list_ports_new, vector=[pos[0], pos[1], pos[2]])

        #6 We empty a part of the 'to_move' lists
        if isinstance(self.list_entities[-1], list):
            a = self.list_entities.pop(-1)
            for entity in a:
                self.list_entities.append(entity)
        else:
            ModelEntity.instances_to_move = None

        if isinstance(self.list_ports[-1], list):
            a = self.list_ports.pop(-1)
            for entity in a:
                self.list_ports.append(entity)
        else:
            Port.instances_to_move = None

        self.cursors.pop(-1)
        return False

    def move_port(func):
        @wraps(func)
        def moved(*args, **kwargs):
            new_args = [args[0], args[1]]  # args[0] = chip, args[1] = name
            compteur = 0
            for i, argument in enumerate(args[2:]):
                if isinstance(argument, str) and (argument in Port.dict_instances):
                    #  if argument is the sting representation of the port
                    new_args.append(Port.dict_instances[argument])
                    compteur+=1
                elif isinstance(argument, Port):
                    #  it the argument is the port itself
                    new_args.append(argument)
                    compteur+=1
                else:
                    new_args.append(argument)
            return func(*new_args, **kwargs)
        return moved

    def port(self, name, pos, ori, widths, subnames, layers, offsets, constraint_port):
        name = check_name(Port, name)
        if constraint_port:
            pos, ori = parse_entry(pos, ori)
            offset=0
            width=50e-6  # 50um
            points = [(0, offset+width/2),
                      (width/3, offset),
                      (0, offset-width/2)]
            self.polyline_2D(points, name='_'+name, layer=layer_PORT, nonmodel=True)
        else:
            pos, ori, widths, offsets = parse_entry(pos, ori, widths, offsets)
            for ii in range(len(widths)):
                width = widths[ii]
                offset = offsets[ii]
                points = [(0, offset+width/2),
                          (width/3, offset),
                          (0, offset-width/2)]
                self.polyline_2D(points, name='_'+name+'_'+subnames[ii], layer=layer_PORT, nonmodel=True)

        return Port(name, pos, ori, widths, subnames, layers, offsets, constraint_port)

    def translate_port(self, ports, vector):
        for port in ports:
            port.pos = port.pos+Vector(vector)

    def rotate_port(self, ports, angle):
        if isinstance(angle, list):
            if len(angle)==2:
                new_angle= np.math.atan2(np.linalg.det([[1,0],angle]),np.dot([1,0],angle))
                new_angle= new_angle/np.pi*180
            else:
                raise Exception("angle should be either a float or a 2-dim array")
        else :
            new_angle=angle
        rad = new_angle/180*np.pi
        rotate_matrix = np.array([[np.cos(rad) ,np.sin(-rad)],[np.sin(rad) ,np.cos(rad)]])
        for port in ports:
            port.ori = rotate_matrix.dot(port.ori)
            posx = port.pos[0]*np.cos(rad)+port.pos[1]*np.sin(-rad)
            posy = port.pos[0]*np.sin(rad)+port.pos[1]*np.cos(rad)
            port.pos = Vector([posx, posy])

    @move_port
    def draw_cable(self, name, *ports, fillet="0.3mm", is_bond=False, to_meander=[[]], meander_length=0, meander_offset=0, is_mesh=False, reverse_adaptor=False):
        '''
        Draws a CPW transmission line along the ports

        Be careful: if the ports are facing eachother and offset by a small distance, sometimes the cable cannot be drawn.

        Inputs:
        -------
        name: (string) base-name of object, draws 'name_adaptor' etc
        iIn: (tuple) input port
        iOut: (tuple) output port
        iMaxfillet: (float), maximum fillet radius
        reverseZ: performs a mirror operation along Z --> useful only when the thickening operation goes in the wrong direction

        '''
        #TODO change the format of the arguments
        meander_length, meander_offset, fillet = parse_entry(meander_length, meander_offset, fillet)

        # to_meander should be a list of list
        # meander_length, meander_offset should be lists
        if not isinstance(to_meander[0], list):
            to_meander = [to_meander]
        if not isinstance(meander_length, list):
            meander_length = [meander_length]*len(to_meander)
        if not isinstance(meander_offset, list):
            meander_offset = [meander_offset]*len(to_meander)

        ports = list(ports)

        indent_level = find_corresponding_list(ports[0], Port.instances_to_move)
        if indent_level is not None:
            if indent_level:  # the found list
                for port in ports:
                    if not port in indent_level:
                        msg = 'Trying to connect ports from different \
                                indentation levels: port %s'%(port.name)
                        raise IndentationError(msg)


        # asserts neither in nor out port are constraint_ports
        if ports[0].constraint_port and ports[-1].constraint_port:
            raise ValueError('At least the first (%s) or last port (%s) \
                             should define the port parameters'%(ports[0].name,
                                                             ports[-1].name))
        elif ports[0].constraint_port:
            ports[0].widths = ports[-1].widths
            ports[0].offsets = ports[-1].offsets
            ports[0].layers = ports[-1].layers
            ports[0].subnames = ports[-1].subnames
            ports[0].N = ports[-1].N
            ports[0].constraint_port = False  # ports[0] is now defined
        elif ports[-1].constraint_port:
            ports[-1].widths = ports[0].widths
            ports[-1].offsets = ports[0].offsets
            ports[-1].layers = ports[0].layers
            ports[-1].subnames = ports[0].subnames
            ports[-1].N = ports[0].N
            ports[-1].constraint_port = False  # ports[-1] is now defined
        else:
            pass
        # recursive approach if there are intermediate non_constraint ports

        cable_portion = 0
        _ports = [ports[0]]
        for port in ports[1:-1]:
            _ports.append(port)
            if not port.constraint_port:
                print(to_meander[cable_portion])
                print(meander_length[cable_portion])
                self.draw_cable(name+'_%d'%cable_portion, *_ports,
                                fillet=fillet, is_bond=is_bond,
                                to_meander=[to_meander[cable_portion]],
                                meander_length=[meander_length[cable_portion]],
                                meander_offset=[meander_offset[cable_portion]],
                                is_mesh=is_mesh,
                                reverse_adaptor=reverse_adaptor)
                cable_portion += 1
                _ports = [port.r]

        if cable_portion != 0:
            name = name+'_%d'%cable_portion
            to_meander = [to_meander[cable_portion]]
            meander_length = [meander_length[cable_portion]]
            meander_offset = [meander_offset[cable_portion]]
        _ports.append(ports[-1])


        # at this stage first and last port are not constraint_port and all
        # intermediate port should be constraint_port

        ports = _ports

        # find and plot adaptor geometry
        if reverse_adaptor:
            points, length_adaptor = ports[-1].compare(ports[0], self.pm)
            index_modified = -1
        else:
            points, length_adaptor = ports[0].compare(ports[-1], self.pm)
            index_modified = 0

        # plot adaptors
        for jj, pts in enumerate(points):
            self.polyline_2D(pts, name=ports[index_modified].name+'_'+ports[index_modified].subnames[jj]+'_adapt', layer=ports[index_modified].layers[jj])

        # define the constraint_port parameters
        for port in ports[1:-1]:
            port.widths = ports[0].widths
            port.offsets = ports[0].offsets
            port.layers = ports[0].layers
            port.subnames = ports[0].subnames
            port.N = ports[0].N

        # find all intermediate paths
        total_path = None
        for ii in range(len(ports)-1):
            path = Path(name, ports[ii], ports[ii+1], fillet)
            if total_path is None:
                total_path = path
            else:
                total_path += path
            ports[ii+1] = ports[ii+1].r # reverse the last port

        total_path.clean()

        # do meandering
        total_path.meander(to_meander[0], meander_length[0], meander_offset[0])

        total_path.clean()
        # plot cable
        self.path_2D(total_path.points, total_path.port_in, total_path.fillet,
                  name=name)

        # if bond plot bonds
        if is_bond:
            self.draw_bond(total_path.to_bond(), *ports[0].bond_params(), name=name+'_wb')

        ports[0].revert()
        ports[-1].revert()

        length = total_path.length() + length_adaptor
        print('Cable "%s" length = %.3f mm'%(name, length*1000))
        return length

    def draw_bond(self, to_bond, ymax, ymin, name='wb', min_dist='0.5mm'):
        # to_bond list of segments
        ymax, ymin, min_dist = parse_entry(ymax, ymin, min_dist)

        min_dist = val(min_dist)
        n_segments = len(to_bond)
        jj=0
        bond_number = 0
        while jj<n_segments:
            elt = to_bond[jj]
            A = elt[0]
            B = elt[1]
            if jj+1 < n_segments and to_bond[jj+1][0] == B:
                B = to_bond[jj+1][1]
                jj+=1
            val_BA = val(B-A)
            ori = way(val_BA)
            length = Vector(val_BA).norm()
            n_bond = int(length/_val(min_dist))+1
            spacing = (B-A).norm()/n_bond
            pos = A+ori*spacing/2
            for ii in range(n_bond):
                self.wirebond_2D(pos, ori, ymax, ymin, layer=layer_Default, name=name+'_%d'%(bond_number))
                bond_number += 1
                pos = pos + ori*spacing
            jj+=1
