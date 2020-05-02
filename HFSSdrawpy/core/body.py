from functools import wraps


from ..utils import Vector, \
                   parse_entry, \
                   check_name, \
                   find_last_list, \
                   find_corresponding_list, \
                   to_move, \
                   _val, \
                   val, \
                   entity_kwargs, equal_float, \
                   way
from .entity import Entity
from ..modeler import Modeler
from ..path_finding.path_finder import Path
from .port import Port

from ..parameters import layer_Default, layer_PORT

class Body(Modeler):

    dict_instances = {}
    def __init__(self, pm=None, name=None, rel_coor=None, ref_name='Global'): #network
        # Note: for now coordinate systems are not reactualized at each run
        if rel_coor is None:
            rel_coor = [[0, 0, 0],  # origin
                        [1, 0, 0],  # new_x
                        [0, 1, 0]]  # new_y
        else:
            rel_coor = parse_entry(rel_coor)

        pm.interface.create_coor_sys(coor_sys=name, rel_coor=rel_coor,
                                     ref_name=ref_name)

        self.pm = pm
        self.name = name
        self.rel_coor = rel_coor
        self.ref_name = ref_name
        self.interface = pm.interface
        self.mode = pm.mode # 'hfss' or 'gds'
        self.dict_instances[name] = self
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
        self.list_entities = to_move(Entity) # save "indentation level"
        self.list_ports = to_move(Port)
        return self

    def __exit__(self, *exc):
        #4 We move the entity that were created by the last function
        list_entities_new = find_last_list(Entity.instances_to_move)
        list_ports_new = find_last_list(Port.instances_to_move)
        pos, angle = self.cursors[-1]

        #5 We move the entities_to_move with the right operation
        if len(list_entities_new)>0:
            self.rotate(list_entities_new, angle=angle)
            self.translate(list_entities_new, vector=[pos[0], pos[1], pos[2]])

        if len(list_ports_new)>0:
            Port.rotate_ports(list_ports_new, angle)
            Port.translate_ports(list_ports_new, vector=[pos[0], pos[1], pos[2]])

        #6 We empty a part of the 'to_move' lists
        if len(self.list_entities) > 0 and isinstance(self.list_entities[-1],
                                                      list):
            a = self.list_entities.pop(-1)
            for entity in a:
                self.list_entities.append(entity)
        else:
            Entity.instances_to_move = None

        if len(self.list_ports) > 0 and isinstance(self.list_ports[-1], list):
            a = self.list_ports.pop(-1)
            for entity in a:
                self.list_ports.append(entity)
        else:
            Port.instances_to_move = None

        self.cursors.pop(-1)
        return False

    def set_body(func):
        """
        Defines a wrapper/decorator which allows the user to always work in the coordinate system of the chosen chip.
        """
        @wraps(func)
        def updated(*args, **kwargs):
            args[0].interface.set_coor_sys(args[0].name)
            return func(*args, **kwargs)
        return updated

    ### Basic drawings

    @set_body
    def box_corner_3D(self, pos, size, **kwargs):
        """
        Draws a 3D box based on the coordinates of its corner.

        Inputs:
        -------
        pos: [x,y,z] the coordinates of the corner of the rectangle in the euclidian basis.
        size: [lx,ly,lz] the dimensions of the rectangle

        **kwargs include layer and name
        Outputs:
        -------
        box: Corresponding 3D Model Entity
        """
        name = self.interface.box_corner_3D(pos, size, **kwargs)
        kwargs = entity_kwargs(kwargs, ['layer', 'nonmodel'])
        return Entity(name, 3, self, **kwargs)

    @set_body
    def box_center_3D(self, pos, size, **kwargs):
        """
        Draws a 3D box based on the coordinates of its center.

        Inputs:
        -------
        pos: [x,y,z] the coordinates of the center of the rectangle in the euclidian basis.
        size: [lx,ly,lz] the dimensions of the rectangle

        **kwargs include layer and name
        Outputs:
        -------
        box: Corresponding 3D Model Entity
        """
        name = self.interface.box_center_3D(pos, size, **kwargs)
        kwargs = entity_kwargs(kwargs, ['layer', 'nonmodel'])
        return Entity(name, 3, self, **kwargs)

    @set_body
    def cylinder_3D(self, pos, radius, height, axis, **kwargs):
        name = self.interface.cylinder_3D(pos, radius, height, axis, **kwargs)
        kwargs = entity_kwargs(kwargs, ['layer', 'nonmodel'])
        return Entity(name, 3, self, **kwargs)


    @set_body
    def disk_2D(self, pos, radius, axis, **kwargs):
        if self.mode=='gds':
            pos = val(pos)
            radius = val(radius)
        name = self.interface.disk_2D(pos, radius, axis, **kwargs)
        kwargs = entity_kwargs(kwargs, ['layer', 'nonmodel'])
        return Entity(name, 2, self, **kwargs)

    @set_body
    def polyline_2D(self, points, closed=True, **kwargs): # among kwargs, name should be given
        i = 0
        while i < len(points[:-1]):
            points_equal = [equal_float(val(p0),val(p1))
                            for p0, p1 in zip(points[i], points[i+1])]
            if all(points_equal):
                points.pop(i)
                print('Warning: Delete two coinciding points on a polyline2D')
            else:
                i+=1
        if self.mode=='gds':
            points = val(points)
        name = self.interface.polyline_2D(points, closed, **kwargs)
        dim = closed + 1
        kwargs = entity_kwargs(kwargs, ['layer', 'nonmodel'])
        return Entity(name, dim, self, **kwargs)

    @set_body
    def path_2D(self, points, port, fillet, **kwargs):
        name = kwargs['name']
        model_entities = []
        if self.mode == 'gds':
            points = val(points)
            fillet = val(fillet)
            _port = port.val()
            names, layers = self.interface.path(points, _port, fillet, name=name)
            for name, layer in zip(names, layers):
                kwargs = {'layer':layer}  # model by default for now
                model_entities.append(Entity(name, 2, self, **kwargs))
        elif self.mode == 'hfss':
            # check that port is at the BEGINNING of the path (hfss only)
            ori = port.ori
            pos = port.pos
            path_entity = self.polyline_2D(points, closed=False,
                                           name=name, layer=layer_Default)
            path_entity.fillet(fillet)

            for ii in range(port.N):
                offset = port.offsets[ii]
                width = port.widths[ii]
                subname = port.subnames[ii]
                layer = port.layers[ii]
                points_starter = [Vector(0, offset+width/2).rot(ori)+pos,
                                  Vector(0, offset-width/2).rot(ori)+pos]
                entity = self.polyline_2D(points_starter, closed=False,
                                          name=name+'_'+subname, layer=layer)
                path_name = name+'_'+subname+'_path'
                current_path_entity = path_entity.copy(new_name=path_name)
                self.interface._sweep_along_path(entity, current_path_entity)
                current_path_entity.delete()
                model_entities.append(entity)

            path_entity.delete()

        return model_entities

    @set_body
    def rect_corner_2D(self, pos, size, **kwargs):
        if self.mode=='gds':
            pos = val(pos)
            size = val(size)
        name = self.interface.rect_corner_2D(pos, size, **kwargs)
        kwargs = entity_kwargs(kwargs, ['layer', 'nonmodel'])
        return Entity(name, 2, self, **kwargs)

    @set_body
    def rect_center_2D(self, pos, size, **kwargs):
        if self.mode=='gds':
            pos = val(pos)
            size = val(size)
        name = self.interface.rect_center_2D(pos, size, **kwargs)
        kwargs = entity_kwargs(kwargs, ['layer', 'nonmodel'])
        return Entity(name, 2, self, **kwargs)

    @set_body
    def wirebond_2D(self, pos, ori, ymax, ymin, **kwargs):
        if self.mode=='gds':
            pos, ori, ymax, ymin = val(pos, ori, ymax, ymin)
            name_a, name_b = self.interface.wirebond_2D(pos, ori, ymax, ymin, **kwargs)
            kwargs = entity_kwargs(kwargs, ['layer', 'nonmodel'])
            return Entity(name_a, 2, self, **kwargs), \
                    Entity(name_b, 2, self, **kwargs)
        else:
            name = self.interface.wirebond_2D(pos, ori, ymax, ymin, **kwargs)
            kwargs = entity_kwargs(kwargs, ['layer', 'nonmodel'])
            return Entity(name, 3, self, **kwargs)

    ### Advanced methods

    def move_port(func):
        """
        This method is supposed to be used as a decorator that helps replace
        the deprecated class 'ConnectElt'.

        Parameters
        ----------
        func : TYPE
            DESCRIPTION.

        Returns
        -------
        TYPE
            DESCRIPTION.

        """
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

    def port(self, name, widths=None, subnames=None, layers=None, offsets=0):
        """
        Creates a port and draws a small triangle for each element of the port

        Parameters
        ----------
        name : str
            Name of the port.
        widths : float, 'VariableString' or list, optional
            Width of the different elements of the port. If None, assumes the
            creation of a constraint_port. The default is None.
        subnames : str or list, optional
            The cable's parts will be name cablename_subname. If None simply
            numbering the cable parts. The default is None.
        layers : int or list, optional
            Each layer is described by an int that is a python constant that one
            should import. If None, layer is the layer_Default The default is
            None.
        offsets : float, 'VariableString' or list, optional
            Describes the offset of the cable part wrt the center of the cable.
            The default is 0.

        Returns
        -------
        'Port'
            Returns a Port object

        """
        # logic for keyword arguments
        if widths is None:
            constraint_port = True
        else:
            constraint_port = False
            if not isinstance(widths, list):
                widths = [widths]

            N = len(widths)

            # default subnames
            if subnames is None:
                subnames = []
                for ii in range(N):
                    subnames.append(str(ii))
            elif not isinstance(subnames, list):
                subnames = [subnames]

            if layers is None:
                layers = [layer_Default]*N
            elif not isinstance(layers, list):
                layers = [layers]*N

            if not isinstance(offsets, list):
                offsets = [offsets]*N

            widths, offsets = parse_entry(widths, offsets)

        # actual drawing and creation

        pos = [0, 0]
        ori = [1, 0]

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

    @move_port
    def draw_cable(self, name, *ports, fillet="0.3mm", is_bond=False, to_meander=[[]], meander_length=0, meander_offset=0, is_mesh=False, reverse_adaptor=False):
        """


        Parameters
        ----------
        name : TYPE
            DESCRIPTION.
        *ports : TYPE
            DESCRIPTION.
        fillet : TYPE, optional
            DESCRIPTION. The default is "0.3mm".
        is_bond : TYPE, optional
            DESCRIPTION. The default is False.
        to_meander : TYPE, optional
            DESCRIPTION. The default is [[]].
        meander_length : TYPE, optional
            DESCRIPTION. The default is 0.
        meander_offset : TYPE, optional
            DESCRIPTION. The default is 0.
        is_mesh : TYPE, optional
            DESCRIPTION. The default is False.
        reverse_adaptor : TYPE, optional
            DESCRIPTION. The default is False.

        Raises
        ------
        IndentationError
            DESCRIPTION.
        ValueError
            DESCRIPTION.

        Returns
        -------
        length : TYPE
            DESCRIPTION.

        """

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
