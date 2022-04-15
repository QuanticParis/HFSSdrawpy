from functools import wraps

import numpy as np
from sympy import posify

from ..parameters import DEFAULT, MASK, MESH, PORT
from ..path_finding.path_finder import Path
from ..path_finding.path_finder_auto import PathAuto
from ..utils import (
    Vector,
    check_name,
    equal_float,
    find_corresponding_list,
    find_last_list,
    find_penultimate_list,
    parse_entry,
    val,
    way,
)
from .entity import Entity
from .modeler import Modeler
from .port import Port


class BodyMover:
    def __init__(self, body):

        self.body = body
        self.id = np.random.rand()

    def __enter__(self):
        # 1 We need to keep track of the entities created during the execution of a function
        if self.body.entities_to_move is None:
            self.body.entities_to_move = []
        else:
            find_last_list(self.body.entities_to_move).append([])

        if self.body.ports_to_move is None:
            self.body.ports_to_move = []
        else:
            find_last_list(self.body.ports_to_move).append([])

    def __exit__(self, *exc):

        # 4 We move the entity that were created by the last function
        list_entities_new = find_last_list(self.body.entities_to_move)
        list_ports_new = find_last_list(self.body.ports_to_move)
        pos, angle = self.body.cursors[-1]

        # 5 We move the entities_to_move with the right operation
        if len(list_entities_new) > 0:
            self.body.rotate(list_entities_new, angle=angle)
            self.body.translate(list_entities_new, vector=[pos[0], pos[1], pos[2]])

        if len(list_ports_new) > 0:
            Port.rotate_ports(list_ports_new, angle)
            Port.translate_ports(list_ports_new, vector=[pos[0], pos[1], pos[2]])

        # 6 We empty a part of the 'to_move' lists
        penultimate_entity_list = find_penultimate_list(self.body.entities_to_move)
        if penultimate_entity_list:
            a = penultimate_entity_list.pop(-1)
            for entity in a:
                penultimate_entity_list.append(entity)
        else:
            self.body.entities_to_move = None

        penultimate_port_list = find_penultimate_list(self.body.ports_to_move)
        if penultimate_port_list:
            a = penultimate_port_list.pop(-1)
            for entity in a:
                penultimate_port_list.append(entity)
        else:
            self.body.ports_to_move = None

        self.body.cursors.pop(-1)
        return False


class Body(Modeler):

    dict_instances = {}

    def __init__(self, pm=None, name=None, rel_coor=None, ref_name="Global"):  # network
        # Note: for now coordinate systems are not reactualized at each run
        if rel_coor is None:
            rel_coor = [[0, 0, 0], [1, 0, 0], [0, 1, 0]]  # origin  # new_x  # new_y
        else:
            rel_coor = parse_entry(rel_coor)

        pm.interface.create_coor_sys(coor_sys=name, rel_coor=rel_coor, ref_name=ref_name)

        self.pm = pm
        self.name = name
        self.rel_coor = rel_coor
        self.ref_name = ref_name
        self.interface = pm.interface
        self.mode = pm.mode  # 'hfss' or 'gds'
        self.dict_instances[name] = self
        self.entities = {DEFAULT: []}  # entities sorted by layer
        self.cursors = []  # tuple to escape list parsing
        self.ports_to_move = None
        self.entities_to_move = None
        self.is_mask = pm.is_mask
        self.gap_mask = parse_entry(pm.gap_mask)
        self.is_litho = pm.is_litho

        pm.bodies.append(self)

    def __call__(self, pos, ori=[1, 0]):
        pos, ori = parse_entry(pos, ori)
        if len(pos) == 2:
            pos.append(0)
        self.cursors.append((pos, ori))
        return BodyMover(self)

    # def __enter__(self):
    #     print("enter(")
    #     #1 We need to keep track of the entities created during the execution of a function
    #     if self.entities_to_move is None:
    #         self.entities_to_move = []
    #     else:
    #         find_last_list(self.entities_to_move).append([])

    #     if self.ports_to_move is None:
    #         self.ports_to_move = []
    #     else:
    #         find_last_list(self.ports_to_move).append([])

    #     return self

    # def __exit__(self, *exc):
    #     print(")exit")
    #     #4 We move the entity that were created by the last function
    #     list_entities_new = find_last_list(self.entities_to_move)
    #     list_ports_new = find_last_list(self.ports_to_move)
    #     pos, angle = self.cursors[-1]

    #     #5 We move the entities_to_move with the right operation
    #     if len(list_entities_new)>0:
    #         self.rotate(list_entities_new, angle=angle)
    #         self.translate(list_entities_new, vector=[pos[0], pos[1], pos[2]])

    #     if len(list_ports_new)>0:
    #         Port.rotate_ports(list_ports_new, angle)
    #         Port.translate_ports(list_ports_new, vector=[pos[0], pos[1], pos[2]])

    #     #6 We empty a part of the 'to_move' lists
    #     penultimate_entity_list = find_penultimate_list(self.entities_to_move)
    #     if penultimate_entity_list:
    #         a = penultimate_entity_list.pop(-1)
    #         for entity in a:
    #             penultimate_entity_list.append(entity)
    #     else:
    #         self.entities_to_move = None

    #     penultimate_port_list = find_penultimate_list(self.ports_to_move)
    #     if penultimate_port_list:
    #         a = penultimate_port_list.pop(-1)
    #         for entity in a:
    #             penultimate_port_list.append(entity)
    #     else:
    #         self.ports_to_move = None

    #     self.cursors.pop(-1)
    #     return False

    def set_body(func):
        """
        Defines a wrapper/decorator which allows the user to always work in the coordinate system of the chosen chip.
        """

        @wraps(func)
        def updated(*args, **kwargs):
            args[0].interface.set_coor_sys(args[0].name)
            if not "layer" in kwargs:
                kwargs["layer"] = DEFAULT
            return func(*args, **kwargs)

        return updated
    
    def current_pos(self):
        """
        Compute the current value of the position of the cursor
        Assume no rotation is done ! (e.g. location of junctions in dosetest)
        TODO do full change of basis
        """
        pos_x = 0
        pos_y = 0
        for coor in self.cursors:
            pos_x += val(coor[0][0])
            pos_y += val(coor[0][1])
        return pos_x, pos_y

    ### Basic drawing

    @set_body
    def box(self, pos, size, name="box_0", **kwargs):
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
        pos, size = parse_entry(pos, size)
        name = check_name(Entity, name)
        kwargs["name"] = name
        self.interface.box(pos, size, **kwargs)
        return Entity(3, self, **kwargs)

    @set_body
    def box_center(self, pos, size, name="box_0", **kwargs):
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
        pos, size = parse_entry(pos, size)
        pos = [p - s / 2 for p, s in zip(pos, size)]
        return self.rect(pos, size, name=name, **kwargs)

    @set_body
    def cylinder(self, pos, radius, height, axis, segments=0, name="cylinder_0", **kwargs):
        pos, radius, height = parse_entry(pos, radius, height)
        name = check_name(Entity, name)
        kwargs["name"] = name
        self.interface.cylinder(pos, radius, height, axis, segments, **kwargs)
        return Entity(3, self, **kwargs)

    @set_body
    def cone(self, pos, radius1, radius2, height, axis, name="cone_0", **kwargs):
        pos, radius1, radius2, height = parse_entry(pos, radius1, radius2, height)
        name = check_name(Entity, name)
        kwargs["name"] = name
        self.interface.cone(pos, radius1, radius2, height, axis, **kwargs)
        return Entity(3, self, **kwargs)

    @set_body
    def sphere(self, pos, radius, name="sphere_0", **kwargs):
        pos, radius = parse_entry(pos, radius)
        name = check_name(Entity, name)
        kwargs["name"] = name
        self.interface.sphere(pos, radius, **kwargs)
        return Entity(3, self, **kwargs)

    @set_body
    def torus(self, pos, majorradius, minorradius, axis, name="torus_0", **kwargs):
        pos, majorradius, minorradius = parse_entry(pos, majorradius, minorradius)
        name = check_name(Entity, name)
        kwargs["name"] = name
        self.interface.torus(pos, majorradius, minorradius, axis, **kwargs)
        return Entity(3, self, **kwargs)

    @set_body
    def disk(self, pos, radius, axis, name="disk_0", **kwargs):
        pos, radius = parse_entry(pos, radius)
        name = check_name(Entity, name)
        kwargs["name"] = name
        if self.mode == "gds":
            pos = val(pos)
            radius = val(radius)
        self.interface.disk(pos, radius, axis, **kwargs)
        return Entity(2, self, **kwargs)

    @set_body
    def polyline(self, points, closed=True, name="polyline_0", **kwargs):
        points = parse_entry(points)
        name = check_name(Entity, name)
        kwargs["name"] = name
        i = 0
        while i < len(points[:-1]):
            points_equal = [
                equal_float(val(p0), val(p1)) for p0, p1 in zip(points[i], points[i + 1])
            ]
            if all(points_equal):
                points.pop(i)
                print("Warning: Delete two coinciding points on a polyline2D")
            else:
                i += 1
        if self.mode == "gds":
            points = val(points)
        self.interface.polyline(points, closed, **kwargs)
        dim = closed + 1
        return Entity(dim, self, **kwargs)

    @set_body
    def rect(self, pos, size, name="rect_0", **kwargs):
        pos, size = parse_entry(pos, size)
        name = check_name(Entity, name)
        kwargs["name"] = name
        if self.mode == "gds":
            pos = val(pos)
            size = val(size)
        self.interface.rect(pos, size, **kwargs)
        return Entity(2, self, **kwargs)

    @set_body
    #####draw arrays of rectangles with dimension (colums x row) with spacing given by a list [x_spacing,y_spacing]
    def rect_array(self, pos, size, columns, rows, spacing, name="rect_array_0", **kwargs):

        pos, size, spacing = parse_entry(Vector(pos), Vector(size), spacing)

        if self.mode == "gds":
            name = check_name(Entity, name)
            kwargs["name"] = name
            pos = val(pos)
            size = val(size)
            spacing = val(spacing)
            self.interface.rect_array(pos, size, columns, rows, spacing, **kwargs)
            return Entity(2, self, **kwargs)
        else:
            columns, rows = parse_entry(columns, rows)
            # rect = self.rect(pos-size/2, size, name, **kwargs)
            rect = self.rect_center(pos, size, name, **kwargs)
            self.duplicate_along_line(
                rect, [0, spacing[1], 0], n=rows)
            self.duplicate_along_line(
                rect, [spacing[0], 0, 0], n=columns)
            return rect

    @set_body
    def rect_center(self, pos, size, name="rect_0", **kwargs):
        pos, size = parse_entry(pos, size)
        pos = [p - s / 2 for p, s in zip(pos, size)]
        return self.rect(pos, size, name=name, **kwargs)

    @set_body
    def wirebond(self, pos, ori, ymax, ymin, name="wb_0", **kwargs):
        pos, ymax, ymin = parse_entry(pos, ymax, ymin)
        pos, ori = Vector(pos), Vector(ori)
        name = check_name(Entity, name)
        kwargs["name"] = name
        if self.mode == "gds":
            pos, ori, ymax, ymin = val(pos, ori, ymax, ymin)
            self.interface.wirebond(pos, ori, ymax, ymin, **kwargs)
            kwargs["name"] = name + "a"
            entity_a = Entity(2, self, **kwargs)
            kwargs["name"] = name + "b"
            entity_b = Entity(2, self, **kwargs)
            return entity_a, entity_b
        else:
            self.interface.wirebond(pos, ori, ymax, ymin, **kwargs)
            return Entity(3, self, **kwargs)

    @set_body
    def text(self, pos, size, text, angle=0, horizontal=True, name="text_0", **kwargs):
        """
        Place text in GDS layout.

        size (number) – Height of the character. The width of a character and the distance between characters are
            this value multiplied by 5 / 9 and 8 / 9, respectively. For vertical text, the distance is multiplied by 11 / 9.
        text (string) – The text to be converted in geometric objects.
        angle (number) – The angle of rotation of the text.
        horizontal (bool) – If True, the text is written from left to right; if False, from top to bottom.
        name (string) – The GDSII layer number for these elements.
        """
        if self.mode == "gds":
            pos, size = parse_entry(pos, size)
            name = check_name(Entity, name)
            kwargs["name"] = name
            pos = val(pos)
            size = val(size)
            self.interface.text(pos, size, text, angle, horizontal, **kwargs)
            return Entity(2, self, **kwargs)
        else:
            pass

    @set_body
    def path(self, points, port, fillet, drop_mask=False, name="path_0", **kwargs):
        # fillet should be either 0 or larger than half of the port width
        name = check_name(Entity, name)
        kwargs["name"] = name
        model_entities = []
        if self.mode == "gds":
            points = val(points)
            fillet = val(fillet)
            _port = port.val(drop_mask=drop_mask)
            

            # TODO, this is a dirty fixe cause of Vector3D

            # due to the 3D vector implementation
            points_2D = []
            for point in points:
                points_2D.append([point[0], point[1]])

            if fillet == 0:
                names, layers = self.interface.path(
                    points_2D, _port, fillet, name=name, corner="natural"
                )
            else:
                names, layers = self.interface.path(
                    points_2D, _port, fillet, name=name, corner="circular bend"
                )

            for name, layer in zip(names, layers):
                kwargs["layer"] = layer
                kwargs["name"] = name
                model_entities.append(Entity(2, self, **kwargs))
        elif self.mode == "hfss":
            # check that port is at the BEGINNING of the path (hfss only)
            ori = Vector(port.ori)
            pos = Vector(port.pos)
            path_entity = self.polyline(points, closed=False, name=name)
            path_entity.fillet(fillet)

            for ii in range(port.N):
                offset = port.offsets[ii]
                width = port.widths[ii]
                subname = port.subnames[ii]
                layer = port.layers[ii]
                points_starter = [
                    Vector(0, offset + width / 2).rot(ori) + pos,
                    Vector(0, offset - width / 2).rot(ori) + pos,
                ]
                entity = self.polyline(
                    points_starter, closed=False, name=name + "_" + subname, layer=layer
                )
                path_name = name + "_" + subname + "_path"
                current_path_entity = path_entity.copy(new_name=path_name)
                self.interface.sweep_along_path(entity, current_path_entity)
                current_path_entity.delete()
                model_entities.append(entity)

            path_entity.delete()

        return model_entities
    
    @set_body
    def duplicate_along_line(self, entity, vec, n=2, new_obj=False, duplicate_assign=False, **kwargs):

        if(self.mode=='hfss'):
            vec, n = parse_entry(vec, n)
            self.interface.duplicate_along_line(
                entity, vec, n=n, new_obj=new_obj, duplicate_assign=duplicate_assign
            )
            return entity
        else:
            pass

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
            new_args = [args[0]]  # args[0] = chip, args[1] = name
            for i, argument in enumerate(args[1:]):
                if isinstance(argument, str) and (argument in Port.dict_instances):
                    #  if argument is the sting representation of the port
                    new_args.append(Port.dict_instances[argument])
                elif isinstance(argument, Port):
                    #  it the argument is the port itself
                    new_args.append(argument)
                else:
                    new_args.append(argument)
            return func(*new_args, **kwargs)

        return moved

    def port(self, widths=None, subnames=None, layers=None, offsets=0, name="port_0"):
        """
        Creates a port and draws a small triangle for each element of the port

        Parameters
        ----------
        name : str, optional
            Name of the port.
        widths : float, 'VariableString' or list, optional
            Width of the different elements of the port. If None, assumes the
            creation of a constraint_port. The default is None.
        subnames : str or list, optional
            The cable's parts will be name cablename_subname. If None simply
            numbering the cable parts. The default is None.
        layers : int or list, optional
            Each layer is described by an int that is a python constant that one
            should import. If None, layer is the DEFAULT The default is
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
                layers = [DEFAULT] * N
            elif not isinstance(layers, list):
                layers = [layers] * N

            if not isinstance(offsets, list):
                offsets = [offsets] * N

            widths, offsets = parse_entry(widths, offsets)

        # actual drawing and creation

        pos = [0, 0]
        ori = [1, 0]

        name = check_name(Port, name)

        if constraint_port:
            pos, ori = parse_entry(pos, ori)
            offset = 0
            width = 50e-6  # 50um
            points = [
                (0, offset + width / 2),
                (width / 3, offset),
                (0, offset - width / 2),
            ]
            self.polyline(points, name="_" + name, layer=PORT, nonmodel=True)
        else:
            pos, ori, widths, offsets = parse_entry(pos, ori, widths, offsets)
            for ii in range(len(widths)):
                width = widths[ii]
                offset = offsets[ii]
                points = [
                    (0, offset + width / 2),
                    (width / 3, offset),
                    (0, offset - width / 2),
                ]
                self.polyline(
                    points,
                    name="_" + name + "_" + subnames[ii],
                    layer=PORT,
                    nonmodel=True,
                )

        result = Port(self, name, pos, ori, widths, subnames, layers, offsets, constraint_port)
        return [result]


    @move_port
    def draw_cable(
        self,
        *ports,
        fillet="0.3mm",
        is_bond=False,
        to_meander=None,
        meander_length=0,
        meander_offset=0,
        reverse_adaptor=False,
        slope=0.5,
        name="cable_0",
        mesh_size=None,
        slanted=False,
        drop_mask=False
    ):
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
            DESCRIPTION. The default is None.
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

        meander_length, meander_offset, fillet, mesh_size = parse_entry(
            meander_length, meander_offset, fillet, mesh_size
        )

        ports = list(ports)

        do_not_beyong = [port.name for port in ports if port.body != self]
        if do_not_beyong:
            raise ValueError("%s ports do not beyond to %s" % (do_not_beyong, self))
        
        # check the indentation level of connected ports
        # you should not try to connect ports different with body(...): levels
        indent_level = find_corresponding_list(ports[0], self.ports_to_move)
        if indent_level is not None:
            if indent_level:  # the found list
                for port in ports:
                    if not port in indent_level:
                        msg = (
                            "Trying to connect ports from different \
                                indentation levels: port %s"
                            % (port.name)
                        )
                        raise IndentationError(msg)

        # asserts neither in nor out port are constraint_ports
        if ports[0].constraint_port and ports[-1].constraint_port:
            raise ValueError(
                "At least the first (%s) or last port (%s) \
                             should define the port parameters"
                % (ports[0].name, ports[-1].name)
            )
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

        if not slanted:  # slanted cables are specific
        
            # first determine the intermediate non constrain ports
            _ports = [[ports[0]]]
            for port in ports[1:-1]:
                if port.constraint_port:
                    _ports[-1].append(port)
                else:
                    _ports[-1].append(port)
                    _ports.append([port])
            _ports[-1].append(ports[-1])
            
            if len(_ports)>1:
                # to_meander should be a list of list: a list for each cable portion
                # meander_length, meander_offset should be lists
                if to_meander is None:
                    to_meander = [[]]*len(_ports)
                if not isinstance(to_meander[0], list):
                    to_meander = [to_meander]
                if not isinstance(meander_length, list):
                    meander_length = [meander_length] * len(_ports)
                if not isinstance(meander_offset, list):
                    meander_offset = [meander_offset] * len(_ports)
                if len(to_meander)!=len(_ports):
                    raise Exception("""to_meander should be a list of list, 
                                    the outer list length being the number of 
                                    cable portions""")
                length = 0
                for ii, port_list in enumerate(_ports):
                    length += self.draw_cable(
                            *port_list,
                            fillet=fillet,
                            is_bond=is_bond,
                            to_meander=to_meander[ii],
                            meander_length=meander_length[ii],
                            meander_offset=meander_offset[ii],
                            reverse_adaptor=bool(ii),
                            mesh_size=mesh_size,
                            slope=slope,
                            name=name + "_%d" % ii
                            )
                print('Total Cable "%s" length = %.3f mm' % (name, 
                                                             length * 1000))
                return length
            else:
                # to_meander should be a list since at this stage we are within
                # one cable portion
                ports = _ports[0]
                if to_meander is None:
                    to_meander = []
                if isinstance(meander_length, list):
                    meander_length = meander_length[0]
                if isinstance(meander_offset, list):
                    meander_offset = meander_offset[0]
                
    
                # find and plot adaptor geometry
                if reverse_adaptor:
                    points, length_adaptor = ports[-1].r.compare(ports[0], self.pm, slope=slope)
                    index_modified = -1
                else:
                    points, length_adaptor = ports[0].compare(ports[-1], self.pm, slope=slope)
                    index_modified = 0
    
                # plot adaptors
                for jj, pts in enumerate(points):
                    adaptor = self.polyline(
                                pts,
                                name=ports[index_modified].name
                                + "_"
                                + ports[index_modified].subnames[jj]
                                + "_adapt",
                                layer=ports[index_modified].layers[jj],
                                )
                    if adaptor.layer==MESH and mesh_size is not None:
                        adaptor.assign_mesh_length(mesh_size)
                        
    
                # define the constraint_port parameters
                for port in ports[1:-1]:
                    port.widths = ports[0].widths
                    port.offsets = ports[0].offsets
                    port.layers = ports[0].layers
                    port.subnames = ports[0].subnames
                    port.N = ports[0].N
    
                # find all intermediate paths
                total_path = None
                for ii in range(len(ports) - 1):
                    path = Path(name, ports[ii], ports[ii + 1], fillet)
                    if total_path is None:
                        total_path = path
                    else:
                        total_path += path
    
                total_path.clean()
    
                # do meandering (not supported for even partially slanted cables)
                if not total_path.is_slanted:
                    total_path.meander(to_meander, meander_length, meander_offset)
    
                total_path.clean()
                # plot cable
                cable = self.path(total_path.points, total_path.port_in, total_path.fillet, name=name, drop_mask=drop_mask)
    
                # assign mesh_size to the mesh layer in the new cable
                if mesh_size is not None:
                    for entity in cable:
                        if entity.layer == MESH:
                            entity.assign_mesh_length(mesh_size)
    
                # if bond plot bonds
                if is_bond:
                    self.draw_bond(total_path.to_bond(), *ports[0].bond_params(), name=name + "_wb")
    
                length = total_path.length() + length_adaptor
                print('Cable "%s" length = %.3f mm' % (name, length * 1000))
                return length
        else:
            if len(ports) > 2:
                raise Exception("Constraint_ports are not supported with slanted cables")
            path = Path(name, ports[0], ports[1], fillet, is_slanted=True)
            cable = self.path(path.points, path.port_in, path.fillet, name=name)
            # assign mesh_size to the mesh layer in the new cable
            if mesh_size is not None:
                for entity in cable:
                    if entity.layer == MESH:
                        entity.assign_mesh_length(mesh_size)
            if is_bond:
                raise Exception("Bonding is not supported with slanted cables")
                
    @move_port
    def draw_cable_auto(
        self,
        *ports,
        fillet="0.3mm",
        is_bond=False,
        mid=0.5,
        meanders=None,
        reverse_adaptor=False,
        slope=0.5,
        name="cable_0",
        mesh_size=None,
        drop_mask=False, 
        target_length=None,
        editable_in_hfss=False # obscur parameter to check that design will adjust fine in HFSS
    ):
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

        fillet, mesh_size = parse_entry(fillet, mesh_size)
        mid, target_length = parse_entry(mid, target_length)
        
        # TODO, make mid a list in the case there are constraints port in the middle
        if meanders is None:
            meanders = []
            
        ports = list(ports)

        do_not_beyong = [port.name for port in ports if port.body != self]
        if do_not_beyong:
            raise ValueError("%s ports do not beyond to %s" % (do_not_beyong, self))
        
        # check the indentation level of connected ports
        # you should not try to connect ports different with body(...): levels
        indent_level = find_corresponding_list(ports[0], self.ports_to_move)
        if indent_level is not None:
            if indent_level:  # the found list
                for port in ports:
                    if not port in indent_level:
                        msg = (
                            "Trying to connect ports from different \
                                indentation levels: port %s"
                            % (port.name)
                        )
                        raise IndentationError(msg)

        # asserts neither in nor out port are constraint_ports
        if ports[0].constraint_port and ports[-1].constraint_port:
            raise ValueError(
                "At least the first (%s) or last port (%s) \
                             should define the port parameters"
                % (ports[0].name, ports[-1].name)
            )
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

        # if not slanted:  # slanted cables are specific
        
        # first determine the intermediate non constrain ports
        _ports = [[ports[0]]]
        for port in ports[1:-1]:
            if port.constraint_port:
                _ports[-1].append(port)
            else:
                _ports[-1].append(port)
                _ports.append([port])
        _ports[-1].append(ports[-1])
        
        if len(_ports)>1:
            raise Exception('Too many ports : \n Please make several cables or use constraint ports')
        else:
            # to_meander should be a list since at this stage we are within
            # one cable portion
            ports = _ports[0]

            # find and plot adaptor geometry
            if reverse_adaptor:
                points, length_adaptor = ports[-1].r.compare(ports[0], self.pm, slope=slope)
                index_modified = -1
            else:
                points, length_adaptor = ports[0].compare(ports[-1], self.pm, slope=slope)
                index_modified = 0

            # plot adaptors
            for jj, pts in enumerate(points):
                adaptor = self.polyline(
                            pts,
                            name=ports[index_modified].name
                            + "_"
                            + ports[index_modified].subnames[jj]
                            + "_adapt",
                            layer=ports[index_modified].layers[jj],
                            )
                if adaptor.layer==MESH and mesh_size is not None:
                    adaptor.assign_mesh_length(mesh_size)
                    

            # define the constraint_port parameters
            for port in ports[1:-1]:
                port.widths = ports[0].widths
                port.offsets = ports[0].offsets
                port.layers = ports[0].layers
                port.subnames = ports[0].subnames
                port.N = ports[0].N

            # find all intermediate paths
            total_path = None
            for ii in range(len(ports) - 1):
                path = PathAuto(name, ports[ii], ports[ii + 1], fillet,
                            mid=mid)
                if total_path is None:
                    total_path = path
                else:
                    total_path += path    
            total_path.analyse_path()     
            
            if target_length is None:
                total_path.meander(meanders, target_length=target_length)
            else:
                total_path.meander(meanders, target_length=target_length-length_adaptor, editable_in_hfss=editable_in_hfss)
            
            total_path.analyse_path()   
            
            # check cable length
            length = val(total_path.length+length_adaptor)
            if target_length is not None:
                if not equal_float(val(total_path.length+length_adaptor), val(target_length)):
                    print('/!\ Cable "%s" length = %.3f (!= %.3f) mm /!\ ' %(name, length * 1000, val(target_length)*1000))
            else:
                print('Cable "%s" length = %.3f mm' % (name, length * 1000))
                
            # plot cable
            cable = self.path(total_path.path, total_path.port_in, total_path.fillet, name=name, drop_mask=drop_mask)

            # assign mesh_size to the mesh layer in the new cable
            if mesh_size is not None:
                for entity in cable:
                    if entity.layer == MESH:
                        entity.assign_mesh_length(mesh_size)
                        
            # if bond, plot bonds
            if is_bond:
                self.draw_bond(total_path.to_bond(), *ports[0].bond_params(), name=name + "_wb")

            return length

    def draw_bond(self, to_bond, ymax, ymin, min_dist="0.5mm", name="wb_0"):
        # to_bond list of segments
        ymax, ymin, min_dist = parse_entry(ymax, ymin, min_dist)

        min_dist = val(min_dist)
        n_segments = len(to_bond)
        jj = 0
        bond_number = 0
        while jj < n_segments:
            elt = to_bond[jj]
            A = elt[0]
            B = elt[1]
            if jj + 1 < n_segments and to_bond[jj + 1][0] == B:
                B = to_bond[jj + 1][1]
                jj += 1
            val_BA = val(B - A)
            ori = way(val_BA)
            length = Vector(val_BA).norm()
            n_bond = int(length / val(min_dist)) + 1
            spacing = (B - A).norm() / n_bond
            pos = A + ori * spacing / 2
            for ii in range(n_bond):
                self.wirebond(pos, ori, ymax, ymin, layer=PORT, 
                              name=name + "_%d" % (bond_number))
                bond_number += 1
                pos = pos + ori * spacing
            jj += 1
    
    def delete_layer(self, layer):
        if layer in self.entities.keys():
            entity_list = self.entities[layer]
            for ii in range(len(self.entities[layer])):
                entity_list[0].delete()
            self.entities.pop(layer)