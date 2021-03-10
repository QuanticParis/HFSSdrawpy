import numpy as np
import math

from ..utils import Vector, parse_entry, check_name, find_last_list, val


class Port:
    dict_instances = {}

    def __init__(
        self,
        body,
        name,
        pos,
        ori,
        widths,
        subnames,
        layers,
        offsets,
        constraint_port,
        key="name",
    ):
        if not (isinstance(key, Port) or key is None):
            name = check_name(self.__class__, name)
        self.name = name
        self.pos = Vector(pos)
        self.ori = Vector(ori)
        self.constraint_port = constraint_port
        self.save = None
        self.body = body
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

        if self.body.ports_to_move is not None:
            find_last_list(self.body.ports_to_move).append(self)
        if key == "name":  # normal initialisation
            self.dict_instances[name] = self

            # create a reversed version of the port that can be called by either
            # port.r or 'port_name.r'
            reversed_ori = -self.ori
            reversed_offsets = None
            if self.offsets is not None:
                reversed_offsets = []
                for ii in range(self.N):
                    reversed_offsets.append(-self.offsets[ii])
            self.r = Port(
                self.body,
                self.name + "_r",
                self.pos,
                reversed_ori,
                self.widths,
                self.subnames,
                self.layers,
                reversed_offsets,
                self.constraint_port,
                key=self,
            )

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
        Port.dict_instances = {}

    @classmethod
    def print_instances(cls):
        for instance_name in cls.dict_instances:
            print(instance_name)  # , cls.dict_instances[instance_name])

    def compare(self, other, pm, slope=0.5):
        points = []

        adapt_dist = pm.set_variable(1e-5, name=self.name + "_adapt")
        max_diff = 0
        for ii in range(self.N):
            if self.layers[ii] != other.layers[ii]:
                raise ValueError(
                    "Tried to connect ports form different \
                                 layers: %s != %s"
                    % (self.layers[ii], other.layers[ii])
                )
            width1 = self.widths[ii]
            width2 = other.widths[ii]

            offset1 = self.offsets[ii]
            offset2 = other.offsets[ii]

            if val(width1) != val(width2) or val(offset1) != val(offset2):
                # need adaptor
                points.append(
                    [
                        Vector(0, offset1 + width1 / 2).rot(self.ori) + self.pos,
                        Vector(adapt_dist, offset2 + width2 / 2).rot(self.ori) + self.pos,
                        Vector(adapt_dist, offset2 - width2 / 2).rot(self.ori) + self.pos,
                        Vector(0, offset1 - width1 / 2).rot(self.ori) + self.pos,
                    ]
                )
            max_diff = max(
                max_diff,
                abs(val(offset1 + width1 / 2 - (offset2 + width2 / 2))),
                abs(val(offset2 - width2 / 2 - (offset1 - width1 / 2))),
            )
        adapt_dist = pm.set_variable(max_diff / slope, name=self.name + "_adapt")

        if len(points) != 0:
            self.save = {
                "pos": self.pos,
                "widths": self.widths,
                "offsets": self.offsets,
            }
            self.pos = self.pos + Vector(max_diff / slope, 0).rot(self.ori)
            self.widths = other.widths
            self.offsets = [-offset for offset in other.offsets]

            self.r.save = {
                "pos": self.r.pos,
                "widths": self.r.widths,
                "offsets": self.r.offsets,
            }
            self.r.pos = self.pos
            self.r.widths = self.widths
            self.r.offsets = other.offsets
        return points, 2 * max_diff

    def val(self):
        _widths = []
        _offsets = []
        for ii in range(self.N):
            width = self.widths[ii]
            offset = self.offsets[ii]
            _widths.append(val(width))
            _offsets.append(val(offset))

        _pos = []
        for coor in self.pos:
            _pos.append(val(coor))
        _pos = Vector(_pos)

        _ori = []
        for coor in self.ori:
            _ori.append(val(coor))
        _ori = Vector(_ori)

        return Port(
            self.body,
            self.name,
            _pos,
            _ori,
            _widths,
            self.subnames,
            self.layers,
            _offsets,
            self.constraint_port,
            key=None,
        )

    def revert(self):
        if self.save is not None:
            self.pos = self.save["pos"]
            self.widths = self.save["widths"]
            self.offsets = self.save["offsets"]

            self.r.pos = self.save["pos"]
            self.r.widths = self.save["widths"]
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
            _y_max_val = val(self.offsets[ii] + self.widths[ii] / 2)
            _y_min_val = val(self.offsets[ii] - self.widths[ii] / 2)
            if _y_max_val > y_max_val:
                y_max = self.offsets[ii] + self.widths[ii] / 2
                y_max_val = _y_max_val
            if _y_min_val < y_min_val:
                y_min = self.offsets[ii] - self.widths[ii] / 2
                y_min_val = _y_min_val
        return y_max, y_min

    @staticmethod
    def translate_ports(ports, vector):
        for port in ports:
            port.pos = port.pos + Vector(vector)

    @staticmethod
    def rotate_ports(ports, angle):
        if isinstance(angle, list):
            if len(angle) == 2:
                new_angle = np.math.atan2(np.linalg.det([[1, 0], angle]), np.dot([1, 0], angle))
                new_angle = new_angle / np.pi * 180
            else:
                raise Exception("angle should be either a float or a 2-dim array")
        else:
            new_angle = angle
        rad = new_angle / 180 * np.pi
        rotate_matrix = np.array([[np.cos(rad), np.sin(-rad)], [np.sin(rad), np.cos(rad)]])
        for port in ports:
            port.ori = rotate_matrix.dot(port.ori[0:2])
            posx = port.pos[0] * math.cos(rad) + port.pos[1] * math.sin(-rad)
            posy = port.pos[0] * math.sin(rad) + port.pos[1] * math.cos(rad)
            port.pos = Vector([posx, posy])
