import numpy as np

from ..parameters import DEFAULT
from ..utils import (
    Vector,
    add_to_corresponding_list,
    check_name,
    find_last_list,
    gen_name,
    general_remove,
    parse_entry,
    val,
)


class Entity:
    # this should be the objects we are handling on the python interface
    # each method of this class should act in return in HFSS/GDS when possible
    dict_instances = {}

    def __init__(
        self, dimension, body, nonmodel=False, layer=DEFAULT, copy=None, name="entity_0", **kwargs
    ):
        name = check_name(self.__class__, name)
        self.name = name
        self.dimension = dimension
        self.body = body
        self.nonmodel = nonmodel
        self.layer = layer

        Entity.dict_instances[name] = self
        if layer in self.body.entities.keys():
            self.body.entities[layer].append(self)
        else:
            self.body.entities[layer] = [self]

        if copy is None:
            if self.body.entities_to_move is not None:
                find_last_list(self.body.entities_to_move).append(self)
            self.is_boolean = False  # did it suffer a bool operation already ?
            self.is_fillet = False  # did it suffer a fillet operation already ?
        else:
            # copy is indeed the original object
            # the new object should be put in the same list indent
            add_to_corresponding_list(copy, self.body.entities_to_move, self)
            self.is_boolean = copy.is_boolean
            self.is_fillet = copy.is_fillet

    def __str__(self):
        return self.name

    def __repr__(self):
        return self.name

    ### General methods

    @staticmethod
    def reset():
        raise NotImplementedError()
        # self.body.entities = {}
        # Entity.dict_instances = {}
        # Entity.instances_to_move = []

    @classmethod
    def print_instances(cls):
        for instance_name in cls.dict_instances:
            print(instance_name)

    ### Modifying methods

    def delete(self):
        # deletes the Entity and its occurences throughout the code
        # it does not delete the entity Python object anymore
        self.body.interface.delete(self)
        self.dict_instances.pop(self.name)
        self.body.entities[self.layer].remove(self)
        if self.body.entities_to_move is not None:
            general_remove(self, self.body.entities_to_move)

    def copy(self, new_name=None, new_layer=None):
        generated_name = gen_name(self.name)
        self.body.interface.copy(self)
        if new_layer is None:
            new_layer = self.layer
        copied = Entity(
            self.dimension,
            self.body,
            nonmodel=self.nonmodel,
            layer=new_layer,
            copy=self,
            name=generated_name,
        )
        if new_name is not None:
            copied.rename(new_name)
        return copied

    def rename(self, new_name):
        self.dict_instances.pop(self.name)
        self.dict_instances[new_name] = self
        self.body.interface.rename(self, new_name)
        self.name = new_name

    def thicken_sheet(self, thickness, bothsides=False):
        self.body.interface.thicken_sheet(self, thickness, bothsides=False)
        self.dimension = 3

    def assign_perfect_E(self, suffix="perfE"):
        self.body.interface.assign_perfect_E(self, self.name + "_" + suffix)

    def assign_waveport(
        self,
        Nmodes=1,
        DoRenorm=False,
        RenormValue="50ohm",
        DoDeembed=False,
        DeembedDist="0mm",
        prefix="port",
    ):
        self.body.interface.assign_waveport(
            self,
            prefix + "_" + self.name,
            Nmodes,
            DoRenorm,
            RenormValue,
            DoDeembed,
            DeembedDist,
        )

    def assign_terminal_auto(self, ground, prefix="port"):
        self.body.interface.assign_terminal_auto(self, prefix + "_" + self.name, ground)

    def connect_faces(self, name, entity1, entity2):
        raise NotImplementedError()

    def duplicate_along_line(self, vec):
        # copy and translate the copy
        vec = Vector(vec)
        copy = self.copy()
        copy.translate(vec)
        return copy

    def find_vertex(self):
        vertices = self.body.interface.get_vertices(self)
        return vertices

    def find_start_vertex(self):
        # finds the lowest vertex in Y in a polygon
        # if there are several, returns the lowest in X
        vertices = self.body.interface.get_vertices(self)
        min_y = vertices[0][1]
        min_x = vertices[0][0]
        indices = [0]
        for ii, vertex in enumerate(vertices[1:]):
            ii += 1
            y = vertex[1]
            if y < min_y:
                indices = [ii]
                min_y = y
                min_x = vertex[0]
            elif y == min_y:
                indices.append(ii)

        result_index = indices[0]
        for index in indices[1:]:
            x = vertices[index][0]
            if x <= min_x:
                min_x = x
                result_index = index
        next_index = (index + 1) % len(vertices)
        prev_index = (index - 1) % len(vertices)
        dx_p = vertices[prev_index][0] - vertices[index][0]
        dx_n = vertices[next_index][0] - vertices[index][0]
        dy_p = vertices[prev_index][1] - vertices[index][1]
        dy_n = vertices[next_index][1] - vertices[index][1]
        angle_p = np.arctan2(dy_p, dx_p)
        angle_n = np.arctan2(dy_n, dx_n)
        is_trigo = angle_p > angle_n
        return result_index, len(vertices), is_trigo

    def fillet(self, radius, vertex_indices=None):

        assert not self.is_fillet, "Cannot fillet an already filleted entity"

        if vertex_indices is None:
            # filleting all vertices
            msg = "Should provide a single radius when filleting all vertices"
            assert not isinstance(radius, list), msg
            if self.body.mode == "gds":
                radius = val(radius)
            self.body.interface.fillet(self, radius)
            self.is_fillet = True
            return None

        if not isinstance(vertex_indices, list):
            vertex_indices = [vertex_indices]
        if isinstance(radius, list):
            # expect vertex_indices also a list of list
            msg = "a vertex_indices list should be given for each radius"
            assert len(radius) == len(vertex_indices), msg
            assert isinstance(vertex_indices[0], list), msg
        else:
            radius = [radius]
            vertex_indices = [vertex_indices]

        if self.is_boolean:
            # should find the lowest/leftest vertex to have consistent behaviour
            index_start, nb_vertices, is_trigo = self.find_start_vertex()
            if not is_trigo:
                vertex_indices = [[-ind for ind in index] for index in vertex_indices]
        else:
            index_start = 0
            if self.body.mode == "gds":
                nb_vertices = len(self.body.interface.get_vertices(self))
            else:
                nb_vertices = len(self.body.interface.get_vertex_ids(self))
        vertex_indices = [
            [(index_start + ind) % nb_vertices for ind in index] for index in vertex_indices
        ]

        radius = parse_entry(radius)
        # check that one index is not present twice
        flat_indices = []
        for indices in vertex_indices:
            flat_indices += indices
        msg = "Vertex index is present more than once in fillet"
        assert len(flat_indices) == len(set(flat_indices)), msg
        if self.body.mode == "gds":
            radius = val(radius)
            self.body.interface.fillet(self, radius, vertex_indices)
        else:
            # manipulate vertex_indices in a good way
            past_indices = []
            for rad, indices in zip(radius, vertex_indices):
                new_indices = []
                append_indices = []
                for index in indices:
                    index += sum(i < index for i in past_indices)
                    new_indices.append(index)
                    if index != 0:
                        append_indices.append(index)
                self.body.interface.fillet(self, rad, new_indices)
                past_indices += append_indices
        self.is_fillet = True
        return None

    #
    #    def fille_edget(self, radius, edge_index):
    #
    #
    #        if vertex_indices is None:
    #            # filleting all vertices
    #            msg = 'Should provide a single radius when filleting all vertices'
    #            assert not isinstance(radius, list), msg
    #            if self.body.mode=='hfss':
    #                self.body.interface._fillet_edges(self, radius, edge_index)
    #                self.is_fillet = True
    #            return None
    #
    #        if not isinstance(edge_index, list):
    #            edge_index = [edge_index]
    #        if isinstance(radius, list):
    #            # expect vertex_indices also a list of list
    #            msg = 'a vertex_indices list should be given for each radius'
    #            assert (len(radius) == len(vertex_indices)), msg
    #            assert (isinstance(vertex_indices[0], list)), msg
    #        else:
    #            radius = [radius]
    #            vertex_indices = [vertex_indices]
    #
    #        if self.is_boolean:
    #            # should find the lowest/leftest vertex to have consistent behaviour
    #            index_start, nb_vertices, is_trigo = self.find_start_vertex()
    #            if not is_trigo:
    #                vertex_indices = [[-ind for ind in index]
    #                                  for index in vertex_indices]
    #        else:
    #            index_start = 0
    #            if self.body.mode == 'gds':
    #                nb_vertices = len(self.body.interface.get_vertices(self))
    #            else:
    #                 nb_vertices = len(self.body.interface.get_vertex_ids(self))
    #        vertex_indices = [[(index_start + ind) % nb_vertices
    #                           for ind in index]
    #                          for index in vertex_indices]
    #
    #        radius = parse_entry(radius)
    #        # check that one index is not present twice
    #        flat_indices = []
    #        for indices in vertex_indices:
    #            flat_indices += indices
    #        msg = 'Vertex index is present more than once in fillet'
    #        assert len(flat_indices)==len(set(flat_indices)), msg
    #        if self.body.mode=='gds':
    #            radius = val(radius)
    #            self.body.interface.fillet(self, radius, vertex_indices)
    #        else:
    #            # manipulate vertex_indices in a good way
    #            past_indices = []
    #            for rad, indices in zip(radius, vertex_indices):
    #                new_indices = []
    #                append_indices = []
    #                for index in indices:
    #                    index += sum(i < index for i in past_indices)
    #                    new_indices.append(index)
    #                    if index != 0:
    #                        append_indices.append(index)
    #                self.body.interface.fillet(self, rad, new_indices)
    #                past_indices += append_indices
    #        self.is_fillet = True
    #        return None

    def assign_material(self, material):
        self.body.interface.assign_material(self, material)

    def assign_impedance(self, ResistanceSq, ReactanceSq, name="impedance"):
        self.body.interface.assign_impedance(self, ResistanceSq, ReactanceSq, name)

    def assign_mesh_length(self, mesh_length):
        mesh_length = parse_entry(mesh_length)
        self.body.interface.assign_mesh_length(self, mesh_length)

    def assign_lumped_RLC(self, points, rlc):

        points = parse_entry(points)
        given_point_0, given_point_1 = Vector(points[0]), Vector(points[1])

        # move the points coordinate in the global coordinate system

        point_0 = given_point_0
        point_1 = given_point_1

        pm = self.body.pm
        current_body = self.body

        while True:

            origin = Vector(current_body.rel_coor[0])
            new_x = Vector(current_body.rel_coor[1])
            new_y = Vector(current_body.rel_coor[2])
            new_z = new_x.cross(new_y)

            change_matrix = np.array([new_x, new_y, new_z])

            point_0 = origin + np.dot(point_0, change_matrix)
            point_1 = origin + np.dot(point_1, change_matrix)

            if current_body.ref_name == "Global":
                break

            for body in pm.bodies:

                if body.name == current_body.ref_name:
                    current_body = body
                    break

        # origin = Vector(current_body.rel_coor[0])
        # new_x = Vector(current_body.rel_coor[1])
        # new_y = Vector(current_body.rel_coor[2])
        # new_z = new_x.cross(new_y)

        # change_matrix = np.array([new_x.as_nda(), new_y.as_nda(), new_z.as_nda()])

        # point_0 = origin + np.dot(point_0.as_nda(), change_matrix)
        # point_1 = origin + np.dot(point_1.as_nda(), change_matrix)

        r, l, c = rlc

        self.body.interface.assign_lumped_rlc(self, r, l, c, point_0, point_1, name="RLC")

    def mirrorZ(self):
        raise NotImplementedError()

    def rotate(self, angle):
        self.body.rotate(self, angle)

    def translate(self, vector):
        self.body.translate(self, vector)

    def subtract(self, tool_entities, keep_originals=False):
        """
        tool_entities: a list of Entity or a Entity
        keep_originals: Boolean, True : the tool entities still exist after
                        boolean operation
        """
        return self.body.subtract([self], tool_entities, keep_originals=keep_originals)

    def unite(self, tool_entities, keep_originals=False, new_name=None):
        """
        tool_entities: a list of Entity or a Entity
        if new_name (str) is provided, the tool_entities + self are kept and
        the union is named new_name
        """
        return self.body.unite(
            tool_entities, main=self, keep_originals=keep_originals, new_name=new_name
        )

    def assign_impedance(self, ResistanceSq, ReactanceSq):
        self.body.interface.assign_impedance(self, ResistanceSq, ReactanceSq)
