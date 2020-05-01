import numpy as np

from .parameters import layer_TRACK, \
                        layer_GAP, \
                        layer_MASK, \
                        layer_Default, \
                        layer_RLC

from .utils import Vector, parse_entry, check_name, find_last_list, \
    add_to_corresponding_list, gen_name, val

class ModelEntity():
    # this should be the objects we are handling on the python interface
    # each method of this class should act in return in HFSS/GDS when possible
    instances_layered = {layer_TRACK:[], layer_GAP:[], layer_MASK:[],
                         layer_Default:[], layer_RLC:[]}
    dict_instances = {}
    instances_to_move = None

    def __init__(self, name, dimension, body, nonmodel=False,
                 layer=layer_Default, coor_sys=None, copy=None):
        name = check_name(self.__class__, name)
        self.name = name
        self.dimension = dimension
        self.body = body
        self.nonmodel = nonmodel
        self.layer = layer
        if coor_sys is None:
            self.coor_sys = body.coor_sys
        else:
            self.coor_sys = coor_sys

        ModelEntity.dict_instances[name] = self
        if layer in self.instances_layered.keys():
            ModelEntity.instances_layered[layer].append(self)
        else:
            ModelEntity.instances_layered[layer]=[self]

        if copy is None:
            if ModelEntity.instances_to_move is not None:
                find_last_list(ModelEntity.instances_to_move).append(self)
            self.is_boolean = False  # did it suffer a bool operation already ?
            self.is_fillet = False  # did it suffer a fillet operation already ?
        else:
            # copy is indeed the original object
            # the new object should be put in the same list indent
            add_to_corresponding_list(copy, ModelEntity.instances_to_move,
                                      self)
            self.is_boolean = copy.is_boolean
            self.is_fillet = copy.is_fillet
    def __str__(self):
        return self.name

    def __repr__(self):
        return self.name

    @staticmethod
    def reset():
        ModelEntity.instances_layered = {}
        ModelEntity.dict_instances = {}
        ModelEntity.instances_to_move = []

    @classmethod
    def print_instances(cls):
        for instance_name in cls.dict_instances:
            print(instance_name)

    def delete(self):
        # deletes the modelentity and its occurences throughout the code
        self.body.interface.delete(self)
        self.dict_instances.pop(self.name)
        self.instances_layered[self.layer].remove(self)
        if self.instances_to_move is not None:
            general_remove(self, self.instances_to_move)
        del self

    def copy(self, new_name=None):
        generated_name = gen_name(self.name)
        self.body.interface.copy(self)
        copied = ModelEntity(generated_name, self.dimension, self.body,
                             nonmodel=self.nonmodel, layer=self.layer,
                             coor_sys=self.coor_sys, copy=self)
        if new_name is not None:
            copied.rename(new_name)
        return copied

    def rename(self, new_name):
        self.dict_instances.pop(self.name)
        self.dict_instances[new_name] = self
        self.body.interface.rename(self, new_name)
        self.name = new_name

    def thicken_sheet(self, thickness, bothsides=False):
        raise NotImplementedError()

    def assign_perfect_E(self, suffix='perfE'):
        self.body.interface.assign_perfect_E(self, self.name+'_'+suffix)

    def connect_faces(self, name, entity1, entity2):
        raise NotImplementedError()

    def duplicate_along_line(self, vec):
        # copy and translate the copy
        vec = Vector(vec)
        copy = self.copy()
        copy.translate(vec)
        return copy

    def find_start_vertex(self):
        # finds the lowest vertex in Y in a polygon
        # if there are several, returns the lowest in X
        vertices = self.body.interface.get_vertices(self)
        min_y = vertices[0][1]
        min_x = vertices[0][0]
        indices = [0]
        for ii, vertex in enumerate(vertices[1:]):
            ii+=1
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
        next_index = (index+1)%len(vertices)
        prev_index = (index-1)%len(vertices)
        dx_p = vertices[prev_index][0] - vertices[index][0]
        dx_n = vertices[next_index][0] - vertices[index][0]
        dy_p = vertices[prev_index][1] - vertices[index][1]
        dy_n = vertices[next_index][1] - vertices[index][1]
        angle_p = np.arctan2(dy_p, dx_p)
        angle_n = np.arctan2(dy_n, dx_n)
        is_trigo = angle_p > angle_n
        return result_index, len(vertices), is_trigo

    def fillet(self, radius, vertex_indices=None):

        assert (not self.is_fillet), 'Cannot fillet an already filleted entity'

        if vertex_indices is None:
            # filleting all vertices
            msg = 'Should provide a single radius when filleting all vertices'
            assert not isinstance(radius, list), msg
            if self.body.mode=='gds':
                radius = val(radius)
            self.body.interface.fillet(self, radius)
            return None

        if not isinstance(vertex_indices, list):
            vertex_indices = [vertex_indices]
        if isinstance(radius, list):
            # expect vertex_indices also a list of list
            msg = 'a vertex_indices list should be given for each radius'
            assert (len(radius) == len(vertex_indices)), msg
            assert (isinstance(vertex_indices[0], list)), msg
        else:
            radius = [radius]
            vertex_indices = [vertex_indices]

        if self.is_boolean:
            # should find the lowest/leftest vertex to have consistent behaviour
            index_start, nb_vertices, is_trigo = self.find_start_vertex()
            if not is_trigo:
                vertex_indices = [[-ind for ind in index]
                                  for index in vertex_indices]
            vertex_indices = [[(index_start + ind) % nb_vertices
                               for ind in index]
                              for index in vertex_indices]

        self.is_fillet = True
        radius = parse_entry(radius)
        # check that one index is not present twice
        flat_indices = []
        for indices in vertex_indices:
            flat_indices += indices
        msg = 'Vertex index is present more than once in fillet'
        assert len(flat_indices)==len(set(flat_indices)), msg
        if self.body.mode=='gds':
            radius = val(radius)
            self.body.interface.fillet(self, radius, vertex_indices)
        else:
            # manipulate vertex_indices in a good way
            past_indices = []
            for rad, indices in zip(radius, vertex_indices):
                new_indices = []
                for index in indices:
                    index += sum(i < index for i in past_indices)
                    new_indices.append(index)
                self.body.interface.fillet(self, rad, new_indices)
                past_indices += indices
        return None

    def assign_material(self, material):
        self.body.interface.assign_material(self, material)

    def assign_mesh_length(self, mesh_length):
        mesh_length = parse_entry(mesh_length)
        self.body.interface.assign_mesh_length(self, mesh_length)

    def assign_lumped_RLC(self, points, rlc):
        points = parse_entry(points)
        # move the points coordinate in the global coordinate system
        if self.body.ref_name != 'Global':
            # TODO do recursive to handle this
            raise NotImplementedError('Do not handle 2nd order relative \
                                      coordinate system yet.')
        origin = self.body.rel_coor[0]
        new_x = self.body.rel_coor[1]
        new_y =self.body.rel_coor[2]
        point_0 = []
        point_1 = []
        for ii in range(3):
            point_0.append(origin[ii] + new_x[ii] * points[0][0] + new_y[ii] * points[0][1])
            point_1.append(origin[ii] + new_x[ii] * points[1][0] + new_y[ii] * points[1][1])

        r, l, c = rlc
        self.body.interface.assign_lumped_rlc(self, r, l, c, point_0,
                                              point_1, name="RLC")

    def mirrorZ(self):
        raise NotImplementedError()

    def rotate(self, angle):
        self.body.rotate(self, angle)

    def translate(self, vector):
        self.body.translate(self, vector)

    def subtract(self, tool_entities, keep_originals=False):
        """
        tool_entities: a list of ModelEntity or a ModelEntity
        keep_originals: Boolean, True : the tool entities still exist after
                        boolean operation
        """
        if not isinstance(tool_entities, list):
            tool_entities = [tool_entities]
        if len(tool_entities)==0:
            pass
        else:
            if not all([entity.dimension==self.dimension
                                                for entity in tool_entities]):
                raise TypeError('All subtracted elements should have the \
                                same dimension')
            else:
                self.body.interface.subtract(self, tool_entities,
                                           keep_originals=True)
                self.is_boolean = True
            if not keep_originals:
                for ii in range(len(tool_entities)):
                    tool_entities[0].delete()

    def unite(self, tool_entities, new_name=None):
        """
        tool_entities: a list of ModelEntity or a ModelEntity
        if new_name (str) is provided, the tool_entities + self are kept and
        the union is named new_name
        """
        return self.body.unite(tool_entities, main=self, new_name=new_name)
