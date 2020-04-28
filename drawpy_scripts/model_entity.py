from .parameters import layer_TRACK, \
                        layer_GAP, \
                        layer_MASK, \
                        layer_Default, \
                        layer_RLC

from .utils import Vector, parse_entry, check_name, find_last_list

def add_to_corresponding_list(elt, nested_list, added_elt):
    # return the last list of a set of nested lists
    if isinstance(nested_list, list):
        if elt in nested_list:
            index = nested_list.index(elt)
            nested_list.insert(index+1, added_elt)
            return True
        else:
            for elt_list in nested_list:
                if isinstance(elt_list, list):
                    if add_to_corresponding_list(elt, elt_list, added_elt):
                        break
            else:
                return False
            return True
    else:
        raise TypeError('Argument is not a list')

def general_remove(elt, nested_list):
    # same as list.remove(elt) but for a nested list
    if isinstance(nested_list, list):
        if elt in nested_list:
            nested_list.remove(elt)
            return True
        else:
            for elt_list in nested_list:
                if isinstance(elt_list, list):
                    success = general_remove(elt, elt_list)
                    if success:
                        break
    else:
        raise TypeError('Argument is not a list')

def gen_name(name):
    # routine to mimic the default naming procedure of HFSS when object
    # already exists
    end = ''
    for ii in name[::-1]:
        if ii.isdigit():
            end+=ii
        else:
            break
    if end=='':
        return name+'1'
    number = int(end[::-1])
    if number==0:
        return name+'1'
    else:
        prefix = name[:-len(str(number))]
        suffix = str(number+1)
        return prefix+suffix

class ModelEntity():
    # this should be the objects we are handling on the python interface
    # each method of this class should act in return in HFSS/GDS when possible
    instances_layered = {layer_TRACK:[], layer_GAP:[], layer_MASK:[],
                         layer_Default:[], layer_RLC:[]}
    dict_instances = {}
    instances_to_move = []
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
            find_last_list(ModelEntity.instances_to_move).append(self)
        else:
            # copy is indeed the original object
            # the new object should be put in the same list indent
            add_to_corresponding_list(copy, ModelEntity.instances_to_move,
                                      self)

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
        # not implemented with the HFSS function for handling the copy better
        # copy and translate the copy
        vec = Vector(vec)
        copy = entity.copy()
        copy.translate(vec)
        return copy

    def fillet(self, radius, vertex_indices):
        # fillet a subset of vertices
        # vertex_indices can be an int or a list of int
        if self.mode=='gds':
            raise NotImplementedError()
        else:
            self.body.interface.fillet(self, radius, vertex_indices)

    def fillets(self, radius):
        # fillet all corner of an entity
        self.body.interface.fillets(self, radius)



    # def make_rlc_boundary(self, corner, size, axis, r, l, c, name="LumpRLC"):
    #     raise NotImplementedError()
#        self.interface.make_rlc_boundary(corner, size, axis, r, l, c, name)

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