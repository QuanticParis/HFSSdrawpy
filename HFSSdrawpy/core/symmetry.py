import copy

from HFSSdrawpy.utils import find_last_list


class BodyMirror:
    def __init__(self, body, angle, magnitude):
        self.normal_vector_polar = (magnitude, angle)
        self.body = body

    def __enter__(self):
        # 1 We need to keep track of the entities created in the symmetry context
        if self.body.entities_to_move is None:
            self.body.entities_to_move = []
        else:
            find_last_list(self.body.entities_to_move).append([])

        if self.body.ports_to_move is None:
            self.body.ports_to_move = []
        else:
            find_last_list(self.body.ports_to_move).append([])

    def __exit__(self, *exception):
        # we have to copy the arrays because when we mirror an element, the symmetric will be appended
        # to the list, so it will be mirrored, and so on. Copying the arrays ensures we only mirror
        # elements that originally were added by the user
        list_entities_new = copy.copy(find_last_list(self.body.entities_to_move))
        list_ports_new = find_last_list(self.body.ports_to_move)

        if len(list_entities_new) > 0:
            symmetric_entities = []
            for entity in list_entities_new:
                symmetric_entity = entity.copy(entity.name + "_symmetric")
                symmetric_entities.append(symmetric_entity)

            self.body.apply_mirror(symmetric_entities, self.normal_vector_polar)
