import numpy as np

from HFSSdrawpy import Port
from HFSSdrawpy.utils import find_last_list, find_penultimate_list


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
        pos, angle = self.body.cursors.pop(-1)

        # 5 We move the entities_to_move with the right operation
        if len(list_entities_new) > 0:
            self.body.rotate(list_entities_new, angle=angle)
            self.body.translate(list_entities_new, vector=[pos[0], pos[1], pos[2]])

        if len(list_ports_new) > 0:
            Port.rotate_ports(list_ports_new, angle)
            Port.translate_ports(list_ports_new, vector=[pos[0], pos[1], pos[2]])

        clean_body_to_move(self.body)
        return False


def clean_body_to_move(body):
    # 6 We empty a part of the 'to_move' lists
    penultimate_entity_list = find_penultimate_list(body.entities_to_move)
    if penultimate_entity_list:
        a = penultimate_entity_list.pop(-1)
        for entity in a:
            penultimate_entity_list.append(entity)
    else:
        body.entities_to_move = None

    penultimate_port_list = find_penultimate_list(body.ports_to_move)
    if penultimate_port_list:
        a = penultimate_port_list.pop(-1)
        for entity in a:
            penultimate_port_list.append(entity)
    else:
        body.ports_to_move = None
