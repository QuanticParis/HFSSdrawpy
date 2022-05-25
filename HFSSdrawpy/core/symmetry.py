import copy

import numpy as np

from HFSSdrawpy import Port
from HFSSdrawpy.core.move import clean_body_to_move
from HFSSdrawpy.parameters import PORT
from HFSSdrawpy.utils import (
    find_last_list,
    find_penultimate_list,
    points_on_line_tangent_to,
    check_name,
)


class BodyMirror:
    def __init__(self, body, angle, magnitude):
        self.normal_vector_polar = (magnitude, angle)
        self.body = body
        # an array that contains in key a port and in value its symmetric
        self.port_symmetry_correspondence = dict()

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
        list_ports_new = copy.copy(find_last_list(self.body.ports_to_move))

        if len(list_entities_new) > 0:
            symmetric_entities = []
            for entity in list_entities_new:
                symmetric_entity = entity.copy(entity.name + "_symmetric")
                symmetric_entities.append(symmetric_entity)

            self.body.apply_mirror(symmetric_entities, self.normal_vector_polar)

        if len(list_ports_new) > 0:
            for port in list_ports_new:
                if port.key != "name":  # if the port is an inverse, we don't copy it
                    continue
                name = check_name(Port, f"{port.name}_symmetric")

                # Compute symmetric position
                position = port.pos[:2]
                orientation = port.ori[:2]
                d_position, d_orientation = compute_translation_rotation(
                    position, orientation, self.normal_vector_polar
                )
                symmetric_position = position + d_position
                symmetric_orientation = orientation + d_orientation

                # draw it
                with self.body(
                    symmetric_position.tolist(), symmetric_orientation.tolist()
                ):
                    port_args = port.val(
                        drop_mask=True
                    )  # the mask will be re-added in `body.port`
                    (new_port,) = self.body.port(
                        port_args.widths or None,
                        port_args.subnames,
                        port_args.layers,
                        port_args.offsets or None,
                        name,
                    )

                self.port_symmetry_correspondence[port] = new_port

        clean_body_to_move(self.body)

    def get_symmetric_of(self, port):
        return self.port_symmetry_correspondence[port]


def compute_translation_rotation(position, orientation, normal_vector_polar):
    position = np.array(position)
    p1, p2 = points_on_line_tangent_to(normal_vector_polar)
    p1, p2 = np.array(p1), np.array(p2)
    n = p1 - p2
    n = n.astype(float)
    n /= np.linalg.norm(n)
    d_position = 2 * (np.eye(2) - np.outer(n, n)) @ (p1 - position)

    # computed symmetric orientation
    orientation = np.array(orientation)
    tangent_vector = np.array(
        [
            np.cos(normal_vector_polar[1]),
            np.sin(normal_vector_polar[1]),
        ]
    )
    d_orientation = -2 * (orientation @ tangent_vector) * tangent_vector

    return d_position, d_orientation
