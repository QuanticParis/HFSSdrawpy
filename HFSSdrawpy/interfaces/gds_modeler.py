import gdspy
import numpy as np
from gdspy import FlexPath

from ..core.symmetry import compute_translation_rotation
from ..utils import Vector, parse_entry, val, points_on_line_tangent_to, check_name

TOLERANCE = 1e-9  # for arcs
print("gdspy_version : ", gdspy.__version__)


class GdsModeler:
    gds_object_instances = {}
    gds_cells = {}
    dict_units = {"km": 1.0e3, "m": 1.0, "cm": 1.0e-2, "mm": 1.0e-3}
    # coor_systems = {'Global':[[0,0,0],[1,0]]}
    # coor_system = coor_systems['Global']

    def __init__(self, unit=1.0e-6, precision=1.0e-9):
        self.unit = unit
        self.precision = precision
        gdspy.current_library = gdspy.GdsLibrary()

    @classmethod
    def print_instances(cls):
        for instance_name in cls.gds_object_instances:
            print(instance_name)

    def reset_cell(self):
        del self.cell

    def create_coor_sys(self, coor_sys="chip", rel_coor=None, ref_name="Global"):
        # this creates a cell, should not care about the rel_coor
        if not (coor_sys in gdspy.current_library.cells.keys()):
            cell = gdspy.Cell(coor_sys)
            self.gds_cells[coor_sys] = cell
        else:
            cell = self.gds_cells[coor_sys]
        # active cell should be the new cell
        self.cell = cell

    def set_coor_sys(self, coor_sys):
        if coor_sys in self.gds_cells.keys():
            self.cell = self.gds_cells[coor_sys]
        else:
            raise ValueError("%s cell do not exist" % coor_sys)

    def copy(self, entity):
        new_polygon = gdspy.copy(self.gds_object_instances[entity.name], 0, 0)
        new_name = check_name(entity.__class__, entity.name)
        self.gds_object_instances[new_name] = new_polygon
        self.cell.add(new_polygon)
        return new_name

    def relayer(self, entity, layer):
        obj = self.gds_object_instances[entity.name]
        obj.layers = [layer] * len(obj.layers)

    def generate_gds(self, file, max_points):
        for instance in self.gds_object_instances.keys():
            obj = self.gds_object_instances[instance]
            if isinstance(obj, gdspy.Polygon) or isinstance(obj, gdspy.PolygonSet):
                self.gds_object_instances[instance] = obj.fracture(
                    max_points=max_points, precision=1e-9
                )
        for cell_name in self.gds_cells.keys():
            filename = file + "_%s.gds" % cell_name
            gdspy.write_gds(filename, cells=[cell_name], unit=1.0, precision=1e-9)

    def import_gds(self, file):
        lib.read_gds(file)
        for instance in self.gds_object_instances.keys():
            obj = self.gds_object_instances[instance]
            if isinstance(obj, gdspy.Polygon) or isinstance(obj, gdspy.PolygonSet):
                self.gds_object_instances[instance] = obj.fracture(
                    max_points=max_points, precision=1e-9
                )
        for cell_name in self.gds_cells.keys():
            filename = file + "_%s.gds" % cell_name
            gdspy.write_gds(filename, cells=[cell_name], unit=1.0, precision=1e-9)

    def get_vertices(self, entity):
        polygon = self.gds_object_instances[entity.name]
        return polygon.polygons[0]

    def set_units(self, units="m"):
        self.unit = self.dict_units[units]

    def box(self, pos, size, **kwargs):
        pass

    def box_center(self, pos, size, **kwargs):
        pass

    def polyline(self, points, closed, **kwargs):
        # TODO sace of open path
        # size is the thickness of the polyline for gds, must be a 2D-list with idential elements
        name = kwargs["name"]
        layer = kwargs["layer"]
        points = parse_entry(points)

        # TODO, this is a dirty fixe cause of Vector3D

        points_2D = []
        for point in points:
            points_2D.append([point[0], point[1]])

        if closed:
            poly1 = gdspy.Polygon(points_2D, layer=layer)
        else:
            poly1 = gdspy.FlexPath(points_2D, 1e-9, layer=layer)

        self.gds_object_instances[name] = poly1
        self.cell.add(poly1)

    def rect(self, pos, size, **kwargs):
        pos, size = parse_entry(pos, size)
        name = kwargs["name"]
        layer = kwargs["layer"]
        # This function neglects the z coordinate
        points = [
            (pos[0], pos[1]),
            (pos[0] + size[0], pos[1] + 0),
            (pos[0] + size[0], pos[1] + size[1]),
            (pos[0], pos[1] + size[1]),
        ]
        poly1 = gdspy.Polygon(points, layer)

        self.gds_object_instances[name] = poly1
        self.cell.add(poly1)
        return poly1

    def text(self, pos, size, text, angle, horizontal, **kwargs):
        pos, size = parse_entry(pos, size)
        name = kwargs["name"]
        layer = kwargs["layer"]

        poly1 = gdspy.Text(
            text, size, pos, horizontal=horizontal, angle=angle, layer=layer
        )

        self.gds_object_instances[name] = poly1
        self.cell.add(poly1)

    def rect_center(self, pos, size, **kwargs):
        pos, size = parse_entry(pos, size)
        corner_pos = [val(p) - val(s) / 2 for p, s in zip(pos, size)]
        self.rect(corner_pos, size, **kwargs)

    def cylinder(self, pos, radius, height, axis, **kwargs):
        pass

    def disk(self, pos, radius, axis, number_of_points=None, **kwargs):
        pos, radius = parse_entry(pos, radius)
        name = kwargs["name"]
        layer = kwargs["layer"]
        assert axis == "Z", "axis must be 'Z' for the gdsModeler"
        round1 = gdspy.Round(
            (pos[0], pos[1]),
            radius,
            layer=layer,
            tolerance=TOLERANCE,
            number_of_points=number_of_points,
        )
        self.gds_object_instances[name] = round1
        self.cell.add(round1)

    def wirebond(
        self, pos, ori, ymax, ymin, height="0.1mm", **kwargs
    ):  # ori should be normed
        bond_diam = "20um"
        bond_pad= "150um"
        pos, ori, ymax, ymin, heigth, bond_diam, bond_pad = parse_entry(
            (pos, ori, ymax, ymin, height, bond_diam, bond_pad)
        )
        bond1 = pos + ori.orth() * (ymax)
        bond2 = pos + ori.orth() * (ymin)

        cable = gdspy.FlexPath(
            [bond1[0:2], bond2[0:2]],
            [bond_diam],
            gdsii_path=False,
            tolerance=TOLERANCE,
            layer=[kwargs["layer_bond"]],
            max_points=0
        )
        polygons = cable.get_polygons()
        names = []
        layers = []
        for ii in range(len(polygons)):
            poly = gdspy.Polygon(polygons[ii])
            poly.layers = [kwargs["layer_bond"]]
            current_name = kwargs["name"] + "connect"
            names.append(current_name)
            layers.append([kwargs["layer_bond"]])
            self.gds_object_instances[current_name] = poly
            self.cell.add(poly)

        self.disk(
            bond1,
            bond_pad / 2,
            "Z",
            layer=kwargs["layer"],
            name=kwargs["name"] + "a",
            number_of_points=None,
        )
        self.disk(
            bond2,
            bond_pad / 2,
            "Z",
            layer=kwargs["layer"],
            name=kwargs["name"] + "b",
            number_of_points=None,
        )

    def path(self, points, port, fillet, name="", corner="circular bend"):

        # TODO, this is a dirty fixe cause of Vector3D

        points_2D = []
        for point in points:
            points_2D.append([point[0], point[1]])

        # use dummy layers to recover the right elements
        layers = [ii for ii in range(len(port.widths))]
        cable = gdspy.FlexPath(
            points_2D,
            port.widths,
            offset=port.offsets,
            corners=corner,
            bend_radius=fillet,
            gdsii_path=False,
            tolerance=TOLERANCE,
            layer=layers,
            max_points=0,
        )  # tolerance (meter) is highly important here should be smaller than the smallest dim typ. 100nm

        polygons = cable.get_polygons()
        names = []
        layers = []
        for ii in range(len(polygons)):
            poly = gdspy.Polygon(polygons[ii])
            poly.layers = [port.layers[ii]]
            current_name = name + "_" + port.subnames[ii]
            names.append(current_name)
            layers.append(port.layers[ii])
            self.gds_object_instances[current_name] = poly
            self.cell.add(poly)
        return names, layers

    def connect_faces(self, entity1, entity2):
        pass

    def delete(self, entity):
        self.cell.polygons.remove(self.gds_object_instances[entity.name])
        self.gds_object_instances.pop(entity.name)

    def rename_entity(self, old_name, new_name):
        polygon = self.gds_object_instances.pop(old_name)
        self.gds_object_instances[new_name] = polygon

    def unite(self, entities, keep_originals=True):

        blank_entity = entities.pop(0)
        blank_polygon = self.gds_object_instances.pop(blank_entity.name)
        self.cell = self.gds_cells[blank_entity.body.name]
        self.cell.polygons.remove(blank_polygon)

        tool_polygons = []
        for tool_entity in entities:
            tool_polygon = self.gds_object_instances[tool_entity.name]
            if isinstance(tool_polygon, gdspy.PolygonSet):
                for polygon in tool_polygon.polygons:
                    tool_polygons.append(polygon)
            else:
                tool_polygons.append(tool_polygon)

        # 2 unite operation
        tool_polygon_set = gdspy.PolygonSet(tool_polygons, layer=blank_entity.layer)
        united = gdspy.boolean(
            blank_polygon,
            tool_polygon_set,
            "or",
            precision=TOLERANCE,
            max_points=0,
            layer=blank_entity.layer,
        )

        self.gds_object_instances[blank_entity.name] = united
        self.cell.add(united)

        return blank_entity

    def intersect(self, entities):
        raise NotImplementedError()

    def scale(self, entities, factor, center):
        if isinstance(entities, list):
            for entity in entities:
                self.scale(
                    self,
                    entity,
                    factor,
                    center
                    )
        else:
            polygon = self.gds_object_instances[entities.name]
            polygon.scale(factor, factor, center)

    def subtract(self, blank_entities, tool_entities, keep_originals=True):
        if isinstance(blank_entities, list):
            for blank_entity in blank_entities:
                self.subtract(
                    blank_entity, tool_entities, keep_originals=keep_originals
                )
        else:
            blank_entity = blank_entities
            # 1 We clear the cell of all elements and create lists to store the polygons
            blank_polygon = self.gds_object_instances.pop(blank_entity.name)
            self.cell = self.gds_cells[
                blank_entity.body.name
            ]  # assumes blank and tool are in same body
            self.cell.polygons.remove(blank_polygon)

            tool_polygons = []
            for tool_entity in tool_entities:
                tool_polygon = self.gds_object_instances[tool_entity.name]
                if isinstance(tool_polygon, gdspy.PolygonSet):
                    for polygon in tool_polygon.polygons:
                        tool_polygons.append(polygon)
                else:
                    tool_polygons.append(tool_polygon)

            # 2 subtract operation
            tool_polygon_set = gdspy.PolygonSet(tool_polygons, layer=blank_entity.layer)
            subtracted = gdspy.boolean(
                blank_polygon,
                tool_polygon_set,
                "not",
                precision=TOLERANCE,
                max_points=0,
                layer=blank_entity.layer,
            )
            if subtracted is not None:
                # 3 At last we update the cell and the gds_object_instance
                self.gds_object_instances[blank_entity.name] = subtracted
                self.cell.add(subtracted)
            else:
                print(
                    "Warning: the entity %s was fully \
                      subtracted"
                    % blank_entity.name
                )
                dummy = gdspy.Polygon([[0, 0]])
                self.gds_object_instances[blank_entity.name] = dummy
                self.cell.add(dummy)
                blank_entity.delete()

    def delete_inside(self, poly_set, mask, keep_originals=False):
        """
        Test if the polygons within the poly_set are in the mask object.
        If so, delete them.
        Parameters
        ----------
        poly_set : Entity
            Typically a hole array.
        mask : Entity
        keep_originals : bool, optional
            Shall we keep the mask element or not. The default is False.
        """
        if isinstance(mask, list):
            if len(mask) > 1:
                mask = mask[0].unite(mask[1:])
            else:
                mask = mask[0]
        poly_set_gds = self.gds_object_instances[poly_set.name]
        mask_gds = self.gds_object_instances[mask.name]
        N = len(poly_set_gds.polygons)
        for ii in range(len(poly_set_gds.polygons)):
            points = poly_set_gds.polygons[N - 1 - ii]
            condition = bool(sum(gdspy.inside(points, mask_gds, precision=TOLERANCE)))
            if condition:
                poly_set_gds.polygons.pop(N - 1 - ii)
        if not keep_originals:
            mask.delete()

    def assign_material(self, *args, **kwargs):
        pass

    def assign_perfect_E(self, entity, name=None):
        pass

    def assign_impedance(self, entities, ResistanceSq, ReactanceSq, name="impedance"):
        pass

    def assign_perfect_E_faces(self, entity):
        pass

    def assign_mesh_length(self, entity, length):  # , suff = '_mesh'):
        pass

    def assign_lumped_rlc(self, entity, r, l, c, start, end, name="RLC"):
        pass

    def assign_lumped_port(self, entity, start, end, name="RLC"):
        pass

    def assign_waveport(self, *args, **kwargs):
        pass

    def assign_terminal_auto(self, *args, **kwargs):
        pass

    def create_object_from_face(self, name):
        pass

    def fillet(self, entity, radius, vertex_indices=None):
        polygon = self.gds_object_instances[entity.name]
        if vertex_indices is None:
            polygon.fillet(radius, max_points=0)
        else:
            vertices_number = len(polygon.polygons[0])
            radii = [0] * vertices_number
            for rad, indices in zip(radius, vertex_indices):
                for index in indices:
                    radii[index] = rad
            polygon.fillet([radii], max_points=0, precision=TOLERANCE)

    def get_vertex_ids(self, entity):
        return None

    def sweep_along_vector(self, names, vector):
        self._modeler.SweepAlongVector(
            self._selections_array(*names),
            [
                "NAME:VectorSweepParameters",
                "DraftAngle:=",
                "0deg",
                "DraftType:=",
                "Round",
                "CheckFaceFaceIntersection:=",
                False,
                "SweepVectorX:=",
                vector[0],
                "SweepVectorY:=",
                vector[1],
                "SweepVectorZ:=",
                vector[2],
            ],
        )

    def thicken_sheet(self, sheet, thickness, bothsides=False):
        self._modeler.ThickenSheet(
            ["NAME:Selections", "Selections:=", sheet, "NewPartsModelFlag:=", "Model"],
            [
                "NAME:SheetThickenParameters",
                "Thickness:=",
                thickness,
                "BothSides:=",
                bothsides,
            ],
        )

    def mirrorZ(self, entity):
        pass

    def translate(self, entities, vector):
        """vector is 3-dimentional but with a z=0 component"""
        if not isinstance(entities, list):
            entities = [entities]
        translation_vector = [vector[0], vector[1]]
        for entity in entities:
            # if entity!=None:
            gds_entity = self.gds_object_instances[entity.name]
            gds_entity.translate(*translation_vector)

    def rotate(self, entities, angle, center=None):
        if center is None:
            center = (0, 0)

        if not isinstance(entities, list):
            entities = [entities]
        for entity in entities:
            # if entity!=None:
            if not entity.esc:
                gds_entity = self.gds_object_instances[entity.name]
                gds_entity.rotate(
                    angle / 360 * 2 * np.pi, center=(val(center[0]), val(center[1]))
                )

    def mirror(self, entities, normal_vector_polar):
        if not isinstance(entities, list):
            entities = [entities]

        p1, p2 = points_on_line_tangent_to(normal_vector_polar)

        mirror_entities = dict()
        for entity in entities:
            mirror_name = entity.name
            print("Register mirror object", mirror_name)
            gds_entity = self.gds_object_instances[mirror_name]
            if isinstance(gds_entity, gdspy.FlexPath):
                points = gds_entity.points
                mirrored_points = []
                for point in points:
                    dp, _ = compute_translation_rotation(
                        point, np.zeros(2), normal_vector_polar
                    )
                    mirrored_point = point + dp
                    mirrored_points.append(mirrored_point)
                mirrored_points = np.array(mirrored_points)
                gds_entity.points = mirrored_points
                mirror = gds_entity
            else:
                mirror = gds_entity.mirror(p1, p2)

            mirror_entities[mirror_name] = mirror
            self.gds_object_instances[mirror_name] = mirror

        return mirror_entities

    def rect_array(self, pos, size, columns, rows, spacing, origin=(0, 0), **kwargs):
        """
        This is a raw function. Better use ab_elt.rect_array
        """
        pos, size, spacing = parse_entry(Vector(pos), Vector(size), spacing)
        name = kwargs["name"]
        layer = kwargs["layer"]

        rect1 = self.rect(pos - size / 2, size, **kwargs)

        cell_to_copy = gdspy.Cell(f"{name}_cell_to_copy")
        # self.gds_cells["cell_to_copy"] = cell_to_copy
        # commenting the previous line allow to ignore the cell at gds generation
        cell_to_copy.add(rect1)

        cell_array = gdspy.CellArray(cell_to_copy, columns, rows, spacing, origin)
        polygon_list = cell_array.get_polygons()
        _rect_array = gdspy.PolygonSet(polygon_list, layer)
        self.cell.add(_rect_array)
        self.gds_object_instances[name] = _rect_array
