import numpy
import sympy
from pint import UnitRegistry
from sympy.parsing import sympy_parser

ureg = UnitRegistry()
Q = ureg.Quantity

LENGTH = "[length]"
INDUCTANCE = "[length] ** 2 * [mass] / [current] ** 2 / [time] ** 2"
CAPACITANCE = "[current] ** 2 * [time] ** 4 / [length] ** 2 / [mass]"
RESISTANCE = "[length] ** 2 * [mass] / [current] ** 2 / [time] ** 3"
FREQUENCY = "1 / [time]"
DIMENSIONLESS = "dimensionless"

LENGTH_UNIT = "meter"
INDUCTANCE_UNIT = "nH"
CAPACITANCE_UNIT = "nF"
RESISTANCE_UNIT = "ohm"
FREQUENCY_UNIT = "GHz"
DIMENSIONLESS_UNIT = ""

### List handling
# Useful function to manipulate to_move entities and ports
def find_last_list(list_entities):
    # return the last list of a set of nested lists

    if isinstance(list_entities, list):
        if len(list_entities) == 0:
            return list_entities
        else:
            if isinstance(list_entities[-1], list):
                return find_last_list(list_entities[-1])
            else:
                return list_entities
    else:
        raise TypeError("There are no list")


def find_penultimate_list(list_entities):
    # return the last list of a set of nested lists
    if isinstance(list_entities, list):
        if len(list_entities) == 0:
            return False
        else:
            if isinstance(list_entities[-1], list):
                if len(list_entities[-1]) == 0:
                    return list_entities
                else:
                    if isinstance(list_entities[-1][-1], list):
                        return find_penultimate_list(list_entities[-1])
                    else:
                        return list_entities
            else:
                return False
    else:
        raise TypeError("There are no list")


def add_to_corresponding_list(elt, nested_list, added_elt):
    # return the last list of a set of nested lists
    if isinstance(nested_list, list):
        if elt in nested_list:
            index = nested_list.index(elt)
            nested_list.insert(index + 1, added_elt)
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
        pass  # raise TypeError('Argument is not a list')


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
        raise TypeError("Argument is not a list")


def find_corresponding_list(elt, nested_list):
    # return the last list of a set of nested lists
    if isinstance(nested_list, list):
        if elt in nested_list:
            return nested_list
        else:
            for elt_list in nested_list:
                if isinstance(elt_list, list):
                    found_list = find_corresponding_list(elt, elt_list)
                    if found_list:
                        break
            else:
                return False
            return found_list
    else:
        return None


### Naming


def gen_name(name):
    # routine to mimic the default naming procedure of HFSS when object
    # already exists
    end = ""
    for ii in name[::-1]:
        if ii.isdigit():
            end += ii
        else:
            break
    if end == "":
        return name + "1"
    number = int(end[::-1])
    if number == 0:
        return name + "1"
    else:
        prefix = name[: -len(str(number))]
        suffix = str(number + 1)
        return prefix + suffix


def check_name(_class, name):
    end = ""
    for ii, char in enumerate(name[::-1]):
        if char.isdigit():
            end += char
        else:
            break
    else:
        ii += 1
    if end == "":
        radical = name
        number = 0
    else:
        radical = name[:-ii]
        number = int(end[::-1])
    new_name = name
    while new_name in _class.dict_instances.keys():
        number += 1
        new_name = radical + str(number)
    if new_name != name:
        print("%s: changed '%s' name into '%s'" % (_class.__name__, name, new_name))
    return new_name


### Litteral Expressions


def equal_float(float1, float2):
    if abs(float1) > 1e-10:
        rel_diff = abs((float1 - float2) / float1)
        if rel_diff < 1e-5:
            return True
        else:
            return False

    elif abs(float2) > 1e-10:
        rel_diff = abs((float1 - float2) / float2)
        if rel_diff < 1e-5:
            return True
        else:
            return False
    else:
        return True


def simplify_arith_expr(expr):
    try:
        out = repr(sympy_parser.parse_expr(str(expr)))
        return out
    except Exception:
        print("Couldn't parse", expr)
        raise


def extract_value_unit(expr, units):
    """
    :type expr: str
    :type units: str
    :return: float
    """
    try:
        return Q(expr).to(units).magnitude
    except Exception:
        try:
            return float(expr)
        except Exception:
            return expr


def extract_value_dim(expr):
    """
    type expr: str
    """
    try:
        return str(Q(expr).dimensionality)
    except Exception:
        return DIMENSIONLESS

def parse_entry(*entries, marker=True):
    # should take a list of tuple of list... of int, float or str...
    parsed = []
    for entry in entries:
        if not isinstance(entry, list) and not isinstance(entry, tuple):
            parsed.append(extract_value_unit(entry, LENGTH_UNIT))
        else:
            if isinstance(entry, list):
                if isinstance(entry, Vector):
                    parsed.append(Vector(parse_entry(*entry, marker=False)))
                else:
                    parsed.append(parse_entry(*entry, marker=False))
            elif isinstance(entry, tuple):
                parsed.append(tuple(parse_entry(*entry, marker=False)))
            else:
                raise TypeError("Not foreseen type: %s" % (type(entry)))
    if len(parsed) == 1 and marker:
        return parsed[0]
    else:
        return parsed


def rem_unit(other):
    try:
        value = extract_value_unit(other, LENGTH_UNIT)
        return value
    except Exception:
        return other


def _val(elt):
    if isinstance(elt, (int, float, numpy.int64, numpy.float64, numpy.int32, numpy.float32)):
        return elt
    elif isinstance(elt, str):
        if LENGTH == extract_value_dim(elt):
            unit = LENGTH_UNIT
        if INDUCTANCE == extract_value_dim(elt):
            unit = INDUCTANCE_UNIT
        if CAPACITANCE == extract_value_dim(elt):
            unit = CAPACITANCE_UNIT
        if RESISTANCE == extract_value_dim(elt):
            unit = RESISTANCE_UNIT
        if FREQUENCY == extract_value_dim(elt):
            unit = FREQUENCY_UNIT
        if DIMENSIONLESS == extract_value_dim(elt):
            unit = DIMENSIONLESS_UNIT
        return extract_value_unit(elt, unit)
    else:
        return float(elt.evalf(subs=variables))


def val(*entries, marker=True):
    # should take a list of tuple of list... of int, float or str...
    parsed = []
    for entry in entries:
        if not isinstance(entry, (list, tuple, Vector)):
            parsed.append(_val(entry))
        else:
            if isinstance(entry, Vector):
                parsed.append(Vector(val(*entry, marker=False)))
            elif isinstance(entry, list):
                parsed.append(val(*entry, marker=False))
            elif isinstance(entry, tuple):
                parsed.append(tuple(val(*entry, marker=False)))
            else:
                raise TypeError("Not foreseen type: %s" % (type(entry)))
    if len(parsed) == 1 and marker:
        return parsed[0]
    else:
        return parsed


def way(vec):
    if vec[1] != 0:
        if abs(vec[0] / vec[1]) < 1e-3:
            if vec[1] > 0:
                return Vector(0, 1)
            elif vec[1] < 0:
                return Vector(0, -1)
    if vec[0] != 0:
        if abs(vec[1] / vec[0]) < 1e-3:
            if vec[0] > 0:
                return Vector(1, 0)
            elif vec[0] < 0:
                return Vector(-1, 0)
    return Vector(0, 0) # diagonal

def way_approx(vec):
    vec_val = val(vec)
    if abs(vec_val[0])>abs(vec_val[1]):
        if vec_val[0]>0:
            return Vector(1, 0)
        else:
            return Vector(-1, 0)
    else:
        if vec_val[1]>0:
            return Vector(0, 1)
        else:
            return Vector(0, -1)

variables = {}


def store_variable(symbol, value):  # put value in SI
    if isinstance(value, str):
        if LENGTH == extract_value_dim(value):
            unit = LENGTH_UNIT
        if INDUCTANCE == extract_value_dim(value):
            unit = INDUCTANCE_UNIT
        if CAPACITANCE == extract_value_dim(value):
            unit = CAPACITANCE_UNIT
        if RESISTANCE == extract_value_dim(value):
            unit = RESISTANCE_UNIT
        if FREQUENCY == extract_value_dim(value):
            unit = FREQUENCY_UNIT
        if DIMENSIONLESS == extract_value_dim(value):
            unit = DIMENSIONLESS_UNIT
        value = extract_value_unit(value, unit)
    variables[symbol] = value


class Vector(numpy.ndarray):

    """
    Vector is a custom 3D vector class, alowing for opperations optimized to
    interface properly with HFSS.
    The class can be instenciate as a 2D vector, how ever, it will effectively
    creat a 3D vector with 0 for z axis.
    """

    def __new__(cls, vec, vec_y=None, vec_z=None):
        """
        Init of the 3D vector:

            If vec_y, and vec_z are None, then vec must a len=2 or len=3 iterable.
            If vec_y is not None, and vec_z is, then creat a vector [vec, vec_y, 0].
            If vec_y and vec_z are not None, then creat a vector [vec, vec_y, vec_z].
        """

        if vec_y is not None:
            vec = [vec, vec_y, 0]
            if vec_z is not None:
                vec[2] = vec_z

        try:
            if not (len(vec) == 2 or len(vec) == 3):
                raise TypeError("vec can only be 2 or 3D, not %iD" % (len(vec)))
        except:
            raise TypeError("vec must be iterable")

        if len(vec) == 2:
            vec = [vec[0], vec[1], 0]

        obj = numpy.asarray(vec).view(cls)
        return obj

    @staticmethod
    def check(elt):

        """
        Utility function to check if an element is compatible with vectors
        opperations. It only requiers to be iterable and of len=3.

        Args:
            elt: The element to be tested

        Returns:
            Boolean, true if elt is compatible with Vector opperations, False
            otherwise.
        """

        try:
            return len(elt) == 3
        except:
            return False

    def __eq__(self, other):
        val_self = val(self)
        val_other = val(other)
        bool_result = (
            equal_float(val_self[0], val_other[0])
            and equal_float(val_self[1], val_other[1])
            and equal_float(val_self[2], val_other[2])
        )
        return bool_result

    def index(self, elt):
        val_self = val(self)
        val_elt = val(elt)
        for ii, item in enumerate(val_self):
            if item == val_elt:
                break
        else:
            return -1
        return ii

    #     def __add__(self, other):
    #         if Vector.check(other):
    #             return Vector([self[0]+other[0], self[1]+other[1], self[2]+other[2]])
    #         else:
    #             try:
    #                 return Vector([self[0]+other, self[1]+other, self[2]+other])
    #             except:
    #                 raise TypeError('Could not perform add operation')

    #     def __radd__(self, other):
    #         return self + other

    #     def __sub__(self, other)    :
    #         if Vector.check(other):
    #             return Vector([self[0]-other[0], self[1]-other[1], self[2]-other[2]])
    #         else:
    #             try:
    #                 Vector([self[0]-other, self[1]-other, self[2]-other])
    #             except:
    #                 raise TypeError('Could not perform sub operation')

    #     def __neg__(self):
    #         return Vector([-self[0], -self[1], -self[2]])

    #     def __rsub__(self, other):
    #         return -self + other

    #     def __mul__(self, other):
    #         if Vector.check(other):
    #             return Vector([self[0]*other[0], self[1]*other[1], self[2]*other[2]])
    #         else:
    #             try:
    #                 return Vector([other*self[0], other*self[1], other*self[2]])
    #             except:
    #                 raise TypeError('Could not perform mul operation')

    #     def __rmul__(self, other):
    #         return self * other

    #     def __truediv__(self, other):
    #         if Vector.check(other):
    #             return Vector([self[0]/other[0], self[1]/other[1], self[2]/other[2]])
    #         else:
    #             try:
    #                 return Vector([self[0]/other, self[1]/other, self[2]/other])
    #             except:
    #                 raise TypeError('Could not perform div operation')

    #     def __rtruediv__(self, other):

    #         self / other

    #     def dot(self, other):
    #         if Vector.check(other):
    #             return self[0]*other[0]+self[1]*other[1]+self[2]*other[2]
    #         else:
    #             raise TypeError('Could not perform dot operation')

    def cross(self, other):

        """
        This function returns the cross product beteween self and other.

        Args:
            other: of type Vector

        Returns:
            type Vector, self x other
        """

        if Vector.check(other) and Vector.check(other):

            return Vector(
                self[1] * other[2] - self[2] * other[1],
                -(self[0] * other[2] - self[2] * other[0]),
                self[0] * other[1] - self[1] * other[0],
            )
        else:
            raise TypeError("Could not perform dot operation")

    def scalar_cross(self, other, ref=None):

        """
        This function is a bit cryptic. It computes the signed magnitude of
        the cross product between self and other, assuming they both are in
        the plan orthogonal to ref.

        Args:
            other: a Vector
            ref: an other Vector, if None, assumed to be [0, 0, 1]

        Returns:
            dot((self x other), ref)
        """

        if ref is None:
            ref = Vector(0, 0, 1)

        if Vector.check(other) and Vector.check(ref):
            return self.cross(other).dot(ref)
        else:
            raise TypeError("Could not perform dot operation")

    def norm(self):
        return (self[0] ** 2 + self[1] ** 2 + self[2] ** 2) ** 0.5

    def abs(self):
        return Vector([abs(self[0]), abs(self[1]), abs(self[2])])

    def unit(self):
        norm = self.norm()
        return Vector([self[0] / norm, self[1] / norm, self[2] / norm])

    def orth(self):
        return Vector([-self[1], self[0]])

    #     def as_nda(self):

    #         return numpy.array([self[0], self[1], self[2]], dtype=object)

    def rot(self, other, ref=None):

        """
        This function is just completly cryptic, Ulysse wrote it a long time ago.
        Here is what it is doing: we assume that self is expressed in x=(100), y=(010), z=(001)
        This function returns the coordinates of self in x'=other,y'=-(other x ref), z'=ref
        In other words, this function computes a 3D change of coordinates.

        Note:
            This function has been writen assuming other and ref are given orthogonal.
            Hence, if not the case, it can have unexpected behaviors.

        Args:
            other: type Vector, the new x reference vector (x')
            ref: type Vector, the new z reference vector (z'), if None, taken to be (0,0,1)

        Returns:
            self expressed in the new coordinate system.
        """

        if ref is None:
            ref = Vector([0, 0, 1])
        else:
            ref = Vector(ref)

        other = Vector(other)

        if Vector.check(other) and Vector.check(ref):

            other = Vector(other).unit()
            ortho = -other.cross(ref)

            return (
                Vector([self.dot(other.refx()), self.dot(other.orth().refy()), 0]) * ref[2]
                + Vector([self.dot(other.orth().refx()), 0, self.dot(other.refz())]) * ref[1]
                + Vector([0, self.dot(other.refy()), self.dot(other.orth().refz())]) * ref[0]
            )
        else:
            raise TypeError("other must be a Vector")

    def px(self):
        return Vector([self[0], 0, 0])

    def py(self):
        return Vector([0, self[1], 0])

    def pz(self):
        return Vector([0, 0, self[2]])

    def refx(self, offset=0):
        return Vector([self[0], -self[1] + 2 * offset, self[2]])

    def refy(self, offset=0):
        return Vector([-self[0] + 2 * offset, self[1], self[2]])

    def refz(self, offset=0):
        return Vector([self[0], self[1], -self[2] + 2 * offset])


# if(__name__ == "__main__"):

#     x = Vector([1, 0, 0])
#     y = Vector([0, -1, 0])

#     print(x.rot(y))


def coor2angle(x, y=None):

    if y is None:
        x, y = x

    norm = (x ** 2 + y ** 2) ** 0.5

    if x != 0 and abs(y / x) < 1:
        angle = numpy.arcsin(y / norm)
        if x < 0:
            angle = numpy.pi - numpy.arcsin(y / norm)
    else:
        angle = numpy.arccos(x / norm)
        if y < 0:
            angle = -numpy.arccos(x / norm) + 2 * numpy.pi

    return angle % (2 * numpy.pi)


# if(__name__=="__main__"):

#     import matplotlib.pyplot as plt

#     plt.figure()

#     for theta in numpy.arange(0, 2*numpy.pi, 0.05):

#         x, y = numpy.cos(theta), numpy.sin(theta)

#         plt.plot(theta, coor2angle(x, y), 'o')

#     plt.show()
