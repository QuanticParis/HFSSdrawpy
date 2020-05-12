# -*- coding: utf-8 -*-
"""
Created on Thu Nov 14 13:46:42 2019

@author: Zaki
"""

from sympy.parsing import sympy_parser
from pint import UnitRegistry
import numpy

ureg = UnitRegistry()
Q = ureg.Quantity

LENGTH = '[length]'
INDUCTANCE = '[length] ** 2 * [mass] / [current] ** 2 / [time] ** 2'
CAPACITANCE = '[current] ** 2 * [time] ** 4 / [length] ** 2 / [mass]'
RESISTANCE = '[length] ** 2 * [mass] / [current] ** 2 / [time] ** 3'

LENGTH_UNIT = 'meter'
INDUCTANCE_UNIT = 'nH'
CAPACITANCE_UNIT = 'fF'
RESISTANCE_UNIT = 'ohm'

### List handling
# Useful function to manipulate to_move entities and ports
def find_last_list(list_entities):
    # return the last list of a set of nested lists

    if isinstance(list_entities, list):
        if len(list_entities)==0:
            return list_entities
        else:
            if isinstance(list_entities[-1], list):
                return find_last_list(list_entities[-1])
            else:
                return list_entities
    else:
        raise TypeError('There are no list')

def find_penultimate_list(list_entities):
    # return the last list of a set of nested lists
    if isinstance(list_entities, list):
        if len(list_entities)==0:
            return False
        else:
            if isinstance(list_entities[-1], list):
                if len(list_entities[-1])==0:
                    return list_entities
                else:
                    if isinstance(list_entities[-1][-1], list):
                        return find_penultimate_list(list_entities[-1])
                    else:
                        return list_entities
            else:
                return False
    else:
        raise TypeError('There are no list')

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
        pass#raise TypeError('Argument is not a list')

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

def check_name(_class, name):
    end = ''
    for ii, char in enumerate(name[::-1]):
        if char.isdigit():
            end+=char
        else:
            break
    else:
        ii += 1
    if end == '':
        radical = name
        number = 0
    else:
        radical = name[:-ii]
        number = int(end[::-1])
    new_name = name
    while(new_name in _class.dict_instances.keys()):
        number+=1
        new_name = radical+str(number)
    if new_name != name:
        print("%s: changed '%s' name into '%s'"%(_class.__name__, name, new_name))
    return new_name

### Litteral Expressions

def equal_float(float1, float2):
    if float1!=0:
        rel_diff = abs((float1-float2)/float1)
        if rel_diff<1e-5:
            return True
        else:
            return False

    elif float2!=0:
        rel_diff = abs((float1-float2)/float2)
        if rel_diff<1e-5:
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
    return str(Q(expr).dimensionality)

def parse_entry(*entries, marker=True):
    #should take a list of tuple of list... of int, float or str...
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
                raise TypeError('Not foreseen type: %s'%(type(entry)))
    if len(parsed)==1 and marker:
        return parsed[0]
    else:
        return parsed

def rem_unit(other):
    try:
        value = extract_value_unit(other, LENGTH_UNIT)
        return value
    except Exception:
        return other

def var(x):
    if isinstance(x, str):
        return VariableString(simplify_arith_expr(x))
    return x

def _val(elt):
    if isinstance(elt, VariableString):
        return elt.value()
    else:
        return elt

def val(*entries, marker=True):
    #should take a list of tuple of list... of int, float or str...
    parsed = []
    for entry in entries:
        if not isinstance(entry, list) and not isinstance(entry, tuple):
            parsed.append(_val(entry))
        else:
            if isinstance(entry, list):
                if isinstance(entry, Vector):
                    parsed.append(Vector(val(*entry, marker=False)))
                else:
                    parsed.append(val(*entry, marker=False))
            elif isinstance(entry, tuple):
                parsed.append(tuple(val(*entry, marker=False)))
            else:
                raise TypeError('Not foreseen type: %s'%(type(entry)))
    if len(parsed)==1 and marker:
        return parsed[0]
    else:
        return parsed

def way(vec):
    if vec[1] != 0:
        if abs(vec[0]/vec[1])<1e-2:
            if vec[1]>0:
                return Vector(0,1)
            elif vec[1]<0:
                return Vector(0,-1)
    if vec[0] != 0 :
        if abs(vec[1]/vec[0])<1e-2:
            if vec[0]>0:
                return Vector(1,0)
            elif vec[0]<0:
                return Vector(-1,0)

class VariableString(str):
    # TODO: What happen with a list (Vector in our case)
    variables = {}
    instances = {}

    def __new__(cls, name, *args, **kwargs):
        # explicitly only pass value to the str constructor
        return super(VariableString, cls).__new__(cls, name)

    def __init__(self, name, value=None):
        # ... and don't even call the str initializer
        if value is not None: # do not result from a computation
            VariableString.store_variable(self, value)
            VariableString.instances[self] = self

    @classmethod
    def store_variable(cls, name, value):  # put value in SI
        if not isinstance(value, VariableString):
            if isinstance(value, (float, int)):
                cls.variables[name] = value
            elif isinstance(value, str):
                if LENGTH == extract_value_dim(value):
                    cls.variables[name] = extract_value_unit(value,
                                                              LENGTH_UNIT)
                if INDUCTANCE == extract_value_dim(value):
                    cls.variables[name] = extract_value_unit(value,
                                                              INDUCTANCE_UNIT)
                if CAPACITANCE == extract_value_dim(value):
                    cls.variables[name] = extract_value_unit(value,
                                                              CAPACITANCE_UNIT)
                if RESISTANCE == extract_value_dim(value):
                    cls.variables[name] = extract_value_unit(value,
                                                          RESISTANCE_UNIT)
            else:
                raise TypeError('Type is not expected')
        else:
            cls.variables[name] = value

    def value(self):
        # print(self.variables)
        try:
            _value = float(eval(str(sympy_parser.parse_expr(str(sympy_parser.parse_expr(self, self.variables)), self.variables))))
        except Exception:
            msg = ('Parsed expression contains a string which does '
                             'not correspond to any design variable')
            raise ValueError(msg)
        return _value


    def __add__(self, other):
        other = rem_unit(other)
        if other=="'":
            return super(VariableString, self, other).__add__()
        return var("(%s) + (%s)" % (self, other))

    def __radd__(self, other):
        other = rem_unit(other)
        if other=="'":
            return super(VariableString, self, other).__radd__()
        return var("(%s) + (%s)" % (other, self))

    def __sub__(self, other):
        other = rem_unit(other)
        return var("(%s) - (%s)" % (self, other))

    def __rsub__(self, other):
        other = rem_unit(other)
        return var("(%s) - (%s)" % (other, self))

    def __mul__(self, other):
        other = rem_unit(other)
        return var("(%s) * (%s)" % (self, other))

    def __rmul__(self, other):
        other = rem_unit(other)
        return var("(%s) * (%s)" % (other, self))

    def __div__(self, other):
        other = rem_unit(other)
        return var("(%s) / (%s)" % (self, other))

    def __rdiv__(self, other):
        other = rem_unit(other)
        return var("(%s) / (%s)" % (other, self))

    def __truediv__(self, other):
        other = rem_unit(other)
        return var("(%s) / (%s)" % (self, other))

    def __rtruediv__(self, other):
        other = rem_unit(other)
        return var("(%s) / (%s)" % (other, self))

    def __pow__(self, other):
        other = rem_unit(other)
        return var("(%s) ** (%s)" % (self, other))

#    def __rpow__(self, other):
#        other = self.rem_unit(other)
#        return var("(%s) ** (%s)" % (other, self))

    def __neg__(self):
        return var("-(%s)" % self)

    def __abs__(self):
        return var("abs(%s)" % self)


class Vector(list):

    """
    Vector is a custom 3D vector class, alowing for opperations optimized to
    interface properly with HFSS.
    The class can be instenciate as a 2D vector, how ever, it will effectively
    creat a 3D vector with 0 for z axis.
    """

    def __init__(self, vec, vec_y=None, vec_z=None):

        """
        Init of the 3D vector:

            If vec_y, and vec_z are None, then vec must a len=2 or len=3 iterable.
            If vec_y is not None, and vec_z is, then creat a vector [vec, vec_y, 0].
            If vec_y and vec_z are not None, then creat a vector [vec, vec_y, vec_z].
        """

        if vec_y is not None:
            vec = [vec, vec_y, 0]
            if(vec_z is not None):
                vec[2] = vec_z

        try:
            if(not (len(vec)==2 or len(vec)==3)):
                raise TypeError('vec can only be 2 or 3D, not %iD' % (len(vec)))
        except:
            raise TypeError('vec must be iterable')

        if(len(vec) == 2):
            vec = [vec[0], vec[1], 0]

        super().__init__(parse_entry(vec))

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
            return len(elt)==3
        except:
            return False

    def __eq__(self, other):
        return (equal_float(val(self[0]), val(other[0])) and
                equal_float(val(self[1]), val(other[1])) and
                equal_float(val(self[2]), val(other[2])))

    def __add__(self, other):
        if Vector.check(other):
            return Vector([self[0]+other[0], self[1]+other[1], self[2]+other[2]])
        else:
            try:
                return Vector([self[0]+other, self[1]+other, self[2]+other])
            except:
                raise TypeError('Could not perform add operation')

    def __radd__(self, other):
        return self + other

    def __sub__(self, other):
        if Vector.check(other):
            return Vector([self[0]-other[0], self[1]-other[1], self[2]-other[2]])
        else:
            try:
                Vector([self[0]-other, self[1]-other, self[2]-other])
            except:
                raise TypeError('Could not perform sub operation')

    def __neg__(self):
        return Vector([-self[0], -self[1], -self[2]])

    def __rsub__(self, other):
        return -self + other

    def __mul__(self, other):
        if Vector.check(other):
            return Vector([self[0]*other[0], self[1]*other[1], self[2]*other[2]])
        else:
            try:
                return Vector([other*self[0], other*self[1], other*self[2]])
            except:
                raise TypeError('Could not perform mul operation')

    def __rmul__(self, other):
        return self * other

    def __truediv__(self, other):
        if Vector.check(other):
            return Vector([self[0]/other[0], self[1]/other[1], self[2]/other[2]])
        else:
            try:
                return Vector([self[0]/other, self[1]/other, self[2]/other])
            except:
                raise TypeError('Could not perform div operation')

    def __rtruediv__(self, other):

        self / other

    def dot(self, other):
        if Vector.check(other):
            return self[0]*other[0]+self[1]*other[1]+self[2]*other[2]
        else:
            raise TypeError('Could not perform dot operation')

    def cross(self, other, ref=None):

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

        if(ref is None):
            ref = Vector(0, 0, 1)

        if(Vector.check(other) and Vector.check(ref)):
            if(ref[0] == 0 and ref[1] == 0):
                return (self[0]*other[1]-self[1]*other[0])*ref[2]
            elif(ref[0] == 0 and ref[2] == 0):
                return -(self[0]*other[2]-self[2]*other[0])*ref[1]
            elif(ref[1] == 0 and ref[2] == 0):
                return (self[1]*other[2]-self[2]*other[1])*ref[0]
            else:
                raise TypeError('ref Vectore must be along x, y or z')
        else:
            raise TypeError('Could not perform dot operation')

    def norm(self):
        return (self[0]**2+self[1]**2+self[2]**2)**0.5

    def abs(self):
        return Vector([abs(self[0]), abs(self[1]), abs(self[2])])

    def unit(self):
        norm = self.norm()
        return Vector([self[0]/norm, self[1]/norm, self[2]/norm])

    def orth(self):
        return Vector([-self[1], self[0]])

    def rot(self, other, ref=None):

        '''
        This function is just completly cryptic, I wrote it a long time ago,
        dont understand it anymore XD.
        '''

        if(ref is None):
            ref = Vector([0, 0, 1])

        if(Vector.check(other) and Vector.check(ref)):

            unitOther = Vector(other).unit()

            if(ref[0]==0 and ref[1]==0):
                return Vector([self.dot(unitOther.refx()), self.dot(unitOther.orth().refy()), 0])
            elif(ref[0]==0 and ref[2]==0):
                return Vector([self.dot(unitOther.orth().refx()), 0, self.dot(unitOther.refz())])
            elif(ref[1]==0 and ref[2]==0):
                return Vector([0, self.dot(unitOther.refy()), self.dot(unitOther.orth().refz())])
            else:
                raise TypeError('ref Vectore must be along x, y or z')
        else:
            raise TypeError('other must be a Vector')

    def px(self):
        return Vector([self[0], 0, 0])

    def py(self):
        return Vector([0, self[1], 0])

    def pz(self):
        return Vector([0, 0, self[2]])

    def refx(self, offset=0):
        return Vector([self[0], -self[1]+2*offset, self[2]])

    def refy(self, offset=0):
        return Vector([-self[0]+2*offset, self[1], self[2]])

    def refz(self, offset=0):
        return Vector([self[0], self[1], -self[2]+2*offset])


if(__name__ == "__main__"):

    x = Vector([1, 0, 0])
    y = Vector([0, -1, 0])

    print(x.rot(y))
