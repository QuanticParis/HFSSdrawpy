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

def entity_kwargs(kwargs, keys):
    entity_kwargs = {}
    for key in keys:
        if key in kwargs.keys():
            entity_kwargs[key] = kwargs[key]
    return entity_kwargs

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

def to_move(cls):
    if cls.instances_to_move is None:
        cls.instances_to_move = []
        return cls.instances_to_move
    else:
        last_list = find_last_list(cls.instances_to_move)
        last_list.append([])
        return last_list

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
    i = 0
    new_name = name
    while(new_name in _class.dict_instances.keys()):
        new_name = name+'_'+str(i)
        i+=1
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

def parse_entry(*entries):
    #should take a list of tuple of list... of int, float or str...
    parsed = []
    for entry in entries:
        if not isinstance(entry, list) and not isinstance(entry, tuple):
            parsed.append(extract_value_unit(entry, LENGTH_UNIT))
        else:
            if isinstance(entry, list):
                if isinstance(entry, Vector):
                    parsed.append(Vector(parse_entry(*entry)))
                else:
                    parsed.append(parse_entry(*entry))
            elif isinstance(entry, tuple):
                parsed.append(tuple(parse_entry(*entry)))
            else:
                raise TypeError('Not foreseen type: %s'%(type(entry)))
    if len(parsed)==1:
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

def val(*entries):
    #should take a list of tuple of list... of int, float or str...
    parsed = []
    for entry in entries:
        if not isinstance(entry, list) and not isinstance(entry, tuple):
            parsed.append(_val(entry))
        else:
            if isinstance(entry, list):
                if isinstance(entry, Vector):
                    parsed.append(Vector(val(*entry)))
                else:
                    parsed.append(val(*entry))
            elif isinstance(entry, tuple):
                parsed.append(tuple(val(*entry)))
            else:
                raise TypeError('Not foreseen type: %s'%(type(entry)))
    if len(parsed)==1:
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
        # print(self.instances)
        # print(self)
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
    def __init__(self, vec, vec_y=None):
        if vec_y is not None:
            vec = [vec, vec_y]
        super().__init__(parse_entry(vec))

    def check(self, elt):
        return isinstance(elt, (list, tuple, numpy.ndarray))

    def check_nb(self, nb):
        return isinstance(nb, float) or isinstance(nb, int) or isinstance(nb, VariableString)

    def __eq__(self, other):
        return (equal_float(val(self[0]), val(other[0])) and
            equal_float(val(self[1]), val(other[1])))

    def __add__(self, other):
        if self.check(other):
            return Vector([self[0]+other[0], self[1]+other[1]])
        else:
            raise TypeError('Could not perform add operation')

    def __radd__(self, other):
        return self + other

    def __sub__(self, other):
        if self.check(other):
            return Vector([self[0]-other[0], self[1]-other[1]])
        else:
            raise TypeError('Could not perform sub operation')

    def __neg__(self):
        return Vector([-self[0], -self[1]])

    def __rsub__(self, other):
        return -self + other

    def __mul__(self, other):
        if self.check(other):
            return Vector([self[0]*other[0], self[1]*other[1]])
        elif self.check_nb(other):
            return Vector([other*self[0], other*self[1]])
        else:
            raise TypeError('Could not perform mul operation')

    def __rmul__(self, other):
        return self * other

    def __truediv__(self, other):
        if self.check(other):
            return Vector([self[0]/other[0], self[1]/other[1]])
        elif self.check_nb(other):
            return Vector([self[0]/other, self[1]/other])
        else:
            raise TypeError('Could not perform div operation')

    def __rtruediv__(self, other):
        if self.check(other):
            return Vector([other[0]/self[0], other[1]/self[1]])
        elif self.check_nb(other):
            return Vector([other/self[0], other/self[1]])
        else:
            raise TypeError('Could not perform rdiv operation')

    def dot(self, other):
        if self.check(other):
            return self[0]*other[0]+self[1]*other[1]
        else:
            raise TypeError('Could not perform dot operation')

    def cross(self, other):
        if self.check(other):
            return self[0]*other[1]-self[1]*other[0]
        else:
            raise TypeError('Could not perform dot operation')

    def norm(self):
        return (self[0]**2+self[1]**2)**0.5

    def abs(self):
        return Vector([abs(self[0]), abs(self[1])])

    def unit(self):
        norm = self.norm()
        return Vector([self[0]/norm, self[1]/norm])

    def orth(self):
        return Vector([-self[1], self[0]])

    def rot(self, other):
        '''
        Inputs:
        -------
        other: vector

        Returns
        -------
        vector: rotation around z of self by an angle given by other w.r.t to x
        '''
        if self.check(other):
            unitOther = Vector(other).unit()
            return Vector([self.dot(unitOther.refx()), self.dot(unitOther.orth().refy())])
        else:
            raise TypeError('Could not perform rdiv operation')

    def px(self):
        return Vector([self[0], 0])

    def py(self):
        return Vector([0, self[1]])

    def refx(self, offset=0):
        return Vector([self[0], -self[1]+2*offset])

    def refy(self, offset=0):
        return Vector([-self[0]+2*offset, self[1]])

