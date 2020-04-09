# -*- coding: utf-8 -*-
"""
Created on Thu Nov 14 13:46:42 2019

@author: Zaki
"""

from sympy.parsing import sympy_parser
from pint import UnitRegistry

import vector

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
                if isinstance(entry, vector.Vector):
                    parsed.append(vector.Vector(parse_entry(*entry)))
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

class VariableString(str):
    # TODO: What happen with a list (Vector in our case)
    variables = {}

    def __new__(cls, name, *args, **kwargs):
        # explicitly only pass value to the str constructor
        return super(VariableString, cls).__new__(cls, name)

    def __init__(self, name, value=None):
        # ... and don't even call the str initializer
        if value is not None:
            VariableString.store_variable(self, value)

    @classmethod
    def store_variable(cls, name, value):  # put value in SI
        if not isinstance(value, VariableString):
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
            cls.variables[name] = value

    def value(self):
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


