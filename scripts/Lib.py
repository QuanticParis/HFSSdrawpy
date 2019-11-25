# -*- coding: utf-8 -*-
"""
Created on Fri Nov 15 11:51:09 2019

@author: Zaki
"""
from functools import wraps

def add_methods_from(*modules):
    print("hihi")
    def decorator(Class):
        for module in modules:
            for method in getattr(module, "__methods__"):
                setattr(Class, method.__name__, method)
        return Class
    return decorator

def register_method(methods):
    @wraps(methods)
    def register_method(method):
        methods.append(method)
        return method # Unchanged
    return register_method
