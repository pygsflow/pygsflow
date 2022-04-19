"""
Class to trick python into installing without numba package support.

"""
from functools import wraps


def jit(**kwargs):
    def wrapper(func):
        @wraps(func)
        def decorated_func(*args, **kwargs):
            return func(args, **kwargs)
        return decorated_func
    return wrapper
