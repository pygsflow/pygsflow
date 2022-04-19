"""
Class to trick python into installing without numba package support.

"""


def jit(func, **kwargs):
    def wrapper(*args, **kwargs):
        func(*args, **kwargs)

    return wrapper
