"""
Class to trick python into installing without ray package support.

"""


def jit(func):
    def wrapper(*args, **kwargs):
        func(*args, **kwargs)

    return wrapper


def get(*args):
    return


def put(*args):
    return
