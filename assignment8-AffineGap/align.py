#!/usr/bin/env python3
########################################################################
# File:align.py
#  executable: align.py
# Purpose:
#
# Author:       Akshar Lohith
# Group(sounding boards):
# History:      AL 11/14/2016 Created
#
########################################################################
import sys

class LocalAligner(object):
    """docstring for LocalAligner."""
    def __init__(self, arg):
        self.arg = arg

    def score(self):
        pass

    def align(self):
        pass

    def traceback(self):
        pass

class GlobalAligner(object):
    """docstring for GlobalAligner."""
    def __init__(self, arg):
        self.arg = arg

    def score(self):
        pass

    def align(self):
        pass

    def traceback(self):
        pass

def main():
    pass

if __name__ == "__main__": # if program is launched alone, this is true and is exececuted. if not, nothing is\
# executedf rom this program and instead objects and variables are made availableto the program that imports this.
    main();
    raise SystemExit
