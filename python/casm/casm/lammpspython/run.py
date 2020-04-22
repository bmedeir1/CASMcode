from __future__ import (absolute_import, division, print_function, unicode_literals)
from builtins import *

import os, shutil, re, subprocess, sys, time, gzip, warnings
from casm.lammpspython.relaxcalc import RelaxCalc


class LAMMPSError(Exception):
    def __init__(self,msg):
        self.msg = msg

    def __str__(self):
        return self.msg


def run(initfile, potentialfile):
    print("Begin lammps run:")
    sys.stdout.flush()
    calc = RelaxCalc(initfile, potentialfile)
    result = calc.run_calc()
    if result is None:
        LAMMPSError("No results from lammps!")

    return result
