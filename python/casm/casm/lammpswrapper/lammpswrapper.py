from __future__ import (absolute_import, division, print_function, unicode_literals)
from builtins import *

import os, shutil, six, re, subprocess, json
import warnings
from casm import lammpspython

class LAMMPSWrapperError(Exception):
    def __init__(self,msg):
        self.msg = msg

    def __str__(self):
        return self.msg

def read_settings(filename):
    try:
        with open(filename, 'rb') as file:
            settings = json.loads(file.read().decode('utf-8'))
    except (IOError, ValueError) as e:
        print("Error reading settings file:", filename)
        raise e

    return settings

def lammps_input_file_names(dirstruc, configname, clex, initfilename="init.mod", potfilename="potential.mod"):
    # Find required input files in CASM project directory tree

    myinfile = dirstruc.settings_path_crawl(initfilename, configname, clex)
    potentialfile = dirstruc.settings_path_crawl(potfilename, configname, clex)
    super_poscarfile = dirstruc.POS(configname)
    speciesfile = dirstruc.settings_path_crawl("SPECIES", configname, clex)

    # Verify that required input files exist
    if myinfile is None:
        raise LAMMPSWrapperError("lammps_input_file_names failed. No file found in CASM project matching: " + initfilename)
    if potentialfile is None:
        raise LAMMPSWrapperError("lammps_input_file_names failed. No file found in CASM project matching: " + potfilename)
    if super_poscarfile is None:
        raise LAMMPSWrapperError("lammps_input_file_names failed. No POS file found for this configuration.")
    if speciesfile is None:
        raise LAMMPSWrapperError("lammps_input_file_names failed. No SPECIES file found in CASM project.")

    return (myinfile, potentialfile, super_poscarfile, speciesfile)