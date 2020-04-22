from __future__ import (absolute_import, division, print_function, unicode_literals)
from builtins import *

import json
import math
import os
import re
import six
import sys
import warnings

from casm.misc import noindent
from casm.project import DirectoryStructure, ProjectSettings
from casm.lammpswrapper.lammpswrapper import LAMMPSWrapperError, lammps_input_file_names, read_settings
from casm import lammpspython


class Relax(object):
    """The Relax class contains functions for setting up, executing, and parsing a Quantum Espresso relaxation.

        The relaxation creates the following directory structure:
        config/
          calctype.name/
              run.0/
              run.1/
              ...
              run.final/

        'run.i' directories are only created when ready.
        'run.final' is a final constant volume run {"calculation":'relax'}

        This automatically looks for Quantum Espresso settings files in .../settings/calctype.name,
        where '...' is the nearest parent directory of 'self.configdir' in the CASM project repository

        Contains:
            self.configdir (.../config)
            self.calcdir   (.../config/calctype.name)

            self.settings = dictionary of settings for job submission and the relaxation, see qewrapper.read_settings

            self.auto = True if using prisms_jobs module's JobDB to manage jobs
            self.sort = True if sorting atoms in POSCAR by type
    """

    def __init__(self, configdir=None):
        """
        Construct a Quantum Espresso relaxation job object.

        Args:
            configdir: path to configuration
            auto: True if using prisms_jobs module's JobDB to manage jobs

        """
        if configdir == None:
            configdir = os.getcwd()

        print("Working on directory " + str(configdir))

        # get the configname from the configdir path
        _res = os.path.split(configdir)
        self.configname = os.path.split(_res[0])[1] + "/" + _res[1]
        print("  Configuration:", self.configname)

        print("Reading CASM settings")
        self.casm_settings = ProjectSettings(configdir)
        if self.casm_settings == None:
            raise LAMMPSWrapperError("Not in a CASM project. The file '.casm' directory was not found.")

        self.casm_directories = DirectoryStructure(configdir)

        print("Constructing a CASM lammps Relax object")
        sys.stdout.flush()

        print("  Setting up directories")
        sys.stdout.flush()

        # store path to .../config, if not existing raise
        self.configdir = os.path.abspath(configdir)
        if not os.path.isdir(self.configdir):
            raise lammpspython.run.LAMMPSError(
                "Error in casm.lammps.Relax: Did not find directory: " + self.configdir)
            sys.stdout.flush()

        # store path to .../config/calctype.name, and create if not existing
        self.calcdir = self.casm_directories.calctype_dir(self.configname, self.casm_settings.default_clex)
        try:
            os.mkdir(self.calcdir)
        except:
            pass

        # read the settings json file
        print("  Reading relax.json settings file")
        sys.stdout.flush()
        setfile = self.casm_directories.settings_path_crawl("relax.json", self.configname,
                                                            self.casm_settings.default_clex)
        if setfile == None:
            raise LAMMPSWrapperError("Could not find \"relax.json\" in an appropriate \"settings\" directory")
            sys.stdout.flush()

        else:
            print("Using " + str(setfile) + " as settings...")

        self.settings = read_settings(setfile)
        self.infiles = None
        self.run_settings = {}

        print("Lammps Relax object constructed\n")
        sys.stdout.flush()
        self.properties =None
        self.setup_done = False
        self.submitted = False
        self.run_done = False


    def run_settings(self):
        self.run_settings['initfile'], self.run_settings['potentialfile'], \
        self.run_settings['poscar'], self.run_settings['species'] = self.infiles
        return self.run_settings


    def setup(self):

        initfilename = self.settings["initfilename"] if "initfilename" in self.settings.keys() else "init.mod"
        potfilename = self.settings["potfilename"] if "potfilename" in self.settings.keys() else "potential.mod"

        self.infiles = lammps_input_file_names(self.casm_directories, self.configname, self.casm_settings.default_clex,
                                     initfilename=initfilename, potfilename=potfilename)
        self.setup_done = True
        return self.infiles

    def submit(self):
        self.run()
        self.submitted = True

    def run(self):
        if not self.setup_done:
            self.setup()
        # construct the Relax object
        relax = lammpspython.Relax(self.calcdir, self.run_settings())
        self.properties = relax.run()
        if self.properties is not None:
            self.finalize(self.properties)
        else:
            LAMMPSWrapperError("Lammps relaxation failed")
        self.run_done = True


    def report_status(self, status, failure_type=None):
        if not self.run_done:
            self.run()
        """Report calculation status to status.json file in configuration directory.

        Args:
            status: string describing calculation status. Currently used values are
                 not_submitted
                 submitted
                 complete
                 failed
             failure_type: optional string describing reason for failure. Currently used values are
                 unknown
                 electronic_convergence
                 run_limit"""

        output = dict()
        output["status"] = status
        if failure_type is not None:
            output["failure_type"] = failure_type

        outputfile = os.path.join(self.calcdir, "status.json")
        with open(outputfile, 'wb') as file:
            file.write(
                six.u(json.dumps(output, cls=noindent.NoIndentEncoder, indent=4, sort_keys=True)).encode('utf-8'))
        print("Wrote " + outputfile)
        sys.stdout.flush()

    def finalize(self, output):
        if not self.run_done:
            self.run()
        # output = self.properties(qedir, outfilename)
        outputfile = os.path.join(self.calcdir, "properties.calc.json")
        with open(outputfile, 'wb') as file:
            file.write(
                six.u(json.dumps(output, cls=noindent.NoIndentEncoder, indent=4, sort_keys=True)).encode('utf-8'))
        print("Wrote " + outputfile)
        sys.stdout.flush()
        self.report_status('complete')

    def is_submitted(self):
        return self.submitted

    def report(self):
        self.finalize()

    @staticmethod
    def properties(calcdir):
        """Report results to properties.calc.json file in configuration directory, after checking for electronic convergence."""
        propfile = os.path.join(calcdir, "properties.calc.json")
        props = open(propfile).read()
        prop_dict = json.JSONDecoder().decode(props)
        return prop_dict


