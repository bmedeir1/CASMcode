from __future__ import (absolute_import, division, print_function, unicode_literals)
from builtins import *

import os, math, sys, shutil, gzip
import casm.quantumespresso
from casm.quantumespresso import qeio

class Relax(object):
    """The Relax class contains functions for setting up, executing, and parsing a Quantum Espresso relaxation.

        The relaxation is initialized in a directory containing Quantum Espresso input files, called 'relaxdir'.
        It then creates the following directory structure:
        .../relaxdir/
            run.0/
            run.1/
            ...
            run.final/



        'run.i' directories are only created when ready.
        'run.final' is a final constant volume run

        Contains:
            self.relaxdir  (.../relax)
            self.rundir    (list of .../relax/run.i)
            self.finaldir (.../relax/run.final)
    """

    def __init__(self, relaxdir=None, settings=None):
        """
        Construct a Quantum Espresso relaxation job object.

        Args:
            relaxdir:  path to quantum espresso relaxation directory
            settings:   dictionary-like object containing settings, or if None, it reads
                        the json file: .../relaxdir/relax.json

                possible settings keys are:
                    used by casm.quantumespresso.run() function:
                        "ncpus": number of ncpus to run mpi on
			"npar" or "ncore": number of ways to parallelize
                        "kpar": number of ways to parallelize k-points
                        "qe_cmd": (default, see quantumespresso.run) shell command to execute quantumespresso, or None to use default mpirun
                        "strict_kpoint": force strict copying of KPOINTS file, otherwise kpoints are scaled based on supercell size
                    used by not_converging():
                        "run_limit": (default 10) maximum number of runs to allow before setting status to "not_converging"

        """

        print("Constructing a Lammps Relax object")
        sys.stdout.flush()

        # store path to .../relaxdir, and create if not existing
        # if relaxdir == None:
        #     relaxdir = os.getcwd()
        # self.relaxdir = os.path.abspath(relaxdir)

        # print("  Relax directory:", self.relaxdir)
        # sys.stdout.flush()

        # find existing .../relaxdir/run.run_index directories, store paths in self.rundir list
        # self.rundir = []
        # self.errdir = []
        # self.update_rundir()
        # self.update_errdir()

        self.finaldir = os.path.join(self.relaxdir, "run.final")

        if settings == None:
            self.settings = dict()
        else:
            self.settings = settings

        self.initfilepath = None
        self.potfilepath = None
        self.setup()
        # set default settings:
        # if not "npar" in self.settings:
        #     self.settings["npar"] = None
        # if not "kpar" in self.settings:
        #     self.settings["kpar"] = None
        # if not "ncore" in self.settings:
        #     self.settings["ncore"] = None
        # if not "qe_cmd" in self.settings:
        #     self.settings["qe_cmd"] = None
        # if not "ncpus" in self.settings:
        #     self.settings["ncpus"] = None
        # if not "run_limit" in self.settings:
        #     self.settings["run_limit"] = 10
        # if not "nrg_convergence" in self.settings:
        #     self.settings["nrg_convergence"] = None
        # if not "compress" in self.settings:
        #     self.settings["compress"] = []
        # if not "err_types" in self.settings:
        #     self.settings["err_types"] = ['SubSpaceMatrixError']
        # ## added these because of key errors
        if not "extra_input_files" in self.settings:
            self.settings["extra_input_files"] = []
        # if not "move" in self.settings:
        #     self.settings["move"] = []
        # if not "copy" in self.settings:
        #     self.settings["copy"] = []
        # if not "remove" in self.settings:
        #     self.settings["remove"] = []
        # if not "compress" in self.settings:
        #     self.settings["compress"] = []
        # if not "backup" in self.settings:
        #     self.settings["backup"] = []
        # if not "initial" in self.settings:
        #     self.settings["initial"] = None
        if not "final" in self.settings:
            self.settings["final"] = None

        print("Lammps Relax object constructed\n")
        sys.stdout.flush()


    def setup(self, initdir, settings):
        """ mv all files and directories (besides initdir) into initdir """

        self.initfilepath = settings['initfile']
        self.potfilepath = settings['potentialfile']

        # print("Moving files into initial run directory:", initdir)
        # initdir = os.path.abspath(initdir)
        # for p in os.listdir(self.relaxdir):
        #     if (p in ([initfilename] + [potfilename] + self.settings["extra_input_files"] )) and (os.path.join(self.relaxdir, p) != initdir):
        #         os.rename(os.path.join(self.relaxdir,p), os.path.join(initdir,p))
        # print("")
        # sys.stdout.flush()
        #
        # # Keep a backup copy of the base Infile
        # shutil.copyfile(os.path.join(initdir,initfilename),os.path.join(self.relaxdir,initfilename + ".base"))

        # If an initial infile is called for, copy it in and set the appropriate flag
        # if (self.settings["initial"] != None) and (os.path.isfile(os.path.join(self.relaxdir,self.settings["initial"]))):
        #     new_values = qeio.Infile(os.path.join(self.relaxdir,self.settings["initial"])).tags
        #     qeio.set_infile_tag(new_values,initfilename,initdir)
        #     print("  Set Infile tags:", new_values, "\n")
        #     sys.stdout.flush()

    def run(self):
        print("Begin lammps relaxation run")
        sys.stdout.flush()
        # initfilename = self.settings["initfilename"]
        # potfilename = self.settings["potfilename"]
        # outfilename = self.settings["outfilename"]
        properties = casm.lammpspython.run(initfile=self.initfilepath, potentialfile=self.potfilepath)
        return properties

