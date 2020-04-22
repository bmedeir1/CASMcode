import os
import sys

import six
from lammps import lammps, PyLammps, IPyLammps
import numpy as np
import json
import re

class LammpsCalc(object):
    def __init__(self, initfilepath, potfilepath):
        # self.lmp = lammps()
        self.initpath = initfilepath
        self.potentialpath = potfilepath
        # self.scel = tokens[0]
        # self.id = tokens[1]
        self.pylmp = PyLammps()
        # self.configurename = configurename
        # self.clex = clex
        self.labels_dict = {}
        self.num_types_cell = None
        self.atoms_num = {}
        self.total_atoms = 0
        self.atoms = None
        self.atom_types = None
        self.cart_coors = None
        self.lattice = None
        self.basis_vec = None
        self.basis_sites = None
        self.cos_alpha = None
        self.cos_beta = None
        self.cos_gamma = None
        self.xlo = None
        self.xhi = None
        self.ylo = None
        self.yhi = None
        self.zlo = None
        self.zhi = None
        self.xy = None
        self.yz = None
        self.xz = None
        self.properties = {}
        self.json_props = None
        self.relaxed_energy = None

    # def load_files(self):
    #     self.lmp.file("init.mod")
    #     self.lmp.file("potential.mod")

    def read_prim(self):
        prim_path = "./prim.json"
        prim_file = open(prim_path).read()
        prim = json.JSONDecoder().decode(prim_file)
        print("prim:\n\n\n\n", prim)
        if prim['coordinate_mode'] == "Fractional":
            self.basis_vec = np.array(prim['lattice_vectors'])
            self.basis_sites = []
            # for basis in prim['']

    def read_structure(self):
        # Dir is defined at CASMcode/python/casm/casm/project/project.py
        # dir.__setattr__("POS", open(self.configurename + "/POS"))
        # pos_file = dir.POS(self.configname)
        pos_file = open(self.configurename + "/POS")
        # assert \
        self.configurename == pos_file.readline().strip()
        a0 = float(pos_file.readline().strip())
        print(a0)
        a1str = pos_file.readline().strip()
        a1 = a1str.split()
        a2str = pos_file.readline().strip()
        a2 = a2str.split()
        a3str = pos_file.readline().strip()
        a3 = a3str.split()

        self.lattice = np.array([a1, a2, a3], float) * a0
        print("lattice")
        print(self.lattice)

        labels = pos_file.readline().strip().split()
        nums = pos_file.readline().strip().split()

        for i in range(len(labels)):
            self.atoms_num[labels[i]] = nums[i]

        self.total_atoms = np.sum(np.array(nums, int))
        self.num_types_cell = len(labels)
        self.atoms = np.zeros((self.total_atoms, 3))
        self.atom_types = np.zeros(self.total_atoms)

        assert 'Direct' == pos_file.readline().strip()

        i = 0
        line = pos_file.readline()
        while line:
            line = line.strip()
            tokens = line.split()
            line = pos_file.readline()
            if len(tokens) <= 0:
                continue
            frac_coor = np.array(tokens[:3], float)
            # cart_coor = frac_coor @ lattice.T # convert fractional coordinates to cartesian
            # print(self.lmp.atoms)
            # self.lmp.atoms[i].position = tuple(cart_coor.reshape(1, -1)[0])
            self.atoms[i] = frac_coor
            # self.lmp.atoms[i].type = int(labels_dict[tokens[3].strip()])
            self.atom_types[i] = int(self.labels_dict[tokens[3].strip()])
            # print(atoms[i], atom_types[i])
            i += 1

    def create_structure(self):
        # Might not need this function anymore
        self.read_structure()
        # self.xlo = round(np.amin(cart_coors[:, 0]), 8)
        # # self.xhi = round(np.amax(cart_coors[:, 0]), 8)
        # self.ylo = round(np.amin(cart_coors[:, 1]), 8)
        # # self.yhi = round(np.amax(cart_coors[:, 1]), 8)
        # self.zlo = round(np.amin(cart_coors[:, 2]), 8)
        # # self.zhi = round(np.amax(cart_coors[:, 2]), 8)
        self.create_region()
        self.read_prim()
        self.create_atoms()


    def create_atoms(self):
        self.cart_coors = self.atoms @ self.lattice

        print('\nAtoms:\n', self.cart_coors[:10], "\n", self.atom_types[:10])

        print("test", self.cart_coors[0][0])
        print("atom types", self.atom_types)
        print("coors: ")
        print(self.cart_coors)
        print(self.pylmp.atoms.natoms)

        # Could optionally use lammps_create_atoms(void *, int, tagint *, int *, double *, double *,
        #                          imageint *, int)
        # But would still need to iterate through list and this method aides debugging.
        for i in range(self.cart_coors.shape[0]):
            # maybe should iterate through possible rounding values
            ### rounding must be consistent across all structure values
            print(
                'create_atoms ' + str(int(self.atom_types[i])) + ' single ' + str(round(self.cart_coors[i][0], 16)) + ' '
                + str(round(self.cart_coors[i][1], 16)) + ' ' + str(round(self.cart_coors[i][2], 16)))
            # self.pylmp.command(
            #     'create_atoms ' + str(int(self.atom_types[i])) + ' single ' + str(round(self.cart_coors[i][0], 16)) + ' '
            #     + str(round(self.cart_coors[i][1], 16)) + ' ' + str(round(self.cart_coors[i][2], 16)))
            self.pylmp.command(
                'create_atoms ' + str(int(self.atom_types[i])) + ' single ' + str(
                    self.cart_coors[i][0]) + ' '
                + str(self.cart_coors[i][1]) + ' ' + str(self.cart_coors[i][2]))

        print(self.cart_coors.shape[0])
        print(self.pylmp.atoms.natoms)



    def get_labels(self):
        raise NotImplementedError()

    def run_calc(self):
        raise NotImplementedError()

    def create_properties(self):
        self.properties['atom_type'] = list(self.labels_dict.keys())
        self.properties['atoms_per_type'] = [self.atoms_num[atom] for atom in self.atoms_num.keys()]
        # enforces same order as atom_type
        self.properties['coord_mode'] = 'cartesian'
        self.properties['is_complete'] = True if self.relaxed_energy is not None else False
        relaxed_basis = []
        relaxed_forces = []
        for i in range(self.total_atoms):
            relaxed_basis.append(list(self.pylmp.atoms[i].position))
            relaxed_forces.append(list(self.pylmp.atoms[i].force))

        self.properties['relaxed_basis'] = relaxed_basis
        self.properties['relaxed_energy'] = self.relaxed_energy
        self.properties['relaxed_forces'] = relaxed_forces
        self.properties['relaxed_lattice'] = [[self.pylmp.system.xhi-self.pylmp.system.xlo, 0, 0],
                                              [self.pylmp.eval('xy'), self.pylmp.system.yhi-self.pylmp.system.ylo, 0],
                                              [self.pylmp.eval('xz'), self.pylmp.eval('yz'), self.pylmp.system.zhi-self.pylmp.system.zlo]]
        j_encoder = json.JSONEncoder()
        self.json_props = j_encoder.encode(self.properties)
        return self.properties
        # self.properties['relaxed_lattice'] = self.lattice

    # def report_results(self):
    #
    #     try:
    #         os.mkdir(self.configurename + '/calctype.lammps')
    #     except:
    #         pass
    #
    #     outputfile = self.configurename + '/calctype.lammps/properties.calc.json'
    #
    #     with open(outputfile, 'wb') as file:
    #         # cls=noindent.NoIndentEncoder,
    #         file.write(six.u(json.dumps(self.properties, indent=4, sort_keys=True)).encode('utf-8'))
    #     print("Wrote " + outputfile)
    #     sys.stdout.flush()
    #     self.report_status()

    # def report_status(self):
    #     output = {'status': 'complete'}
    #     outputfile = self.configurename + '/calctype.lammps/status.json'
    #
    #     with open(outputfile, 'wb') as file:
    #         # cls=noindent.NoIndentEncoder,
    #         file.write(six.u(json.dumps(output, indent=4, sort_keys=True)).encode('utf-8'))
    #     print("Wrote " + outputfile)
    #     sys.stdout.flush()

    def create_region(self):
        # In the future will have to read prim and take max or min vector for each dimension
        # Need to add space for next position b/c otherwise periodicity will make atoms at border overlap
        # basis_vec = np.array([0.666666666667, 0.333333333334, 0.500000000000])
        # lat_vec = np.array([[3.184000000000, 0.000000000000, 0.000000000000],
        #                     [-1.592000000000, 2.757424885650, 0.000000000000],
        #                     [0.000000000000, 0.000000000000, 5.249000000000]])
        # dist = basis_vec @ lat_vec
        # self.xhi += dist[0]
        # self.yhi += dist[1]
        # self.zhi += dist[2]
        # a = 3.2094
        # self.lmp.lattice("hcp ", a)
        # self.xy = self.lattice[1, 0]
        # self.yz = self.lattice[2, 1]
        # self.xz = self.lattice[2, 0]
        self.xlo = 0
        self.ylo = 0
        self.zlo = 0
        self.cos_alpha = self.cos_angle(self.lattice[1], self.lattice[2])
        self.cos_beta = self.cos_angle(self.lattice[0], self.lattice[2])
        self.cos_gamma = self.cos_angle(self.lattice[0], self.lattice[1])
        self.xy = self.norm(self.lattice[1]) * self.cos_gamma
        self.xz = self.norm(self.lattice[2]) * self.cos_beta
        self.xhi = self.norm(self.lattice[0])
        self.yhi = np.sqrt(self.norm(self.lattice[1])**2 - self.xy**2)
        self.yz = ((self.norm(self.lattice[1]) * self.norm(self.lattice[2]) * self.cos_alpha) - (self.xy * self.xz)) / self.yhi
        self.zhi = np.sqrt((self.norm(self.lattice[2])**2) - (self.xz**2) - (self.yz**2))
        # self.xhi = self.lattice[0, 0] - self.lattice[1, 0] - self.lattice[2, 0]
        # self.yhi = self.lattice[1, 1] - self.lattice[0, 1] - self.lattice[2, 1]
        # self.zhi = self.lattice[2, 2] - self.lattice[0, 2] - self.lattice[1, 2]

        self.pylmp.atom_style('atomic')
        print('prism ', self.xlo, ' ', self.xhi, ' ', self.ylo, ' ', self.yhi, ' ', self.zlo, ' ', self.zhi, ' ', self.xy, ' ', self.xz, ' ', self.yz)
        self.pylmp.region('box prism ', self.xlo, ' ', self.xhi, ' ', self.ylo, ' ', self.yhi, ' ', self.zlo, ' ', self.zhi, ' ', self.xy, ' ', self.xz, ' ', self.yz)
        self.pylmp.command("boundary	p p p")
        # self.lmp.region("box block", self.xlo, self.xhi, self.ylo, self.yhi, self.zlo, self.zhi)
        # self.lmp.command("create_box {0:d} box".format(num_types))
        self.pylmp.command("create_box {0:d} box".format(len(self.labels_dict.keys())))
        print("bounding box", self.pylmp.system.xlo, self.pylmp.system.xhi, self.pylmp.system.ylo, self.pylmp.system.yhi, self.pylmp.system.zlo, self.pylmp.system.zhi)

        # self.lmp.create_atoms("1 box")
        for i in self.labels_dict.values():
            self.pylmp.command("mass {0:d} 1.0e-20".format(i))
            # self.pylmp.command("mass 2 1.0e-20")

    def norm(self, vec):
        return np.sqrt(np.sum(vec**2))

    def dot_prod(self, vec1, vec2):
        return np.sum(vec1 * vec2)

    def cos_angle(self, vec1, vec2):
        cos_theta = self.dot_prod(vec1, vec2) / (self.norm(vec1) * self.norm(vec2))
        assert cos_theta >= -1 and cos_theta <= 1
        if cos_theta < -1:
            cos_theta = -1
        elif cos_theta > 1:
            cos_theta = 1

        return cos_theta

    def visualize(self):
        try:
            print("creating", "./dumpfiles/" + self.scel)
            os.mkdir("./dumpfiles/" + self.scel)
        except:
            pass

        try:
            print("creating", "./dumpfiles/" + self.configurename)
            os.mkdir("./dumpfiles/" + self.configurename)
        except:
            pass

        self.pylmp.dump("dump1", "all", "atom", 1, "./dumpfiles/" + self.configurename + "/dump.atom") # , "zoom 2.0"
