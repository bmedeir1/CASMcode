from casm.lammpspython import lammpscalc

class RelaxCalc(lammpscalc.LammpsCalc):
    def __init__(self, initfilepath, potfilepath):
        super().__init__(initfilepath, potfilepath)
        self.pylmp.file(initfilepath)
        self.get_labels()
        self.create_structure()
        # self.lmp.command("read_data         data.lmp")
        self.pylmp.file(potfilepath)
        # self.visualize()
        # self.run_calc()
        # self.create_properties()
        # self.report_results()
        # self.report_status()

    def get_labels(self):
        with open("./potential.mod") as fp:
            pair_coeff = None
            line = fp.readline()
            while line:
                tokens = line.split()
                if len(tokens) > 0 and "pair_coeff" == tokens[0]:
                    pair_coeff = tokens
                    break
                line = fp.readline()

            if pair_coeff is None:
                raise Exception("pair_coeff command not found in potential.mod or potential.mod does not exist!")

            potential_filename = None
            nums = []
            atoms = []
            i = 1
            last_word = tokens[0]
            for word in pair_coeff[1:]:
                if word.isnumeric():
                    nums.append(int(word))
                elif word == '*':
                    nums.append(i)
                    i += 1
                elif last_word.isnumeric() or last_word == '*':
                    potential_filename = word
                elif potential_filename is not None:
                    atoms.append(word)

                last_word = word

        self.labels_dict = {atoms[i]: nums[i] for i in range(min(len(nums), len(atoms)))}
        print(self.labels_dict)

    def run_calc(self):
        self.pylmp.variable("Ti          equal   550.0")
        self.pylmp.variable("Tfinal      equal   800.0")
        self.pylmp.variable("timestep    equal   0.004")
        # self.lmp.command("timestep          ${timestep}")
        self.pylmp.timestep("${timestep}")
        self.pylmp.variable("nequil      equal   5e4")
        self.pylmp.variable("nprod       equal   5e4")
        self.pylmp.variable(" P           equal   0.0")
        self.pylmp.variable("Tdamp       equal   50*dt")
        self.pylmp.variable("Pdamp       equal   500*dt")
        self.pylmp.variable("thermofreq  equal   500")
        self.pylmp.variable("dumpfreq    equal   20*${thermofreq}")
        self.pylmp.variable("seed        equal   1808071")
        self.pylmp.variable("nevery      equal   ${thermofreq}  ")
        self.pylmp.variable("nrepeat     equal   ${nprod}/${thermofreq}")
        self.pylmp.variable("nfreq       equal   ${nprod}")
        self.pylmp.group("grpA        type 1")
        self.pylmp.group("grpB        type 2")
        print(self.pylmp.groups)
        self.pylmp.variable("N1          equal count(grpA)")
        self.pylmp.variable("N2          equal count(grpB)")
        print(self.pylmp.variables['N1'].value, self.pylmp.variables['N2'].value)
        # print(self.lmp.eval('temp'))
        self.pylmp.variable("rhon equal atoms/vol")
        self.pylmp.variable("boxlx  equal lx")
        self.pylmp.variable("boxly  equal ly")
        self.pylmp.variable("boxlz  equal lz")
        self.pylmp.variable("step equal step")
        # self.lmp.variable("xcm  equal  xcm(all,x)")
        # self.lmp.variable("ycm  equal  xcm(all,y)")
        # self.lmp.variable("zcm  equal  xcm(all,z)")
        self.pylmp.compute("c1 all  temp/com")
        self.pylmp.variable("petot equal pe")

        # self.lmp.variable("dumppath", "string", "dumpfiles")
        # print(self.lmp.variables['dumppath'].value)
        # maybe fix later to check for existence of dumpfiles
        # process = subprocess.run(['mkdir', 'dumpfiles'],
        #                          stdout=subprocess.PIPE,
        #                          universal_newlines=True)
        # print(process)
        self.pylmp.compute("xu all property/atom xu")
        self.pylmp.compute("yu all property/atom yu")
        self.pylmp.compute("zu all property/atom zu")
        self.pylmp.compute("dummy1 all reduce max c_xu")
        self.pylmp.compute("dummy2 all reduce max c_yu")
        self.pylmp.compute("dummy3 all reduce max c_zu")
        self.pylmp.variable("dummy equal sum(c_xu)+sum(c_yu)+sum(c_zu)")
        self.pylmp.variable("dummy equal c_dummy1+c_dummy2+c_dummy3")
        self.pylmp.variable("s_pe   format   petot  %10.10f")
        # self.lmp.variable("variable  s_pe   format petot   %10.10f")
        self.pylmp.command("thermo ${thermofreq}")
        self.pylmp.command("thermo_style  custom step temp ke pe press pxx pyy pzz pxy pxz pyz lx ly lz vol enthalpy v_dummy")
        self.pylmp.command("thermo_modify  norm no flush yes")
        self.pylmp.fix("f3 all box/relax aniso $P")
        self.pylmp.command("minimize ${etol} ${ftol} ${maxiter} ${maxeval}")
        self.pylmp.unfix("f3")
        # print(self.lmp.eval('c1'))
        print("potential energy:", self.pylmp.eval('pe'), self.pylmp.variables['s_pe'].value)
        self.relaxed_energy = self.pylmp.eval('pe')
        print("bounding box", self.pylmp.system.xlo, self.pylmp.system.xhi, self.pylmp.system.ylo, self.pylmp.system.yhi, self.pylmp.system.zlo, self.pylmp.system.zhi)
        # print(self.pylmp.tilt)
        print(self.pylmp.eval('xy'))
        return self.create_properties()
    # def make_IPylammps(self):
    #     super().make_IPylammps()