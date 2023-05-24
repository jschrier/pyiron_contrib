"""
I/O module for DFT code (in this case ADF).
"""

# Imports
import re
from ase.atoms import Atoms as ASEAtoms, Atom as ASEAtom
import ase.units as units
import numpy as np


# Functions
def write_input(structure, total_charge=0, total_spin=0, functional="GGA PBE", basis="Type TZP", 
                coretype="None", opt=False, geo_iteration=100,
                e_iteration=300,  electronic_conv=1e-6, unrestricted=True,
                filePath='input.adf'):
    """write_input write input file for DFT (ADF) code.

    Parameters
    ----------
    structure : ase.atoms.Atoms
        ASE atoms of structure you want to use.
    total_charge : int, optional
        total charge on the system you are calculating, by default 0
    total_spin : int, optional
        spin (number of unpaired electrons) on the system, by default 0
    functional : str, optional
        dft xc functional to use, by default "GGA PBE"
    basis : str, optional
       basis set to use, by default "Type TZP"
    coretype : str, optional
        How large of frozen-core to use, by default "Large", "Large","Small", or "None"
    opt : bool, optional
        perform Geo optimization, by default False
    geo_iteration : int, optional
        how many steps to take in geometry optimization, default 100
    e_iteration : int, optional
       max number of electronic (SCF) iterations per step, by default 300
    electronic_conv : float, optional
        electronic convergence tolerance, by default 1e-6
    unrestricted : bool, optional
        perform calculation as unrestricted?, by default True
    filePath : str, optional
        path to input file to write, by default 'input.adf'
    """
    # Variables
    file_lines = []

    # Format ATOMS section
    atom_lines = [
        "{} {} {} {}\n".format(s, pos[0], pos[1], pos[2])
        for s, pos in zip(structure.get_chemical_symbols(), structure.positions)
    ]

    # Setup input lines
    file_lines.append("System\n")
    file_lines.append('Charge {}\n'.format(total_charge))
    file_lines.append("Atoms\n")

    for atom in atom_lines:
        file_lines.append(atom)

    file_lines.append('End\n')
    file_lines.append('End\n')

    if opt:
        file_lines.append('Task GeometryOptimization\n')
        file_lines.append('GeometryOptimization\n')
        file_lines.append('\tMaxIterations {}\n'.format(geo_iteration))
        file_lines.append('End\n')
    else:
        file_lines.append('Properties\n')
        if 'hybrid' not in functional.lower():
            file_lines.append('\tGradients Yes\n')
        file_lines.append('End\n')
        file_lines.append('Task SinglePoint\n')
    file_lines.append('Engine ADF\n')
    file_lines.append('Basis\n')
    file_lines.append('\t{}\n'.format(basis))
    file_lines.append('\tCore {}\n'.format(coretype))
    file_lines.append('End\n')
    file_lines.append('XC\n')
    file_lines.append('\t{}\n'.format(functional))
    file_lines.append('End\n')
    file_lines.append('SAVE TAPE21\n')
    file_lines.append('SCF\n')
    electronic_conv_2nd_str = "%.e" % 1e-3
    electronic_conv_str = "%.e" % electronic_conv
    file_lines.append('\titerations {}\n'.format(e_iteration))
    file_lines.append('\tconverge ' + electronic_conv_str + ' ' + electronic_conv_2nd_str + '\n')
    file_lines.append('End\n')
    file_lines.append('BeckeGrid\n')
    file_lines.append('End\n')
    if unrestricted:
        file_lines.append('Unrestricted Yes\n')
        file_lines.append('SpinPolarization {}\n'.format(total_spin))
    else:
        file_lines.append('Unrestricted No\n')
    file_lines.append('EndEngine\n')

    # Write input file
    with open(filePath, 'w') as outFile:
        outFile.writelines(file_lines)


class DFTOutput:
    """ Output parser
    """
    def __init__(self, filePath):
        """__init__ execute the following block when Class is called.

        Parameters
        ----------
        filePath : str
            path to output file.
        """
        with open(filePath,'r') as file1:
            self.fileLines = file1.readlines()
        with open(filePath,'r') as file1:
            self.fileStr = file1.read()
        self.read_properties()
        del self.fileStr
        del self.fileLines

    def read_adf_version(self):
        """read_adf_version get the adf version
        """
        for line in self.fileLines:
            sline = line.strip().split()
            if ('RunTime' in line) and ((sline[0] == 'AMS') or (sline[0] == 'ADF')):
                self.version = sline[1]
                break

    def read_energy(self):
        """
        Read energ(ies) from output file.
        """

        # Variables
        energy = "None"
        # Find energy in a.u.
        for line in reversed(self.fileLines):
            # Find energy
            if 'Energy (hartree)' in line:
                energy = float(line.split()[-1]) * units.Hartree
                break

        self.energy = energy

    def read_geo_steps_energies(self):
        """read_geo_steps read in the number of geometry steps taken in output file and energies
        """
        steps = 0
        energies = []
        try:
            for line in self.fileLines:
                if "current energy" in line:
                    steps += 1
                    sline = line.strip().split()
                    if 'current' == sline[0]:
                        energies.append(float(sline[2])*units.Hartree)
            self.steps = steps
            if len(energies) > 0:
                self.energies = energies
                self.energy = energies[-1]
            else:
                self.read_energy()
                self.energies = [self.energy]
        except:
            self.steps = 0
            self.energies = []
            self.read_energy()

    def read_coordinates(self):
        """
        Read structure from geometry optimization output file.
        """

        # Variables
        atomsList = []
        atomIdx = -1

        read = False
        # Get atoms
        tmp_atomlist = []
        for line in self.fileLines[atomIdx + 1:]:
            # Format line
            line = line.strip().split()

            if read and (len(line) > 0):
                # End if "Index" encountered
                if (line[0] == "Index"):
                    continue

                # Gather atoms
                symbol = line[1]
                tmp_atomlist.append([symbol, float(line[2]), float(line[3]), float(line[4])])
            else:
                if read:
                    read = False
                    atomsList.append(tmp_atomlist.copy())
                    tmp_atomlist = []
            if (len(line) > 0):
                if line[0] == 'Atoms':
                    read = True

        if len(atomsList) > 1:
            return atomsList[:-1] # Last one is repeat of final geometry
        elif len(atomsList) == 1:
            return atomsList
        else:
            return None

    def read_forces(self):
        """
        Read structure from geometry optimization output file.
        """
        # Variables
        forcesList = []
        atomIdx = -1
        read = False
        # Get atoms
        tmp_atomlist = []
        
        for line in self.fileLines[atomIdx + 1:]:
            # Format line
            line = line.strip().split()

            if read and (len(line) > 0) and (skip == 0):
                # Skip lines "---------"
                if (len(line) != 5):
                    continue

                # Gather atoms
                symbol = line[1]
                try: #Convert gradients to forces!
                    tmp_atomlist.append([symbol, -float(line[2]), -float(line[3]), -float(line[4])])
                except ValueError:
                    tmp_atomlist.append([symbol,np.inf,np.inf,np.inf])
            elif read and (skip > 0):
                skip -= 1
            else:
                if read:
                    read = False
                    forcesList.append(tmp_atomlist.copy())
                    tmp_atomlist = []
            if len(line) > 0:
                if (line[0] == 'Energy') and (line[1] == 'gradients'):
                    read = True
                    skip=3
        if len(forcesList) > 1:
            return forcesList[:-1] # Last one is repeat of final geometry
        elif len(forcesList) == 1:
            return forcesList
        else:
            return None

    def get_rxyz_string(self,atoms,forces,energy,skip_forces=False):
        """dump the ase_atoms to rxyz file string"""
        natom = len(atoms)
        ss = ''
        # write an xyz file first
        try:
            ss += "%d\n\n" % natom
            for atom in atoms:
                symb = atom[0]
                coords = atom[1:]
                ss += "%2s %12.6f %12.6f %12.6f\n" % (symb, coords[0], coords[1], coords[2])
            ss += "\n"
            # then force
            if not skip_forces:
                ss += "FORCES\n"
                for force in forces:
                    symb = force[0]
                    atom_force = force[1:]
                    ss += "%3s %22.14e %22.14e %22.14e\n" % (symb, atom_force[0], atom_force[1], atom_force[2])
                ss += '\n'

            # then different properties
            ss += "ENERGY %22.14f\n" % (energy/units.Hartree) # rxyz has hartree
            ss += 'SOURCE %s\n' % ('ADF' + self.version)
            return ss
        except:
            return None


    def read_trajectory_rxyz(self):
        """read_trajectory_rxyz read optimization trajectory into rxyz files with energy/forces
        """
        atomsList = self.read_coordinates()
        forcesList = self.read_forces()
        energies = self.energies
        skip_forces = False
        if forcesList is None:
            skip_forces = True 
            forcesList = [[]]*len(atomsList)
        self.rxyz_trajectory = []
        for i,tatoms in enumerate(atomsList):
            if (i+1 > len(forcesList)) or (i+1 > len(energies)):
                break
            rxyz = self.get_rxyz_string(atoms=tatoms,forces=forcesList[i],energy=energies[i],skip_forces=skip_forces)
            self.rxyz_trajectory.append(rxyz)
        

    def read_homo_lumo_gap(self):
        """
        Read homo-lumo gap information.
        """
        try:
            # Variables
            hlDict = {'spin_1': {'homo': 0.0, 'lumo': 0.0}, 'spin_2': {'homo': 0.0, 'lumo': 0.0}}
            reKey_homo = "\s*HOMO.*\n"
            reKey_lumo = "\s*LUMO.*\n"

            # Apply regex
            result_homo = (re.findall(reKey_homo, self.fileStr))[-2:]
            result_lumo = (re.findall(reKey_lumo, self.fileStr))[-2:]

            # Format result into dictionary
            hlDict['spin_1']['homo'] = float(result_homo[0].strip().split()[4])*units.Hartree
            hlDict['spin_1']['lumo'] = float(result_lumo[0].strip().split()[4])*units.Hartree
            hlDict['spin_2']['homo'] = float(result_homo[1].strip().split()[4])*units.Hartree
            hlDict['spin_2']['lumo'] = float(result_lumo[1].strip().split()[4])*units.Hartree

            self.homo_lumo = hlDict
        except:
            self.homo_lumo = {}


    def get_final_geo(self):
        """convert_rxyz_ase take in rxyz string and convert to ase atoms
        """
        try:
            rxyzstr = self.trajectory[-1]
            lines = rxyzstr.split('\n')
            atoms = []
            for line in lines:
                sline = line.split()
                if 'FORCES' in line:
                    break
                elif len(sline) == 4:
                    symbol = sline[0]
                    coords = (float(sline[1]),float(sline[2]),float(sline[3]))
                    if len(atoms) == 0:
                        atoms = ASEAtoms([ASEAtom(symbol,coords)])
                    else:
                        atoms.append(ASEAtom(symbol,coords))
            self.final_structure = atoms
        except:
            self.final_structure = None


    def _check_error(self, errorKey):
        # Look for error key
        for line in self.fileLines:
            if (line == errorKey):
                return True
        return False

    def error_scf_incomplete(self):
        """
        Identify if there was an error with SCF convergence.
        """

        # Variables
        errorKey = "SCF did not converge\n"

        return self._check_error(errorKey=errorKey)


    def error_geomopt_incomplete(self):
        """
        Identify if there was an error with the geometry optimization.
        """

        # Variables
        errorKey = "Geometry optimization did NOT converge\n"

        return self._check_error(errorKey=errorKey)


    def error_job_quit(self):
        """
        Identify if there was an error with the job running made by the scheduler.
        """

        # Variables
        errorKey = "Process received SIGTERM\n"

        return self._check_error(errorKey=errorKey)

    def error_thermo_incomplete(self):
        """
        Identify if there was an error with the vibrational analysis
        """

        errorKey = "" ######## NEED TO FILL

        return self._check_error(errorKey=errorKey)


    def error_check(self):
        """
        Get errors from output file.
        """

        # Variables
        errorcodeList = []

        # Get errors
        errorDict = {"SCF_ERROR": self.error_scf_incomplete(),
                    "GEOMOPT_ERROR": self.error_geomopt_incomplete(),
                    "JOB_ERROR": self.error_job_quit(),
                    "THERMO_ERROR":self.error_thermo_incomplete()}

        # Populate errorcode list
        for key in errorDict:
            if errorDict[key] is True:
                errorcodeList.append(key)

        self.errorcodeList = errorcodeList



    def read_properties(self):
        """
        Read properties from ADF output file.
        """

        # Adf version
        self.read_adf_version()

        # HOMO-LUMO
        self.read_homo_lumo_gap()

        # Steps
        self.read_geo_steps_energies()

        # Get trajectory:
        self.read_trajectory_rxyz()

        # Get final geometry as ase_atoms
        self.get_final_geo()

        # Check Errors
        self.error_check()


# Main
if __name__ == '__main__':
    pass
