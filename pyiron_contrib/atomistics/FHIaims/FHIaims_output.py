# TODO: implement output parsing; this is just the generic DFT stub from Michael Taylor 

"""
Output parsing module for FHI-aims
"""

# Imports
import re
import os
import numpy as np
#from ase.atoms import Atoms as ASEAtoms, Atom as ASEAtom
#import ase.units as units

FHI_OUT_KEYWORD_AIMS_UUID_TAG = "aims_uuid"
EXPORT_AIMS_UUID = FHI_OUT_KEYWORD_AIMS_UUID_TAG
EXPORT_FHI_AIMS_VERSION = "fhi-aims-version"
EXPORT_FHI_AIMS_PARAMETERS = "fhi-aims-parameters"
EXPORT_TOTAL_TIME = "total_time"


class InitialGeometryStreamParser:
    def __init__(self):
        self._input_geometry_block_flag = False
        self._inp_geom_unit_cell_flag = False

        # accumulator lists for current atomic structure
        self.lattice_vectors_lst = []
        self.positions_lst = []
        self.species_lst = []

        self._stop_processing = False

    def process_line(self, line):
        if self._stop_processing:
            return

        line = line.strip(" \t\n")

        if not line.startswith("|") and self._input_geometry_block_flag:
            self._input_geometry_block_flag = False
            self._stop_processing = True

        if line.startswith("Input geometry:"):
            self._input_geometry_block_flag = True

        if self._input_geometry_block_flag and line.startswith("|"):
            line = line.strip(" \t\n|")

            if len(self.lattice_vectors_lst) >= 3:
                self._inp_geom_unit_cell_flag = False

            if self._inp_geom_unit_cell_flag:
                self.lattice_vectors_lst.append([float(s) for s in line.split()[:3]])

            if line.startswith("Unit cell:"):
                self._inp_geom_unit_cell_flag = True

            if "Species" in line:
                line_tags = line.split()
                atom_positions = [float(s) for s in line_tags[3:6]]
                atom_species = line_tags[2]

                self.positions_lst.append(atom_positions)
                self.species_lst.append(atom_species)


class UpdateGeometryStreamParser:
    def __init__(self):
        self._atomic_structure_block_flag = False

        # accumulator lists for all atoomic structures
        self.lattice_vectors_traj = []
        self.positions_traj = []
        self.species_traj = []

        # accumulator lists for current atomic structure
        self._lattice_vectors_lst = []
        self._positions_lst = []
        self._species_lst = []

    def process_line(self, line):

        line = line.strip(" \t\n")

        if line.startswith("Updated atomic structure:"):  # or line.startswith("Final atomic structure:"):
            self._atomic_structure_block_flag = True

        if line.startswith("-------------") and self._atomic_structure_block_flag:
            self._atomic_structure_block_flag = False

            # save collected structure
            if len(self._lattice_vectors_lst) > 0:
                self.lattice_vectors_traj.append(self._lattice_vectors_lst)
            self.positions_traj.append(self._positions_lst)
            self.species_traj.append(self._species_lst)

            # reset accumulator lists
            self._lattice_vectors_lst = []
            self._positions_lst = []
            self._species_lst = []

        if self._atomic_structure_block_flag:

            if "lattice_vector" in line:
                self._lattice_vectors_lst.append([float(s) for s in line.split()[1:4]])

            if line.startswith("atom "):
                line_tags = line.split()
                atom_positions = [float(s) for s in line_tags[1:4]]
                atom_species = line_tags[4]

                self._positions_lst.append(atom_positions)
                self._species_lst.append(atom_species)


class EnergyForcesStressesStreamParser:
    def __init__(self):
        self.free_energies_list = []
        self.energies_corrected_list = []
        self.energies_uncorrected_list = []

        self.forces_lst = []
        self.stresses_lst = []

        self.block_flag = False
        self.force_block_flag = False
        self.stress_block_flag = False

        self.stress_line_counter = 0
        self.current_forces = []
        self.current_stresses = []

    def process_line(self, line):
        if "Energy and forces in a compact form:" in line:
            self.block_flag = True

        if self.block_flag and "------------------------------------" in line:
            self.block_flag = False
            self.force_block_flag = False
            self.forces_lst.append(self.current_forces)

        if self.block_flag and 'Total energy corrected        :' in line:
            E0 = float(line.split()[5])
            self.energies_corrected_list.append(E0)
        elif self.block_flag and 'Electronic free energy        :' in line:
            F = float(line.split()[5])
            self.free_energies_list.append(F)
        elif self.block_flag and 'Total energy uncorrected      :' in line:
            E_uncorr = float(line.split()[5])
            self.energies_uncorrected_list.append(E_uncorr)

        if self.block_flag and "Total atomic forces" in line:
            self.force_block_flag = True
            self.current_forces = []

        if self.force_block_flag and line.strip().startswith("|"):
            self.current_forces.append([float(f) for f in line.split()[-3:]])

        if "|              Analytical stress tensor" in line or "Numerical stress tensor" in line:
            self.stress_block_flag = True
            self.current_stresses = []
            self.stress_line_counter = 0

        if self.stress_block_flag:
            self.stress_line_counter += 1

        if self.stress_block_flag and self.stress_line_counter in [6, 7, 8]:
            sline = [float(f) for f in line.split()[2:5]]
            self.current_stresses.append(sline)

        if self.stress_line_counter > 8:
            self.stress_block_flag = False
            self.stress_line_counter = 0
            self.stresses_lst.append(self.current_stresses)


class MetaInfoStreamParser:
    _total_time_pattern = re.compile("Total time\s*:\s*([0-9.]*)\s*s")

    FHI_OUT_KEYWORD_TOTAL_TIME = "Total time"
    FHI_OUT_KEYWORD_AIMS_UUID_TAG = "aims_uuid"
    FHI_OUT_KEYWORD_VERSION_TAG = "Version"

    def __init__(self):
        self.version = None
        self.aims_uuid = None
        self.total_time = None

    def process_line(self, line):

        line = line.strip(" \t\n")

        if line.startswith(self.FHI_OUT_KEYWORD_VERSION_TAG):
            self.version = " ".join(line.split()[1:])
        elif line.startswith(self.FHI_OUT_KEYWORD_AIMS_UUID_TAG):
            self.aims_uuid = " ".join(line.split(":")[1:]).strip()

        if (self.FHI_OUT_KEYWORD_TOTAL_TIME in line) and (self.total_time is None):
            line = line.strip()
            res = self._total_time_pattern.findall(line)
            if len(res) > 0:
                self.total_time = res[0]

def collect_output(working_directory="", FHI_output_file="aims.log"):
    FHI_output_file = os.path.join(working_directory, FHI_output_file)

    init_geom_stream_parser = InitialGeometryStreamParser()
    upd_geom_stream_parser = UpdateGeometryStreamParser()
    efs_stream_parser = EnergyForcesStressesStreamParser()
    metainfo_stream_parser = MetaInfoStreamParser()

    with open(FHI_output_file, 'r') as f:
        for line in f:
            line = line.strip(" \t\n")
            if line.startswith("#"):
                continue

            init_geom_stream_parser.process_line(line)
            upd_geom_stream_parser.process_line(line)
            efs_stream_parser.process_line(line)
            metainfo_stream_parser.process_line(line)

    if len(efs_stream_parser.free_energies_list) == 0 or len(efs_stream_parser.forces_lst) == 0:
        raise ValueError("No free electronic energies found. Calculation could be broken")

    if len(init_geom_stream_parser.lattice_vectors_lst) > 0:
        lattice_vectors_traj = [
                                   init_geom_stream_parser.lattice_vectors_lst] + upd_geom_stream_parser.lattice_vectors_traj
    else:
        lattice_vectors_traj = upd_geom_stream_parser.lattice_vectors_traj

    positions_traj = [init_geom_stream_parser.positions_lst] + upd_geom_stream_parser.positions_traj

    if len(positions_traj) == 0:
        raise ValueError("No cells or positions info found. Calculation could be broken")

    output_generic_dict = {
        'cells': np.array(lattice_vectors_traj),
        'energy_pot': np.array(efs_stream_parser.free_energies_list),
        'energy_tot': np.array(efs_stream_parser.free_energies_list),
        'forces': np.array(efs_stream_parser.forces_lst),
        'positions': np.array(positions_traj),
        # 'steps'
        # 'temperature'
        # 'computation_time'
        # 'unwrapped_positions'
        # 'indices'
    }

    output_dft_dict = {
        'free_energy': np.array(efs_stream_parser.free_energies_list),
        'energy_corrected': np.array(efs_stream_parser.energies_corrected_list),
        'energy_uncorrected': np.array(efs_stream_parser.energies_uncorrected_list),
    }

    if len(efs_stream_parser.stresses_lst) > 0:
        output_generic_dict["stresses"] = efs_stream_parser.stresses_lst
        stresses = output_generic_dict["stresses"]
        output_generic_dict["pressures"] = np.array([-np.trace(stress) / 3.0 for stress in stresses])

    if len(output_generic_dict["cells"]) > 0:
        cells = output_generic_dict["cells"]
        output_generic_dict["volume"] = np.array([np.linalg.det(cell) for cell in cells])

    meta_info_dict = {}
    if metainfo_stream_parser.version is not None:
        meta_info_dict[EXPORT_FHI_AIMS_VERSION] = metainfo_stream_parser.version
    if metainfo_stream_parser.aims_uuid is not None:
        meta_info_dict[EXPORT_AIMS_UUID] = metainfo_stream_parser.aims_uuid
    if metainfo_stream_parser.total_time is not None:
        meta_info_dict[EXPORT_TOTAL_TIME] = metainfo_stream_parser.total_time

    return output_generic_dict, output_dft_dict, meta_info_dict