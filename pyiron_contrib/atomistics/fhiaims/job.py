import os
from pyiron_atomistics.atomistics.job.atomistic import AtomisticGenericJob
from pyiron_atomistics.atomistics.structure.atoms import ase_to_pyiron
from pyiron_base import DataContainer
from .dft_io import (write_input, DFTOutput)


input_dict = {
    "total_spin": 0,
    "total_charge": 0,
    "functional": "GGA PBE",
    "basis": "Type TZP",
    "coretype": "Large",
    "opt": False,
    "geo_iteration":100,
    "e_iteration": 300,
    "electronic_conv": 1e-6,
    "unrestricted":True
}

class DFTjob(AtomisticGenericJob):
    def __init__(self, project, job_name):
        super(DFTjob, self).__init__(project, job_name)
        self.__version__ = "0.1"
        self.__name__ = "DFT"
        self.input = DataContainer(input_dict, table_name="control")
        self._executable_activate()
        self._compress_by_default = True

    def write_input(self):
        write_input(
            structure=self.structure,
            total_spin=int(self.input["total_spin"]),
            total_charge=int(self.input["total_charge"]),
            functional=self.input["functional"],
            basis=self.input["basis"],
            dispersion=self.input["dispersion"],
            opt=self.input["opt"],
            geo_iteration=int(self.input["geo_iteration"]),
            e_iteration=int(self.input["e_iteration"]),
            electronic_conv=float(self.input["electronic_conv"]),
            unrestricted=self.input['unrestricted'],
            filePath=os.path.join(self.working_directory, "input.ams")
        )

    def collect_output(self):
        output_file = os.path.join(self.working_directory, 'output.ams')
        if os.path.exists(output_file):
            output = DFTOutput(filePath=output_file)
            output_dict = output.__dict__
            final_structure = ase_to_pyiron(output.final_structure)
            with self.project_hdf5.open("output/generic") as h5out:
                h5out["energy_tot"] = output_dict["energy"]
            with self.project_hdf5.open("output") as h5:
                final_structure.to_hdf(hdf=h5)
                del output_dict['final_structure'] # Remove ase atoms (incompatible with hd5)
                del output_dict['init_structure'] # Remove ase atoms (incompatible with hd5)
                h5["dft"] = output_dict

    def drop_status_to_aborted(self):
        """
        Workaround to prevent MPI error
        """
        self.status.finished = True

    def to_hdf(self, hdf=None, group_name=None):
        super(DFTjob, self).to_hdf(hdf=hdf, group_name=group_name)
        self._structure_to_hdf()
        with self.project_hdf5.open("input") as h5in:
            self.input.to_hdf(h5in)

    def from_hdf(self, hdf=None, group_name=None):
        super(DFTjob, self).from_hdf(hdf=hdf, group_name=group_name)
        self._structure_from_hdf()
        with self.project_hdf5.open("input") as h5in:
            self.input.from_hdf(h5in)
