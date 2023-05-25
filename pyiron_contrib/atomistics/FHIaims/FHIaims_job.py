import os
from pyiron_atomistics.atomistics.job.atomistic import AtomisticGenericJob
from pyiron_atomistics.atomistics.structure.atoms import ase_to_pyiron
from pyiron_base import DataContainer
from .FHIaims_input import write_input
from .FHIaims_output import FHIaimsOutput

# everything from the Basics of Running FHI-aims tutorials "https://fhi-aims-club.gitlab.io/tutorials/basics-of-running-fhi-aims/1-Molecules/"

# these should be the default settings.  
# TODO: think more clearly about what behavior we actually want...I'm defaulting to molecular for now
# TODO:  "species" was the best name that I could find for this, but I don't know...

input_dict = {
    "xc": "pbe",
    "relativistic": "atomic_zora scalar",
    "relax_geometry": None,
    "relax_unit_cell": None,
    "spin": "none",  # "none" is defined string keyword, *not* python None
    "species": "light",
    "k_grid": None,
    "k_grid_density": None
}

class FHIaims(AtomisticGenericJob):
    def __init__(self, project, job_name):
        super(FHIaims, self).__init__(project, job_name)
        self.__version__ = "0.1"
        self.__name__ = "FHIaims"
        self.input = DataContainer(input_dict, table_name="control")
        self._executable_activate()
        self._compress_by_default = True

    def write_input(self):
        write_input(
            structure=self.structure,
            xc=self.input["xc"],
            relativistic=self.input["relativistic"],
            relax_geometry=self.input["relax_geometry"],
            relax_unit_cell=self.input["relax_unit_cell"],
            spin=self.input["spin"],
            species=self.input["species"],
            k_grid=self.input["k_grid"],
            k_grid_density=self.input["k_grid_density"],
            control_filepath=os.path.join(self.working_directory, "control.in"),
            geometry_filepath=os.path.join(self.working_directory, "geometry.in")
        )

# TODO: implement output

    def collect_output(self):
        output_file = os.path.join(self.working_directory, 'aims.log')
        if os.path.exists(output_file):
            output = FHIaimsOutput(filePath=output_file)
            output_dict = output.__dict__
            final_structure = ase_to_pyiron(output.final_structure)
            with self.project_hdf5.open("output/generic") as h5out:
                h5out["energy_tot"] = output_dict["energy"]
            with self.project_hdf5.open("output") as h5:
                final_structure.to_hdf(hdf=h5)
                del output_dict['final_structure'] # Remove ase atoms (incompatible with hd5)
                del output_dict['init_structure'] # Remove ase atoms (incompatible with hd5)
                h5["FHIaims"] = output_dict

# TODO: think about MPI and HDF

    def drop_status_to_aborted(self):
        """
        Workaround to prevent MPI error
        """
        self.status.finished = True

    def to_hdf(self, hdf=None, group_name=None):
        super(FHIaims, self).to_hdf(hdf=hdf, group_name=group_name)
        self._structure_to_hdf()
        with self.project_hdf5.open("input") as h5in:
            self.input.to_hdf(h5in)

    def from_hdf(self, hdf=None, group_name=None):
        super(FHIaims, self).from_hdf(hdf=hdf, group_name=group_name)
        self._structure_from_hdf()
        with self.project_hdf5.open("input") as h5in:
            self.input.from_hdf(h5in)
