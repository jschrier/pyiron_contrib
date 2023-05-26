import os
from pyiron_atomistics.atomistics.job.atomistic import AtomisticGenericJob
from pyiron_atomistics.atomistics.structure.atoms import ase_to_pyiron
from pyiron_base import DataContainer, state
from .FHIaims_input import write_input
from .FHIaims_output import collect_output

# everything from the Basics of Running FHI-aims tutorials "https://fhi-aims-club.gitlab.io/tutorials/basics-of-running-fhi-aims/1-Molecules/"

# these should be the default settings.  
# TODO: think more clearly about what behavior we actually want...I'm defaulting to molecular for now

# setup some defaults
input_dict = {
    "xc": "pbe", 
    "relativistic": "atomic_zora scalar",
    "relax_geometry": "bfgs 5e-3",
    "spin": "none",
    "precision": "light"
}

class FHIaims(AtomisticGenericJob):
    def __init__(self, project, job_name):
        super(FHIaims, self).__init__(project, job_name)
        self.__version__ = "0.1"
        self.__name__ = "FHIaims"
        self.input = DataContainer(input_dict, table_name="control")
        self._executable_activate()
        self._compress_by_default = True

        project_config = state.settings
        # TODO: for now we assume the first one is right; this is typically the user home version rather than conda...

        # I'm getting an error that project doesn't have this
        self.resource_directory = state.settings.resource_paths[0]

        self.species_directory = os.path.join(self.resource_directory, 'FHIaims',  'species_defaults', 'defaults_2020')



# TODO: Implement convenience functions and appropriate keywords
    def calc_minimize(self, keywords):
        print("TODO:  calc_minimize is not yet implemented; use calc_generic")

    def calc_static(self, keywords):
        print("TODO: calc_static is not yet implemented; use calc generic")

    def calc_magnetic(self, keywords):
        print("TODO: calc_magnetic is not yet implemented; use calc generic")

    def calc_magnetic(self, keywords):
        print("TODO: calc_magnetic is not yet implemented; use calc generic")


    def write_input(self):
        keywords = self.input.to_builtin()
        species = keywords["precision"]
        del keywords["precision"]

        write_input(
            structure=self.structure,
            keyword_definitions = keywords,
            control_filepath=os.path.join(self.working_directory, "control.in"),
            geometry_filepath=os.path.join(self.working_directory, "geometry.in"),
            species_directory = os.path.join(self.species_directory, species))


# TODO: implement output

    def collect_output(self):
        output_dict, output_dft_dict, meta_info_dict = collect_output      (
             working_directory=self.working_directory,
            FHI_output_file='aims.log')

        with self.project_hdf5.open("output") as hdf5_output:
            with hdf5_output.open("generic") as hdf5_generic:
                for k, v in output_dict.items():
                    hdf5_generic[k] = v
            with hdf5_output.open("dft") as hdf5_dft:
                for k, v in output_dft_dict.items():
                    hdf5_dft[k] = v
            hdf5_output["meta_info"] = meta_info_dict


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
