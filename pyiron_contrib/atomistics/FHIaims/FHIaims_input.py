"""
Input generation for FHI-aims
"""

# Imports
import re
import os
from ase.atoms import Atoms as ASEAtoms, Atom as ASEAtom
import ase.units as units



# Functions

def write_input(structure, 
                keyword_definitions,
                control_filepath = "control.in",
                geometry_filepath = "geometry.in",
                species_directory = "~/pyiron/resources/species_defaults/defaults_2020/light"):
    """write_input write input file for FHI-aims DFT code.

    Parameters
    ----------
    structure : ase.atoms.Atoms
        ASE atoms of structure you want to use.
    keyword_definitions : dictionary
        Keys define the FHIaims keyword, values define settings
    species : string (default = "light")
        defines the type of species definition used 
    control_filepath : string; default `control.in`
    geometry_filepath : string; default `geometry.in`
    species_filepath : string; default `~/pyiron/resources/FHIAims/species_defaults/defaults_2020/light`

    """
    # Variables
    control_filelines = []
    geometry_filelines = []
    unique_species_filenames = {}

    # define the first part of control.in
    for key, value in keyword_definitions.items():
        control_filelines.append("{}\t{}\n".format(key, value))

    # Find the unique elements and identify the relevant filenames
    unique_species_filenames = set( 
        ["{:02d}_{}_default".format(Z, abbrev) 
            for Z, abbrev 
            in zip(structure.get_atomic_numbers(),
                   structure.get_chemical_symbols())
        ]
    )

    #these should go into special subdirectories; make sure they exist
    os.makedirs(os.path.dirname(control_filepath), 
                exist_ok=True)

    # write control.in file, along with the species info
    with open(control_filepath, 'w') as outFile:
        outFile.writelines(control_filelines)
        # append the species info.
        # safer to use python than to assume that cat exists
        for f in unique_species_filenames:
            species_file = os.path.join(species_directory, f)
            with open(species_file, 'r') as inFile:
                outFile.writelines(inFile.readlines())
        

    # Format geometry.in file 
    # TODO: Add functionality for charge guesses, moments, etc.
    geometry_filelines = [
        "atom {} {} {} {}\n".format(pos[0], pos[1], pos[2], s)
        for s, pos in zip(structure.get_chemical_symbols(), structure.positions)
    ]    

    os.makedirs(os.path.dirname(geometry_filepath), 
                exist_ok=True)
    # Write geometry.in file
    with open(geometry_filepath, 'w') as outFile:
        outFile.writelines(geometry_filelines)



# Main
if __name__ == '__main__':
    pass
