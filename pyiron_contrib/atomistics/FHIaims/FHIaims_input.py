"""
Input generation for FHI-aims
"""

# Imports
import re
from ase.atoms import Atoms as ASEAtoms, Atom as ASEAtom
import ase.units as units


# Functions
def write_input(structure, 
                xc = "pbe", 
                relativistic = "atomic_zora scalar", 
                relax_geometry = None,
                relax_unit_cell = None,
                spin = "none",
                species = "light",
                k_grid = None,
                k_grid_density = None,
                control_filepath = "control.in",
                geometry_filepath = "geometry.in",
                species_filepath = "species_defaults/defaults_2020/"):
    """write_input write input file for FHI-aims DFT code.

    Parameters
    ----------
    structure : ase.atoms.Atoms
        ASE atoms of structure you want to use.
    xc : str, optional
        dft xc functional to use, by default "pbe"
    relativistic: str, optional
        relativistic treatment to use, by default "atomic_zora scalar"
    relax_unit_cell: Boolean, optional
        should unit cell be relaxed in periodic systems, by default None
    ...
    # TODO:  write some more, once we decide on correct default choices
    """
    # Variables
    control_filelines = []
    geometry_filelines = []
    unique_species_filenames = {}

    control_filelines.append("xc {}\n".format(xc))
    control_filelines.append("relativistic {}\n".format(relativistic))
    control_filelines.append("spin {}\n".format(spin))

    # TODO: smarter checking that supplied keywords are valid

    if (relax_geometry != None):
        control_filelines.append("relax_geometry {}\n".format(relax_geometry))

    if (relax_unit_cell != None):
        control_filelines.append("relax_unit_cell {}\n".format(relax_unit_cell))

    # TODO: smarter checking between mutually exclusive options

    if (k_grid != None):
        control_filelines.append("k_grid {}\n".format(k_grid))

    if (k_grid_density != None):
        control_filelines.append("k_grid_density {}\n".format(k_grid_density))


    # Find the elements and append the potentials
    # TODO:  for now we only treat the defaults...
    # TODO:  It's dodgy to assume unix-style filepaths joining... 

    unique_species_filenames = set( 
        ["{}/{}/{:02d}_{}_default".format(
            species_filepath, species, Z, abbrev) 
            for Z, abbrev 
            in zip(structure.get_atomic_numbers(),
                   structure.get_chemical_symbols())
        ]
    )

    # write control.in file, along with the species info

    with open(control_filepath, 'w') as outFile:
        outFile.writelines(control_filelines)
        # append the species info.
        # safer to use python than to assume that cat exists
        for f in unique_species_filenames:
            with open(f, 'r') as inFile:
                outFile.writelines(inFile.readlines())
        

    # Format geometry.in file 
    # TODO: Add functionality for charge guesses, moments, etc.
    geometry_filelines = [
        "atom {} {} {} {}\n".format(pos[0], pos[1], pos[2], s)
        for s, pos in zip(structure.get_chemical_symbols(), structure.positions)
    ]    
    # Write geometry.in file
    with open(geometry_filepath, 'w') as outFile:
        outFile.writelines(geometry_filelines)



# Main
if __name__ == '__main__':
    pass
