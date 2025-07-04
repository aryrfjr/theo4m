"""

Loads information from extended .xyz files.

See also:
https://realpython.com/introduction-to-python-generators/#example-1-reading-large-files
https://wiki.fysik.dtu.dk/ase/ase/io/formatoptions.html?highlight=ase%20io%20extxyz#ase.io.extxyz.read_extxyz

Don't know how to use generators yet.
See: https://realpython.com/introduction-to-python-generators

from ase.io.extxyz import read_extxyz

"""
from numpy import array
from numpy import dot
from ase import Atom
from ase import Atoms
from ase.io.extxyz import key_val_str_to_dict # to read extra params

__all__ = [
        "read_extended_xyz",
        ]

def read_extended_xyz(xyzfile, frac = False):
    fileobj = open(xyzfile)
    lines = fileobj.readlines()
    fileobj.close()
    il = 0
    frames_atoms = []
    frames_params = []
    while il < len(lines):
        nat = int(lines[il])
        frame_params = key_val_str_to_dict(lines[il+1])
        cell = frame_params["Lattice"]
        if frac:
            geom_start = il + 1
            geom_stop = geom_start + nat
            species = [line.split()[0] for line in lines[geom_start:geom_stop]]
            geom = dot(array([[float(col) for col in line.split()[1:4]]
                for line in lines[geom_start:geom_stop]]), cell)
        else:
            geom_start = il + 2
            geom_stop = geom_start + nat
            species = [line.split()[0] for line in lines[geom_start:geom_stop]]
            geom = array([[float(col) for col in line.split()[1:4]]
                for line in lines[geom_start:geom_stop]])
        frame_atoms = Atoms(
        cell=cell,
        pbc=True,
        positions=geom,
        symbols=''.join(species))
        frames_atoms.append(frame_atoms)
        frames_params.append(frame_params)
        il += nat + 2
    return frames_atoms, frames_params

