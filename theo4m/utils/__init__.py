"""

Theo4M (Theory for Materials).

Miscellaneous utils.

"""

__all__ = [
          "get_clusters_from_cell", 
          "cryst2angs", 
          "angs2cryst", 
          "dump_supercell", 
          "adjust_pbc_sa",
          "", 
          "write_ext_xyz"
          ]

def get_clusters_from_cell(atoms, out_cell, centered = True):
    from theo4m.utils import get_clusters_from_cell
    return get_clusters_from_cell(atoms, out_cell, centered = True)

def adjust_pbc(atoms):
    from theo4m.utils.supercell import adjust_pbc
    return adjust_pbc(atoms)

def adjust_pbc_sa(iat, atoms):
    from theo4m.utils.supercell import adjust_pbc_sa
    return adjust_pbc_sa(iat, atoms)

def cryst2angs(atoms):
    from theo4m.utils.supercell import cryst2angs
    return cryst2angs(atoms)

def angs2cryst(atoms):
    from theo4m.utils.supercell import angs2cryst
    return angs2cryst2(atoms)

def dump_supercell(fdump, adump, spc_symbs = [], output="xyz", coord_type="angs"):
    from theo4m.utils.supercell import dump_supercell

def get_item_ext_xyz(atoms, params, frac, 
                     extra_propstr = "", extra_prop = []):
    from theo4m.utils.supercell import get_item_ext_xyz
    return get_item_ext_xyz(atoms, params, frac, extra_propstr, extra_prop)

def write_ext_xyz(dir2w, file_name, atoms, params):
    from theo4m.utils.supercell import write_ext_xyz
    return write_ext_xyz(dir2w, file_name, atoms, params)


