"""

Theo4M (Theory for Materials).

Machine learning package (MacLe4M).

"""

__all__ = [
           "topological_index", 
           "thermodynamic_index", 
           "move_atom", 
           "unmove_atom", 
           "swap_atoms",
           "write_ASAP_ext_xyz",
           "average",
           "best_match",
           "min_max",
           "get_dist"
        ]

def topological_index(atoms, tbonds, tcoords, spc_symbs, short_bond_thr):
    from theo4m.macle4m.dbutils import topological_index
    return topological_index(atoms, tbonds, tcoords, spc_symbs, short_bond_thr)

def thermodynamic_index(atoms):
    from theo4m.macle4m.dbutils import thermodynamic_index
    return thermodynamic_index(atoms)

def move_atom(iat, atoms, dd):
    from theo4m.macle4m.dbutils import move_atom
    return move_atom(iat, atoms, dd)

def unmove_atom(iat, atoms, ux, uy, uz):
    from theo4m.macle4m.dbutils import unmove_atom
    return unmove_atom(iat, atoms, ux, uy, uz)

def swap_atoms(iata, iatb, atoms):
    from theo4m.macle4m.dbutils import swap_atoms
    return swap_atom(iata, iatb, atoms)

def write_ASAP_ext_xyz(dir2w, file_name, atoms_set, params_set, lmax, nmax, pcs_SOAPs = True):
    from theo4m.macle4m.dbutils import write_ASAP_ext_xyz
    return write_ASAP_ext_xyz(dir2w, file_name, atoms_set, params_set, lmax, nmax, pcs_SOAPs = True)

def min_max(atoms_A, atoms_B, S_A, S_B, wtd = "MIN", sort = False):
    from theo4m.macle4m.kernelutils import min_max
    return min_max(atoms_A, atoms_B, S_A, S_B, wtd, sort)

def average(structures, soapsl):
    from theo4m.macle4m.kernelutils import average
    return average(structures, soapsl)

def best_match(structures, soapsl):
    from theo4m.macle4m.kernelutils import best_match
    return best_match(structures, soapsl)

def boxplot(structures, soapsl):
    from theo4m.macle4m.kernelutils import boxplot
    return boxplot(structures, soapsl)

def get_dist(atomA, atomB, cell):
    from theo4m.macle4m.dbutils import get_dist
    return get_dist(atomA, atomB, cell)

def write_QE_pwi(dir2w, optid, atoms, kpts, ecut, calculation):
    from theo4m.macle4m.dbutils import write_QE_pwi
    return write_QE_pwi(dir2w, optid, atoms, kpts, ecut, calculation)

