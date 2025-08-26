from Bio.PDB.Residue import Residue

def refine_occupancy(residue: Residue, delta_occ_percent: int):
    for atom in residue.get_atoms():
        if atom != '' and atom.is_disordered() != 0:
            conformation_a = atom.disordered_get('A')
            conformation_b = atom.disordered_get('B')
            conformation_a.set_occupancy(round(conformation_a.occupancy + delta_occ_percent / 100, 2))
            conformation_b.set_occupancy(round(conformation_b.occupancy - delta_occ_percent / 100, 2))

def change_occupancy(residue: Residue, occ_percent: int):
    for atom in residue.get_atoms():
        if atom != '' and atom.is_disordered() != 0:
            conformation_a = atom.disordered_get('A')
            conformation_b = atom.disordered_get('B')
            conformation_a.set_occupancy(round(occ_percent / 100, 2))
            conformation_b.set_occupancy(round(occ_percent / 100, 2))
