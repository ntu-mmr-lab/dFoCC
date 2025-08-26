# The following script is mostly contributed by Prof. LOE
# 
# It automatically rotates all the dihedral angle within the residue

from numpy import pi
from Bio.PDB import Select, vectors, rotaxis
from Bio.PDB.vectors import calc_angle, calc_dihedral

from .dict_of_atom_vec_from_residue import dict_of_atom_vec_from_residue

from Bio.PDB.Structure import Structure
from Bio.PDB.Residue import Residue
from Bio.PDB.Atom import Atom
from Bio.PDB import Vector

def get_dihedrals(residue: Residue, alt_conf_id=None) -> tuple[list[float], list[list[str]]]:
    res = dict_of_atom_vec_from_residue(residue, alt_conf_id)
    dihedral_angles: list[float] = []
    dihedral_atoms: list[list[str]] = []

    # the atom combinations of dihedral angles:
    # [0]--[1]
    #         \ rot axis
    #          [2]--[3]
    # the actual atom that will be rotated is [3] (and the remaining side chain)
    # the bond between [1] and [2] is the rotation axis, and either of them can be considered as a fixed point that the rotating atoms base on
    dihedral_atoms_selections = [
        # the first dihedral angle
        ['N','CA','CB','CG'],['N','CA','CB','CG1'],['N','CA','CB','OG1'],['N','CA','CB','OG'],['N','CA','CB','SG'],
        # the second dihedral angle
        ['CA','CB','CG','CD'],['CA','CB','CG1','CD1'],['CA','CB','CG','CD1'],['CA','CB','CG','OD1'],['CA','CB','CG','ND1'],['CA','CB','CG','SD'],
        # the third dihedral angle
        ['CB','CG','CD','OE1'],['CB','CG','CD','CE'],['CB','CG','CD','NE'],['CB','CG','SD','CE'],
        # the fourth dihedral angle
        ['CG','CD','CE','NZ'],['CG','CD','NE','CZ'],
    ]

    for atoms in dihedral_atoms_selections:
        # for chi in chi_combis: for atoms in chi:
        try:
            dihedral_angles.append(calc_dihedral(res[atoms[0]], res[atoms[1]], res[atoms[2]], res[atoms[3]]) * 180 / pi)
            dihedral_atoms.append(atoms)
        except:
            pass

    return dihedral_angles, dihedral_atoms

def rotate_dihedrals(residue, *rxn_coords, alt_conf_id=None):
    """
    @param rxn_coords: the tenth degrees of changes of the rotations of the dihedral angles (at most 4 reaction coordinates)
    """
    # rxn_coords = [chi1, chi2, chi3, chi4]
    
    # 1. get the dihedral atoms
    vec = dict_of_atom_vec_from_residue(residue, alt_conf_id)
    dihedral_angles, dihedral_atoms = get_dihedrals(residue, alt_conf_id)
    if len(dihedral_atoms) < 1: # skip the boring residues that can barely rotate, like glycine and alanine
        return
    atom_list = list(vec.keys())

    # 2. first remove the atoms on the main chain, which will never be rotated
    atom_list = [elem for elem in atom_list if elem not in ['N','C','O']]

    # 3. then for each group of dihedral atoms, drop the 3 atoms that will not be rotated
    #    and also drop those previously fixed atoms
    rot_groups = []
    for chi in dihedral_atoms:
        atom_list = [elem for elem in atom_list if elem not in chi[0:3]]

        # BUG: without this fix, CG2 of isoleucine may rotate with the second dihedral angle rotation
        # Although it seems like Ile is the only counter example, this issue should be fixed generally
        # (e.g. to deal with hydrogen in the future)
        if 'CG1' not in atom_list and 'CG2' in atom_list:
            atom_list.remove('CG2')
        
        rot_groups.append(atom_list)

    # 4. rotate the remaining atoms for each group of dihedral atoms
    for index in range(len(rot_groups)):
        vec = dict_of_atom_vec_from_residue(residue, alt_conf_id) # refresh the residue with the new coordinate
        if rxn_coords[index] != 0.0:
            #             [1]
            # Recall that    \    is the rotation axis
            #                 [2]
            CBtoCAv = vec[dihedral_atoms[index][2]] - vec[dihedral_atoms[index][1]]
            rot1 = rotaxis(rxn_coords[index] * pi / 1800, CBtoCAv)

            # for the remaining atoms, rotate them base on the [1] atom
            for rot_atom_name in rot_groups[index]:
                atom_to_move = vec[rot_atom_name] - vec[dihedral_atoms[index][1]]
                move_to_ref = atom_to_move.left_multiply(rot1)
                if (rot_atom := residue[rot_atom_name]).is_disordered() != 0:
                    rot_atom = rot_atom.disordered_get(alt_conf_id)
                rot_atom.set_coord(move_to_ref+vec[dihedral_atoms[index][1]])

    return dihedral_angles
