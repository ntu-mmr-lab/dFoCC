from typing import Literal
from numpy import pi, result_type
from Bio.PDB import Select, vectors, rotaxis
from Bio.PDB.vectors import calc_angle, calc_dihedral

from .dict_of_atom_vec_from_residue import dict_of_atom_vec_from_residue

from Bio.PDB.Structure import Structure
from Bio.PDB.Chain import Chain
from Bio.PDB.Residue import Residue
from Bio.PDB.Atom import Atom, DisorderedAtom
from Bio.PDB import Vector

ResseqAtomPair = tuple[int, Literal['N']|Literal['CA']|Literal['C']]

def calculate_main_chain_general(chain: Chain, res_from: ResseqAtomPair, alt_conf_id: str | None = None) -> tuple[str, float]:
    """
    Compare the rotation angle of a main chain rotation reaction coordinate
    by calulating the dihedral angle near the first fixed atom
    Note that by inputing a residue i with atom N, CA and C
    the result is identical to ω(i-1), φ(i) and ψ(i) respectively

    @param chain: the chain of modified residues
    @param res_from: tuple (resseq, atom_name) of the first fixed atom
    """
    # get the actual index of the residue
    # as there are HETATM residues that have strange convention of index numbers
    resseq_list: list[int] = [full_id[1] for full_id in chain.child_dict.keys()]
    index = resseq_list.index(res_from[0])

    res_prev = dict_of_atom_vec_from_residue(chain.child_list[index-1])
    res_curr = dict_of_atom_vec_from_residue(chain.child_list[index])
    res_next = dict_of_atom_vec_from_residue(chain.child_list[index+1])
    match res_from[1]:
        case 'N': # i.e. ω(i-1)
            return f'omega({res_from[0]-1})', calc_dihedral(
                res_prev['CA'],
                res_prev['C'],
                res_curr['N'],
                res_curr['CA'],
            ) * 180 / pi
        case 'CA': # i.e. φ(i)
            return f'phi({res_from[0]})', calc_dihedral(
                res_prev['C'],
                res_curr['N'],
                res_curr['CA'],
                res_curr['C'],
            ) * 180 / pi
        case 'C': # i.e. ψ(i)
            return f'psi({res_from[0]})', calc_dihedral(
                res_curr['N'],
                res_curr['CA'],
                res_curr['C'],
                res_next['N'],
            ) * 180 / pi
        # case 'O':
        #     return calc_dihedral(
        #         res_curr['CA'],
        #         res_curr['C'],
        #         res_curr['O'],
        #         res_next['N'],
        #     ) * 180 / pi
        case _:
            assert res_from[1] in ('N', 'CA', 'C'), 'Atom is not in the main chain or is not valid'
def rotate_main_chain_general(chain: Chain, res_from: ResseqAtomPair, res_to: ResseqAtomPair, angle: float, alt_conf_id: str|None = None):
    """
    To rotate the main chain by fixing two atoms and rotate the atoms between them
    Note that this function assumes that there is no hydrogen in the structure model
    and a side chain or a carboxyl oxygen should not be moved if linked to a fixed atom

    @param chain: the chain of residues to modify
    @param res_from: tuple (resseq, atom_name) of the first fixed atom
    @param res_to: tuple (resseq, atom_name) of the second fixed atom
    """
    # get the list of all residue involved
    # the number is inclusive, i.e. if `from_res, to_res = 67, 69`, res 67, res 68 and res 69 are selected
    resseq_list: list[int] = [full_id[1] for full_id in chain.child_dict.keys()]
    index_from, index_to = resseq_list.index(res_from[0]), resseq_list.index(res_to[0])
    all_res: list[Residue] = chain.child_list[index_from:index_to+1]
    
    # make a list of all atoms in the specified range of residue, except the first and the last one
    all_rot_atom: list[Atom] = [atom for res in all_res[1:-1] for atom in res]

    # append the atoms from the first residue
    match res_from[1]:
        case 'N':
            for atom in all_res[0]:
                if atom.name != 'N':
                    all_rot_atom.append(atom)
        case 'CA':
            # for atom in all_res[0]:
            #     if atom.name not in ('CA', 'N'):
            #         all_rot_atom.append(atom)
            all_rot_atom.append(all_res[0]['C'])
            all_rot_atom.append(all_res[0]['O'])
        case 'C':
            # all_rot_atom.append(all_res[0]['O'])
            pass
        # case 'O':
        #     pass
        case _:
            assert res_from[1] in ('N', 'CA', 'C'), 'Atom is not in the main chain or is not valid'
    fixed_from = all_res[0][res_from[1]]
    if fixed_from.is_disordered() != 0:
        fixed_from = fixed_from.disordered_get(alt_conf_id)

    # append the atoms from the last residue
    match res_to[1]:
        case 'N':
            pass
        case 'CA':
            all_rot_atom.append(all_res[-1]['N'])
        case 'C':
            for atom in all_res[0]:
                if atom.name not in ('O', 'C'):
                    all_rot_atom.append(atom)
        # case 'O':
        #     for atom in all_res[0]:
        #         if atom.name != 'O':
        #             all_rot_atom.append(atom)
        case _:
            assert res_from[1] in ('N', 'CA', 'C'), 'Atom is not in the main chain or is not valid'
    fixed_to = all_res[-1][res_to[1]]
    if fixed_to.is_disordered() != 0:
        fixed_to = fixed_to.disordered_get(alt_conf_id)

    rotation_axis = fixed_to.get_vector() - fixed_from.get_vector()
    rotation_matrix = rotaxis(angle *pi/1800, rotation_axis)
    def rotate(atom: Atom):
        old_vec = atom.get_vector() - fixed_from.get_vector()
        move_vec = old_vec.left_multiply(rotation_matrix)
        new_atom_pos = move_vec + fixed_from.get_vector()
        atom.set_coord(new_atom_pos)

    for rot_atom in all_rot_atom:
        if rot_atom.is_disordered() == 0:
            rotate(rot_atom)
        else:
            rotate(rot_atom.disordered_get(alt_conf_id)) # type: ignore
