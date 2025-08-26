import typing

from numpy import float64, pi
from numpy import typing as npt
from Bio.PDB import Select, vectors, rotaxis
from Bio.PDB.vectors import calc_angle, calc_dihedral

from .dict_of_atom_vec_from_residue import dict_of_atom_vec_from_residue

from Bio.PDB.Structure import Structure
from Bio.PDB.Residue import Residue
from Bio.PDB.Atom import Atom, DisorderedAtom
from Bio.PDB import Vector

def rotate_atom(
    rot_matrix: npt.NDArray[float64],
    rot_atom: Atom,
    ref_atom_pos: Vector,
):
    old_vec = rot_atom.get_vector() - ref_atom_pos
    move_vec = old_vec.left_multiply(rot_matrix)
    new_atom_pos = move_vec + ref_atom_pos
    rot_atom.set_coord(new_atom_pos)

def get_alt_from_atom(atom: Atom, alt_conf_id: str | None = None):
    if atom.is_disordered() == 0:
        return atom
    else:
        return typing.cast(DisorderedAtom, atom).disordered_get(alt_conf_id)

def rotate_chi_angle_by_atom(
    angle: int,
    axis_atoms: tuple[Atom, Atom],
    rot_atoms: list[Atom],
):
    rotation_axis = axis_atoms[1].get_vector() - axis_atoms[0].get_vector()
    rotation_matrix = rotaxis(angle * pi / 1800, rotation_axis)

    for rot_atom in rot_atoms:
        rotate_atom(rotation_matrix, rot_atom, axis_atoms[1].get_vector())

def bend_bond_angle_by_atom(
    angle: int,
    ref_atoms: tuple[Atom, Atom],
    rot_atoms: list[Atom],
):
    vecs_ref = [atom.get_vector() for atom in ref_atoms]
    rotation_axis = (vecs_ref[1] - vecs_ref[0]) ** (vecs_ref[2] - vecs_ref[1])
    rotation_matrix = rotaxis(angle * pi / 1800, rotation_axis)

    for rot_atom in rot_atoms:
        rotate_atom(rotation_matrix, rot_atom, vecs_ref[1])

def translate_atoms_by_atom(
    x: int, y: int, z: int,
    tr_atoms: list[Atom],
):
    translation = Vector(x / 1000, y / 1000, z / 1000)

    for tr_atom in tr_atoms:
        tr_atom.set_coord(tr_atom.get_vector() + translation)

def rotate_chi_angle(
        residue: Residue,
        angle: int,
        axis_atoms: tuple[str, str],
        rot_atoms: list[str],
        alt_conf_id=None
):
    def get_alt(atom_name: str) -> Atom:
        if residue[atom_name].is_disordered() == 0:
            return residue[atom_name]
        else:
            return residue[atom_name].disordered_get(alt_conf_id)

    vec_axis_0, vec_axis_1 = get_alt(axis_atoms[0]).get_vector(), get_alt(axis_atoms[1]).get_vector()
    rotation_axis = vec_axis_1 - vec_axis_0
    rotation_matrix = rotaxis(angle * pi / 1800, rotation_axis)

    def rotate(atom: Atom):
        old_vec = atom.get_vector() - vec_axis_1
        move_vec = old_vec.left_multiply(rotation_matrix)
        new_atom_pos = move_vec + vec_axis_1
        atom.set_coord(new_atom_pos)

    for rot_atom_name in rot_atoms:
        rotate(get_alt(rot_atom_name))

def bend_bond_angle(
        residue: Residue,
        angle: int,
        ref_atoms: tuple[str, str, str],
        rot_atoms: list[str],
        alt_conf_id=None,
):
    def get_alt(atom_name: str) -> Atom:
        if residue[atom_name].is_disordered() == 0:
            return residue[atom_name]
        else:
            return residue[atom_name].disordered_get(alt_conf_id)

    vecs_ref = [get_alt(atom_name).get_vector() for atom_name in ref_atoms]
    rotation_axis = (vecs_ref[1] - vecs_ref[0]) ** (vecs_ref[2] - vecs_ref[1])
    rotation_matrix = rotaxis(angle * pi / 1800, rotation_axis)

    def rotate(atom: Atom):
        old_vec = atom.get_vector() - vecs_ref[1]
        move_vec = old_vec.left_multiply(rotation_matrix)
        new_atom_pos = move_vec + vecs_ref[1]
        atom.set_coord(new_atom_pos)

    for rot_atom_name in rot_atoms:
        rotate(get_alt(rot_atom_name))

def translate_atoms(
    residue: Residue,
    x: int, y: int, z: int,
    tr_atoms: list[str],
    alt_conf_id=None,
):
    def get_alt(atom_name: str) -> Atom:
        if residue[atom_name].is_disordered() == 0:
            return residue[atom_name]
        else:
            return residue[atom_name].disordered_get(alt_conf_id)

    translation = Vector(x / 1000, y / 1000, z / 1000)
    
    for tr_atom_name in tr_atoms:
        tr_atom = get_alt(tr_atom_name)
        tr_atom.set_coord(tr_atom.get_vector() + translation)

def translate_atoms_by_vector(
    magnitude: int,
    vector: Vector,
    tr_atoms: list[Atom],
):
    translation = vector.normalized() ** magnitude # scalar multiplication
    for tr_atom in tr_atoms:
        tr_atom.set_coord(tr_atom.get_vector() + translation)

def rotate_atoms_by_vector(
    angle: int,
    rot_vector: Vector,
    fixed_position: Vector,
    rot_atoms: list[Atom],
):
    rotation_matrix = rotaxis(angle * pi / 1800, rot_vector)

    def rotate(atom: Atom):
        old_vec = atom.get_vector() - fixed_position
        move_vec = old_vec.left_multiply(rotation_matrix)
        new_atom_pos = move_vec + fixed_position
        atom.set_coord(new_atom_pos)

    for rot_atom in rot_atoms:
        rotate(rot_atom)
