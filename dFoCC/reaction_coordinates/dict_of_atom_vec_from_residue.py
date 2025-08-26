# The `dict_of_atom_vec_from_residue` function below is introduced by LOE
# 
# It reduces the need to `get_vector` every atom in the residue

from Bio.PDB.Residue import Residue
from Bio.PDB import Vector

def dict_of_atom_vec_from_residue(res: Residue, alt_conf_id=None) -> dict[str, Vector]:
    """
    Return a dict of vectors of the position of the atoms
    We can make use of this dict instead of storing all the vectors in variables one by one
    This function is scripted by LOE

    .. code-block:: python
    # before
    CA = residue['CA']
    CB = residue['CB']
    CAv = CA.get_vector()
    CBv = CB.get_vector()

    # after
    vector_dict = dict_of_atom_vec_from_residue(residue)
    CAv = vector_dict['CA']
    CBv = vector_dict['CB']
    """
    residue = res.get_atoms()
    residue_dict = {}
    for x in residue:
        if x != '':
            if x.is_disordered() != 0:
                x = x.disordered_get(alt_conf_id)
            residue_dict[x.get_name()] = x.get_vector()
    return residue_dict
