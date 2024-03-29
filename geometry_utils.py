import numpy as np
from rdkit import Chem

def boxcal(mymol):
    icn = 0
    tlist, x_coords, y_coords, z_coords = [], [], [], []
    # sum_x_coords, sum_y_coords, sum_z_coords, sum_denominator = 0, 0, 0, 0  # accumulators
    # pd = Chem.GetPeriodicTable()  # get the periodic table for atomic mass

    for atom in mymol.GetAtoms():
        # get the 3d coordinates of each atom
        pos = mymol.GetConformer().GetAtomPosition(icn)
        tlist.append(pos)
        icn += 1
    coords = np.array(tlist)
    x_coords, y_coords, z_coords = coords.T

    sim_center = np.array([np.mean(x_coords), np.mean(y_coords), np.mean(z_coords)])

    max_ext = np.array([np.max(x_coords), np.max(y_coords), np.max(z_coords)])
    min_ext = np.array([np.min(x_coords), np.min(y_coords), np.min(z_coords)])
    box_size = max_ext - min_ext

    return sim_center, box_size

def round_up_to_even(num):
    return np.ceil(num / 2) * 2

def center_from_pdbqt(pdbqt_file):
    with open(pdbqt_file) as tempfile:
        target_ligand = [' '.join(line.split()).split(' ') for line in tempfile]

    target_coords = np.array([[float(y) for y in x[-7:-4:1]] for x in target_ligand])

    x_coords, y_coords, z_coords = target_coords.T

    target_center = np.array([np.mean(x_coords), np.mean(y_coords), np.mean(z_coords)])

    return target_center