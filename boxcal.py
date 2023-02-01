from rdkit import Chem
import numpy as np
import os
from pathlib import Path

def boxcal(mymol):

    '''For a molecule, return its name, max extent, min extent
        single centre of gravity, weighted centre of gravity, and box size'''
    icn = 0
    tlist, x_coords, y_coords, z_coords = [], [], [], []
    # sum_x_coords, sum_y_coords, sum_z_coords, sum_denominator = 0, 0, 0, 0  # accumulators
    # pd = Chem.GetPeriodicTable()  # get the periodic table for atomic mass

    for atom in mymol.GetAtoms():
        # get the 3d coordinates of each atom
        pos = mymol.GetConformer().GetAtomPosition(icn)
        # Used to calculate weighted centre
        # sum_x_coords += pos[0] * pd.GetAtomicWeight(atom.GetAtomicNum())
        # sum_y_coords += pos[1] * pd.GetAtomicWeight(atom.GetAtomicNum())
        # sum_z_coords += pos[2] * pd.GetAtomicWeight(atom.GetAtomicNum())
        # sum_denominator += pd.GetAtomicWeight(atom.GetAtomicNum())
        tlist.append(pos)
        icn += 1
    coords = np.array(tlist)
    x_coords, y_coords, z_coords = coords.T

    sim_center = [round(sum(x_coords)/len(x_coords), 5), # simple centre of gravity
                  round(sum(y_coords)/len(y_coords), 5),
                  round(sum(z_coords)/len(z_coords), 5)]

    # weighted_center = np.array([round((sum_x_coords / sum_denominator), 5), # weighted centre of gravity
    #                             round((sum_y_coords / sum_denominator), 5),
    #                             round((sum_z_coords / sum_denominator), 5)])

    max_ext = np.array([max(x_coords), max(y_coords), max(z_coords)]) # max extent
    min_ext = np.array([min(x_coords), min(y_coords), min(z_coords)]) # min extent
    box_size = (max_ext - min_ext).round(5) # box_size = max - min

    return sim_center, box_size

def round_up_to_even(num):
    return np.ceil(num/2) * 2

def center_from_pdbqt(pdbqt_file):
    with open(pdbqt_file) as tempfile:

        target_ligand = [' '.join(line.split()).split(' ') for line in tempfile]

        target_coords = np.array([[float(y) for y in x[-7:-4:1]] for x in target_ligand])

        x_coords, y_coords, z_coords = target_coords.T

        target_center = [round(sum(x_coords) / len(x_coords), 5),  # simple centre of gravity
                         round(sum(y_coords) / len(y_coords), 5),
                         round(sum(z_coords) / len(z_coords), 5)]
    return target_center

if __name__ == "__main__":

    lib_path = "/root/autodl-tmp/lbvs/all_dude/pur2/docking_test/lib/"
    lig_path = "/root/autodl-tmp/lbvs/all_dude/pur2/docking_test/lig/"

    for files in os.listdir(lig_path):
        if files.endswith(".pdbqt"):
            pdbqt_path = os.path.join(lig_path, files)
            ref_path = pdbqt_path.replace(".pdbqt", ".sdf")
            print(ref_path)
            os.system(f"obabel -ipdbqt {pdbqt_path} -osdf -O {ref_path}")
        # elif files.endswith(".sdf"):
        #     ref_path = os.path.join(lig_path, files)

    ref_mol = Chem.MolFromMolFile(ref_path)
    center = center_from_pdbqt(pdbqt_path)

    mols = [ref_mol]
    for files in os.listdir(lib_path):
        if files.endswith(".sdf"):
            sdf_path = os.path.join(lib_path, files)
            mols.append(Chem.MolFromMolFile(sdf_path))

    filtered_mols = [x for x in mols if x is not None]
    size_lis = np.array([boxcal(i)[1] for i in filtered_mols])
    x, y, z = size_lis.T
    largest_npts = [int(round_up_to_even(max(x)/0.375)),
                    int(round_up_to_even(max(y)/0.375)),
                    int(round_up_to_even(max(z)/0.375))]

    with open(f"{lig_path}/{Path(pdbqt_path).stem}.grid.txt", "w") as gridfile:
        gridfile.write(f"""\
{Path(pdbqt_path).stem}
spacing    0.375
npts       {largest_npts[0]} {largest_npts[1]} {largest_npts[-1]}
center    {center[0]:.3f} {center[1]:.3f} {center[2]:.3f}\
""")