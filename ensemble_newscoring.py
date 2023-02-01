#!/root/miniconda3/envs/ncflbvs/bin python
from gridmap import gridmap_information
from pathlib import Path
from create_mapfld import *
from ensemble_mapping import *
from boxcal import boxcal, center_from_pdbqt, round_up_to_even
from create_folder import createFolder
from rdkit import Chem
from modify_pdbqt import modify_pdbqt
import numpy as np
from subsample import create_subsample
from roc_auc import roc_auc_ef
import os, subprocess, time, argparse
import functools
from extract_best_pose import extract_best_pos
from mpire import WorkerPool
#======================================================================================================================
# parser = argparse.ArgumentParser(description='dude target')
# parser.add_argument('target')
# args = parser.parse_args()
# dude_target = args.target

targets = set(os.listdir("/root/autodl-tmp/lbvs/test_ensemble_full/"))
# targets = ['ampc']
for dude_target in targets:
#======================================================================================================================
#======================================================================================================================
# define all the paths needed
    os.chdir("/root/autodl-tmp/lbvs/new_scoring/")
    lib_path = f'/root/autodl-tmp/lbvs/test_ensemble_full/{dude_target}/docking_test/lib/'
    lig_path = f'/root/autodl-tmp/lbvs/test_ensemble_full/{dude_target}/docking_test/lig/'
    dlg_original_path = "/root/autodl-tmp/lbvs/new_scoring/"
    for files in os.listdir(lig_path):
        if files.endswith(".pdbqt"):
            pdbqt_path = os.path.join(lig_path, files)
            grid_file_path = pdbqt_path.replace(".pdbqt", ".grid.txt")
    modify_pdbqt(pdbqt_path) # make the reference pdbqt file rigid
    #======================================================================================================================
    # actives_to_be_deleted = create_subsample(f"{lib_path}/actives_final.sdf.gz")[0]
    # decoys_to_be_deleted = create_subsample(f"{lib_path}/decoys_final.sdf.gz")[0]
    #======================================================================================================================
    subprocess.check_call(f"obabel -isdf {lib_path}actives_final.sdf.gz -osdf -O *.sdf --split", shell=True, cwd=lib_path)
    time.sleep(5)
    subprocess.check_call(f"obabel -isdf {lib_path}decoys_final.sdf.gz -osdf -O *.sdf --split", shell=True, cwd=lib_path)
    time.sleep(5)
    subprocess.check_call(f"rm {lib_path}actives_final.sdf.gz", shell=True)
    subprocess.check_call(f"rm {lib_path}actives_final.sdf", shell=True)
    subprocess.check_call(f"rm {lib_path}decoys_final.sdf.gz", shell=True)
    time.sleep(15)
    #======================================================================================================================
    # def multiprocess_subsample(files):
    #     subprocess.check_call(f"rm {lib_path}{files}.sdf", shell=True)
    # with WorkerPool(n_jobs=96) as pool:
    #     pool.map(multiprocess_subsample, [*actives_to_be_deleted, *decoys_to_be_deleted], progress_bar=True)
    #======================================================================================================================
    for files in os.listdir(lig_path):
        if files.endswith(".pdbqt"):
            pdbqt_path = os.path.join(lig_path, files)
            ref_path = pdbqt_path.replace(".pdbqt", ".sdf")
            subprocess.check_call(f"obabel -ipdbqt {pdbqt_path} -osdf -O {ref_path}", shell=True)
            

    ref_mol = Chem.MolFromMolFile(ref_path)
    center = center_from_pdbqt(pdbqt_path)

    mols = [ref_mol]
    for files in os.listdir(lib_path):
        if files.endswith(".sdf"):
            sdf_path = os.path.join(lib_path, files)
            mols.append(Chem.MolFromMolFile(sdf_path))

    filtered_mols = [x for x in mols if x is not None]
    size_lis = np.array([boxcal(i)[1] for i in filtered_mols])
    ref_mol_size = size_lis[0]

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
    print("Finished creating 'grid.txt'")
    #======================================================================================================================
    #======================================================================================================================
    supported_atom_types = ["HD", "C", "A", "N", "NA", "OA", "F", "P", "SA", "S",
                            "Cl", "Br", "I", "Mg", "Ca", "Mn", "Fe", "Zn", "d"]

    pur2_test = gridmap_information.from_gridmap(grid_file_path, pdbqt_path)
    atomtype_not_in_rec = list(set(supported_atom_types)-set(pur2_test.target_atomtype))

    elec_map(folder=lig_path, spacing=pur2_test.spacing, npts=pur2_test.npts,
                gridcenter=pur2_test.gridcenter, gridmap_coords=pur2_test.target_gridmap_coords,
                ligand_coords=pur2_test.target_coords, atomic_partial_charge=pur2_test.target_charge,
                filename=Path(pdbqt_path).stem)

    affinity_mapping_partial = functools.partial(affinity_mapping, pur2_test.target_ligand, lig_path, pur2_test.spacing,
                                                    pur2_test.npts, pur2_test.gridcenter, pur2_test.target_gridmap_coords,
                                                    pur2_test.target_coords, Path(pdbqt_path).stem)
    general_map_partial = functools.partial(general_map, pur2_test.target_ligand, lig_path, pur2_test.spacing,
                                            pur2_test.npts, pur2_test.gridcenter, pur2_test.target_gridmap_coords,
                                            pur2_test.target_coords, Path(pdbqt_path).stem)



    with WorkerPool(n_jobs=4) as pool:
        pool.map(affinity_mapping_partial, pur2_test.target_atomtype)
        pool.map(general_map_partial, atomtype_not_in_rec)

    create_mapfld(folder=lig_path, filename=Path(pdbqt_path).stem, spacing=pur2_test.spacing,
                    npts=pur2_test.npts, center=pur2_test.gridcenter)
    print("Finished creating map files")
    #======================================================================================================================
    #======================================================================================================================
    createFolder(libpath=lib_path)
    for files in os.listdir(lig_path):
        if files.endswith(".maps.fld"):
            mapfld_path = os.path.join(lig_path, files)
    with open(f"{lig_path}/docking_file.lst", "w") as lstfile:
        lstfile.write(f"{mapfld_path}\n")

    def multiprocess_libpath(folder):
        abs_path = os.path.join(lib_path, folder)
        if os.path.isdir(abs_path):
            for filename in os.listdir(abs_path):
                file_path = os.path.join(abs_path, filename)

                if os.path.isfile(file_path) and ".sdf" in filename:
                    pdb_filepath = file_path.replace(".sdf", ".pdb")
                    pdbqt_filename = filename.replace(".sdf", ".pdbqt")
                    pdbqt_filepath = file_path.replace(".sdf", ".pdbqt")

                    # run openbabel to convert sdf file into pdb file
                    os.system("obabel %s -isdf -opdb -O %s" % (file_path, pdb_filepath))
                    # run script convert pdb to pdbqt
                    os.system("python ligandprep.py -l %s -o %s" % (pdb_filepath, pdbqt_filepath))
                    with open(f"{lig_path}/docking_file.lst", "a") as tempfile:
                        tempfile.write(f"{pdbqt_filepath}\n{Path(pdbqt_filepath).stem}\n")
                        # tempfile.write(f"{Path(pdbqt_filepath).stem}\n")

    with WorkerPool(n_jobs=192) as pool:
        pool.map(multiprocess_libpath, os.listdir(lib_path))

    env = os.environ.copy().update(OMP_NUM_THREADS=384)
    subprocess.run(f"/root/AutoDock-GPU/bin/autodock_gpu_128wi -B {lig_path}/docking_file.lst -D all --nrun 100 --xmloutput 0", shell=True, env=env)

    def multiprocess_move_dlg(files):
        if files.endswith(".dlg"):
            dlg_path = os.path.join(dlg_original_path, files)
            subprocess.check_call(f"mv {dlg_path} {lib_path}{Path(dlg_path).stem}/{Path(dlg_path).stem}.dlg", shell=True)
    with WorkerPool(n_jobs=192) as pool:
        pool.map(multiprocess_move_dlg, os.listdir(dlg_original_path), progress_bar=True)
    # for files in os.listdir(dlg_original_path):
    #     if files.endswith(".dlg"):
    #         dlg_path = os.path.join(dlg_original_path, files)
    #         subprocess.check_call(f"mv {dlg_path} {lib_path}{Path(dlg_path).stem}/{Path(dlg_path).stem}.dlg", shell=True)
    #======================================================================================================
    #======================================================================================================
    name_lis, score_lis, smile_lis = [], [], []
    actives_name = []

    for folder in os.listdir(lib_path):
        abs_path = os.path.join(lib_path, folder)
        if os.path.isdir(abs_path):
            for filename in os.listdir(abs_path):
                if filename.endswith("dlg"):
                    file_path = os.path.join(abs_path, filename)
                    name = Path(file_path).stem
                    name_lis.append(name)
                    with open(file_path, "r") as dlgfile:
                        lines = dlgfile.readlines()
                        index = lines.index("    CLUSTERING HISTOGRAM\n")
                        score = float(lines[index+10].split("|")[1].strip())
                        score_lis.append(score)

    for name in name_lis:
        if name.startswith("CHEMBL"):
            actives_name.append(name)

    # ref_ind = name_lis.index(f"{Path(ref_path).stem}")
    # ref_score = score_lis[ref_ind]
    # del name_lis[ref_ind]
    # del score_lis[ref_ind]

    # normalized_score_lis = -1/(1-np.array(score_lis)/ref_score)
    sorted_index = np.argsort(score_lis)
    sorted_name = [name_lis[i] for i in sorted_index]
    sorted_score = [score_lis[i] for i in sorted_index]

    with open(f"{lib_path.replace('lib/', '')}ranked.txt", "w") as tempfile:
        [tempfile.write(f"{name}\n") for name in sorted_name]
        tempfile.truncate(tempfile.tell()-1)

    with open(f"{lib_path.replace('lib/', '')}actives.txt", "w") as tempfile:
        [tempfile.write(f"{name}\n") for name in set(actives_name)]
        tempfile.truncate(tempfile.tell()-1)

    with open(f"{lib_path.replace('lib/', '')}scores.txt", "w") as tempfile:
        [tempfile.write(f"{score}\n") for score in sorted_score]
        tempfile.truncate(tempfile.tell()-1)

    oldwd = os.getcwd()
    os.chdir(f"{lib_path.replace('lib/', '')}")

    try:
        roc_auc_ef(ranked_file="ranked.txt", 
                    actives_file="actives.txt",
                    candidate_name=dude_target)
    except:
        pass

    os.chdir(oldwd)

    def multiprocess_extraction(subfolder):
        candi_path = os.path.join(lib_path, subfolder)
        if os.path.isdir(candi_path):
            for file in os.listdir(candi_path):
                if file.endswith(".dlg"):
                    dlg_file_path = os.path.join(candi_path, file)
                    # print(dlg_file_path)
                    with open(f"{candi_path}/{Path(dlg_file_path).stem}_best_position.pdbqt", "w") as tmpfile:
                        for what in extract_best_pos(dlg_file_path):
                            tmpfile.write(what.replace("DOCKED: ", ""))
    with WorkerPool(n_jobs=96) as pool:
        pool.map(multiprocess_extraction, os.listdir(lib_path), progress_bar=True)

    subprocess.check_call(f"find . -name '*.dlg' -type f -delete", shell=True, cwd=lib_path)
    subprocess.check_call(f"find . -name '*.pdb' -type f -delete", shell=True, cwd=lib_path)
    subprocess.check_call(f"find . -name '*.sdf' -type f -delete", shell=True, cwd=lib_path)
