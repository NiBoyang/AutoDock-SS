import os, subprocess
import argparse
import functools
from pathlib import Path
from rdkit import Chem
import numpy as np
from importlib.resources import path
import os
from mpire import WorkerPool
from grid_map_utils import *
from geometry_utils import *
from file_utils import *
from roc_auc import *
import glob
import psutil

# Common Variables
# Change them to yours
autodock_gpu_exec_path = "/root/AutoDock-GPU/bin/autodock_gpu_128wi" # directory of Autodock-GPU exec file
lib_path = f'/root/autodl-tmp/cache4/enamine_2/lib/' # directory of your VS library
lig_path = f'/root/autodl-tmp/cache4/enamine_2/lig/' # directory of your reference ligand
path_of_scripts = "/root/autodl-tmp/adss-scripts/" # directory of Autodock-SS scripts
os.chdir(path_of_scripts)
dlg_original_path = lig_path # directory for storing dlg files
sdfgz_files = glob.glob(f"{lib_path}*.sdf.gz") # the VS library file should have the extension of .sdf.gz
supported_atom_types = ["HD", "C", "A", "N", "NA", "OA", "F", "P", "SA", "S",
                        "Cl", "Br", "I", "Mg", "Ca", "Mn", "Fe", "Zn", "d"]
env = os.environ.copy().update(OMP_NUM_THREADS=psutil.cpu_count(logical=True))
n_jobs = psutil.cpu_count(logical=False)

# Global variables
def get_files_with_extension(directory, extension):
    return [os.path.join(directory, file) for file in os.listdir(directory) if file.endswith(extension)]

pdbqt_files = get_files_with_extension(lig_path, ".pdbqt") # the reference ligand file should have .pdbqt extension
if not pdbqt_files:
    raise ValueError("No .pdbqt files found in the specified ligand path.")
pdbqt_path = pdbqt_files[0]
grid_file_path = pdbqt_path.replace(".pdbqt", ".grid.txt")
ref_path = pdbqt_path.replace(".pdbqt", ".sdf")

# Functions
def calculate_grid_dimensions(lig_path, lib_path):

    """
    This function calculates the grid dimensions for the docking process. It uses the reference ligand and the library of molecules
    to determine the size of the grid box that will be used in the docking simulations. The grid dimensions are calculated based on 
    the size of the molecules and their spatial arrangement. The function uses subprocess calls to run external programs for file 
    conversion and manipulation. The results are written to a grid.txt file which is used in the subsequent docking process.
    """

    subprocess.check_call(f"obabel -ipdbqt {pdbqt_path} -osdf -O {ref_path}", shell=True)
    subprocess.check_call(f"cp {ref_path} {lib_path}", shell=True)
    
    ref_mol = Chem.MolFromMolFile(ref_path)
    center = center_from_pdbqt(pdbqt_path)

    mols = [ref_mol] + [Chem.MolFromMolFile(sdf_path) for sdf_path in get_files_with_extension(lib_path, ".sdf")]
    filtered_mols = [mol for mol in mols if mol is not None]

    size_lis = np.array([boxcal(mol)[1] for mol in filtered_mols])
    x, y, z = size_lis.T

    largest_npts = [int(round_up_to_even(max(dim)/0.375)) for dim in [x, y, z]]

    with open(f"{lig_path}/{Path(pdbqt_path).stem}.grid.txt", "w") as gridfile:
        gridfile.write(f"""\
{Path(pdbqt_path).stem}
spacing    0.375
npts       {largest_npts[0]} {largest_npts[1]} {largest_npts[-1]}
center    {center[0]:.3f} {center[1]:.3f} {center[2]:.3f}\
""")
    print("Finished creating 'grid.txt'")

def create_map_files(lig_path, grid_file_path, pdbqt_path, n_jobs):
    
    """
    This function creates map files for the docking process. It uses the reference ligand and the grid dimensions calculated 
    in the previous step to generate map files for each atom type. The function uses multiprocessing to speed up the map file 
    generation process.
    """
    
    supported_atom_types = ["HD", "C", "A", "N", "NA", "OA", "F", "P", "SA", "S",
                            "Cl", "Br", "I", "Mg", "Ca", "Mn", "Fe", "Zn", "d"]

    query_ligmap = GridMapInformation.from_gridmap(grid_file_path, pdbqt_path)
    atomtype_not_in_rec = list(set(supported_atom_types) - set(query_ligmap.target_atomtype))

    elec_map(folder=lig_path, 
             spacing=query_ligmap.spacing, 
             npts=query_ligmap.npts,
             gridcenter=query_ligmap.grid_center, 
             gridmap_coords=query_ligmap.target_gridmap_coords,
             ligand_coords=query_ligmap.target_coords, 
             atomic_partial_charge=query_ligmap.target_charge,
             filename=Path(pdbqt_path).stem)

    mapping_args = (query_ligmap.target_ligand, 
                    lig_path, 
                    query_ligmap.spacing,
                    query_ligmap.npts, 
                    query_ligmap.grid_center, 
                    query_ligmap.target_gridmap_coords,
                    query_ligmap.target_coords, 
                    Path(pdbqt_path).stem)

    with WorkerPool(n_jobs=n_jobs) as pool:
        pool.map(functools.partial(affinity_map, *mapping_args), query_ligmap.target_atomtype)
        pool.map(functools.partial(general_map, *mapping_args), atomtype_not_in_rec)

    create_mapfld(folder=lig_path, 
                  filename=Path(pdbqt_path).stem, 
                  spacing=query_ligmap.spacing,
                  npts=query_ligmap.npts, 
                  center=query_ligmap.grid_center)
    
    print("Finished creating map files")


def create_folder_and_initialize_docking_file(lib_path, lig_path):
    """Creates a new folder and initializes a docking file with path to .maps.fld files."""
    create_folder(libpath=lib_path)
    for files in os.listdir(lig_path):
        if files.endswith(".maps.fld"):
            mapfld_path = os.path.join(lig_path, files)

    with open(f"{lig_path}/docking_file.lst", "w") as lstfile:
        lstfile.write(f"{mapfld_path}\n")

def run_ligand_prep_commands(file_path, pdb_filepath, pdbqt_filepath):
    """Runs openbabel and another script to prepare the file."""
    os.system(f"obabel {file_path} -isdf -opdb -O {pdb_filepath}")
    os.system(f"python ligandprep.py -l {pdb_filepath} -o {pdbqt_filepath}")

def write_to_docking_file(lig_path, pdbqt_filepath):
    """Writes to the docking file."""
    with open(f"{lig_path}/docking_file.lst", "a") as tempfile:
        tempfile.write(f"{pdbqt_filepath}\n{Path(pdbqt_filepath).stem}\n")

def process_folder_and_files(folder):
    """Processes all the .sdf files in a given folder, converting and preparing them."""
    abs_path = os.path.join(lib_path, folder)
    if os.path.isdir(abs_path):
        for filename in os.listdir(abs_path):
            file_path = os.path.join(abs_path, filename)

            if os.path.isfile(file_path) and ".sdf" in filename:
                pdb_filepath = file_path.replace(".sdf", ".pdb")
                pdbqt_filepath = file_path.replace(".sdf", ".pdbqt")
                run_ligand_prep_commands(file_path, pdb_filepath, pdbqt_filepath)
                write_to_docking_file(lig_path, pdbqt_filepath)

def run_autodock(autodock_gpu_exec_path, lig_path, env):
    """Runs AutoDock-GPU with the prepared input."""
    env = os.environ.copy().update(OMP_NUM_THREADS=192)
    subprocess.run(f"{autodock_gpu_exec_path} -B {lig_path}/docking_file.lst -D all --nrun 100 --xmloutput 0", shell=True, env=env)

def move_dlg(files):
    if files.endswith(".dlg"):
        dlg_path = os.path.join(dlg_original_path, files)
        subprocess.check_call(f"mv {dlg_path} {lib_path}{Path(dlg_path).stem}/{Path(dlg_path).stem}.dlg", shell=True)
        

def extract_scores_from_file(file_path):
    name = Path(file_path).stem
    with open(file_path, "r") as dlgfile:
        lines = dlgfile.readlines()
        index = lines.index("    CLUSTERING HISTOGRAM\n")
        score = float(lines[index+10].split("|")[1].strip())
        return name, score

def read_dlg_file_and_extract_scores(lib_path):
    name_lis, score_lis = [], []
    paths = []

    for folder in os.listdir(lib_path):
        abs_path = os.path.join(lib_path, folder)
        if os.path.isdir(abs_path):
            for filename in os.listdir(abs_path):
                if filename.endswith("dlg"):
                    file_path = os.path.join(abs_path, filename)
                    paths.append(file_path)

    with WorkerPool(n_jobs=n_jobs) as pool:
        for name, score in pool.map(extract_scores_from_file, paths, progress_bar=True):
            name_lis.append(name)
            score_lis.append(score)
            
    return name_lis, score_lis

def extract_actives_names(name_lis):
    actives_name = []
    for name in name_lis:
        if name.startswith("CHEMBL"):
            actives_name.append(name)
    return actives_name

def normalize_and_sort_scores(name_lis, score_lis, ref_path):
    """Normalizes and sorts scores."""
    ref_ind = name_lis.index(f"{Path(ref_path).stem}")
    ref_score = score_lis[ref_ind]
    del name_lis[ref_ind]
    del score_lis[ref_ind]

    normalized_score_lis = np.array(score_lis)/ref_score
    sorted_index = np.argsort(normalized_score_lis)[::-1]
    sorted_name = [name_lis[i] for i in sorted_index]
    sorted_score = [normalized_score_lis[i] for i in sorted_index]

    return sorted_name, sorted_score

def write_to_files(lib_path, sorted_name, actives_name, sorted_score):
    """Writes sorted names, active names and sorted scores to their respective files."""
    base_path = lib_path.replace('lib/', '')
    
    with open(f"{base_path}ranked.txt", "w") as tempfile:
        tempfile.write('\n'.join(map(str, sorted_name)))
        
    if actives_name:    
        with open(f"{base_path}actives.txt", "w") as tempfile:
            tempfile.write('\n'.join(map(str, actives_name)))
        
    with open(f"{base_path}scores.txt", "w") as tempfile:
        tempfile.write('\n'.join(map(str, sorted_score)))


def compute_roc_auc_ef(base_path, dude_target):
    """Computes ROC AUC EF."""
    oldwd = os.getcwd()
    os.chdir(base_path)
    roc_auc_ef(ranked_file="ranked.txt",
               actives_file="actives.txt",
               candidate_name=dude_target) # change reference ligand name to yours
    os.chdir(oldwd)

def multiprocess_extraction(subfolder):
    """Extracts best positions using multiple processes."""
    candi_path = os.path.join(lib_path, subfolder)
    if os.path.isdir(candi_path):
        for file in os.listdir(candi_path):
            if file.endswith(".dlg"):
                dlg_file_path = os.path.join(candi_path, file)
                with open(f"{candi_path}/{Path(dlg_file_path).stem}_best_position.pdbqt", "w") as tmpfile:
                    for what in extract_best_pos(dlg_file_path):
                        tmpfile.write(what.replace("DOCKED: ", ""))

def remove_files(lib_path, extensions=['dlg', 'pdb', 'sdf']):
    """Removes files with specific extensions."""
    for ext in extensions:
        subprocess.run(["find", ".", "-name", f'*.{ext}', "-type", "f", "-delete"], cwd=lib_path)

def main():
    
    print(f"Confirmation: The PDBQT file you input is : {pdbqt_path}")
    
    modify_pdbqt(pdbqt_path) # make the reference pdbqt file rigid

    for sdfgz_file in sdfgz_files:
    # Convert each .sdf.gz file to .sdf
        subprocess.check_call(f"obabel -isdf {sdfgz_file} -osdf -O *.sdf --split", shell=True, cwd=lib_path)
        # Try to remove the .sdf.gz file as it was split
        try:
            os.remove(sdfgz_file)
            os.remove(sdfgz_file.replace(".gz", ""))
        except OSError as e:
            print(f"Error: {sdfgz_file} : {e.strerror}")
    
    calculate_grid_dimensions(lig_path, lib_path)

    create_map_files(lig_path, grid_file_path, pdbqt_path, n_jobs)

    create_folder_and_initialize_docking_file(lib_path, lig_path)

    with WorkerPool(n_jobs=n_jobs) as pool:
        pool.map(process_folder_and_files, os.listdir(lib_path))    

    run_autodock(autodock_gpu_exec_path, lig_path, env)

    with WorkerPool (n_jobs=n_jobs) as pool:
        pool.map(move_dlg, os.listdir(dlg_original_path), progress_bar=True)
    
    name_lis, score_lis = read_dlg_file_and_extract_scores(lib_path)
    
    actives_name = extract_actives_names(name_lis)

    sorted_name, sorted_score = normalize_and_sort_scores(name_lis, score_lis, ref_path)

    write_to_files(lib_path, sorted_name, actives_name, sorted_score)

    # If you are NOT doing benchmarking, please mute the next line
    # compute_roc_auc_ef(base_path=f"{lib_path.replace('lib/', '')}", dude_target=lib_path)

    with WorkerPool(n_jobs=n_jobs) as pool:
        pool.map(multiprocess_extraction, os.listdir(lib_path), progress_bar=True)

    remove_files(lib_path, extensions=['dlg', 'pdb', 'sdf'])


if __name__ == "__main__":
    main()
