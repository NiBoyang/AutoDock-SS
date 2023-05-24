import os, subprocess, time
import functools
from pathlib import Path
from rdkit import Chem
import numpy as np
from importlib.resources import path
import os
from mpire import WorkerPool

# from gridmap import gridmap_information
class gridmap_information():

    def __init__(self, target_gridmap_coords, target_atomtype, target_charge,
                 target_coords, target_ligand, target_center, target_size,
                 gridcenter, spacing, npts):

        self.target_gridmap_coords = target_gridmap_coords
        self.target_charge = target_charge
        self.target_atomtype = target_atomtype
        self.target_coords = target_coords
        self.target_ligand = target_ligand
        self.target_center = target_center
        self.target_size = target_size
        self.gridcenter = gridcenter
        self.spacing = spacing
        self.npts = npts

    @classmethod
    def from_gridmap(cls, gridfile, pdbqt):

        with open(gridfile) as tempfile:

            for line in tempfile:
                if 'center' in line:
                    gridcenter = [float(x) for x in ' '.join(line.split()).split(' ')[1:4]]
                elif 'spacing' in line:
                    spacing = float(' '.join(line.split()).split(' ')[1])
                elif 'npts' in line:
                    npts = [int(x) for x in ' '.join(line.split()).split(' ')[1:4]]
        
            xmin = gridcenter[0] - (npts[0] / 2) * spacing
            xcoord = [round(xmin + i * spacing, 4) for i in range(0, npts[0] + 1)]
            ymin = gridcenter[1] - (npts[1] / 2) * spacing
            ycoord = [round(ymin + i * spacing, 4) for i in range(0, npts[1] + 1)]
            zmin = gridcenter[2] - (npts[2] / 2) * spacing
            zcoord = [round(zmin + i * spacing, 4) for i in range(0, npts[2] + 1)]

            target_gridmap_coords = [[xcoord[k], ycoord[j], zcoord[i]]
                    for i in range(0, npts[2] + 1) for j in range(0, npts[1] + 1) for k in range(0, npts[0] + 1)]
        
        with open(pdbqt) as tempfile:

            target_ligand = [' '.join(line.split()).split(' ') for line in tempfile]
            # extract partial charge from pdbqt
            target_charge = [float(x[-2]) for x in target_ligand]
            # extract target ligand's coordinates from pdbqt
            target_coords = np.array([[float(y) for y in x[-7:-4:1]] for x in target_ligand])
            # extract atom types of target ligand from pdbqt
            target_atomtype = list(set([z[-1] for z in target_ligand]))

            x_coords, y_coords, z_coords = target_coords.T

            target_center = [round(sum(x_coords) / len(x_coords), 5),  # simple centre of gravity
                             round(sum(y_coords) / len(y_coords), 5),
                             round(sum(z_coords) / len(z_coords), 5)]

            max_ext = np.array([max(x_coords), max(y_coords), max(z_coords)])  # max extent
            min_ext = np.array([min(x_coords), min(y_coords), min(z_coords)])  # min extent
            target_size = (max_ext - min_ext).round(5)

        return cls(target_gridmap_coords, target_atomtype, target_charge,
                   target_coords, target_ligand, target_center, target_size,
                   gridcenter, spacing, npts)
    
# from create_mapfld import *

def create_mapfld(folder, filename, spacing, npts, center):

    with open(f"{folder}{filename}.maps.fld", "w") as mapfld:

        mapfld.write(f"""
# AVS field file
#
# AutoDock Atomic Affinity and Electrostatic Grids
#
# Created by ./../autogrid4.
#
#SPACING {spacing}
#NELEMENTS {npts[0]} {npts[1]} {npts[-1]}
#CENTER {center[0]:.3f}  {center[1]:.3f} {center[-1]:.3f}
#MACROMOLECULE ZINC000000002212_receptor.pdbqt
#GRID_PARAMETER_FILE D:\AutoDock_Workplace\script test\..\scripttest.gpf
#
ndim=3			# number of dimensions in the field
dim1={npts[0]+1}			# number of x-elements
dim2={npts[1]+1}			# number of y-elements
dim3={npts[-1]+1}			# number of z-elements
nspace=3		# number of physical coordinates per point
veclen=9		# number of affinity values at each point
data=float		# data type (byte, integer, float, double)
field=uniform		# field type (uniform, rectilinear, irregular)
coord 1 file=ligand.maps.xyz filetype=ascii offset=0
coord 2 file=ligand.maps.xyz filetype=ascii offset=2
coord 3 file=ligand.maps.xyz filetype=ascii offset=4
label=A-affinity	# component label for variable 1
label=C-affinity	# component label for variable 2
label=HD-affinity	# component label for variable 3
label=N-affinity	# component label for variable 4
label=NA-affinity	# component label for variable 5
label=OA-affinity	# component label for variable 6
label=F-affinity
label=P-affinity
label=SA-affinity
label=S-affinity
label=Cl-affinity
label=Br-affinity
label=I-affinity
label=Electrostatics
label=Desolvation
#
# location of affinity grid files and how to read them
#
variable 1 file={filename}.A.map filetype=ascii skip=6
variable 2 file={filename}.C.map filetype=ascii skip=6
variable 3 file={filename}.HD.map filetype=ascii skip=6
variable 4 file={filename}.N.map filetype=ascii skip=6
variable 5 file={filename}.NA.map filetype=ascii skip=6
variable 6 file={filename}.OA.map filetype=ascii skip=6
variable 7 file={filename}.F.map filetype=ascii skip=6
variable 8 file={filename}.P.map filetype=ascii skip=6
variable 9 file={filename}.SA.map filetype=ascii skip=6
variable 10 file={filename}.S.map filetype=ascii skip=6
variable 11 file={filename}.Cl.map filetype=ascii skip=6
variable 12 file={filename}.Br.map filetype=ascii skip=6
variable 13 file={filename}.I.map filetype=ascii skip=6
variable 14 file={filename}.e.map filetype=ascii skip=6
variable 15 file={filename}.d.map filetype=ascii skip=6
""")
        
# from mapping import *
def affinity_mapping(ligand, folder, spacing, npts, gridcenter, gridcoord, ligandcoord, filename, at):
    import numpy as np
    # halogens = ["F", "Cl", "Br", "I"]
    # aromatic_atoms = ["A", "NA", "OA", "SA"]
    # if at in halogens:
    #     compensation = set(halogens)-set(at)
    # elif at in aromatic_atoms:
    #     compensation = set(aromatic_atoms)-set(at)
    # else:
    #     compensation = [None]
    # extract data of an atom type
    # extract the coordiante and convert from string to float
    at_pdbqt = []
    # compensation_pdbqt = []

    for i in range(0, len(ligand)):
        if ligand[i][-1] == at:
            at_pdbqt.append(ligandcoord[i])
        # elif ligand[i][-1] in compensation:
        #     compensation_pdbqt.append(ligandcoord[i])

    # generate map file of carbon
    with open(f"{folder}{filename}.{at}.map", "w") as gridmap:

        # write the document information into each map file
        gridmap.write("GRID_PARAMETER_FILE " + folder + "\ligand.gpf\n")
        gridmap.write(f"GRID_DATA_FILE {filename}.maps.fld\n")
        gridmap.write(f"MACROMOLECULE {filename}.pdbqt\n")
        gridmap.write("SPACING " + str(spacing) + "\n")
        gridmap.write("NELEMENTS " + str(npts[0]) + " " + str(npts[1]) + " " + str(npts[2]) + "\n")
        gridmap.write("CENTER " + str(gridcenter[0]) + " " + str(gridcenter[1]) + " " + str(gridcenter[2]) + "\n")

        # if the carbon atom is within 5 angstrom of a grid point
        # energy will be counted into that grid point based on the distance
        for i in range(0, len(gridcoord)):
            energy = 0
            dist_lis = []
            for k in range(0, len(at_pdbqt)):
            #     if (gridcoord[i][0] - 1.54) <= at_pdbqt[k][0] <= (gridcoord[i][0] + 1.54)\
            #         and (gridcoord[i][1] - 1.54) <= at_pdbqt[k][1] <= (gridcoord[i][1] + 1.54)\
            #             and (gridcoord[i][2] - 1.54) <= at_pdbqt[k][2] <= (gridcoord[i][2] + 1.54):
                if (gridcoord[i][0]-at_pdbqt[k][0])**2+ \
                    (gridcoord[i][1]-at_pdbqt[k][1])**2+ \
                    (gridcoord[i][2]-at_pdbqt[k][2])**2 <= (1*1.54)**2:

                    dist = np.linalg.norm(np.array(gridcoord[i]) - np.array(at_pdbqt[k]))
                    # energy = 0.3125*(energy - 10) / (dist**0.5)
                    
                    dist_lis.append(dist)

            if len(dist_lis) != 0:
                dist_lis = np.array(dist_lis)
                sorted_ind = np.argsort(dist_lis)
                sorted_dist = dist_lis[sorted_ind[::-1]]
                for d in sorted_dist:
                    energy = 0.225*(energy - 10) / (d**0.5)
    
            # add positive value to the grid point with 0 energy
            # if energy == 0:
            #     dist1 = [np.linalg.norm(np.array(gridcoord[i]) - np.array(at_pdbqt[j])) for j in range(0, len(at_pdbqt))]
            #     energy = min(dist1)

            # if 0 <= energy <= 50:
            #     energy = 0.029 * energy ** 2.5 + 1

            gridmap.write(f"{energy:.3f}\n")



def general_map(ligand, folder, spacing, npts, gridcenter, gridcoord, ligandcoord, filename, at):
    with open(f"{folder}{filename}.{at}.map", "w") as gridmap:
        # write the document information into each map file
        gridmap.write("GRID_PARAMETER_FILE " + folder + "\ligand.gpf\n")
        gridmap.write(f"GRID_DATA_FILE {filename}.maps.fld\n")
        gridmap.write(f"MACROMOLECULE {filename}.pdbqt\n")
        gridmap.write("SPACING " + str(spacing) + "\n")
        gridmap.write("NELEMENTS " + str(npts[0]) + " " + str(npts[1]) + " " + str(npts[2]) + "\n")
        gridmap.write("CENTER " + str(gridcenter[0]) + " " + str(gridcenter[1]) + " " + str(gridcenter[2]) + "\n")
        for i in range(0, len(gridcoord)):
            gridmap.write("0" + "\n")



def elec_map(folder, spacing, npts, gridcenter, gridmap_coords, ligand_coords, atomic_partial_charge, filename):
    import numpy as np
    with open(f"{folder}{filename}.e.map", "w") as e_gridmap:
        e_gridmap.write("GRID_PARAMETER_FILE " + folder + "\ligand.gpf\n")
        e_gridmap.write(f"GRID_DATA_FILE {filename}.maps.fld\n")
        e_gridmap.write(f"MACROMOLECULE {filename}.pdbqt\n")
        e_gridmap.write("SPACING " + str(spacing) + "\n")
        e_gridmap.write("NELEMENTS " + str(npts[0]) + " " + str(npts[1]) + " " + str(npts[2]) + "\n")
        e_gridmap.write("CENTER " + str(gridcenter[0]) + " " + str(gridcenter[1]) + " " + str(gridcenter[2]) + "\n")


        for i in gridmap_coords:
            potential, count = 0, 0
            for j in ligand_coords:
                r = np.linalg.norm(np.array(i)-np.array(j))
                loop_potential = atomic_partial_charge[count] / (-0.1465*r)
                potential += loop_potential
                count += 1
            # if potential <= 0:
            #     potential = 0
            # potential = (potential) * 5
            e_gridmap.write(f"{potential:.3f}\n")

def des_map(at, ligand, folder, spacing, npts, gridcenter, gridcoord, ligandcoord, filename):
    import numpy as np
    # extract data of an atom type
    # extract the coordiante and convert from string to float
    at_pdbqt = []
    for i in range(0, len(ligand)):
        if ligand[i][-1] == at:
            at_pdbqt.append(ligandcoord[i])

    # generate map file of carbon
    with open(f"{folder}{filename}.{at}.map", "w") as gridmap:

        # write the document information into each map file
        gridmap.write("GRID_PARAMETER_FILE " + folder + "\ligand.gpf\n")
        gridmap.write(f"GRID_DATA_FILE {filename}.maps.fld\n")
        gridmap.write(f"MACROMOLECULE {filename}.pdbqt\n")
        gridmap.write("SPACING " + str(spacing) + "\n")
        gridmap.write("NELEMENTS " + str(npts[0]) + " " + str(npts[1]) + " " + str(npts[2]) + "\n")
        gridmap.write("CENTER " + str(gridcenter[0]) + " " + str(gridcenter[1]) + " " + str(gridcenter[2]) + "\n")

        # if the carbon atom is within 5 angstrom of a grid point
        # energy will be counted into that grid point based on the distance
        for i in range(0, len(gridcoord)):
            energy = 0
            for k in range(0, len(at_pdbqt)):
                if (gridcoord[i][0] - 1.54) <= at_pdbqt[k][0] <= (gridcoord[i][0] + 1.54) \
                        and (gridcoord[i][1] - 1.54) <= at_pdbqt[k][1] <= (gridcoord[i][1] + 1.54) \
                        and (gridcoord[i][2] - 1.54) <= at_pdbqt[k][2] <= (gridcoord[i][2] + 1.54):
                    # grid points collect energy based on distance to atoms
                    x_diff = gridcoord[i][0] - at_pdbqt[k][0]
                    y_diff = gridcoord[i][1] - at_pdbqt[k][1]
                    z_diff = gridcoord[i][2] - at_pdbqt[k][2]
                    dist = np.sqrt(x_diff ** 2 + y_diff ** 2 + z_diff ** 2)
                    energy = round(energy - 10 / dist, 8)

            # add positive value to the grid point with 0 energy
            gridmap.write(str(energy) + "\n")

# from boxcal import boxcal, center_from_pdbqt, round_up_to_even
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

# from create_folder import createFolder
def createFolder(libpath):
    import os, shutil
    dir_name = libpath
    files = os.listdir(dir_name)
    for i in files:
        os.mkdir(os.path.join(dir_name , i.split(".")[0]))
        shutil.copy(os.path.join(dir_name , i), os.path.join(dir_name , i.split(".")[0]))
        os.remove(os.path.join(dir_name, i))

# from modify_pdbqt import modify_pdbqt
def modify_pdbqt(pdbqt_file):
    # list to store file lines
    lines = []
    keywords = ("REMARK", "ROOT", "ENDROOT", "BRANCH", "ENDBRANCH", "TORSDOF")
    # read file
    with open(pdbqt_file, 'r') as fp:
        # read an store all lines into list
        lines = fp.readlines()

    # Write file
    with open(pdbqt_file, 'w') as fp:
        for line in lines:
            if line.startswith(keywords):
                continue
            else:
                fp.write(line)

# from roc_auc import roc_auc_ef

def roc_auc_ef(ranked_file, actives_file, candidate_name):

    from sklearn import metrics
    import matplotlib.pyplot as plt
    import numpy as np

    with open(actives_file, "r") as tempfile:
        actives = [line.strip("\n") for line in tempfile.readlines()]

    with open(ranked_file, "r") as tempfile:
        total = [line.strip("\n") for line in tempfile.readlines()]

    num_of_actives = len(actives)
    num_of_decoys = len(total) - num_of_actives

    act, dcy = 0, 0
    fpr, tpr = [], []
    for candidate in total:
        if candidate in actives:
            act += 1
        else:
            dcy += 1
        fpr.append(dcy/num_of_decoys)
        tpr.append(act/num_of_actives)

    with plt.style.context(['science', "high-vis", "grid"]):
        pparam = dict(xlabel='False Positive Rate', ylabel='True Positive Rate')
        fig, ax = plt.subplots()
        ax.set_title(f"ROC of {candidate_name}")
        ax.plot(fpr, tpr, label="NCF-LBVS")
        ax.plot(fpr, fpr, label="Random")
        ax.legend(title=f"AUC = {metrics.auc(fpr, tpr):.5f}")
        ax.autoscale(tight=True)
        ax.set(**pparam)
        ax.set_ylim(bottom=-0.03, top=1.03)
        ax.set_xlim(left=-0.03, right=1.03)
        fig.savefig(f"ROC_{candidate_name}.svg")

    data_set_count = len(total)
    active_ratio = np.divide(num_of_actives, data_set_count)
    Subsets = [[int(round(data_set_count * 0.01)), 1 ],
               [int(round(data_set_count * 0.05)), 5 ],
               [int(round(data_set_count * 0.10)), 10],
               [int(round(data_set_count * 0.15)), 15],
               [int(round(data_set_count * 0.20)), 20],
               [int(round(data_set_count * 0.25)), 25],
               [int(round(data_set_count * 0.50)), 50],
               [int(round(data_set_count * 0.99)),100]]

    EF_active_found  = 0.0
    active_found     = 0.0
    active_not_found = 0.0
    Active_Found     = []

    with open(f"EF_{candidate_name}.txt", "w") as tempfile:
        tempfile.write(f"""\
Number of total ligands: {data_set_count}\n
Number of known actives: {num_of_actives}\n
Number of known decoys: {num_of_decoys}\n
AUC : {metrics.auc(fpr, tpr):.5f}\n
""")
        for data_set_index, compound in enumerate(total):
            Name = compound.split()

            if data_set_index <= num_of_actives and Name[0] in actives:
                active_found += 1.0
            elif data_set_index > num_of_actives and Name[0] in actives:
                active_not_found += 1.0

            ## For enrichment factor calculation
            if Name[0] in actives:
                EF_active_found += 1.0
                Active_Found.append(Name[0])

        # Calculate the Enrichment factor at certain subset ratio
            for Subset in Subsets:
                if data_set_index == Subset[0]:
                    EF = float(EF_active_found / data_set_index) / active_ratio
                    tempfile.write(f"EF at {Subset[1]} % : {EF:.2f}\n")

# from extract_best_pose import extract_best_pos
def extract_best_pos(dlgfile):
    with open(dlgfile, "r") as tmpfile:
        lines = tmpfile.readlines()
        index = lines.index("    CLUSTERING HISTOGRAM\n")
        try:
            best_run_ind = int(lines[index + 10].split()[4])
        except:
            best_run_ind = int(lines[index + 10].split()[3])
        best_run_start = [lines.index(l) for l in lines if l.startswith(f"Run:   {best_run_ind} /")][0] + 3
        if best_run_ind == 100:
            return lines[best_run_start:]
        else:
            best_run_stop = [lines.index(l) for l in lines if l.startswith(f"Run:   {best_run_ind+1} /")][0] - 7
            return lines[best_run_start:best_run_stop]

to_be_screened = set(os.listdir("/root/autodl-tmp/dude-single/"))-set(['ampc'])
# to_be_screened = ["comt03"]
for dude_target in to_be_screened:
    #======================================================================================================================
    ##
    # Define directories of VS library and reference ligand
    ##
    #======================================================================================================================
    # define all the paths needed
    os.chdir("/root/autodl-tmp/adss-scripts/") # directory of this script
    lib_path = f'/root/autodl-tmp/dude-single/{dude_target}/docking_test/lib/' # directory of your VS library
    lig_path = f'/root/autodl-tmp/dude-single/{dude_target}/docking_test/lig/' # directory of your reference ligand
    # vs_lib_name = f'' # name of your VS library file including extension
    dlg_original_path = os.getcwd()
    for files in os.listdir(lig_path):
        if files.endswith(".pdbqt"):
            pdbqt_path = os.path.join(lig_path, files)
            grid_file_path = pdbqt_path.replace(".pdbqt", ".grid.txt")
    modify_pdbqt(pdbqt_path) # make the reference pdbqt file rigid
    #======================================================================================================================
    ##
    # Decompress the VS library and split it into single molecules (used for DUD-E benchmark)
    ##
    #======================================================================================================================
    subprocess.check_call(f"obabel -isdf {lib_path}actives_final.sdf.gz -osdf -O *.sdf --split", shell=True, cwd=lib_path)
    subprocess.check_call(f"obabel -isdf {lib_path}decoys_final.sdf.gz -osdf -O *.sdf --split", shell=True, cwd=lib_path)
    time.sleep(3)
    subprocess.check_call(f"rm {lib_path}actives_final.sdf.gz", shell=True)
    subprocess.check_call(f"rm {lib_path}actives_final.sdf", shell=True)
    subprocess.check_call(f"rm {lib_path}decoys_final.sdf.gz", shell=True)
    time.sleep(15)
    
    # subprocess.check_call(f"obabel -isdf {lib_path}{vs_lib_name} -osdf -O *.sdf --split", shell=True, cwd=lib_path)
    # time.sleep(3)
    # try:
    #     subprocess.call(f"rm {vs_lib_name.replace('.sdf.gz', '')}*", shell=True, cwd=lib_path)
    # except:
    #     pass
    #======================================================================================================================
    ##
    # Calculate the dimension of the grid box
    ##
    #======================================================================================================================
    for files in os.listdir(lig_path):
        if files.endswith(".pdbqt"):
            pdbqt_path = os.path.join(lig_path, files)
            ref_path = pdbqt_path.replace(".pdbqt", ".sdf")
            subprocess.check_call(f"obabel -ipdbqt {pdbqt_path} -osdf -O {ref_path}", shell=True)
            subprocess.check_call(f"cp {ref_path} {lib_path}", shell=True)

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
   
    supported_atom_types = ["HD", "C", "A", "N", "NA", "OA", "F", "P", "SA", "S",
                            "Cl", "Br", "I", "Mg", "Ca", "Mn", "Fe", "Zn", "d"]

    query_ligmap = gridmap_information.from_gridmap(grid_file_path, pdbqt_path)
    atomtype_not_in_rec = list(set(supported_atom_types) - set(query_ligmap.target_atomtype))

    elec_map(folder=lig_path, spacing=query_ligmap.spacing, npts=query_ligmap.npts,
             gridcenter=query_ligmap.gridcenter, gridmap_coords=query_ligmap.target_gridmap_coords,
             ligand_coords=query_ligmap.target_coords, atomic_partial_charge=query_ligmap.target_charge,
             filename=Path(pdbqt_path).stem)

    affinity_mapping_partial = functools.partial(affinity_mapping, query_ligmap.target_ligand, lig_path, query_ligmap.spacing,
                                                 query_ligmap.npts, query_ligmap.gridcenter, query_ligmap.target_gridmap_coords,
                                                 query_ligmap.target_coords, Path(pdbqt_path).stem)
    
    general_map_partial = functools.partial(general_map, query_ligmap.target_ligand, lig_path, query_ligmap.spacing,
                                            query_ligmap.npts, query_ligmap.gridcenter, query_ligmap.target_gridmap_coords,
                                            query_ligmap.target_coords, Path(pdbqt_path).stem)

    with WorkerPool(n_jobs=4) as pool:
        pool.map(affinity_mapping_partial, query_ligmap.target_atomtype)
        pool.map(general_map_partial, atomtype_not_in_rec)

    create_mapfld(folder=lig_path, filename=Path(pdbqt_path).stem, spacing=query_ligmap.spacing,
                  npts=query_ligmap.npts, center=query_ligmap.gridcenter)
    print("Finished creating map files")
   
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

    with WorkerPool(n_jobs=192) as pool:
        pool.map(multiprocess_libpath, os.listdir(lib_path))

    env = os.environ.copy().update(OMP_NUM_THREADS=192)
    # !!! Change the autodock-gpu path to yours
    subprocess.run(f"/root/AutoDock-GPU/bin/autodock_gpu_128wi -B {lig_path}/docking_file.lst -D all --nrun 100 --xmloutput 0", shell=True, env=env)

    def multiprocess_move_dlg(files):
        if files.endswith(".dlg"):
            dlg_path = os.path.join(dlg_original_path, files)
            subprocess.check_call(f"mv {dlg_path} {lib_path}{Path(dlg_path).stem}/{Path(dlg_path).stem}.dlg", shell=True)
    with WorkerPool(n_jobs=192) as pool:
        pool.map(multiprocess_move_dlg, os.listdir(dlg_original_path), progress_bar=True)
   
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

    ref_ind = name_lis.index(f"{Path(ref_path).stem}")
    ref_score = score_lis[ref_ind]
    del name_lis[ref_ind]
    del score_lis[ref_ind]

    normalized_score_lis = np.array(score_lis)/ref_score
    sorted_index = np.argsort(normalized_score_lis)[::-1]
    sorted_name = [name_lis[i] for i in sorted_index]
    sorted_score = [normalized_score_lis[i] for i in sorted_index]

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

    roc_auc_ef(ranked_file="ranked.txt",
               actives_file="actives.txt",
               candidate_name=f"{dude_target}") # change reference ligand name to yours
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
    with WorkerPool(n_jobs=192) as pool:
        pool.map(multiprocess_extraction, os.listdir(lib_path), progress_bar=True)

    subprocess.check_call(f"find . -name '*.dlg' -type f -delete", shell=True, cwd=lib_path)
    subprocess.check_call(f"find . -name '*.pdb' -type f -delete", shell=True, cwd=lib_path)
    subprocess.check_call(f"find . -name '*.sdf' -type f -delete", shell=True, cwd=lib_path)
