import subprocess, os, argparse
from pathlib import Path
import numpy as np
from roc_auc import roc_auc_ef

parser = argparse.ArgumentParser(description='dude target')
parser.add_argument('target')
args = parser.parse_args()
dude_target = args.target

os.chdir("/root/autodl-tmp/lbvs/improved_code/")
lib_path = f'/root/autodl-tmp/lbvs/all_dude/{dude_target}/docking_test/lib/'
lig_path = f'/root/autodl-tmp/lbvs/all_dude/{dude_target}/docking_test/lig/'
dlg_original_path = "/root/autodl-tmp/lbvs/improved_code/"

subprocess.check_call(f"/root/AutoDock-GPU/bin/autodock_gpu_64wi -B {lig_path}/docking_file.lst -D all --nrun 100 --xmloutput 0", shell=True)

for files in os.listdir(dlg_original_path):
    if files.endswith(".dlg"):
        dlg_path = os.path.join(dlg_original_path, files)
        subprocess.check_call(f"mv {dlg_path} {lib_path}{Path(dlg_path).stem}/{Path(dlg_path).stem}.dlg", shell=True)
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

sorted_index = np.argsort(np.array(score_lis))
sorted_name = [name_lis[i] for i in sorted_index]

with open(f"{lib_path.replace('lib/', '')}ranked.txt", "w") as tempfile:
    [tempfile.write(f"{name}\n") for name in sorted_name]
    tempfile.truncate(tempfile.tell()-1)

with open(f"{lib_path.replace('lib/', '')}actives.txt", "w") as tempfile:
    [tempfile.write(f"{name}\n") for name in set(actives_name)]
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