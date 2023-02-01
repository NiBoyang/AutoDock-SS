import os, shutil, re, subprocess
from pathlib import Path
import numpy as np
from rdkit import Chem

lib_path = "/root/autodl-tmp/lbvs/all_dude/ampc/docking_test/lib/"

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
                    score = float(lines[index+10].split()[2])
                    score_lis.append(score)
                # if filename.startswith("CHEMBL"):
                #     actives_name.append(name)
for names in name_lis:
    if names.startswith("ZINC") == False:
        actives_name.append(names)

sorted_index = np.argsort(np.array(score_lis))
sorted_name = [name_lis[i] for i in sorted_index]
#==============================================================
# actives_path = "/root/autodl-tmp/lbvs/dude/actives_final.sdf"
# mol_parse = [x for x in Chem.SDMolSupplier(actives_path)]
# actives_smiles = [Chem.MolToSmiles(y) for y in mol_parse]
# actives_name = [mymol.GetProp("_Name") for mymol in mol_parse]
# print(actives_smiles)

with open("/root/autodl-tmp/lbvs/all_dude/ampc/docking_test/ranked.txt", "w") as tempfile:
    [tempfile.write(f"{name}\n") for name in sorted_name]
    tempfile.truncate(tempfile.tell()-1)
    
with open("/root/autodl-tmp/lbvs/all_dude/ampc/docking_test/actives.txt", "w") as tempfile:
    [tempfile.write(f"{name}\n") for name in set(actives_name)]
    tempfile.truncate(tempfile.tell()-1)