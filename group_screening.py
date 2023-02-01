import os, subprocess
from pathlib import Path

lib_path = '/root/autodl-tmp/lbvs/all_dude/pur2/docking_test/lib/'
lig_path = '/root/autodl-tmp/lbvs/all_dude/pur2/docking_test/lig/'
for files in os.listdir(lig_path):
    if files.endswith(".maps.fld"):
        mapfld_path = os.path.join(lig_path, files)
# mapfld_path = "/root/autodl-tmp/lbvs/all_dude/abl1/docking_test/lig/lig.maps.fld"

with open(f"{lig_path}/docking_file.lst", "w") as lstfile:
    lstfile.write(f"{mapfld_path}\n")

    for folder in os.listdir(lib_path):
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
                    os.system("python ligandprep.py -l %s -v -o %s" % (pdb_filepath, pdbqt_filepath))

                    lstfile.write(f"{pdbqt_filepath}\n")
                    lstfile.write(f"{Path(pdbqt_filepath).stem}\n")

subprocess.call("cd /root/autodl-tmp/lbvs/dlg-storage/", shell=True)
subprocess.call(f"/root/AutoDock-GPU/bin/autodock_gpu_64wi -B {lig_path}/docking_file.lst -D all --xmloutput 0", shell=True)