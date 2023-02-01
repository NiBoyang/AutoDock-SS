import os, subprocess

lib_path = '/root/autodl-tmp/lbvs/all_dude/vgfr2/docking_test/lib/'

for folders in os.listdir(lib_path):
    mol_folders = os.path.join(lib_path, folders)
    for files in os.listdir(mol_folders):
        if files.endswith(".dlg"):
            dlg_files = os.path.join(mol_folders, files)
            subprocess.check_call(f"rm {dlg_files}", shell=True)