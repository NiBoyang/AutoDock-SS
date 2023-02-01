import sys, os, subprocess


os.chdir("/root/autodl-tmp/lbvs/new_scoring/")
lib_path = f'/root/autodl-tmp/lbvs/all_dude/aa2ar/docking_test/lib/'
lig_path = f'/root/autodl-tmp/lbvs/all_dude/aa2ar/docking_test/lig/'
dlg_original_path = "/root/autodl-tmp/lbvs/new_scoring/"

a = int(subprocess.check_output(f"find {lig_path} -name '*.pdbqt' -type f | wc -l", shell=True).decode(sys.stdout.encoding).strip())
print(a)