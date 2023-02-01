from pathlib import Path
import os
import subprocess

dlg_original_path = "/root/autodl-tmp/lbvs/code/"
to_path = "/root/autodl-tmp/lbvs/all_dude/pur2/docking_test/lib/"

for files in os.listdir(dlg_original_path):
    if files.endswith(".dlg"):
        dlg_path = os.path.join(dlg_original_path, files)
        subprocess.check_call(f"mv {dlg_path} {to_path}{Path(dlg_path).stem}/{Path(dlg_path).stem}.dlg", shell=True)
        # Path(dlg_path).rename(f"{to_path}{Path(dlg_path).stem}.dlg")
