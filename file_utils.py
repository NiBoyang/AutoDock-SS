import os
import shutil

def create_folder(libpath):
    """
    Create a folder for each file in the given directory and move the file into its respective folder.
    """
    files = os.listdir(libpath)
    for file in files:
        folder_path = os.path.join(libpath, file.split(".")[0])
        os.mkdir(folder_path)
        shutil.move(os.path.join(libpath, file), folder_path)

# from modify_pdbqt import modify_pdbqt
def modify_pdbqt(pdbqt_file):
    """
    Modify the given pdbqt file by removing lines that start with specific keywords.
    """
    keywords = ("REMARK", "ROOT", "ENDROOT", "BRANCH", "ENDBRANCH", "TORSDOF")
    
    with open(pdbqt_file, 'r') as fp:
        lines = [line for line in fp if not line.startswith(keywords)]

    with open(pdbqt_file, 'w') as fp:
        fp.writelines(lines)

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