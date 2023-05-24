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
