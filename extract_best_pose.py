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



if __name__ == "__main__": 
    import os
    from pathlib import Path
    from mpire import WorkerPool

    folders = os.listdir("/root/autodl-tmp/lbvs/all_dude/")
    for i in folders:
        lib_path = f"/root/autodl-tmp/lbvs/all_dude/{i}/docking_test/lib/"
    # for subfolder in os.listdir(lib_path):
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
        with WorkerPool(n_jobs=96) as pool:
            pool.map(multiprocess_extraction, os.listdir(lib_path), progress_bar=True)
    
