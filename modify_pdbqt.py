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