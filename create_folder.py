def createFolder(libpath):
    import os, shutil
    dir_name = libpath
    files = os.listdir(dir_name)
    for i in files:
        os.mkdir(os.path.join(dir_name , i.split(".")[0]))
        shutil.copy(os.path.join(dir_name , i), os.path.join(dir_name , i.split(".")[0]))
        os.remove(os.path.join(dir_name, i))