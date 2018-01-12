import os

def make_directory(dir_name):
    path = os.path.join(os.getcwd(), dir_name)
    if not os.path.isdir(path):
        os.mkdir(path)

    return path
