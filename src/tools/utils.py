import os

def checkdirs(directories):
    missing_directories = []
    for directory in directories:
        if not os.path.isdir(directory):
            print("No such directory: %s" % directory)
            missing_directories += [directory]
    if len(missing_directories) > 0:
        raise
    return