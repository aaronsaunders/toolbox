# http://stackoverflow.com/questions/2186525/use-a-glob-to-find-files-recursively-in-python

import os.path.walk, fnmatch

def find_files(treeroot, pattern):
    results = []
    for base, dirs, files in os.path.walk(treeroot):
        goodfiles = fnmatch.filter(files, pattern)
        results.extend(os.path.join(base, f) for f in goodfiles)
    return results
