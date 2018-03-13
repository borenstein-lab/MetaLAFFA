#!/net/borenstein/vol1/PROGRAMS/python2/bin/python
#
# Author: Adrian Verster
# Date: 11/30/2017

import sys
sys.path = ['', '/net/borenstein/vol1/PROGRAMS/python2/lib/python27.zip', '/net/borenstein/vol1/PROGRAMS/python2/lib/python2.7', '/net/borenstein/vol1/PROGRAMS/python2/lib/python2.7/plat-linux2', '/net/borenstein/vol1/PROGRAMS/python2/lib/python2.7/lib-tk', '/net/borenstein/vol1/PROGRAMS/python2/lib/python2.7/lib-old', '/net/borenstein/vol1/PROGRAMS/python2/lib/python2.7/lib-dynload', '/net/borenstein/vol1/PROGRAMS/python2/lib/python2.7/site-packages']
import pandas as pd

if __name__ == "__main__":
    infiles = sys.argv[1:]
    DfFinal = pd.DataFrame(columns = ["KO"])
    for inf in infiles:
        Df = pd.read_csv( inf, sep = "\t", header = 0, compression = "gzip")
        DfFinal = DfFinal.merge(Df, how = "outer", on = ["KO"])
    DfFinal.fillna(0, inplace = True)
    DfFinal.to_csv(sys.stdout, sep = "\t", index = False)
