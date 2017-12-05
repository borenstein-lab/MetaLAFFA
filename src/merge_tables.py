#!/bin/env python
#
# Author: Adrian Verster
# Date: 11/30/2017

import sys
import pandas as pd

if __name__ == "__main__":
    infiles = sys.argv[1:]
    DfFinal = pd.DataFrame(columns = ["KO"])
    for inf in infiles:
        Df = pd.read_csv( inf, sep = "\t", header = 0, compression = "gzip")
        DfFinal = DfFinal.merge(Df, how = "outer", on = ["KO"])
    DfFinal.fillna(0, inplace = True)
    DfFinal.to_csv(sys.stdout, sep = "\t", index = False)
