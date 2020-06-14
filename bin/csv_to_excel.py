#!/usr/bin/env python

import pandas as pd
import sys

path_in = sys.argv[1]
assert path_in.endswith(".csv"), f"wrong input expected file to end with .csv but got {path_in[-len('.csv'):]}"
#path_out = path_in[:-len("csv")]+"xlsx"
path_out = sys.argv[2]
df = pd.read_csv(path_in,sep=",")
df.to_excel(path_out)