### DYNAmical Multi-planet Injection TEster (DYNAMITE) ###
### HPC Multi-Core Merging Capabilities ###
### Jeremy Dietrich ###
### jdietrich1@email.arizona.edu ###
### 2022 January 10 ###
### Version 2.0 ###
### Dietrich & Apai (2020), AJ, 160, 107D ###
### Dietrich & Apai (2021), AJ, 161, 17D ###
### Dietrich, Apai, & Malhotra (2022), accepted to AJ ###
### https://iopscience.iop.org/article/10.3847/1538-3881/aba61d ###

import os
import sys
import glob
import numpy as np
from dynamite_v2 import dynamite

glb = len(glob.glob1(".","saved_data_*.npz"))

for i in range(glb):
    with np.load("saved_data_" + str(i+1) + ".npz", allow_pickle=True) as data:
        if i == 0:
            data1 = data["data"]

        else:
            for j in range(len(data["data"])):
                if isinstance(data1[j], list):
                    if isinstance(data["data"][j], list):
                        try:
                            if data1[j] != data["data"][j]:
                                data1[j] = data1[j] + data["data"][j]

                        except:
                            print(data1[j], data["data"][j])

                    else:
                        if data1[j] != list(data["data"][j]):
                            data1[j].append(data["data"][j])

                elif isinstance(data1[j], np.ndarray):
                    if list(data1[j]) != list(data["data"][j]):
                        np.append(data1[j], data["data"][j])

with open("rejected_values.txt", "w") as f:
    for i in range(glb):
        with open("rejected_values_" + str(i+1) + ".txt") as fi:
            f.write(fi.read())

with open("unstable_times.txt", "w") as f:
    for i in range(glb):
        with open("unstable_times_" + str(i+1) + ".txt") as fi:
            f.write(fi.read())

if 'data1' in locals():
    np.savez("saved_data.npz", data=data1)
    dynamite(merged_data=data1)

    if len(sys.argv) == 1:
        for i in range(glb):
            os.remove("saved_data_" + str(i+1) + ".npz")
            os.remove("rejected_values_" + str(i+1) + ".txt")
            os.remove("unstable_times_" + str(i+1) + ".txt")
