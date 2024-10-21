### DYNAmical Multi-planet Injection TEster (DYNAMITE) ###
### HPC Multi-Core Merging Capabilities ###
### Jamie Dietrich ###
### jdietrich@asu.edu ###
### 2024 October 21 ###
### Version 3.0 ###
### Dietrich & Apai (2020), AJ, 160, 107D ###
### https://iopscience.iop.org/article/10.3847/1538-3881/aba61d ###
### Dietrich & Apai (2021), AJ, 161, 17D ###
### Dietrich, Apai, & Malhotra (2022), AJ, 163, 88D ###
### Basant, Dietrich, & Apai (2022), AJ, 164, 12B ###
### Basant, Dietrich, & Apai (2022), RNAAS, 6, 213 ###
### Dietrich (2024), AJ, 168, 119 ###

import os
import sys
import glob
import numpy as np
from dynamite import dynamite

glb = len(glob.glob1(".","saved_data_*.npz"))
print(len(glb))
for i in range(glb):
    with np.load("saved_data_" + str(i+1) + ".npz", allow_pickle=True) as data:
        if i == 0:
            data1 = data["data"]
            print("Data1!")
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

if 'data1' in locals():
    np.savez("saved_data.npz", data=data1)
    dynamite(merged_data=data1)

    if len(sys.argv) == 1:
        for i in range(glb):
            os.remove("saved_data_" + str(i+1) + ".npz")
            os.remove("rejected_values_" + str(i+1) + ".txt")
            os.remove("unstable_times_" + str(i+1) + ".txt")
