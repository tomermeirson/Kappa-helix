# Kappa-helix

# This is an example to assign the secondary structure k-helix (commonly known as polyproline type II helix) 

# Setting the path to STRIDE is required to run the function
# STRIDE (stride_WIN32.exe) can be downloaded from the following URL:
# http://webclu.bio.wzw.tum.de/stride/

path.stride = getwd()

source('Assign_khelix.R')

assign.ss.khelix.stride(pdb = '3D06',chain.id = 'A',path.stride=path.stride)
