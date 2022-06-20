#-----------------------------/Running STAGATE/--------------------------------#
#------------------------------------------------------------------------------#
# Running STAGATE on Simulated data sets for the purpose of benchmarking
#------------------------------------------------------------------------------#
import os
os.environ["OMP_NUM_THREADS"] = "1" # export OMP_NUM_THREADS=1
os.environ["OPENBLAS_NUM_THREADS"] = "1" # export OPENBLAS_NUM_THREADS=1
os.environ["MKL_NUM_THREADS"] = "1" # export MKL_NUM_THREADS=1
os.environ["VECLIB_MAXIMUM_THREADS"] = "1" # export VECLIB_MAXIMUM_THREADS=1
os.environ["NUMEXPR_NUM_THREADS"] = "1" # export NUMEXPR_NUM_THREADS=1
import re
import pandas as pd
import numpy as np
import scanpy as sc
#import matplotlib.pyplot as plt

import sys
import STAGATE
import timeit
import time


import tensorflow as tf

session_conf = tf.compat.v1.ConfigProto(
      intra_op_parallelism_threads=1,
      inter_op_parallelism_threads=1)
sess = tf.compat.v1.Session(config=session_conf)

from igraph import compare_communities

#------------------------------------------------------------------------------#
# NOTE: Simulated data sets are produced by the simulate.R file and
# saved as csv files.
#------------------------------------------------------------------------------#


# Create output directory
input = '/home/pcnmartin/group/slide_seqV2/'
output = os.path.join(input,'stagateSim/')
if not os.path.exists(output):
    os.mkdir(output)

# create time file
time = os.path.join(input,'stagateSim/time.txt')
if not os.path.exists(time):
    ft = open(time,'x').close()

# create performance file
perf = os.path.join(input,'stagateSim/performance.txt')
if not os.path.exists(perf):
    fp = open(perf,'x').close()

# Listing simulation files
simInput = "/home/pcnmartin/Vesalius/Simulation"
simFiles = os.listdir(simInput)
simFiles = [i for i in simFiles if i.endswith(".csv")]
fileTag = [re.sub(".csv","",i) for i in simFiles]
simFiles = [os.path.join(simInput,i) for i in simFiles]

# Loading counts
counts_file =  os.path.join(input,'Puck_200115_08.digital_expression.txt.gz')
counts = pd.read_csv(counts_file, sep='\t', index_col=0)



# Run benchmarking
for i in range(0,len(simFiles)):
    sim = pd.read_csv(simFiles[i],index_col=6)
    #sim.rows = list(sim["simBarcode"])
    subCounts = counts.loc[:,list(sim["barcodes"])]
    subCounts.columns = list(sim.index)
    # if re.search("dot",fileTag[i]):
    #     q = 6
    # else:
    #     q = 3
    s = timeit.default_timer()
    adata = sc.AnnData(subCounts.T,dtype = np.float32)
    adata.var_names_make_unique()
    coord_df = sim.loc[:, ['x', 'y']]
    adata.obsm["spatial"] = coord_df.to_numpy()
    sc.pp.calculate_qc_metrics(adata, inplace=True)
    sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=2000)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    STAGATE.Cal_Spatial_Net(adata, rad_cutoff=50)
    STAGATE.Stats_Spatial_Net(adata)

    adata = STAGATE.train_STAGATE(adata, alpha=0)

    sc.pp.neighbors(adata, use_rep='STAGATE')
    sc.tl.umap(adata)
    sc.tl.louvain(adata, resolution=0.5)
    export = pd.concat([sim,adata.obs["louvain"]], axis =1)

    e = timeit.default_timer() - s
    ftime = fileTag[i] + "," + str(e) + "," + "sec" + "\n"
    ft = open(time,'a')
    ft.write(ftime)
    ft.close()
    ### TO adjust!!!!!!!!
    ### Should I run this in R? Make sure that the same package is used
    ### in all cases


    #fperf = fileTag[0] + "," + str(ari) + "," + str(vi) + "\n"
    #fp = open(perf,'a')
    #fp.write(fperf)
    #fp.close()

    filename = 'STAGATE_' + fileTag[i] + ".csv"
    filename = os.path.join(output,filename)
    export.to_csv(filename)
