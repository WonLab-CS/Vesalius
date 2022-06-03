#-----------------------------/Running SpaGCN/---------------------------------#
#------------------------------------------------------------------------------#
# Running SpaGCN on Simulated data sets for the purpose of benchmarking
#------------------------------------------------------------------------------#
import os,csv,re,sys,timeit,time
import pandas as pd
import numpy as np
import scanpy as sc
import math
import SpaGCN as spg
from scipy.sparse import issparse
import random, torch
import warnings
from igraph import compare_communities
#------------------------------------------------------------------------------#
# NOTE: Simulated data sets are produced by the simulate.R file and
# saved as csv files.
#------------------------------------------------------------------------------#


# Create output directory
input = '/home/pcnmartin/group/slide_seqV2/'
output = os.path.join(input,'spaGCNSim/')
if !path.exists(output):
    os.mkdir(output)

# create time file
time = os.path.join(input,'spaGCNSim/time.txt')
if !path.exists(time):
    ft = open(time,'x').close()

# create performance file
perf = os.path.join(input,'spaGCNSim/performance.txt')
if !path.exists(perf):
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
    adata = sc.AnnData(X :subCounts)
    adata.var_names_make_unique()
    ##### read data
    # adata.obs["x1"]=spatial[1]
    # adata.obs["x2"]=spatial[2]
    # adata.obs["x3"]=spatial[3]
    # adata.obs["x4"]=spatial[4]
    # adata.obs["x5"]=spatial[5]
    # adata.obs["x_array"]=adata.obs["x2"]
    # adata.obs["y_array"]=adata.obs["x3"]
    # adata.obs["x_pixel"]=adata.obs["x4"]
    # adata.obs["y_pixel"]=adata.obs["x5"]
    # adata=adata[adata.obs["x1"]==1]
    # adata.var_names=[i.upper() for i in list(adata.var_names)]
    # adata.var["genename"]=adata.var.index.astype("str")
    # adata.write_h5ad(f"{dir_output}/sample_data.h5ad")

    #Calculate adjacent matrix
    b=49
    a=1
    # Check for no image approach!!!!
    adj=calculate_adj_matrix(x=x_pixel,y=y_pixel, histology=False)
    np.savetxt(f'{dir_output}/adj.csv', adj, delimiter=',')


    ##### Spatial domain detection using SpaGCN
    spg.prefilter_genes(adata, min_cells=3) # avoiding all genes are zeros
    spg.prefilter_specialgenes(adata)
    #Normalize and take log for UMI
    sc.pp.normalize_per_cell(adata)
    sc.pp.log1p(adata)


    ### 4.2 Set hyper-parameters

    p=0.5
    spg.test_l(adj,[1, 10, 100, 500, 1000])
    l=spg.find_l(p=p,adj=adj,start=100, end=500,sep=1, tol=0.01)

    r_seed=t_seed=n_seed=100
    res=spg.search_res(adata, adj, l, q, start=0.7, step=0.1, tol=5e-3, lr=0.05, max_epochs=20, r_seed=r_seed,
    t_seed=t_seed, n_seed=n_seed)

    ### 4.3 Run SpaGCN
    clf=spg.SpaGCN()
    clf.set_l(l)

    #Set seed
    random.seed(r_seed)
    torch.manual_seed(t_seed)
    np.random.seed(n_seed)
    #Run
    clf.train(adata,adj,init_spa=True,init="louvain",res=res, tol=5e-3, lr=0.05, max_epochs=200)
    y_pred, prob=clf.predict()
    adata.obs["pred"]= y_pred
    adata.obs["pred"]=adata.obs["pred"].astype('category')
    #Do cluster refinement(optional)
    adj_2d=spg.calculate_adj_matrix(x=x_array,y=y_array, histology=False)
    refined_pred=spg.refine(sample_id=adata.obs.index.tolist(), pred=adata.obs["pred"].tolist(), dis=adj_2d, shape="hexagon")
    adata.obs["refined_pred"]=refined_pred
    adata.obs["refined_pred"]=adata.obs["refined_pred"].astype('category')
    e = timeit.default_timer() - s

    ftime = fileTag[0] + "," + str(e) + "," + "sec" + "\n"
    ft = open(time,'a')
    ft.write(ftime)
    ft.close()
    ### TO adjust!!!!!!!!
    ### Should I run this in R? Make sure that the same package is used
    ### in all cases
    #ari = compare_communities(comm1,comm2,method = "adjusted_rand")
    #vi = compare_communities(comm1,comm2,method = "vi")
    fperf = fileTag[0] + "," + str(ari) + "," + str(vi) + "\n"
    fp = open(perf,'a')
    fp.write(fperf)
    fp.close()
    #Save results
    adata.write_h5ad(f"{dir_output}/results.h5ad")

    adata.obs.to_csv(f'{dir_output}/metadata.tsv', sep='\t')
