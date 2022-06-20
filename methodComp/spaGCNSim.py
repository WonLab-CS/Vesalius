#-----------------------------/Running SpaGCN/---------------------------------#
#------------------------------------------------------------------------------#
# Running SpaGCN on Simulated data sets for the purpose of benchmarking
#------------------------------------------------------------------------------#
import os,csv,re
import pandas as pd
import numpy as np
import scanpy as sc
import math
import SpaGCN as spg
from scipy.sparse import issparse
import random, torch
import warnings
import timeit
import time
warnings.filterwarnings("ignore")
torch.set_num_threads(1)

#------------------------------------------------------------------------------#
# NOTE: Simulated data sets are produced by the simulate.R file and
# saved as csv files.
#------------------------------------------------------------------------------#


# Create output directory
input = '/home/pcnmartin/group/slide_seqV2/'
output = os.path.join(input,'spagcnSim/')
if not os.path.exists(output):
    os.mkdir(output)

# create time file
time = os.path.join(input,'spagcnSim/time.txt')
if not os.path.exists(time):
    ft = open(time,'x').close()

# create performance file
perf = os.path.join(input,'spagcnSim/performance.txt')
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


for i in range(0,len(simFiles)):
    sim = pd.read_csv(simFiles[i],index_col=6)
    #sim.rows = list(sim["simBarcode"])
    subCounts = counts.loc[:,list(sim["barcodes"])]
    subCounts.columns = list(sim.index)
    if re.search("dot",fileTag[i]):
         q = 6
    else:
         q = 3
    s = timeit.default_timer()

    adata = sc.AnnData(subCounts.T)
    adata.var_names_make_unique()
    ##### read data

    adata.obs["x_array"]=sim["x"]
    adata.obs["y_array"]=sim["y"]
    adata.obs["x_pixel"]=sim["x"]
    adata.obs["y_pixel"]=sim["y"]
    x_array=adata.obs["x_array"].tolist()
    y_array=adata.obs["y_array"].tolist()
    x_pixel=adata.obs["x_pixel"].tolist()
    y_pixel=adata.obs["y_pixel"].tolist()

    #Calculate adjacent matrix
    b=49
    a=1
    # Check for no image approach!!!!
    adj=spg.calculate_adj_matrix(x=x_pixel,y=y_pixel, histology=False)



    ##### Spatial domain detection using SpaGCN
    spg.prefilter_genes(adata, min_cells=3) # avoiding all genes are zeros
    spg.prefilter_specialgenes(adata)
    #Normalize and take log for UMI
    sc.pp.normalize_per_cell(adata)
    sc.pp.log1p(adata)


    ### 4.2 Set hyper-parameters

    p=0.5
    spg.test_l(adj,[1, 10, 100, 500, 1000])

    l=spg.find_l(p=p,adj=adj,start=0, end=2500,sep=0.5, tol=0.01)
    print(type(l))
    if l is not np.float64:
        l=np.float64(10)

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
    export = pd.concat([sim,adata.obs["pred"],adata.obs["refined_pred"]], axis =1)
    e = timeit.default_timer() - s

    ftime = fileTag[i] + "," + str(e) + "," + "sec" + "\n"
    ft = open(time,'a')
    ft.write(ftime)
    ft.close()

    filename = 'SpaGCN_' + fileTag[i] + ".csv"
    filename = os.path.join(output,filename)
    export.to_csv(filename)
