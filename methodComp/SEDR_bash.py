

#-----------------------------/ Running SEDR /---------------------------------#
#------------------------------------------------------------------------------#
# Running SEDR on Simulated data sets for the purpose of benchmarking
#------------------------------------------------------------------------------#
import os
os.environ["OMP_NUM_THREADS"] = "1" # export OMP_NUM_THREADS=1
os.environ["OPENBLAS_NUM_THREADS"] = "1" # export OPENBLAS_NUM_THREADS=1
os.environ["MKL_NUM_THREADS"] = "1" # export MKL_NUM_THREADS=1
os.environ["VECLIB_MAXIMUM_THREADS"] = "1" # export VECLIB_MAXIMUM_THREADS=1
os.environ["NUMEXPR_NUM_THREADS"] = "1" # export NUMEXPR_NUM_THREADS=1

import torch
device = torch.device('cpu')
torch.set_num_threads(1)
import argparse
import warnings
import numpy as np
import pandas as pd
from src.graph_func import graph_construction
from src.utils_func import mk_dir, adata_preprocess, load_ST_file
import anndata
from src.SEDR_train import SEDR_Train
from sklearn import metrics
import scanpy as sc
import timeit
import time
import sys,re
import scipy.sparse as sp
from scipy.spatial import distance


#torch.cuda.cudnn_enabled = False

np.random.seed(0)
torch.manual_seed(0)





# ################ Parameter setting

parser = argparse.ArgumentParser()
parser.add_argument('--file', type=str)
parser.add_argument('--k', type=int, default=10, help='parameter k in spatial graph')
parser.add_argument('--knn_distanceType', type=str, default='euclidean',
                    help='graph distance type: euclidean/cosine/correlation')
parser.add_argument('--epochs', type=int, default=200, help='Number of epochs to train.')
parser.add_argument('--cell_feat_dim', type=int, default=200, help='Dim of PCA')
parser.add_argument('--feat_hidden1', type=int, default=100, help='Dim of DNN hidden 1-layer.')
parser.add_argument('--feat_hidden2', type=int, default=20, help='Dim of DNN hidden 2-layer.')
parser.add_argument('--gcn_hidden1', type=int, default=32, help='Dim of GCN hidden 1-layer.')
parser.add_argument('--gcn_hidden2', type=int, default=8, help='Dim of GCN hidden 2-layer.')
parser.add_argument('--p_drop', type=float, default=0.2, help='Dropout rate.')
parser.add_argument('--using_dec', type=bool, default=True, help='Using DEC loss.')
parser.add_argument('--using_mask', type=bool, default=False, help='Using mask for multi-dataset.')
parser.add_argument('--feat_w', type=float, default=10, help='Weight of DNN loss.')
parser.add_argument('--gcn_w', type=float, default=0.1, help='Weight of GCN loss.')
parser.add_argument('--dec_kl_w', type=float, default=10, help='Weight of DEC loss.')
parser.add_argument('--gcn_lr', type=float, default=0.01, help='Initial GNN learning rate.')
parser.add_argument('--gcn_decay', type=float, default=0.01, help='Initial decay rate.')
parser.add_argument('--dec_cluster_n', type=int, default=10, help='DEC cluster number.')
parser.add_argument('--dec_interval', type=int, default=20, help='DEC interval nnumber.')
parser.add_argument('--dec_tol', type=float, default=0.00, help='DEC tol.')
# ______________ Eval clustering Setting ______________
parser.add_argument('--eval_resolution', type=int, default=1, help='Eval cluster number.')
parser.add_argument('--eval_graph_n', type=int, default=20, help='Eval graph kN tol.')

params = parser.parse_args()
params.device = device

warnings.filterwarnings("ignore")

#------------------------------------------------------------------------------#
# NOTE: Simulated data sets are produced by the simulate.R file and
# saved as csv files.
#------------------------------------------------------------------------------#


# Create output directory
input = '/home/pcnmartin/group/slide_seqV2/'
output = os.path.join(input,'sedrSim/')
if not os.path.exists(output):
    os.mkdir(output)

# create time file
time = os.path.join(input,'sedrSim/time.txt')
if not os.path.exists(time):
    ft = open(time,'x').close()

# create performance file
perf = os.path.join(input,'sedrSim/performance.txt')
if not os.path.exists(perf):
    fp = open(perf,'x').close()


# Loading counts
counts_file =  os.path.join(input,'Puck_200115_08.digital_expression.txt.gz')
counts = pd.read_csv(counts_file, sep='\t', index_col=0)





def res_search_fixed_clus(adata, fixed_clus_count, increment=0.02):
    '''
        arg1(adata)[AnnData matrix]
        arg2(fixed_clus_count)[int]

        return:
            resolution[int]
    '''
    for res in sorted(list(np.arange(0.2, 2.5, increment)), reverse=True):
        sc.tl.leiden(adata, random_state=0, resolution=res)
        count_unique_leiden = len(pd.DataFrame(adata.obs['leiden']).leiden.unique())
        if count_unique_leiden == fixed_clus_count:
            break
    return res



simFiles = os.path.join(params.file)
fileTag = re.sub(".csv","",simFiles)
fileTag = re.sub("/home/pcnmartin/Vesalius/Simulation/","",fileTag)
sim = pd.read_csv(simFiles,index_col=6)

subCounts = counts.loc[:,list(sim["barcodes"])]
subCounts.columns = list(sim.index)

s = timeit.default_timer()
adata_h5 = sc.AnnData(subCounts.T,dtype = np.float32)
adata_h5.var_names_make_unique()
coord_df = sim.loc[:, ['x', 'y']]
adata_h5.obsm["spatial"] = coord_df.to_numpy()
adata_h5.uns["spatial"] = coord_df.to_numpy()

adata_X = adata_preprocess(adata_h5, min_cells=5, pca_n_comps=params.cell_feat_dim)
graph_dict = graph_construction(adata_h5.obsm['spatial'], adata_h5.shape[0], params)


params.cell_num = adata_h5.shape[0]
print('==== Graph Construction Finished')



sedr_net = SEDR_Train(adata_X, graph_dict, params)
if params.using_dec:
    sedr_net.train_with_dec()
else:
    sedr_net.train_without_dec()
sedr_feat, _, _, _ = sedr_net.process()



adata_sedr = anndata.AnnData(sedr_feat)
adata_sedr.uns['spatial'] = adata_h5.uns['spatial']
adata_sedr.obsm['spatial'] = adata_h5.obsm['spatial']

sc.pp.neighbors(adata_sedr, n_neighbors=params.eval_graph_n)
sc.tl.umap(adata_sedr)


if re.search("dot",simFiles) is not None:
    n_clusters = 6
else:
    n_clusters = 3
eval_resolution = res_search_fixed_clus(adata_sedr, n_clusters)

sc.tl.leiden(adata_sedr, key_added="SEDR_leiden", resolution=eval_resolution)
tmp = adata_sedr.obs["SEDR_leiden"]
tmp = tmp.set_axis(list(sim.index))
export = pd.concat([sim,tmp], axis =1)
e = timeit.default_timer() - s
ftime = fileTag + "," + str(e) + "," + "sec" + "\n"
ft = open(time,'a')
ft.write(ftime)
ft.close()


filename = 'SEDR_' + fileTag + ".csv"
filename = os.path.join(output,filename)
export.to_csv(filename)
