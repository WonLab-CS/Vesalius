#-----------------------------/Running Giotto/---------------------------------#
#------------------------------------------------------------------------------#
# This file contains all code related to Giotto runs for the purpose of
# performance comaprison. We were not able to run Giotto on Slide-seq V2 this
# file only contains code related to Visium data.
# Giotto on slide seq V2 lead to segmentation faults
# NOTE The majority of this code is taken from
# https://github.com/JinmiaoChenLab/SEDR_analyses
#------------------------------------------------------------------------------#
library(Giotto)
library(Seurat)
library(ggplot2)
library(patchwork)
#------------------------------------------------------------------------------#
# SSV2
# Only running on two data sets
# DOES NOT WORK! Uses over 400 gb of ram and crashes !!!!!!
# No idea why 
#------------------------------------------------------------------------------#
slideTag <- c("Puck_200115_08","Puck_190926_03")

slideBeads <-c("~/group/slide_seqV2/Puck_200115_08_bead_locations.csv",
               "~/group/slide_seqV2/Puck_190926_03_bead_locations.csv")

slideCounts <- c("~/group/slide_seqV2/Puck_200115_08.digital_expression.txt.gz",
                 "~/group/slide_seqV2/Puck_190926_03.digital_expression.txt.gz")


input <- "/home/pcnmartin/group/slide_seqV2"

#slideBeads <-paste0(input,list.files(input, pattern ="locations.csv"))
#slideCounts <- paste0(input,list.files(input, pattern ="digital_expression.txt.gz"))

#slideTag <- sapply(strsplit(slideBeads,"V2/"),"[[",2)
#slideTag <- gsub("_bead_locations.csv","",slideTag)


output <- paste0(input,"/GiottoBenchMarking")
if(!dir.exists(output)){
    dir.create(output)
}

## removing all the intermeidate files
## set to out put dir
setwd(output)
time <-vector("list",length(slideBeads))


for(i in seq_along(slideTag)){

  counts <- read.table(slideCounts[i], header = TRUE )
  rownames(counts) <- counts[,1]
  counts <- counts[,-1]
  coord <- read.csv(slideBeads[i], header=T)
  rownames(coord) <- coord$barcodes
  s <- Sys.time()
  instruc <- createGiottoInstructions(save_plot=F, show_plot=F, save_dir=output)
  giotto <- createGiottoObject(raw_exprs = counts,
                                     spatial_locs = coord[,c("xcoord","ycoord")],
                                     instructions = instruc,
                                     cell_metadata = coord[,c("xcoord","ycoord")])

  giotto <- filterGiotto(gobject = giotto,
                         expression_threshold = 1,
                         expression_values = c('raw'))

  giotto <- normalizeGiotto(gobject = giotto, scalefactor = 6000)

  giotto <- addStatistics(gobject = giotto)



  giotto <- calculateHVG(gobject = giotto)

  gene_metadata <- fDataDT(giotto)
  featgenes <- gene_metadata[hvg == 'yes' & perc_cells > 3 & mean_expr_det > 0.4]$gene_ID

  giotto <- runPCA(gobject = giotto,
                   genes_to_use = featgenes,
                   scale_unit = F,
                   center=T,
                   method="irlba")


  giotto <- runUMAP(giotto, dimensions_to_use = 1:30)

  giotto <- createSpatialNetwork(gobject=giotto,
                                 method='kNN',
                                 k=15,
                                 maximum_distance_knn=500,
                                 name='spatial_network')
  spatial_genes <- silhouetteRank(giotto)
  ext_spatial_genes <-spatial_genes[1:2000,]$gene
  #----------------------------------------------------------------------------#
  # using this thing again to make it work - i havce no idea why this needs to
  # done again
  #----------------------------------------------------------------------------#
  gobject <- giotto
  expression_values  <- 'scaled'
  subset_genes <- ext_spatial_genes
  spatial_network_name  <- 'spatial_network'
  b <- 0


  select_spatialNetwork <- function(gobject,
                                    name = NULL,
                                    return_network_Obj = FALSE) {

    if (!is.element(name, names(gobject@spatial_network))){
      message = sprintf("spatial network %s has not been created. Returning NULL.
                        check which spatial networks exist with showNetworks() \n", name)
      warning(message)
      return(NULL)
    }else{
      networkObj = gobject@spatial_network[[name]]
      networkDT = networkObj$networkDT
    }

    if (return_network_Obj == TRUE){
      return(networkObj)
    }else{
      return(networkDT)
    }
  }


  spatial_network = select_spatialNetwork(gobject,name = spatial_network_name,return_network_Obj = FALSE)

  select_expression_values <- function(gobject, values) {

    if(values == 'scaled' & is.null(gobject@norm_scaled_expr)) {
      stop('run first scaling step')
    } else if(values == 'scaled') {
      expr_values = gobject@norm_scaled_expr
    } else if(values == 'normalized' & is.null(gobject@norm_expr)) {
      stop('run first normalization step')
    } else if(values == 'normalized') {
      expr_values = gobject@norm_expr
    } else if(values == 'custom' & is.null(gobject@custom_expr)) {
      stop('first add custom expression matrix')
    } else if(values == 'custom') {
      expr_values = gobject@custom_expr
    } else if(values == 'raw') {
      expr_values = gobject@raw_exprs
    }

    return(expr_values)

  }



  # get expression matrix
  values = match.arg(expression_values, c('normalized', 'scaled', 'custom'))
  expr_values = select_expression_values(gobject = gobject, values = values)

  if(!is.null(subset_genes)) {
    expr_values = expr_values[rownames(expr_values) %in% subset_genes,]
  }

  # data.table variables
  gene_ID = value = NULL

  # merge spatial network with expression data
  expr_values_dt = data.table::as.data.table(expr_values); expr_values_dt[, gene_ID := rownames(expr_values)]
  expr_values_dt_m = data.table::melt.data.table(expr_values_dt, id.vars = 'gene_ID', variable.name = 'cell_ID')

  convert_to_full_spatial_network =  function(reduced_spatial_network_DT) {

    # data.table variables
    distance = rank_int = NULL

    # find location coordinates
    coordinates = grep('sdim', colnames(reduced_spatial_network_DT), value = T)

    begin_coordinates = grep('begin', coordinates, value = T)
    new_begin_coordinates = gsub(x = begin_coordinates, pattern = '_begin', replacement = '')
    new_begin_coordinates = gsub(x = new_begin_coordinates, pattern = 'sdim', replacement = 'source_')

    end_coordinates = grep('end', coordinates, value = T)
    new_end_coordinates = gsub(x = end_coordinates, pattern = '_end', replacement = '')
    new_end_coordinates = gsub(x = new_end_coordinates, pattern = 'sdim', replacement = 'target_')

    # create normal source --> target
    part1 = data.table::copy(reduced_spatial_network_DT)
    part1 = part1[, c('from', 'to', begin_coordinates, end_coordinates, 'distance', 'weight'), with = F]
    colnames(part1) = c('source', 'target', new_begin_coordinates, new_end_coordinates, 'distance', 'weight')

    # revert order target (now source) --> source (now target)
    part2 = data.table::copy(reduced_spatial_network_DT[, c('to', 'from', end_coordinates, begin_coordinates, 'distance', 'weight'), with = F])
    colnames(part2) = c('source', 'target', new_begin_coordinates, new_end_coordinates, 'distance', 'weight')

    # combine and remove duplicates
    full_spatial_network_DT = rbind(part1, part2)
    full_spatial_network_DT = unique(full_spatial_network_DT)

    # create ranking of interactions
    data.table::setorder(full_spatial_network_DT, source, distance)
    full_spatial_network_DT[, rank_int := 1:.N, by = 'source']

    # create unified column
    full_spatial_network_DT = sort_combine_two_DT_columns(full_spatial_network_DT, 'source', 'target', 'rnk_src_trgt')

    return(full_spatial_network_DT)

  }


  sort_combine_two_DT_columns = function(DT,
                                         column1,
                                         column2,
                                         myname = 'unif_gene_gene') {

    # data.table variables
    values_1_num = values_2_num = scolumn_1 = scolumn_2 = unif_sort_column = NULL

    # maybe faster with converting to factors??

    # make sure columns are character
    selected_columns = c(column1, column2)
    DT[,(selected_columns):= lapply(.SD, as.character), .SDcols = selected_columns]

    # convert characters into numeric values
    uniq_values = sort(unique(c(DT[[column1]], DT[[column2]])))
    uniq_values_num = 1:length(uniq_values)
    names(uniq_values_num) = uniq_values


    DT[,values_1_num := uniq_values_num[get(column1)]]
    DT[,values_2_num := uniq_values_num[get(column2)]]


    DT[, scolumn_1 := ifelse(values_1_num < values_2_num, get(column1), get(column2))]
    DT[, scolumn_2 := ifelse(values_1_num < values_2_num, get(column2), get(column1))]

    DT[, unif_sort_column := paste0(scolumn_1,'--',scolumn_2)]
    DT[, c('values_1_num', 'values_2_num', 'scolumn_1', 'scolumn_2') := NULL]
    data.table::setnames(DT, 'unif_sort_column', myname)

    return(DT)
  }





  spatial_network = convert_to_full_spatial_network(spatial_network)


  spatial_network_ext = data.table::merge.data.table(spatial_network, expr_values_dt_m, by.x = 'target', by.y = 'cell_ID', allow.cartesian = T)


  spatial_network_ext_smooth = spatial_network_ext[, mean(value), by = c('source', 'gene_ID')]


  dt_to_matrix <- function(x) {
    rownames = as.character(x[[1]])
    mat = methods::as(as.matrix(x[,-1]), 'Matrix')
    rownames(mat) = rownames
    return(mat)
  }



  spatial_smooth_dc = data.table::dcast.data.table(data = spatial_network_ext_smooth, formula = gene_ID~source, value.var = 'V1')
  spatial_smooth_matrix = dt_to_matrix(spatial_smooth_dc)


  all_cells = colnames(expr_values)
  smoothed_cells = colnames(spatial_smooth_matrix)
  missing_cells = all_cells[!all_cells %in% smoothed_cells]



  metadata = pDataDT(giotto)
  subset_cell_IDs = subset(metadata, !(cell_ID %in% missing_cells))$cell_ID
  giotto = subsetGiotto(giotto, cell_ids = subset_cell_IDs)


  #####################################


  spat_cor_netw_DT=detectSpatialCorGenes(giotto,
                                         expression_values = 'scaled',
                                         method='network',
                                         spatial_network_name='spatial_network',
                                         subset_genes=ext_spatial_genes,
                                         network_smoothing=0)
  # cluster spatial genes
  spat_cor_netw_DT=clusterSpatialCorGenes(spat_cor_netw_DT,
                                          name='spat_netw_clus', k=15)
                 #heatmap_legend_param=list(title=NULL))


  sample_rate=2
  target=500
  tot=0
  num_cluster=15
  gene_list=list()
  clust=spat_cor_netw_DT$cor_clusters$spat_netw_clus
  for(i in seq(1, num_cluster)){
    gene_list[[i]]=colnames(t(clust[which(clust==i)]))
  }
  for(i in seq(1, num_cluster)){
    num_g=length(gene_list[[i]])
    tot=tot+num_g/(num_g^(1/sample_rate))
  }
  factor=target/tot
  num_sample=c()
  for(i in seq(1, num_cluster)){
    num_g=length(gene_list[[i]])
    num_sample[i]=round(num_g/(num_g^(1/sample_rate)) * factor)
  }
  set.seed(10)
  samples=list()
  union_genes=c()
  for(i in seq(1, num_cluster)){
    if(length(gene_list[[i]])<num_sample[i]){
      samples[[i]]=gene_list[[i]]
    }else{
      samples[[i]]=sample(gene_list[[i]], num_sample[i])
    }
    union_genes=union(union_genes, samples[[i]])
  }
  union_genes=unique(union_genes)

  my_spatial_genes <- union_genes

  HMRF_spatial_genes <- doHMRF(gobject=giotto,
                               expression_values='scaled',
                               spatial_genes=my_spatial_genes,
                               k=15,
                               spatial_network_name="spatial_network",
                               betas=c(0, 10, 5),
                               output_folder=output)

  giotto <- addHMRF(gobject=giotto,
                    HMRFoutput=HMRF_spatial_genes,
                    k=15,
                    betas_to_add=c(0, 10, 20, 30, 40),
                    hmrf_name='HMRF')
  giotto <- giotto@cell_metadata
  fileOut <- paste0(output,"/",slideTag[i],"_SSV2_BM.Rda")
  save(giotto,file =fileOut,sep =",",quote=F)
  time[[i]] <- Sys.time() -s
  # Removing all the stuff giotto outputs...
  rm(giotto);gc()
  frem <- list.files(pattern =".txt")
  for(i in frem){file.remove(i)}
  unlink("result.spatial.zscore", recursive =TRUE)
}

save(time, file = paste0(output,"/Giotto_SSV2_time.Rda"))






#------------------------------------------------------------------------------#
# DPFLC
#------------------------------------------------------------------------------#
# set.seed(1)
#
# input <- list.dirs("~/group/visium/DLPFC_globus", recursive =F)
# time <-vector("list",length(input))
# count <- 1
# n_cluster <-c(7,7,7,7,5,5,5,5,7,7,7,7)
# for(k in input){
#
#   s<- Sys.time()
#   sec <- Read10X(k)
#   spatial_locations <- data.table::fread(paste0(k,"/tissue_positions_list.csv"))
#   spatial_locations <- spatial_locations[match(colnames(sec), V1)]
#   colnames(spatial_locations) = c('barcode', 'in_tissue', 'array_row', 'array_col', 'col_pxl', 'row_pxl')
#
#   output <- paste0(k,"/Giotto")
#   if(!dir.exists(output)){
#       dir.create(output)
#   }
#
#   python_path <- NULL
#   if(is.null(python_path)) {
#     installGiottoEnvironment()
#   }
#
#
# myinst=createGiottoInstructions(save_plot=T, show_plot=F, save_dir=output)
#
# visium_brain <- createGiottoObject(raw_exprs = sec,
#                                    spatial_locs = spatial_locations[,.(row_pxl,-col_pxl)],
#                                    instructions = myinst,
#                                    cell_metadata = spatial_locations[,.(in_tissue, array_row, array_col)])
#
# metadata = pDataDT(visium_brain)
# in_tissue_barcodes = metadata[in_tissue == 1]$cell_ID
# visium_brain = subsetGiotto(visium_brain, cell_ids = in_tissue_barcodes)
#
# visium_brain <- filterGiotto(gobject = visium_brain,
#                              expression_threshold = 1,
#
#                              expression_values = c('raw'),
#                              verbose = T)
#
# visium_brain <- normalizeGiotto(gobject = visium_brain, scalefactor = 6000, verbose = T)
#
# visium_brain <- addStatistics(gobject = visium_brain)
#
#
#
# visium_brain <- calculateHVG(gobject = visium_brain)
#
# gene_metadata = fDataDT(visium_brain)
# featgenes = gene_metadata[hvg == 'yes' & perc_cells > 3 & mean_expr_det > 0.4]$gene_ID
#
# visium_brain <- runPCA(gobject = visium_brain,
#                        genes_to_use = featgenes,
#                        scale_unit = F,
#                        center=T,
#                        method="irlba")
#
#
# visium_brain <- runUMAP(visium_brain, dimensions_to_use = 1:20)
# visium_brain <- runtSNE(visium_brain, dimensions_to_use = 1:20)
#
#
#
# visium_brain <- createSpatialNetwork(gobject=visium_brain,
#                                      method='kNN',
#                                      k=5,
#                                      maximum_distance_knn=400,
#                                      name='spatial_network')
#
#
# spatial_genes <- silhouetteRank(visium_brain)
#
#
# ext_spatial_genes=spatial_genes[1:1500,]$gene
#
#
#
# ####################################
# ## Custom function from https://github.com/JinmiaoChenLab/SEDR_analyses
# ## Without these functions Giotto will not run to completion
# ####################################
# gobject <- visium_brain
# expression_values  <- 'scaled'
# subset_genes <- ext_spatial_genes
# spatial_network_name  <- 'spatial_network'
# b <- 0
#
#
# select_spatialNetwork <- function(gobject,
#                                   name = NULL,
#                                   return_network_Obj = FALSE) {
#
#   if (!is.element(name, names(gobject@spatial_network))){
#     message = sprintf("spatial network %s has not been created. Returning NULL.
#                       check which spatial networks exist with showNetworks() \n", name)
#     warning(message)
#     return(NULL)
#   }else{
#     networkObj = gobject@spatial_network[[name]]
#     networkDT = networkObj$networkDT
#   }
#
#   if (return_network_Obj == TRUE){
#     return(networkObj)
#   }else{
#     return(networkDT)
#   }
# }
#
#
# spatial_network = select_spatialNetwork(gobject,name = spatial_network_name,return_network_Obj = FALSE)
#
# select_expression_values <- function(gobject, values) {
#
#   if(values == 'scaled' & is.null(gobject@norm_scaled_expr)) {
#     stop('run first scaling step')
#   } else if(values == 'scaled') {
#     expr_values = gobject@norm_scaled_expr
#   } else if(values == 'normalized' & is.null(gobject@norm_expr)) {
#     stop('run first normalization step')
#   } else if(values == 'normalized') {
#     expr_values = gobject@norm_expr
#   } else if(values == 'custom' & is.null(gobject@custom_expr)) {
#     stop('first add custom expression matrix')
#   } else if(values == 'custom') {
#     expr_values = gobject@custom_expr
#   } else if(values == 'raw') {
#     expr_values = gobject@raw_exprs
#   }
#
#   return(expr_values)
#
# }
#
#
#
# # get expression matrix
# values = match.arg(expression_values, c('normalized', 'scaled', 'custom'))
# expr_values = select_expression_values(gobject = gobject, values = values)
#
# if(!is.null(subset_genes)) {
#   expr_values = expr_values[rownames(expr_values) %in% subset_genes,]
# }
#
# # data.table variables
# gene_ID = value = NULL
#
# # merge spatial network with expression data
# expr_values_dt = data.table::as.data.table(expr_values); expr_values_dt[, gene_ID := rownames(expr_values)]
# expr_values_dt_m = data.table::melt.data.table(expr_values_dt, id.vars = 'gene_ID', variable.name = 'cell_ID')
#
# convert_to_full_spatial_network =  function(reduced_spatial_network_DT) {
#
#   # data.table variables
#   distance = rank_int = NULL
#
#   # find location coordinates
#   coordinates = grep('sdim', colnames(reduced_spatial_network_DT), value = T)
#
#   begin_coordinates = grep('begin', coordinates, value = T)
#   new_begin_coordinates = gsub(x = begin_coordinates, pattern = '_begin', replacement = '')
#   new_begin_coordinates = gsub(x = new_begin_coordinates, pattern = 'sdim', replacement = 'source_')
#
#   end_coordinates = grep('end', coordinates, value = T)
#   new_end_coordinates = gsub(x = end_coordinates, pattern = '_end', replacement = '')
#   new_end_coordinates = gsub(x = new_end_coordinates, pattern = 'sdim', replacement = 'target_')
#
#   # create normal source --> target
#   part1 = data.table::copy(reduced_spatial_network_DT)
#   part1 = part1[, c('from', 'to', begin_coordinates, end_coordinates, 'distance', 'weight'), with = F]
#   colnames(part1) = c('source', 'target', new_begin_coordinates, new_end_coordinates, 'distance', 'weight')
#
#   # revert order target (now source) --> source (now target)
#   part2 = data.table::copy(reduced_spatial_network_DT[, c('to', 'from', end_coordinates, begin_coordinates, 'distance', 'weight'), with = F])
#   colnames(part2) = c('source', 'target', new_begin_coordinates, new_end_coordinates, 'distance', 'weight')
#
#   # combine and remove duplicates
#   full_spatial_network_DT = rbind(part1, part2)
#   full_spatial_network_DT = unique(full_spatial_network_DT)
#
#   # create ranking of interactions
#   data.table::setorder(full_spatial_network_DT, source, distance)
#   full_spatial_network_DT[, rank_int := 1:.N, by = 'source']
#
#   # create unified column
#   full_spatial_network_DT = sort_combine_two_DT_columns(full_spatial_network_DT, 'source', 'target', 'rnk_src_trgt')
#
#   return(full_spatial_network_DT)
#
# }
#
#
# sort_combine_two_DT_columns = function(DT,
#                                        column1,
#                                        column2,
#                                        myname = 'unif_gene_gene') {
#
#   # data.table variables
#   values_1_num = values_2_num = scolumn_1 = scolumn_2 = unif_sort_column = NULL
#
#   # maybe faster with converting to factors??
#
#   # make sure columns are character
#   selected_columns = c(column1, column2)
#   DT[,(selected_columns):= lapply(.SD, as.character), .SDcols = selected_columns]
#
#   # convert characters into numeric values
#   uniq_values = sort(unique(c(DT[[column1]], DT[[column2]])))
#   uniq_values_num = 1:length(uniq_values)
#   names(uniq_values_num) = uniq_values
#
#
#   DT[,values_1_num := uniq_values_num[get(column1)]]
#   DT[,values_2_num := uniq_values_num[get(column2)]]
#
#
#   DT[, scolumn_1 := ifelse(values_1_num < values_2_num, get(column1), get(column2))]
#   DT[, scolumn_2 := ifelse(values_1_num < values_2_num, get(column2), get(column1))]
#
#   DT[, unif_sort_column := paste0(scolumn_1,'--',scolumn_2)]
#   DT[, c('values_1_num', 'values_2_num', 'scolumn_1', 'scolumn_2') := NULL]
#   data.table::setnames(DT, 'unif_sort_column', myname)
#
#   return(DT)
# }
#
#
#
#
#
# spatial_network = convert_to_full_spatial_network(spatial_network)
#
#
# spatial_network_ext = data.table::merge.data.table(spatial_network, expr_values_dt_m, by.x = 'target', by.y = 'cell_ID', allow.cartesian = T)
#
#
# spatial_network_ext_smooth = spatial_network_ext[, mean(value), by = c('source', 'gene_ID')]
#
#
# dt_to_matrix <- function(x) {
#   rownames = as.character(x[[1]])
#   mat = methods::as(as.matrix(x[,-1]), 'Matrix')
#   rownames(mat) = rownames
#   return(mat)
# }
#
#
#
# spatial_smooth_dc = data.table::dcast.data.table(data = spatial_network_ext_smooth, formula = gene_ID~source, value.var = 'V1')
# spatial_smooth_matrix = dt_to_matrix(spatial_smooth_dc)
#
#
# all_cells = colnames(expr_values)
# smoothed_cells = colnames(spatial_smooth_matrix)
# missing_cells = all_cells[!all_cells %in% smoothed_cells]
#
#
#
# metadata = pDataDT(visium_brain)
# subset_cell_IDs = subset(metadata, !(cell_ID %in% missing_cells))$cell_ID
# visium_brain = subsetGiotto(visium_brain, cell_ids = subset_cell_IDs)
#
#
# #####################################
#
#
# spat_cor_netw_DT=detectSpatialCorGenes(visium_brain,
#                                        expression_values = 'scaled',
#                                        method='network',
#                                        spatial_network_name='spatial_network',
#                                        subset_genes=ext_spatial_genes,
#                                        network_smoothing=0)
# # cluster spatial genes
# spat_cor_netw_DT=clusterSpatialCorGenes(spat_cor_netw_DT,
#                                         name='spat_netw_clus', k=15)
#                heatmap_legend_param=list(title=NULL))
#
#
# sample_rate=2
# target=500
# tot=0
# num_cluster=n_cluster[count]
# gene_list=list()
# clust=spat_cor_netw_DT$cor_clusters$spat_netw_clus
# for(i in seq(1, num_cluster)){
#   gene_list[[i]]=colnames(t(clust[which(clust==i)]))
# }
# for(i in seq(1, num_cluster)){
#   num_g=length(gene_list[[i]])
#   tot=tot+num_g/(num_g^(1/sample_rate))
# }
# factor=target/tot
# num_sample=c()
# for(i in seq(1, num_cluster)){
#   num_g=length(gene_list[[i]])
#   num_sample[i]=round(num_g/(num_g^(1/sample_rate)) * factor)
# }
# set.seed(10)
# samples=list()
# union_genes=c()
# for(i in seq(1, num_cluster)){
#   if(length(gene_list[[i]])<num_sample[i]){
#     samples[[i]]=gene_list[[i]]
#   }else{
#     samples[[i]]=sample(gene_list[[i]], num_sample[i])
#   }
#   union_genes=union(union_genes, samples[[i]])
# }
# union_genes=unique(union_genes)
#
#
#
# my_spatial_genes <- union_genes
#
# HMRF_spatial_genes=doHMRF(gobject=visium_brain,
#                           expression_values='scaled',
#                           spatial_genes=my_spatial_genes,
#                           k=n_cluster[count],
#                           spatial_network_name="spatial_network",
#                           betas=c(0, 10, 5),
#                           output_folder=output)
#
#
# ### Visualize HMRF result
# visium_brain=addHMRF(gobject=visium_brain,
#                      HMRFoutput=HMRF_spatial_genes,
#                      k=n_cluster[count],
#                      betas_to_add=c(0, 10, 20, 30, 40),
#                      hmrf_name='HMRF')
#
#
#
# tag <- strsplit(k, "/")[[1]]
# tag <- tag[length(tag)]
# file <- paste0(k,"/Giotto/",tag,"_giotto.rda")
# save(visium_brain,file = file)
# e<- Sys.time()
# time[[count]] <- e - s
# count <- count +1
# }
# save(time,file = "~/group/visium/DLPFC_globus/Giotto_time.Rda")
