# 4_cosmx_spatial_cell_com.R
## This 

## s_interactive --mem=150GB --partition=caslake 
## srun --time=04:00:00 --mem=180G -p tier2q  --pty bash -I 
## conda activate /gpfs/data/biocore-workshop/spatial_transcriptomics_bruker_2026_workshop2/conda_envs/seurat_spatial
## R --vanilla

## Set correct library path. 
.libPaths("/gpfs/data/biocore-workshop/spatial_transcriptomics_bruker_2026_workshop2/conda_envs/seurat_spatial/lib/R/library")

print(.libPaths())

## Suppress progress bars and verbose cli output for clean SLURM logs
progressr::handlers("void")
options(progressr.enable = FALSE)
options(cli.default_handler = function(...) invisible(NULL))

## Increase future globals limit (default 500MB is too small for large datasets)
options(future.globals.maxSize = 2 * 1024^3)  # 2GB

## LIBRARIES
library(Seurat)
library(SpatialCellChat)
library(Matrix)
library(patchwork)

## PATHS
w_dir = "/gpfs/data/biocore-workshop/spatial_transcriptomics_bruker_2026_workshop2"
in_dir =  file.path(w_dir, 'CosMx/results/BreastCancer/Subset_MM')
out_dir = file.path(w_dir, 'CosMx/results/BreastCancer/SpatialCellChat')
dir.create(out_dir, showWarnings = F)

## Print starting time. 
start_time = Sys.time()
print(paste("Starting script at:", Sys.time()))

############
## Create seurat object with Subset data. 
############

# 1. Load counts and create base object
## Read mtx file directly 
counts = readMM(file.path(in_dir, "matrix.mtx"))
barcodes = read.table(file.path(in_dir, "barcodes.tsv"), header = FALSE)$V1
features = read.table(file.path(in_dir, "features.tsv"), header = FALSE, sep = "\t")

rownames(counts) = features$V2  # gene names (col 2), or V1 for IDs
colnames(counts) = barcodes

seurat_obj = CreateSeuratObject(counts = counts)

# 2. Add cell metadata (with annotations)
metadata = read.csv(file.path(in_dir, "metadata_annotation.csv"), row.names = 1)
seurat_obj = AddMetaData(seurat_obj, metadata)

## SUBSET - TESTING SCRIPT ONLY.
seurat_obj = subset(seurat_obj, CenterX_global_px > 35000  & CenterY_global_px < 98000)

saveRDS(seurat_obj, file = file.path(out_dir, "seurat_obj_subset_SpatialCellChat.rds"))

## 3. Normalize data
seurat_obj = NormalizeData(seurat_obj)

print(seurat_obj)

############
## Create CellChat object and identify over-expressed genes and interactions
############

print(paste("Starting CellChat analysis at:", Sys.time()))

## 1. Spatial coordinates. 
conversion.factor = 0.12028 ## Conversion factor for CosMx. 

## Spatial coordinates - CenterX_global_px and CenterY_global_px columns in metadata.
spatial_coords = dplyr::select(seurat_obj@meta.data, CenterX_global_px, CenterY_global_px) %>%
    rename(x = CenterX_global_px, y = CenterY_global_px) %>% as.matrix()

d = dist(spatial_coords, method = 'euclidean') ## Instead of computeCellDistance, function is not working well. 
## Calculate approximate spot size based on the minimum distance between spots.
spot.size = min(d)*conversion.factor # converting the distance in Pixels to Micrometers
spatial.factors = data.frame(ratio = conversion.factor, tol = spot.size/2)

## 2. Create Spatial CellChat obj. 
chat = createSpatialCellChat(
    object = seurat_obj, 
    group.by = "final_annotation", 
    assay = "RNA", 
    datatype = "spatial",
    coordinates = spatial_coords, 
    spatial.factors = spatial.factors
    )

## Set database. Here we use the built-in CellChatDB.human, (alt, mouse.)
CellChatDB = CellChatDB.human

## Set the database to use for this analysis: FILTERING STEP - For this case, only cell-cell contact. 
CellChatDB_subset = subsetDB(CellChatDB, search = list(c("Cell-Cell Contact")), key = c("annotation"))
chat@DB = CellChatDB_subset

## 3. Preprocess the expression data for cell-cell communication analysis.
chat = subsetData(chat) # This step is necessary even if using the whole database
chat = preProcessing(chat) # the function now requires only 1 argument, which is the object itself

print(paste("Pre-processing completed at:", Sys.time()))

## 4. Identify over-expressed genes and interactions.

# Identify over-expressed genes and interactions EXPENSIVE STEP.
chat = identifyOverExpressedGenes(
  chat,
  ## Methods: 'wilcox', 'meringue' & "moransi"
  selection.method = "meringue", # method for selecting (spatially) variable features
  do.grid = FALSE # if true, do "grid" operation to speed up computation
)

print(paste("Over-expressed genes identified at:", Sys.time()))

## Identify over-expressed interactions. Only require that either ligand or receptor from one pair is over-expressed.
chat = identifyOverExpressedInteractions(chat, variable.both = F)

## 5. Compute communication probability and infer cellular communication network.
## Play with the parameters to find best scale.distance value, or set to NULL. 
## Ranges and scale is base on the units of the spatial coordinates, which for CosMx is micrometers.
chat = computeCommunProb(
  chat,
  distance.use = TRUE, 
  scale.distance = 1,  ## 1-1 Since CosMx is already adjusted to micrometers. 
  #contact.dependent = TRUE,
  interaction.range = 250, ## Based on expected micrometers range - paracrine signaling. 
  contact.range = 10 ## Based on expected micrometers range - juxtacrine/contatc-dependent signaling.
)

print(paste("Communication probability computed at:", Sys.time()))

## Filter out to keep only significant interactions. 
chat = filterProbability(chat)

## Filter based on number of cells supporting the interaction.
chat = filterCommunication(
  chat, min.cells = NULL,
  min.links = 10,
  min.cells.sr = 10
)

print(paste("Filtered communication network at:", Sys.time()))

saveRDS(chat, file = file.path(out_dir, "SpatialCellChat.rds"))

## Store cell-cell Communication df. 
chat@LR$LRsig %>% write.csv(file.path(out_dir, "cell_cell_communication.csv"), row.names = F)


## Aggregate the cell-cell communication network.
chat = computeCommunProbPathway(chat)

saveRDS(chat, file = file.path(out_dir, "SpatialCellChat_communication.rds"))

quit()

## Add pathway information. 
chat = aggregateNet(chat)

print(paste("Aggregated communication network at:", Sys.time()))


## Print total time: 
print(paste("Total time for script:", Sys.time() - start_time))

saveRDS(chat, file = file.path(out_dir, "SpatialCellChat_expanded.rds"))
