ls()
rm(list=ls())
ls()
memory.size() ### Checking your memory size
memory.limit() ## Checking the set limit
memory.limit(size=56000) ### expanding your memory
memory.limit()
library(NMF)
library(ggalluvial)
library(SeuratDisk)
library(SeuratData)
library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)
library(CellChat)

## CellChat requires 2 user inputs: (1) normalized gene expression data, (2) cell labels ##
## Prepare counts csv file and metadata csv file that will be used to make Cell Chat object ##
## Can make both files after loading an H5 Seurat file or and R object ##
## After the files are made, read the csv files back in and make the Cell Chat object ##

## You Can Read In As An RDS File ##
Male_Seurat <- readRDS("MaleLow.Robj")
Male_Seurat

## Or, Read In As An H5 File ##
Male_Seurat <- LoadH5Seurat("MaleLow.h5Seurat")
Male_Seurat

## Get metadata from the Seurat object and store in variable called meta ##
meta <- Male_Seurat@meta.data
dim(meta) ## get dimentions of metadata table ##
head(meta) ## read first 6 lines of metadata table ##
meta2 <- as.data.frame(meta) ## convert table to a data frame ##
write.table(meta2, sep=",", file="MetadataCellChat.csv", row.names=TRUE, col.names=NA, quote=FALSE) ## export as a csv file ##

## Get The Normalized Counts from the Seurat object ##
counts <- Male_Seurat@assays$RNA@data
head(counts)
OligoCounts <- as.data.frame(counts)
write.table(OligoCounts, sep=",", file="CountsCellChat.csv", row.names=TRUE, col.names=NA, quote=FALSE)

## Read in the counts data stored as a csv file ##
## Make matrix of gene count information ##
data.input <-read.table(file.choose(), sep = ",", header = TRUE, row.names = 1) ## read in csv file ##
class(data.input) ## look at class of the object - should be a data frame ##
data.input <- as.matrix(data.input) ## convert it to a matrix ##
class(data.input) ## check class of the object - should be a "matrix" "array"
data.input <- as(data.input, "dgCMatrix") ## convert to dgCMatrix object ##
head(data.input) ## look at first 6 lines of counts matrix ##
class(data.input) ## confirm this is a matrix object ##
is.atomic(data.input) ## is it atomic - should be FALSE ##
dim(data.input) ## look at matrix dimensions ##
nrow(data.input) ## look at number of rows ##
ncol(data.input) ## look at number of columns ##
colnames(data.input) ## look at column names ##

## read in csv file with metadata and store in variable called meta ##
meta <- read.table(file.choose(), sep = ",", header = TRUE, row.names = 1) 
head(meta) ## look at first 6 lines of metadata table ##
rownames(meta) ## look at rownames of table ##
class(meta) ## look at object class - should be a data frame ##
dim(meta) ## look at table dimensions ##
nrow(meta) ## look rows in the table ##
all(rownames(meta)==colnames(data.input)) ## this must be TRUE to make CellChat object ##

             ### Prepare To Make CellChat Object ###
## extract the cell names for the desired group - e.g. control vs drug-treated ##
## we want to look at a group called Male_Low ##
cell.use = rownames(meta)[meta$type == "Male_Low"] 

## Prepare input data for CelChat analysis ##
data.input = data.input[, cell.use] ## create data frame with desired counts info ##
head(data.input) ## check first 6 lines of counts info ##
meta = meta[cell.use, ] ## create data frame with desired meta info ##
head(meta) ## check first 6 lines of metadata info ##
unique(meta$monaco2) # check the cell labels

## Create CellChat Object ##
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "monaco2")
cellchat

## Save CellChat Object for use at a later time ##
saveRDS(cellchat,"MaleLow.Robj")

## Read In Cell Chat Object ##
cellchat <- readRDS("MaleLow.Robj")
cellchat ## check object

## Choose CellChat Ligand Receptor Database (this is for mouse) ##
CellChatDB <- CellChatDB.mouse
## Get a visual (pie chart) of the CellChat mouse database ##
showDatabaseCategory(CellChatDB)
## Show the structure of the database ##
dplyr::glimpse(CellChatDB$interaction)

## Optional: use a subset of CellChatDB for cell-cell communication analysis ##
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")

## Preferred: use entire CellChatDB for cell-cell communication analysis ##
CellChatDB.use <- CellChatDB

## Now, set the database to be used in the object ##
cellchat@DB <- CellChatDB.use

## Pre-Processing Steps - Three steps that must be run ##
# subset the expression data of signaling genes for saving computation cost
cellchat <- subsetData(cellchat) ## Step 1 is necessary even if using the whole database ##
cellchat <- identifyOverExpressedGenes(cellchat) ## Step 2 ##
cellchat <- identifyOverExpressedInteractions(cellchat) ## Step 3 ##

## Compute the communication probability and infer the cellular communication network ##
cellchat <- computeCommunProb(cellchat) ## this step will take a few minutes ##

# Filter out the cell-cell communication if there are only few number of cells in certain cell groups ##
cellchat <- filterCommunication(cellchat, min.cells = 10) ## this is optional ##

## Now access the inferred cell-to-cell communications of interest ##
## Set slot.name = "netP" to access the the inferred communications at the level of signaling pathways ##
## Other options: 
## df.net <- subsetCommunication(cellchat, sources.use = c(1,2), targets.use = c(4,5)) gives the inferred ##
## cell-cell communications sending from cell groups 1 and 2 to cell groups 4 and 5. ##
## df.net <- subsetCommunication(cellchat, signaling = c("WNT", "TGFb")) gives the inferred cell-cell communications ##
## mediated by signaling WNT and TGFb. ##
df.net <- subsetCommunication(cellchat, slot.name = "netP") ## go with "netP" ##

## export info as a csv file ##
write.table(df.net, sep=",", file="MaleHighCellChat.csv", row.names=TRUE, col.names=NA, quote=FALSE) 

## Communication Probabilities Of A Signaling Pathway Level ##
cellchat <- computeCommunProbPathway(cellchat)

## Calculate Aggregated Cell to Cell Communication Network ##
cellchat <- aggregateNet(cellchat)
## Visualize Aggregated Cell To Cell Communication Network ##
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

## Examine Signaling Sent From Each Group ##
mat <- cellchat@net$weight
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}

## Take a look at Signaling Pathways Showing Significant Communications ##
cellchat@netP$pathways

dev.off() ## reset par to default

## Visualize each signaling pathway using hierarchy plot, chord diagram or circle plot ##
pathways.show <- c("NRXN") 
# Hierarchy plot
# Here we define `vertex.receiver` so we look at signaling between specific cell types ## 
vertex.receiver = seq(1,4) # a numeric vector. 
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver)
# Circle plot
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
# Chord diagram
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")
# Heatmap
par(mfrow=c(1,1))
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")

#Compute the contribution of each ligand-receptor pair to the overall signaling pathway and ##
## visualize cell-cell communication mediated by a single ligand-receptor pair ##
netAnalysis_contribution(cellchat, signaling = pathways.show)

## Visualize the cell-cell communication mediated by a single ligand-receptor pair. ##
## Use a function extractEnrichedLR to extract all the significant interactions (L-R pairs) ##
## and related signaling genes for a given signaling pathway ##
pairLR.NRXN <- extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.NRXN[1,] # show one ligand-receptor pair
## Hierarchy plot ##
vertex.receiver = seq(1,4) # a numeric vector
netVisual_individual(cellchat, signaling = pathways.show,  pairLR.use = LR.show, vertex.receiver = vertex.receiver)
# Circle plot
netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")
# Chord diagram
netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "chord")

## Automatically save the plots of the all inferred network for quick exploration ##
## Access all the signaling pathways showing significant communications ##
pathways.show.all <- cellchat@netP$pathways
## check the order of cell identity to set suitable vertex.receiver ##
levels(cellchat@idents)
vertex.receiver = seq(1,4)
for (i in 1:length(pathways.show.all)) {
  ## Visualize communication network associated with both signaling pathway and individual L-R pairs ##
  netVisual(cellchat, signaling = pathways.show.all[i], vertex.receiver = vertex.receiver, layout = "hierarchy")
  ## Compute and visualize the contribution of each ligand-receptor pair to the overall signaling pathway ##
  gg <- netAnalysis_contribution(cellchat, signaling = pathways.show.all[i])
  ggsave(filename=paste0(pathways.show.all[i], "_L-R_contribution.pdf"), plot=gg, width = 3, height = 2, units = 'in', dpi = 300)
}

## Visualize cell-cell communication mediated by multiple ligand-receptors or signaling pathways ##
## Bubble Plot ##
## show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use') ##
## to other cell groups (defined by 'targets.use') ##
netVisual_bubble(cellchat, sources.use = c(1:12), targets.use = c(6,9), remove.isolate = FALSE)

## show all the significant interactions (L-R pairs) associated with certain signaling pathways ##
netVisual_bubble(cellchat, sources.use = c(1:12), targets.use = c(6,9), signaling = c("VTN", "SEMA6", "LAMININ", "CNTN",
                                                                                      "BMP", "WNT"), remove.isolate = FALSE)
## show all the significant interactions (L-R pairs) associated with certain signaling pathways ##
netVisual_bubble(cellchat, sources.use = c(1:12), targets.use = c(8), remove.isolate = FALSE)
## show all the significant interactions (L-R pairs) associated with certain signaling pathways##
netVisual_bubble(cellchat, sources.use = c(1:12), targets.use = c(8), signaling = c("IGF", "NCAM", "CX3C",
                                                                                    "GAS", "TGFb", "SPP1",
                                                                                    "EPHA", "JAM", "CSF", "NRXN"),
                 remove.isolate = FALSE)
## Too Microglia - show all the significant interactions (L-R pairs) associated with certain signaling pathways ##
netVisual_bubble(cellchat, sources.use = c(1:12), return.data = TRUE, targets.use = c(8), signaling = c("NCAM", "CX3C",
                                                                                                        "GAS", "TGFb", "SPP1",
                                                                                                        "EPHA", "NRXN"),
                 remove.isolate = FALSE)

## From Microglia - show all the significant interactions (L-R pairs) associated with certain signaling pathways ##
netVisual_bubble(cellchat, sources.use = c(8), targets.use = c(1:7, 9:12), 
                 remove.isolate = FALSE)
## From Microglia - show all the significant interactions (L-R pairs) associated with certain signaling pathways ##
netVisual_bubble(cellchat, sources.use = c(8), targets.use = c(1:7, 9:12), 
                 remove.isolate = FALSE, signaling = c("NRXN", "NRG", "NCAM", "EPHA", "LAMININ", "EPHB",
                                                       "COLLAGEN", "PSAP"))

## show all the significant interactions (L-R pairs) based on user's input (defined by `pairLR.use`) ##
pairLR.use <- extractEnrichedLR(cellchat, signaling = c("CCL","CXCL","FGF"))
netVisual_bubble(cellchat, sources.use = c(3,4), targets.use = c(5:8), pairLR.use = pairLR.use, remove.isolate = TRUE)

## Plot the signaling gene expression distribution using violin/dot plot ##
plotGeneExpression(cellchat, signaling = "TGFb", type = "violin")

## Compute the network centrality scores ##
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular ##
## communication network of signaling pathways ##
## Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups ##
netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)

## Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways ##
## What signals contribute most to outgoing and incoming signaling ##
ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")
ht1 + ht2

##  Signaling role analysis on the cell-cell communication networks of interest ##
ht <- netAnalysis_signalingRole_heatmap(cellchat, signaling = c("PDGF", "SEMA3", "SEMA4", "SEMA6", "MPZ"), pattern = "incoming")
ht

###  DO ThIS  ###
## Visualize dominant senders and receivers in a 2D space ##
## Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways#
gg1 <- netAnalysis_signalingRole_scatter(cellchat)
## Signaling role analysis on the cell-cell communication networks of interest ##
gg2 <- netAnalysis_signalingRole_scatter(cellchat, signaling = c("PDGF"))
gg1 + gg2

## save CellChat object ##
saveRDS(cellchat, file = "cellchat_MaleLow.rds")
## save table with pathway info ##
rf.net <- pathways.show.all <- cellchat@netP
rf.net
write.table(rf.net, sep=",", file="FemaleHighCellChat.csv", row.names=TRUE, col.names=NA, quote=FALSE)

## look at communication patterns of secreting cells ##
selectK(cellchat, pattern = "outgoing") ## will give you a graph ##
nPatterns = 3
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns) ## gives a heatmap ##

## river plot ##
netAnalysis_river(cellchat, pattern = "outgoing", sources.use = c(8))

## dot plot ##
netAnalysis_dot(cellchat, pattern = "outgoing")

## incoming communication patterns of target cells ##
## Incoming patterns show how the target cells (i.e. cells as signal receivers) coordinate with each other ##
## as well as how they coordinate with certain signaling pathways to respond to incoming signals. ##
selectK(cellchat, pattern = "incoming")
nPatterns = 4
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming", k = nPatterns)
## river plot ##
netAnalysis_river(cellchat, pattern = "incoming")

## dot plot ##
netAnalysis_dot(cellchat, pattern = "incoming")


