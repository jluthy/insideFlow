# # InsideFlow R package for FlowJo!
# #
# library(flowCore)
# library(FlowSOM)
# library(data.table)
# # This is an that shows how to load a data matrix into an insideFlow object.
# # Start of user inputs
# runID <- "aloha_aina"
# fjXDim <- 10
# fjYDim <- 10
# EventNumberDP <- "EventNumberDP"
# CellIdDP <- "CellId"
# fjPopName <- "Ungated"
# fjDims <- 2
# fjSmooth <- 0 # good range -3 (sharper) to 3 (smoother)
# fjK <- 25 # good range 5 to 100
# fjAdjust <- 0 / 100.0 # good range 0 to 3
# fjSeed <- 1 # any
# fjUmap.Dist <- 0.5 # good range 0.1 to 2
# fjXdim <- 20 # it's sometimes useful to asymmetric SOM (Xdim != Ydim) between 5 to 64
# fjYdim <- 20 # it's sometimes useful to asymmetric SOM (Xdim != Ydim) between 5 to 64
# fjRLEN <- 15 # good range here is between 5 to 50
# fjMapSize <- 5000 # 1000 is good in almost all cases. Good range 50 to 5000
# fjSuperWeight <- 0.2 #values can range from 0.0 to 1.0
# fjEmbedSom.Map <- "FALSE"
# fjFlowSOM.Map <- "FALSE"
# fjUMAP.Map <- "FALSE"
# fjTSNE.Map <- "FALSE"
# fjSuper.Map <- "FALSE"
# fjRefMapped <- "FALSE"
# # path <- getwd()
# # path
# source("/Users/joshualuthy/Documents/code/insideFlow/R/insideFlow.R")
# # LOAD the data
# load(file="data/exMatrix.rda")
#
# # create an insideFlow object
# myFJobj <- new("insideFlow", runID = "runID")
#
# # Load the data into the new object
# myFJobj <- loadFlowJoDataTable(myFJobj, exMatrix)
#
# # Display tfhe parameter names from loaded data
# getParNames(myFJobj)
# # TODO do a better job of filtering only parameter names and stain names
#
# # convert csv to flowFrame
# myFJobj <- insideFlow::flowFrameFromFlowJoCSV(myFJobj)
#
# # run FlowSOM clustering
# # numClust <- 10
#
# EventNumberDP <- "EventNumberDP"
# CellIdDP <- "CellId"
# addCellIdsToResult <- FALSE
# runningInFlowJo <- !addCellIdsToResult
# myFJobj <- calcFlowSOM(myFJobj, applyOnMap = "None", 8, 8, 12)
#
# # plot MST
# plotBgd <- "Rainbow"
# radarColor <- "Rainbow"
# outs <- outPutFolder
# plotAllAsStars <- TRUE
# bgMetaClusters <- TRUE
# plotView <- "MST"
# plotViewOptions <- c("MST", "grid")
# plotPar <- "All as piecharts"
# outputFolder <- "/Users/joshualuthy/Desktop/outs/"
#
#
# fsMap <- plotFlowSOMResults(object=myFJobj, plotAllAsStars = TRUE, plotPar=plotPar, radarColor=radarColor, plotBgd=plotBgd, numClust=12)
# # png(filename = paste0(outputFolder, "fSOM12", ".png" ))
# # plot(fsMap, main=paste0("MST"))
# # dev.off() # This will cause the plot to be saved, then close the 'png device'
# ggplot2::ggsave(filename=paste0(outputFolder, "SOM_8x15", runID, "_", plotView, ".png"), plot = fsMap, width = 8, height = 15, units = "in", device = "png", bg="grey")
