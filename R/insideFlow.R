#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Class definition
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' @title The InsideFlow Class
#' @description
#' The insideFlow object is a representation of a gated population of cells that
#' are exported from a FlowJo plugin process. These files loaded into the object
#' can be either .fcs files or .csv files containing expression data from a flow
#' cytometery experiment. The object stores metadata as well as representations
#' of low-dimensional data from clustering or dimensionality reductions methods.
#'
#' Each insideFlow object has a number of slots which store information. Key
#' slots to access are listed below.
#'
#' @slot runID The unique id for plugin process
#' @slot inputCSV The expression matrix exported from FlowJo
#' @slot parameterNames The metadata for parameter and stain names
#' @slot flowFrameFJ The expression matrix stored as flowFrame
#' @slot BatchCorrect The expression matrix both uncorrected and batch corrected
#' @slot PeacoQCPops The results from performing QC on flow data with PeacoQC
#' @slot FlowAIPops The results from performing QC on flow data with FlowAI
#' @slot Rphenograph The results from clustering with Phenograph
#' @slot FastPhenograph The results from the FastPG algorithm
#' @slot FlowSOMClusters The results from FlowSOM clustering
#' @slot DimReduction The low-dimensional embedding result and map
#' @slot ClustRCheck The results from cluster evaluation with ClustRCheck
#' @slot Taylor The results from cluster evaluation via Taylor Index
#' @docType class
#' @name insideFlow-class
#'
insideFlow <- setClass("insideFlow",
         slots = c(
           runID = "character",
           inputCSV = "list",
           parameterNames = "list",
           flowFrameFJ = "list",
           BatchCorrect = "list",
           PeacoQCPops = "list",
           FlowAIPops = "list",
           Rphenograph = "list",
           FastPhenograph = "list",
           FlowSOMClusters = "list",
           DimReduction = "list",
           ClustRCheck = "list",
           Taylor = "list"
         ),
         prototype = list(
           runID = "runID",
           inputCSV = list(X=matrix(),
                           cellIds=matrix(),
                           eventNumbers=matrix()
           ),
           parameterNames = list(fjParNames=NULL,
                                 fjStainNames=NULL

           ),
           flowFrameFJ = list(ffFJ=list()
           ),
           BatchCorrect = list(uncorrected=matrix(),
                               corrected.asinh.transformed=matrix(),
                               corrected.reverse.transformed=matrix(),
                               markersToNorm="character",
                               labels=matrix(),
                               id=matrix(),
                               colNames=matrix()

           ),
           PeacoQCPops = list(results=matrix()
           ),
           FlowAIPops = list(goodEvents=matrix(),
                             badEvents=matrix()
           ),
           Rphenograph = list(clusters=matrix(),
                              clustNumber=matrix(),
                              modularity=matrix()
           ),
           FastPhenograph = list(cluster = matrix(),
                                 clustNumber=matrix(),
                                 modularity=matrix()
           ),
           FlowSOMClusters = list(clusters=matrix(),
                                  clustNumber=matrix(),
                                  som=list()
           ),
           DimReduction = list(embeddingMat=matrix(),
                                    embeddingType="character",
                                    map="ANY"
           ),
           ClustRCheck = list(values=matrix()
           ),
           Taylor = list(values=matrix()
           )
         ),
)
#######################################################################
##         BASIC SLOT QUERIES     ##
#######################################################################
#' @title View the unique run identifier
#' @return Returns the runID generated from the plugin process.
#' @param object A insideFlow object
#' @rdname getRunID
#' @examples \dontrun{
#' getRunID(object)}
#' @docType methods
#' @export
setGeneric("getRunID", function(object) {
  standardGeneric("getRunID")
})
#' @rdname getRunID
setMethod("getRunID",
          "insideFlow",
          function(object){
            # cat("RunID:", object@runID, "\n")
            runid <- object@runID
            return(runid)
          }
)
#'
#' @title View parameter names
#' @return Returns the parameter names from the exported population
#' @param object A insideFlow object
#' @rdname getParNames
#' @examples \dontrun{
#' getParNames(object)}
#' @docType methods
#' @export
setGeneric("getParNames", function(object) {
  standardGeneric("getParNames")
})
#' @rdname getParNames
setMethod("getParNames",
          "insideFlow",
          function(object){
            if(is.null(object@parameterNames$fjParNames)) {
              print("Parameter Names: None detected")
              return(NULL)
            } else {
              parNames <- object@parameterNames$fjParNames
              # cat("Parameter Names: ", parNames, "\n", "\n")
              return(parNames)
            }
          })
#'
#' @title View stain names
#' @return Returns the stain names if present from the exported population
#' @param object A insideFlow object
#' @rdname getStainNames
#' @examples \dontrun{
#' getStainNames(object)}
#' @docType methods
#' @export
setGeneric("getStainNames", function(object) {
  standardGeneric("getStainNames")
})
#' @rdname getStainNames
setMethod("getStainNames",
          "insideFlow",
          function(object){
            if(is.null(object@parameterNames$fjStainNames)) {
              print("Stain Names: None detected")
              return(NULL)
            } else {
              stNames <- object@parameterNames$fjStainNames
              # cat("Stain Names: ", stNames, "\n", "\n")
              return(stNames)
            }
          })
#######################################################################
##         Load a file from single csv file     ##
#######################################################################
#' @title Load data into insideFlow object
#' @return Returns the newly created object
#' @param object A insideFlow object
#' @param csvPath The path to the exported csv file from FlowJo
#' @importFrom data.table fread
#' @examples \dontrun{
#' inputCSVpath <- "/Documents/ExtNode.csv"
#' myobject <- loadFlowJoCSV(object, inputCSVpath)}
#' @docType methods
#' @export
setGeneric("loadFlowJoCSV", function(object, csvPath) {
  standardGeneric("loadFlowJoCSV")
})
#' @rdname loadFlowJoCSV
setMethod("loadFlowJoCSV",
          "insideFlow",
          function(object, csvPath){
            data <- fread(csvPath, check.names=FALSE)
            cellIdsColumn <- rep(NA, nrow(data))
            cellIds <- data.frame(CellIdDP = cellIdsColumn)
            names(cellIds) <- CellIdDP
            if (CellIdDP %in% names(data)) {
              cellIds <- data[, CellIdDP]
              data <- data[, -which(names(data) %in% c(CellIdDP))]
            }

            EventNumberColumn <- rep(NA, nrow(data))
            eventNumbers <- data.frame(EventNumberDP = EventNumberColumn)
            names(eventNumbers) <- EventNumberDP
            if (EventNumberDP %in% names(data)) {
              eventNumbers <- data[, EventNumberDP]
              data <- data[, -which(names(data) %in% c(EventNumberDP))]
            }

            ## Again splitting by " :: " but this time taking the second part (if present) as description
            colNames <- colnames(data)
            npar <- length(colNames)
            parNames <- unlist(lapply(1:npar, function(i) {
                       unlist(strsplit(colNames[i], " :: "))[1]
                   }))
            stNames <- unlist(lapply(1:npar, function(i) {
              s <- unlist(strsplit(colNames[i], " :: "))
              if (length(s) >= 2) s[2]
              else NA
            }))

            # Trim the leading and tailing whitespaces
            colnames(data) <- gsub("^\\s+|\\s+$", "", colnames(data))
            object@parameterNames <- list(fjParNames=parNames,
                                          fjStainNames=stNames)

            object@inputCSV <- list(X=as.matrix(data),cellIds=as.matrix(cellIds),eventNumbers=as.matrix(eventNumbers))
            colnames(object@inputCSV$cellIds) <- "CellId"
            colnames(object@inputCSV$eventNumbers) <- "EventNumberDP"
            return(object)
          }

)
#######################################################################
##         Load a data frame or Table     ##
#######################################################################
#' @title Load a data.table into insideFlow object
#' @return Returns the newly created object
#' @param object A insideFlow object
#' @param data The data object to load into insideFlow object
#' @importFrom data.table fread
#' @examples \dontrun{
#' exMatrix <- fread("/Documents/ExtNode.csv")
#' myobject <- loadFlowJoDataTable(object, exMatrix)}
#' @docType methods
#' @export
setGeneric("loadFlowJoDataTable", function(object, data) {
  standardGeneric("loadFlowJoDataTable")
})
#' @rdname loadFlowJoDataTable
setMethod("loadFlowJoDataTable",
          "insideFlow",
          function(object, data){
            cellIdsColumn <- rep(NA, nrow(data))
            cellIds <- data.frame(CellIdDP = cellIdsColumn)
            names(cellIds) <- CellIdDP
            if (CellIdDP %in% names(data)) {
              cellIds <- data[, CellIdDP]
              data <- data[, -which(names(data) %in% c(CellIdDP))]
            }

            EventNumberColumn <- rep(NA, nrow(data))
            eventNumbers <- data.frame(EventNumberDP = EventNumberColumn)
            names(eventNumbers) <- EventNumberDP
            if (EventNumberDP %in% names(data)) {
              eventNumbers <- data[, EventNumberDP]
              data <- data[, -which(names(data) %in% c(EventNumberDP))]
            }

            ## Again splitting by " :: " but this time taking the second part (if present) as description
            colNames <- colnames(data)
            npar <- length(colNames)
            parNames <- unlist(lapply(1:npar, function(i) {
              unlist(strsplit(colNames[i], " :: "))[1]
            }))
            stNames <- unlist(lapply(1:npar, function(i) {
              s <- unlist(strsplit(colNames[i], " :: "))
              if (length(s) >= 2) s[2]
              else NA
            }))

            # Trim the leading and tailing whitespaces
            colnames(data) <- gsub("^\\s+|\\s+$", "", colnames(data))
            object@parameterNames <- list(fjParNames=parNames,
                                          fjStainNames=stNames)

            object@inputCSV <- list(X=as.matrix(data),cellIds=as.matrix(cellIds),eventNumbers=as.matrix(eventNumbers))
            colnames(object@inputCSV$cellIds) <- "CellId"
            colnames(object@inputCSV$eventNumbers) <- "EventNumberDP"
            return(object)
          }

)
################################################################
#    Prepare for batch corrections cyCombine #
################################################################
#' @title Transform data before batch correction step
#' @return Returns the object with asinh-transformed values to batch correct
#' @param object A insideFlow object
#' @param fjMarkersList The selected list of markers to normalize
#' @param fjCofactor The selected Cofactor value to apply
#' @importFrom cyCombine transform_asinh
#' @examples \dontrun{
#' myobject <- prepareBEC(object, markers, cofact)}
#' @docType methods
#' @export
setGeneric("prepareBEC", function(object, fjMarkersList, fjCofactor){
  standardGeneric("prepareBEC")
})
#' @rdname prepareBEC
setMethod("prepareBEC",
          "insideFlow",
          function(object, fjMarkersList, fjCofactor){
            object@BatchCorrect$markersToNorm <- fjMarkersList

            object@BatchCorrect$uncorrected <- transform_asinh(
              df = as.data.frame(object@inputCSV$X),
              markers = object@BatchCorrect$markersToNorm,
              cofactor = fjCofactor,
              .keep = TRUE # Lets you keep all columns, in case they are useful to you
            )
            return(object)
          }
)
################################################################
#    Perform Batch Correction cyCombine #
################################################################
#' @title Performs batch correction with cyCombine
#' @return Returns the batch corrected object
#' @param object A insideFlow object
#' @param batchID The selected list of markers to normalize
#' @param normMethod The selected Cofactor value to apply
#' @param xdim The first dimension of SOM used
#' @param ydim The second dimension of SOM used
#' @importFrom cyCombine normalize
#' @importFrom cyCombine create_som
#' @examples \dontrun{
#' myobject <- runBatch_correct(object, bID, nrmMethod, xdim, ydim)}
#' @docType methods
#' @export
setGeneric("runBatch_correct", function(object, batchID, normMethod, xdim, ydim){
  standardGeneric("runBatch_correct")
})
#' @rdname runBatch_correct
setMethod("runBatch_correct",
          "insideFlow",
          function(object, batchID, normMethod, xdim, ydim){
            object@BatchCorrect$uncorrected$batch <- as.factor(object@BatchCorrect$uncorrected[,batchID])
            # Run batch correction
            object@BatchCorrect$labels <- object@BatchCorrect$uncorrected %>%
              normalize(markers = object@BatchCorrect$markersToNorm,
                        norm_method = normMethod) %>%
              create_som(markers = object@BatchCorrect$markersToNorm,
                         rlen = 10,
                         xdim = xdim,
                         ydim = ydim)

            # 'covar' is the optional second keyword present in the exported data.
            object@BatchCorrect$corrected.asinh.transformed <- object@BatchCorrect$uncorrected %>%
              correct_data(label = object@BatchCorrect$labels,
                           markers = object@BatchCorrect$markersToNorm)#,
            # covar = "status")

            return(object)
          }
)
################################################################
##    Returns Plots Post Batch Correction With cyCombine     ##
################################################################
#' @title Return plots comparing markers after batch normalization
#' @return Returns the plots from cyCombine showing normalized histograms
#' @param object A insideFlow object
#' @param fjPopName The name of selected population in FlowJo
#' @param fjCofactor The selected Cofactor value to apply
#' @param outPutFolder The path to the output folder
#' @importFrom cyCombine plot_density
#' @examples \dontrun{
#' normalizedPlots(object, popName, coFact, outsFolder)}
#' @docType methods
#' @export
setGeneric("normalizedPlots", function(object, fjPopName, fjCofactor, outPutFolder){
  standardGeneric("normalizedPlots")
})
#' @rdname normalizedPlots
setMethod("normalizedPlots",
          "insideFlow",
          function(object, fjPopName, fjCofactor, outPutFolder){
            if(is.null(myFJobj@BatchCorrect$uncorrected$label)){
              print("Adding cluster labels")
              myFJobj@BatchCorrect$uncorrected$label <- myFJobj@BatchCorrect$corrected.asinh.transformed$label
            }else{print("Cluster labels exist")}
            if(is.null(myFJobj@BatchCorrect$uncorrected$id)){
              print("Adding cell ids")
              myFJobj@BatchCorrect$uncorrected$id <- myFJobj@BatchCorrect$corrected.asinh.transformed$id
            }else{print("Cell ids exist")}

            # Plots can be made like so:
            # png("NUL")
            # pdf("NUL")
            plot_density(uncorrected = myFJobj@BatchCorrect$uncorrected,
                         corrected = myFJobj@BatchCorrect$corrected.asinh.transformed,
                         markers = myFJobj@BatchCorrect$markersToNorm,
                         filename = paste0(outPutFolder, "Density_plot_", fjPopName, "_cofactor_" , fjCofactor,"_.png" ))

            return(object)
          })
#######################################################################
##         Reverse transform before exporting FCS post cyCombine     ##
#######################################################################
#' @title Reverse transforms data matix before exporting to FlowJo
#' @return Returns reversed transformed data matrix before export to FlowJo
#' @param object A insideFlow object
#' @param fjCofactor The selected Cofactor value to apply
#' @importFrom cyCombine transform_asinh
#' @examples \dontrun{
#' myobject <- reverseTransform(object, popName, coFact, outsFolder)}
#' @docType methods
#' @export
setGeneric("reverseTransform", function(object, fjCofactor){
  standardGeneric("reverseTransform")
})
#' @rdname reverseTransform
setMethod("reverseTransform",
          "insideFlow",
          function(object, fjCofactor){
            # Transform back to scaled values before export
            object@inputCSV$X <- transform_asinh(
              df = object@BatchCorrect$corrected.asinh.transformed,
              markers = object@BatchCorrect$markersToNorm,
              cofactor = fjCofactor,
              reverse = TRUE,
              .keep = TRUE # Lets you keep all columns, in case they are useful to you
            )
            return(object)
          }
)
#######################################################################
##         Cleanup Object before exporting FCS post cyCombine        ##
#######################################################################
#' @title Cleanup Object
#' @return Returns the light-weight object cleaned of transformed data slots
#' from batch correction processes
#' @param object A insideFlow object
#' @examples \dontrun{
#' myobject <- cleanupFJobj(object)}
#' @docType methods
#' @export
setGeneric("cleanupFJobj", function(object){
  standardGeneric("cleanupFJobj")
})
#' @rdname cleanupFJobj
setMethod("cleanupFJobj",
          "insideFlow",
          function(object){
            if(!is.null(object@inputCSV$X)){
              # remove added cols from cyCombine
              if("id" %in% colnames(object@inputCSV$X)){
                object@inputCSV$X[["id"]] <- NULL
              }
              if("label" %in% colnames(object@inputCSV$X)){
                object@inputCSV$X[["label"]] <- NULL
              }
              if("batch" %in% colnames(object@inputCSV$X)){
                object@inputCSV$X[["batch"]] <- NULL
              }
              # get the original ordering before exporting
              paramOrderIdx <- match(object@parameterNames$fjParNames, colnames(object@inputCSV$X))
              # Replace original expression matrix with batch corrected one
              object@inputCSV$X <- object@inputCSV$X[ ,paramOrderIdx]

              # Clean up to make the object light
              object@BatchCorrect$corrected.asinh.transformed <- as.matrix(0)
              object@BatchCorrect$corrected.reverse.transformed <- as.matrix(0)
              print("cleaned FlowJo R object")
            }
            return(object)
          }
)
################################################################
#    Build a FlowFrame  #
################################################################
#' @title Build a Flowframe from CSV file expression matrix
#' @return Returns reversed transformed data matrix before export to FlowJo
#' @param object A insideFlow object
#' @param dropCompPrefix Boolean to allow for paramName cleanup
#' @importFrom methods new
#' @examples \dontrun{
#' myobject <- flowFrameFromFlowJoCSV(object, popName, coFact, outsFolder)}
#' @docType methods
#' @export
setGeneric("flowFrameFromFlowJoCSV", function(object, dropCompPrefix=FALSE){
  standardGeneric("flowFrameFromFlowJoCSV")
})
#' @rdname flowFrameFromFlowJoCSV
setMethod("flowFrameFromFlowJoCSV",
          "insideFlow",
          function(object, dropCompPrefix=FALSE){
            data <- object@inputCSV$X
            colNames <- colnames(data)
            npar <- length(colNames)
            id <- paste("$P",1:npar,sep="")

            maxRange <- unlist(lapply(1:npar, function(x) {
              max(data[,x])
            }))
            range <- maxRange + 1
            minRange <- unlist(lapply(1:npar, function(x) {
              min(0, min(data[,x]))
            }))

            colNames <- gsub("^\\s+|\\s+$", "", colNames) 	# Trim the leading and tailing whitespaces

            ## FlowJo splits parameter names and description by " :: ", e.g., "Comp-PE-Cy5-A :: CD19#PE-Cy5"
            names <- unlist(lapply(1:npar, function(i) {
              unlist(strsplit(colNames[i], " :: "))[1]
            }))
            ## FlowJo indicates that a patameter is compensated by the Comp- prefix in its name, we can drop that if undesirable
            if (dropCompPrefix) {
              names <- gsub("^\\Comp-", "", names)
            }
            ## Again splitting by " :: " but this time taking the second part (if present) as description
            desc <- unlist(lapply(1:npar, function(i) {
              s <- unlist(strsplit(colNames[i], " :: "))
              if (length(s) >= 2) s[2]
              else NA
            }))

            pars <- new("AnnotatedDataFrame",
                        data=data.frame(
                          row.names=I(id), name=I(names), desc=I(desc),
                          range=range, minRange=minRange, maxRange=maxRange),
                        varMetadata=data.frame(
                          row.names=I(c("name","desc","range", "minRange", "maxRange")),
                          labelDescription=I(c("Name of Parameter", "Description of Parameter", "Range of Parameter", "Minimum Parameter Value after Transforamtion", "Maximum Parameter Value after Transformation"))
                        )
            )

            ## Build the descripton
            txt <- list()
            txt["FCSversion"] <- "3"
            txt["$PAR"] <- as.character(npar)
            txt["$TOT"] <- as.character(dim(data)[1])
            txt["$BYTEORD"] <- "4,3,2,1"
            txt["$DATATYPE"] <- "F"
            txt["$MODE"] <- "L"
            txt["$NEXTDATA"] <- "0"

            txt["$BEGINANALYSIS"] <- "0"
            txt["$ENDANALYSIS"] <- "0"
            txt["$BEGINSTEXT"] <- "0"
            txt["$ENDSTEXT"] <- "0"
            dataStart <- 1000 ## This doesn't really matter, it's not an FCS file anyway, but let's provide something reasonable
            txt["$BEGINDATA"] <- as.character(dataStart)
            txt["$ENDDATA"] <- as.character(dataStart + npar * dim(data)[1] * 4 - 1) ## $DATATYPE is F and PnB is 32, so here we do: start + rows x cols x 4 bytes per value - 1

            txt[paste(id,"N",sep="")] <- names
            txt[paste(id,"S",sep="")] <- desc
            txt[paste(id,"R",sep="")] <- as.character(range)
            txt[paste(id,"B",sep="")] <- "32"
            txt[paste(id,"E",sep="")] <- "0,0"
            txt[paste(id,"G",sep="")] <- "1"
            txt[paste("flowCore_",id,"Rmin",sep="")] <- as.character(minRange)
            txt[paste("flowCore_",id,"Rmax",sep="")] <- as.character(maxRange)

            txt["$FIL"] <- object@runID
            txt["FILENAME"] <- object@runID
            txt["$SYS"] <- paste(Sys.info()['sysname'], Sys.info()['release'])
            txt["CREATOR"] <- R.version.string
            txt["transformation"] <- "applied"
            txt["GUID"] <- object@runID  #generateGUID(data)
            txt["ORIGINALGUID"] <- object@runID #generateGUID(data)

            txt["$INST"] <- "Unknown institution"
            txt["$OP"] <- "Unknown operator"
            txt["$CYT"] <- "Unknown cytometer"
            txt["$SRC"] <- "Unknown source"

            curtime <- Sys.time()
            txt["$DATE"] <- toupper(format(curtime, "%d-%b-%Y")) ## Current data in the format of "05-OCT-2017"
            timeDelta <- 120 ## picking 2 minutes as cytometer run time
            if ("time" %in% tolower(names)) {
              txt["$TIMESTEP"] <- "0.01" ## We don't know that, so just something that would be reasonable-ish
              timeIndex <- which(tolower(names) == "time")
              timeDelta <- ceiling((max(data[,timeIndex]) - min(data[,timeIndex])) * 0.01)
            }
            if (toupper(format(curtime - timeDelta, "%d-%b-%Y")) == toupper(format(curtime, "%d-%b-%Y"))) {
              ## Subtracting timeDelta is still the same day, so let's do it that way.
              txt["$BTIM"] <- format(curtime - timeDelta, "%H:%M:%S")
              txt["$ETIM"] <- format(curtime, "%H:%M:%S")
            } else {
              ## Subtracting timeDelta gets us to yesterday, that's not good, so let's just say we started at midnight and add delta to it as ETIM...
              ## (all this is arbitrary anyway, we just want resonable values)
              txt["$BTIM"] <- "00:00:00"
              txt["$ETIM"] <- sprintf("%02i:%02i:%02i", floor(timeDelta / 3600), floor((timeDelta%%3600) / 60), floor(timeDelta%%60))
            }
            colnames(data) <- names
            object@flowFrameFJ <- list(fcs=new("flowFrame",
                                               exprs=as.matrix(data),
                                               parameters=pars,
                                               description=txt),
                                       cellIds=object@inputCSV$eventNumbers)

            return(object)
          }
)
################################################################
#    Perform FlowSOM Clustering  #
################################################################
#' @title Perform FlowSOM Clustering
#' @return Returns object with FlowSOM clustering results
#' @param object An insideFlow object
#' @param applyOnMap Boolean to apply to trained map or not
#' @param trainedMap Object which contains trained map
#' @param fsXdim Dimension 1 for grid size used for clustering
#' @param fsYdim Dimension 2 for grid size used for clustering
#' @param numClust Number of clusters to return from FlowSOM
#' @importFrom FlowSOM UpdateFlowSOM FlowSOM NewData
#' @examples \dontrun{
#' object <- calcFlowSOM(object, applyOnMap=FALSE, fsXdim=10, fsYdim=10, numClust=8)}
#' @docType methods
#' @export
setGeneric("calcFlowSOM", function(object, applyOnMap, trainedMap, fsXdim, fsYdim, numClust){
  standardGeneric("calcFlowSOM")
})
#' @rdname calcFlowSOM
setMethod("calcFlowSOM",
          "insideFlow",
          function(object, applyOnMap, trainedMap, fsXdim, fsYdim, numClust) {
            # parNames <- c(colnames(myFJobj@flowFrameFJ$fcs))
            # parIndices <- match(parNames, colnames(object@flowFrameFJ$fcs))
            # if (length(parNames) == 0 || length(parIndices) == 0)
            #   stop("Something seems wrong, it's like the input FCS file does not contain the provided input parameters.", call.=FALSE)

            if (applyOnMap) {
              ## This is FS2.0 method with NewData() wrapped in UpdateFlowSOM() in case a user loads older FlowSOM object##
              object@FlowSOMClusters[["som"]] <- UpdateFlowSOM(NewData(fsom = trainedMap@FlowSOMClusters[["som"]], input = object@flowFrameFJ$fcs, madAllowed = 4))
            } else {
              object@FlowSOMClusters[["som"]] <- FlowSOM(
                object@flowFrameFJ$fcs,
                compensate = FALSE,
                transform = FALSE,
                scale = FALSE,
                colsToUse = NULL,
                xdim = fsXdim,
                ydim = fsYdim,
                nClus = numClust,
                seed = 23) #FJ_SEED)
            }

            return(object)
          })
################################################################
#    Plot FlowSOM Result Images  #
################################################################
#' @title Plot FlowSOM Images
#' @return Returns reversed transformed data matrix before export to FlowJo
#' @param object A insideFlow object
#' @param plotAllAsStars Boolean to plot MST by selected param or with piecharts
#' @param plotPar Selected parameter to plot as heatmapped MST
#' @param radarColor The selected color palette for minimum spanning tree plot
#' @param plotBgd Color palette for mst nodes
#' @param numClust Number of metaclusters requested by user
#' @importFrom FlowSOM AddStars AddStarsPies AddLabels
#' @importFrom FlowSOM UpdateFlowSOM GetChannels PlotFlowSOM FlowSOM_colors
#' @importFrom grDevices hcl.colors
#' @importFrom grDevices rainbow
#' @importFrom ggplot2 aes coord_fixed geom_segment ggplot
#' @importFrom ggplot2 theme theme_void labs xlim
#' @importFrom dplyr count
#' @importFrom ggnewscale new_scale
#' @importFrom grDevices colorRampPalette
#' @importFrom ggpubr as_ggplot get_legend ggarrange
#' @examples \dontrun{
#' myobject <- plotFlowSOMResults(object, popName, coFact, outsFolder)}
#' @docType methods
#' @export
setGeneric("plotFlowSOMResults", function(object, plotAllAsStars, plotPar, radarColor, plotBgd, numClust){
  standardGeneric("plotFlowSOMResults")
})
#' @rdname plotFlowSOMResults
setMethod("plotFlowSOMResults",
          "insideFlow",
          function(object, plotAllAsStars, plotPar, radarColor, plotBgd, numClust) {

            ## FlowSOM functions so we can modify text size of radar legend ##
            ParseArcs <- function(x, y, arcValues, arcHeights){
              markers <- names(arcValues)
              arcValues <- c(0, (arcValues / sum(arcValues)) * (2 * pi))
              arcValues <- cumsum(arcValues)
              resDf <- data.frame(Markers = markers, x0 = x, y0 = y,
                                  start = arcValues[-length(arcValues)],
                                  end = arcValues[-1], value = arcHeights)
              return(resDf)
            }
            # AddScale
            AddScale <- function(p,
                                 values = NULL,
                                 colors = NULL,
                                 limits = NULL,
                                 showLegend = TRUE,
                                 labelLegend = "",
                                 type = "fill"){
              p <- p + ggnewscale::new_scale(type)
              if(is.character(values) | is.logical(values)) values <- factor(values)
              if (!is.null(limits) &&  is.null(colors)){
                args <- list()
                args[[type]] <- limits
                p <- p + do.call(ggplot2::lims, args)
              }
              if(!is.null(colors)){
                # Discrete values
                if (is.factor(values)){
                  # If a function: make into the right amount of values
                  if(is.function(colors)){
                    colors <- colors(length(levels(values)))
                  }
                  # Check color names vs backgroundValues
                  colorNames <- names(colors)
                  levels <- levels(values)
                  # If no names available, use levels
                  if(is.null(colorNames)){
                    colorNames <- levels[seq_along(colors)]
                  }
                  # If any names missing, show warning
                  notAvailable <- ! levels %in% colorNames
                  if(any(notAvailable)){
                    warning("You did not provide a color for ",
                            paste(levels[notAvailable], collapse =", "))
                    colors <- c(colors,
                                rep("#FFFFFF", sum(notAvailable)))
                    names(colors) <- c(colorNames,
                                       levels[notAvailable])
                  }
                  scale_function <- eval(parse(text = paste0("ggplot2::scale_",
                                                             type,
                                                             "_manual")))
                  if(showLegend) {
                    p <- p + scale_function(values = colors)
                  } else {
                    p <- p + scale_function(values = colors,
                                            guide = FALSE)
                  }
                } else if (is.numeric(values)){ # Continuous values
                  if (is.function(colors)){
                    colors <- colors(100)
                  } else {
                    colors <- grDevices::colorRampPalette(colors)(100)
                  }
                  scale_function <- eval(parse(text = paste0("ggplot2::scale_",
                                                             type,
                                                             "_gradientn")))
                  if(showLegend) {
                    p <- p + scale_function(colors = colors,
                                            limits = limits)
                  } else {
                    p <- p + scale_function(colors = colors,
                                            limits = limits,
                                            guide = FALSE)
                  }
                }
              }
              p <- p + ggplot2::labs(fill = labelLegend)
              return(p)
            }
            ## Legends
            PlotStarLegend.FJ <- function(markers,
                                          colors,
                                          starHeight = 1){
              requireNamespace("ggplot2")
              markers <- factor(markers, levels = markers)
              nMarkers <- length(markers)
              circularCoords <- seq(from = 2 * pi / (nMarkers * 2),
                                    by = 2 * pi / nMarkers,
                                    length.out = nMarkers)
              dfSegments <- data.frame(Markers = markers,
                                       x = sin(circularCoords),
                                       y = cos(circularCoords),
                                       xend = 1.1 *  ifelse(sin(circularCoords) >= 0, 1, -1),
                                       yend = NA)
              nLeftRight <- dplyr::count(dfSegments, .data$x >= 0)
              if (nrow(nLeftRight) != 1){
                by <- ifelse(nMarkers <= 8, 1, 0.65)
                left <- seq(from = 0, by = by, length.out = nLeftRight[1, 2])
                right <- -seq(from = 0, by = by, length.out = nLeftRight[2, 2])
                dfSegments[which(dfSegments$x < 0), ]$yend <- left - mean(left)
                dfSegments[which(dfSegments$x >= 0), ]$yend <- right - mean(right)
              } else {
                dfSegments$yend <- -2
              }
              horizontalLines <- data.frame(Markers = dfSegments$Markers,
                                            x = dfSegments$xend, y = dfSegments$yend,
                                            xend = dfSegments$xend +
                                              ifelse(dfSegments$xend >= 0, 0.5, -0.5),
                                            yend = dfSegments$yend,
                                            stringsAsFactors = FALSE)
              dfSegments <- rbind(dfSegments, horizontalLines)
              dfLabels <- data.frame(x = horizontalLines$xend +
                                       ifelse(horizontalLines$xend >= 0, 0.3, -0.3),
                                     y = horizontalLines$yend)
              markers_tmp <- rep(1, nMarkers)
              names(markers_tmp) <- markers
              dfStar <- ParseArcs(0, 0, markers_tmp, starHeight)
              dfStar$Markers <- factor(dfStar$Marker, levels = markers)
              l <- ggplot2::ggplot() +
                ggplot2::xlim(c(-5, 6)) +
                ggplot2::coord_fixed(clip = "off") +
                ggplot2::theme_void() +
                ggplot2::theme(legend.position = "none")
              l <- AddStarsPies(l,
                                dfStar,
                                colors)
              l <- AddScale(p = l,
                            values = markers,
                            colors = colors,
                            showLegend = FALSE,
                            type = "color")
              l <- l +
                ggplot2::geom_segment(data = dfSegments,
                                      ggplot2::aes(x = .data$x, y = .data$y,
                                                   xend = .data$xend,
                                                   yend = .data$yend,
                                                   color = .data$Markers),
                                      size = 0.6)

              l <- AddLabels(l,
                             labels = markers,
                             layout = dfLabels,
                             hjust = ifelse(horizontalLines$xend >= 0, 0, 1),
                             textSize = 3) # This is where we can control text size of radar legends
              return(l)
            }
            # Stars
            PlotStars.FJ <- function(fsom,
                                     markers = fsom$map$colsUsed,
                                     colorPalette = FlowSOM_colors,
                                     list_insteadof_ggarrange = FALSE,
                                     ...){
              fsom <- UpdateFlowSOM(fsom)
              channels <- GetChannels(fsom, markers)
              p <- PlotFlowSOM(fsom = fsom,...)

              if(!is.null(names(colorPalette))) {
                names(colorPalette) <- GetChannels(fsom, names(colorPalette))
              }
              p <- AddStars(p = p,
                            fsom = fsom,
                            markers = channels,
                            colorPalette = colorPalette)
              if(!is.null(names(colorPalette))) {
                names(colorPalette) <- fsom$prettyColnames[GetChannels(fsom,names(colorPalette))]
              }
              l1 <- PlotStarLegend.FJ(fsom$prettyColnames[channels],colorPalette)
              l2 <- ggpubr::get_legend(p)
              if (list_insteadof_ggarrange){
                p <- p + ggplot2::theme(legend.position = "none")
                l2 <- ggpubr::as_ggplot(l2)
                return(list(tree = p,
                            starLegend = l1,
                            backgroundLegend = l2))
              } else {
                p <- ggpubr::ggarrange(p,
                                       ggpubr::ggarrange(l1, l2, ncol = 1),
                                       NULL,
                                       ncol = 3, widths = c(3, 1.5,0.3),
                                       legend = "none")
                return(p)
              }
            }

            if(plotBgd == "Rainbow"){
              nodeBackgroundPalette <- grDevices::rainbow(n=numClust, alpha = 0.3 )
            } else {
              nodeBackgroundPalette <- hcl.colors(n=numClust, palette = plotBgd, alpha = 0.6 )
            }

            n_radarSlices <-length(object@FlowSOMClusters$som$map$colsUsed)

            if(radarColor == "Rainbow"){
              radarPalette <- grDevices::rainbow(n=n_radarSlices, alpha = 0.3 )
            } else {
              radarPalette <- hcl.colors(n=n_radarSlices, palette = radarColor, alpha = 0.6 )
            }

            # if (!plotAllAsStars) {
            #   ## If not plotting all then make sure plotPar is among parameters in the FCS file
            #   if (!(plotPar %in% colnames(object@flowFrameFJ$fcs))) {
            #     plotPar2 <- gsub("[", "<", plotPar, fixed=TRUE)
            #     plotPar2 <- gsub("]", ">", plotPar2, fixed=TRUE)
            #     if (plotPar2 %in% colnames(object@flowFrameFJ$fcs)) {
            #       plotPar <- plotPar2
            #     } else {
            #       plotPar2 <- gsub("_", "/", plotPar, fixed=TRUE)
            #       if (plotPar2 %in% colnames(object@flowFrameFJ$fcs)) {
            #         plotPar <- plotPar2
            #       } else {
            #         plotPar3 <- gsub("[", "<", plotPar2, fixed=TRUE)
            #         plotPar3 <- gsub("]", ">", plotPar3, fixed=TRUE)
            #         if (plotPar3 %in% colnames(object@flowFrameFJ$fcs)) {
            #           plotPar <- plotPar3
            #         } else {
            #           stop(paste("The input FCS file does not contain ", plotPar), call.=FALSE)
            #         }
            #       }
            #     }
            #   }
            # }

            if (plotAllAsStars) {
              ## Use PlotStars
              if (bgMetaClusters) {
                p2 <- PlotStars.FJ(
                  object@FlowSOMClusters$som,
                  view=plotView,
                  colorPalette = radarPalette,
                  backgroundColors = nodeBackgroundPalette,
                  maxNodeSize = (1*0.01*200),
                  equalNodeSize = FALSE,
                  backgroundValues = as.factor(as.numeric(object@FlowSOMClusters$som$metaclustering))) ## Used in FlowJo as -1 to make it go 0..n-1 instead 1..n to match FlowJo's populations
              } else {
                p2 <- PlotStars.FJ(
                  object@FlowSOMClusters$som,
                  view=plotView,
                  colorPalette = radarPalette,
                  maxNodeSize = (1*0.01*200),
                  equalNodeSize = FALSE)
              }
            } else {
              ## Use PlotMarker
              if (bgMetaClusters) {
                p2 <- PlotMarker(
                  object@FlowSOMClusters$som,
                  plotPar,
                  view=plotView,
                  colorPalette = nodeBackgroundPalette,
                  maxNodeSize = (1*0.01*200),
                  equalNodeSize = FALSE,
                  backgroundValues = as.factor(as.numeric(object@FlowSOMClusters$som$metaclustering))) ## Used in FlowJo as -1 to make it go 0..n-1 instead 1..n to match FlowJo's populations
              } else {
                p2 <- PlotMarker(
                  object@FlowSOMClusters$som,
                  plotPar,
                  view=plotView,
                  colorPalette = nodeBackgroundPalette,
                  maxNodeSize = (1*0.01*200),
                  equalNodeSize = FALSE,
                )
              }
            }

            object@FlowSOMClusters$som$prettyColnames <- gsub("(<)(.+?)(>)","", object@FlowSOMClusters$som$prettyColnames)

            return(p2)


          })
################################################################
#    Perform Phenograph Clustering  #
################################################################
