runID <- "aloha_aina"
fjXDim <- 10
fjYDim <- 10
EventNumberDP <- "EventNumberDP"
CellIdDP <- "CellId"
fjPopName <- "CD3+"
fjDims <- 2
fjSmooth <- 0 # good range -3 (sharper) to 3 (smoother)
fjK <- 25 # good range 5 to 100
fjAdjust <- 0 / 100.0 # good range 0 to 3
fjSeed <- 1 # any
fjUmap.Dist <- 0.5 # good range 0.1 to 2
fjXdim <- 20 # it's sometimes useful to asymmetric SOM (Xdim != Ydim) between 5 to 64
fjYdim <- 20 # it's sometimes useful to asymmetric SOM (Xdim != Ydim) between 5 to 64
fjRLEN <- 15 # good range here is between 5 to 50
fjMapSize <- 5000 # 1000 is good in almost all cases. Good range 50 to 5000
fjSuperWeight <- 0.2 #values can range from 0.0 to 1.0
fjEmbedSom.Map <- "FALSE"
fjFlowSOM.Map <- "FALSE"
fjUMAP.Map <- "FALSE"
fjTSNE.Map <- "FALSE"
fjSuper.Map <- "FALSE"

load(file="data/exMatrix.rda")

testthat::test_that("Create a new insideFlow object", {
  myFJobj <- new("insideFlow", runID = "runID")
})

testthat::test_that("Can load new dataset", {
  myFJobj <- loadFlowJoDataTable(myFJobj, exMatrix)
})

testthat::test_that("Gets the parameter names" {
  getParNames(myFJobj)
})
