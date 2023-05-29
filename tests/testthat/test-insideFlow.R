# set up test
# globalVariables(c('runID','EventNumberDP','CellIdDP', 'exMatrix'))
# runID <- "aloha_aina"
EventNumberDP <- "EventNumberDP"
CellIdDP <- "CellId"
# exMatrix <- readRDS(file= system.file('tests/testthat/', 'exMatrix.rds', package = 'insideFlow'))

testthat::test_that("Create a new insideFlow object", {
  myFJobj <- new("insideFlow", runID = "aloha_aina")
  getRunID(myFJobj)
})

# testthat::test_that("Can load new dataset", {
#   myFJobj <- new("insideFlow", runID = "testID")
#   myFJobj <- loadFlowJoDataTable(myFJobj, exMatrix)
#   getParNames(myFJobj)
# })

# testthat::test_that("Gets the parameter names", {
#   getParNames(myFJobj)
# })
# saveRDS(exMatrix, file="./tests/testthat/exMatrix.rds")
