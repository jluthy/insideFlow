test_that("Data is correctly loaded from CSV", {

  EventNumberDP <- "Event #" # This is the name of the column in the CSV file that contains the event number
  data_path <- system.file("extdata", "exportPar_and_Stain_Names_Live.csv", package = "insideFlow")
  obj <- new("insideFlow", runID = "aloha_aina")
  obj <- loadFlowJoCSV(obj, data_path, eventNumber = EventNumberDP)

  # Check that the data is correctly loaded as Matrix
  expect_true(is.matrix(obj@inputCSV$X))
  # Check that parameter names are correctly extracted
  expect_true(length(getParNames(obj)) > 0)
  # Check that stain names are correctly extracted
  expect_true(length(getStainNames(obj)) > 0)
})
