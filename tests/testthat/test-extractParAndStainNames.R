test_that("Fallback to parameter names when stain names are missing", {

  EventNumberDP <- "Event #" # This is the name of the column in the CSV file that contains the event number
  data_path <- system.file("extdata", "exportParameterNames_Live.csv", package = "insideFlow")
  obj <- new("insideFlow", runID = "aloha_aina")
  obj <- loadFlowJoCSV(obj, data_path, eventNumber = EventNumberDP)

  # Check that parameter names are used when stain names are missing
  expect_true(length(getParNames(obj)) > 0)
  # Check that the Stain Names contain NA
  expect_true(all(is.na(getStainNames(obj))))
  # expect_true(length(getStainNames(obj)) > 0)
})
