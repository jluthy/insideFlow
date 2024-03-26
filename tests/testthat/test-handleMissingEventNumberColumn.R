test_that("Function handles missing event number column gracefully", {

  data_path <- system.file("extdata", "sample_data_without_event_num.csv", package = "yourPackageName")
  obj <- new("insideFlow", runID = "aloha_aina")

  expect_warning(
    obj <- loadFlowJoCSV(obj, data_path, eventNumber = "NonExistentColumn"),
    "Column NonExistentColumn not found in the provided CSV. Filling eventNumbers with NA."
  )

  # Assuming the function fills with NA or similar fallback
  expect_true(all(is.na(obj@inputCSV$eventNumbers)))
})
