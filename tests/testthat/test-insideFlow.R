test_that("Create an insideFlow object", {
  obj <- new("insideFlow", runID = "aloha_aina")
  # Check that the object is an S4 insideFlow object
  testthat::expect_s4_class(obj, "insideFlow")
  # Check that the runID is set correctly
  expect_equal(obj@runID, "aloha_aina")
})
