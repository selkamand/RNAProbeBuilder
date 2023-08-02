test_that("alignment QC works", {
  input = c("A", "A","A","-","A","C",'-','-','A',"-")
  expected = c(1, 2, 3, NA, 4, 5, NA, NA, 6, NA)
  observed <- calculate_position(input)

  expect_equal(observed, expected)
})
