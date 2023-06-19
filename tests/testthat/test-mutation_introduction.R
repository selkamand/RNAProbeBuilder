test_that("inroduction of dup mutations work", {
  wt = "ACTGA"
  expected_mut = "ACTTGA"
  position = 3
  observed_mut= mutate_seq_dup_scalar(sequence = wt, position = position)
  expect_equal(observed_mut, expected_mut)

  wt = "ACTGA"
  expected_mut = "ACTGAA"
  position = 5
  observed_mut= mutate_seq_dup_scalar(sequence = wt, position = position)
  expect_equal(observed_mut, expected_mut)

  wt = "ACTGA"
  expected_mut = "AACTGA"
  position = 1
  observed_mut= mutate_seq_dup_scalar(sequence = wt, position = position)
  expect_equal(observed_mut, expected_mut)
})

test_that("inroduction of SNV mutations work", {
  wt = "ACTGA"
  expected_mut = "ACCGA"
  position = 3
  observed_mut= mutate_seq_snv(sequence = wt,position = position, ref = "T", alt = "C")
  expect_equal(observed_mut, expected_mut)
})
