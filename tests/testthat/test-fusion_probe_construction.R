test_that("fusion sequence construction", {

  fusion_string = "CACTG|AAAAT"


  # Even Probe Sizes
  expect_equal(
    probes_construct_fusion_sequence(fusion = fusion_string, probe_size = 2),
    "GA"
    )

  expect_equal(
    probes_construct_fusion_sequence(fusion = fusion_string, probe_size = 4),
    "TGAA"
  )

  # Odd Probe Sizes (centered base should be first base of downstream fusion partner)
  expect_equal(
    probes_construct_fusion_sequence(fusion = fusion_string, probe_size = 5),
    "TGAAA"
  )

  expect_error(
    probes_construct_fusion_sequence(fusion = fusion_string, probe_size = 11),
    "long"
    )


})
