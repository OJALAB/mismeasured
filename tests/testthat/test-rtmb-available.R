# CI guard: the onestep estimator's tests are gated on RTMB. Without this
# check a failed RTMB install on CI silently skips all of them and the
# estimator ships untested. CI sets CI_REQUIRE_RTMB=true.

test_that("RTMB is installed when the CI requires it", {
  skip_if(!nzchar(Sys.getenv("CI_REQUIRE_RTMB")),
          "CI_REQUIRE_RTMB not set")
  expect_true(requireNamespace("RTMB", quietly = TRUE))
})
