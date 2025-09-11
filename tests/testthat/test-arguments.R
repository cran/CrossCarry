test_that("Invalid arguments produce errors", {
  data("Arterial", package = "CrossCarry")
  cc <- createCarry(Arterial, "Treatment", "Period", "Subject", TRUE)

  # Invalid Mv
  expect_error(CrossGEE("Pressure", "Treatment", "Period", "Subject",
                        carry = cc$carryover, data = cc$data,
                        correlation = "AR-M", Mv = 0))

  # Nonexistent id variable
  expect_error(CrossGEE("Pressure", "Treatment", "Period", "WrongID",
                        carry = cc$carryover, data = cc$data))

  # Nonexistent time variable in Kron
  expect_error(CrossGEEKron("Pressure", "Treatment", "Period", "Subject",
                            time = "WrongTime",
                            carry = cc$carryover, data = cc$data,
                            correlation = "AR-M", Mv = 1))

  # Invalid nodes argument
  expect_error(CrossGEESP("Pressure", "Treatment", "Period", "Subject",
                          time = "Time",
                          carry = cc$carryover, data = cc$data,
                          correlation = "AR-M", Mv = 1,
                          nodes = -2))
})
