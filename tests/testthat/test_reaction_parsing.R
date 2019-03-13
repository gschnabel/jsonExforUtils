context("test parsing of REACTION field")
library(testthat)

test_that("reactions are parsed properly", {

  reacStr <- "(100-FM-257(N,ABS),,RI,,LIM,RECOM) ignored free text"
  reacStruc <- parseReacStr(reacStr)
  expect_equal(reacStruc$type, "reac")
  expect_equal(reacStruc$reac, "(100-FM-257(N,ABS),,RI,,LIM,RECOM)")
  expect_equal(reacStruc$target$Z, 100)
  expect_equal(reacStruc$target$A, 257)
  expect_equal(reacStruc$target$sym, "FM")
  expect_equal(reacStruc$projectile, "N")
  expect_equal(reacStruc$process, "ABS")
  expect_null(reacStruc$residual)
  expect_equal(reacStruc$quantspec, ",RI,,LIM,RECOM")

  # something with residual
  reacStr <- "(98-CF-252(0,F)0-NN-1,PR/PAR,KE)"
  reacStruc <- parseReacStr(reacStr)
  expect_equal(reacStruc$residual$Z, 0)
  expect_equal(reacStruc$residual$A, 1)

  # something with metastable state
  reacStr <- "(78-PT-198(N,G)78-PT-199-M,,SIG,,MXW,RECOM)"
  reacStruc <- parseReacStr(reacStr)
  expect_equal(reacStruc$quantspec, ",SIG,,MXW,RECOM")
  expect_equal(reacStruc$residual$meta, "M")
})


test_that("reaction expressions are parsed properly", {

  reacStr <- "((64-GD-158(N,G)64-GD-159,,SIG,,MXW,RECOM)/ (64-GD-158(N,G)64-GD-159,,SIG,,,RECOM))"
  reacStruc <- parseReacExpr(reacStr)
  expect_equal(reacStruc$type, "expr")
  expect_equal(reacStruc$ops, "/")
  expect_true(length(reacStruc$nodes)==2)
  expect_equal(reacStruc$nodes[[1]]$reac, "(64-GD-158(N,G)64-GD-159,,SIG,,MXW,RECOM)")
  expect_equal(reacStruc$nodes[[2]]$reac, "(64-GD-158(N,G)64-GD-159,,SIG,,,RECOM)")
})





