library(Rcpp)

source(file.path('R functions', 'nNLS functions.R'))

test.summary <- data.frame('p value'=c(0.01, 0.02, 0.2, 0.5), check.names=FALSE)
test.result <- nparcomp.adjust(test.summary)
test.expected <- c(0.03940399, 0.058808, 0.36, 0.5)

stopifnot(isTRUE(all.equal(test.result$'p adjusted', test.expected, tolerance=1e-08)))

test.empty <- data.frame('p value'=numeric(0), check.names=FALSE)
test.empty.result <- nparcomp.adjust(test.empty)

stopifnot(nrow(test.empty.result)==0)
stopifnot(identical(test.empty.result$'p adjusted', numeric(0)))
stopifnot(inherits(try(nparcomp.adjust(NULL), silent=TRUE), 'try-error'))
stopifnot(inherits(try(nparcomp.adjust(data.frame(p=0.5)), silent=TRUE), 'try-error'))
stopifnot(inherits(try(nparcomp.adjust(data.frame('p value'=NA_real_, check.names=FALSE)), silent=TRUE), 'try-error'))

test.large <- data.frame('p value'=seq(0, 1, length.out=100000), check.names=FALSE)
test.large.result <- nparcomp.adjust(test.large)

stopifnot(nrow(test.large.result)==100000)
stopifnot(all(is.finite(test.large.result$'p adjusted')))
