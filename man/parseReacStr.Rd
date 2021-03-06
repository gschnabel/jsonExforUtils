% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{parseReacStr}
\alias{parseReacStr}
\title{Parse Reaction String}
\usage{
parseReacStr(reacStr, pos = 1, quiet = TRUE)
}
\arguments{
\item{reacStr}{the reaction string}

\item{pos}{position in \code{reacStr} where parsing should start}
}
\value{
A list with the following components:
        \describe{
          \item{start}{start position of reaction string}
          \item{end}{end position of reaction string}
          \item{reac}{the matched reaction string}
          \item{type}{contains the string "reac"}
          \item{projectile}{projectile}
          \item{process}{usually gives ejectile information}
          \item{target}{list describing target with elements \code{A, Z, sym}}
          \item{residual}{list describing residual nucleus with elements \code{A, Z, sym, meta}}
          \item{quantspec}{string with additional flags indicating among other things observable type, e.g. DA or SIG}
        }
}
\description{
Parse Reaction String
}
