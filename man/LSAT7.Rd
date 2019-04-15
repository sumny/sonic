\name{LSAT7}
\alias{LSAT7}

\title{Law School Admissions Test Section 7 (Data)}

\description{
  Data from Bock & Lieberman (1970). Contains responses of 1000 persons to 5
  dichotomously scored items (0, 1), obtained from the Law School Admissions
  Test, Section 7.
}

\usage{data("LSAT7")}

\format{
A numeric matrix containing responses of 1000 persons to 5 items.
}

\references{
  Bock RD, Lieberman M (1970).
    Fitting a Response Model for n Dichotomously Scored Items.
    \emph{Psychometrika}, \bold{35}(2), 179--197.
}

\examples{
data("LSAT7", package = "sonic")
head(LSAT7)
}

\keyword{datasets}
