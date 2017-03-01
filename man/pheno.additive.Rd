\name{pheno.additive}
\alias{pheno.additive}

\title{ Additive Phenotype Simulations }

\description{ .... }
\usage{ pheno.additive(x, beta, sparse = TRUE) }

\arguments{
  \item{x}{ A \code{\link{bed.matrix}} }
}

\details{ bla bla

NOTE: monorphic SNPs are skipped
}

\value{
blu blu
}

% \seealso{  }

\examples{
# Load data
data(AGT)
x <- as.bed.matrix(AGT.gen, AGT.fam, AGT.bim)  # 361 SNPs

# define (sparse) vector of effects
beta <- runif(361)
beta[ sample.int(361, 250) ] <- 0

# three different ways to compute the phenotypes
standardize(x) <- "p"
phenos.1 <- as.vector(x \%*\% beta)
phenos.2 <- pheno.additive(x, beta, sparse = FALSE)
phenos.3 <- pheno.additive(x, beta, sparse = TRUE)

# all three should be equal up to rounding errors
range( phenos.1 - phenos.2 )
range( phenos.1 - phenos.3 )

}

% \keyword{ Burden }
