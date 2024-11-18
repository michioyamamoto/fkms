\name{fkms}
\alias{fkms}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	Functional K-means Clustering for Regularly and Irregularly Observed Longitudinal Data
}
\description{
	This function conducts the functional k-means clustering proposed	by Yamamoto and Terada (2024, submitted).
}
\usage{
	fkms(fdata, dformat=1, N.clust=2, N.bs=NULL, N.random=1, N.rep=100, ini.cls=NA, show.random.ite=FALSE, basis="B", lambda=0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
	\item{fdata}{
		Specifies the data. It must be in either list format or a long-format data frame. In list format, each element of the list is a matrix of \eqn{N_i\times2} with the number of time points \eqn{N_i}, where the first column represents measurement points and the second column represents measurement values. In long-format data frame, the first column should contain the subject identifier, the second column the measurement times, and the third column the measurement values. Optionally, the fourth column can specify weights for each subject’s measurements. If the fourth column is not provided, weights are generated automatically.
	}
	\item{dformat}{
		Data format; 1 for list format, 2 for data frame format.
	}
	\item{N.clust}{
		Number of clusters.
	}
	\item{N.bs}{
		Number of basis functions. When using Fourier basis, this must be an odd number. If an even number is provided, the next odd number will be used as the value of \code{N.bs}.
	}
	\item{N.random}{
		Number of random initial starts.
	}
	\item{N.rep}{
		Maximum number of algorithm iterations.
	}
	\item{ini.cls}{
		Initial cluster vector. If not specified, random assignment is used.
	}
	\item{show.random.ite}{
		If \code{TRUE}, displays the number of the random start currently being executed.
	}
	\item{basis}{
		Specifies the basis functions to use. \code{"B"}=B-spline basis, \code{"F"}=Fourier basis.
	}
	\item{lambda}{
		Specifies the smoothing parameter. Currently, this is only applicable when \code{basis}=\code{"B"}.
	}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
	\item{cluster}{
		A vector containing the estimated cluster assignments for each subject.
	}	
	\item{centers}{
		A matrix consisting of the estimated cluster center functions.
	}	
	\item{x.centers}{
		The measurement points at which the functions stored in centers are evaluated.
	}	
	\item{res.fit}{
		The estimated regression coefficients for the cluster center functions.
	}	
	\item{n.opt}{
		The index of the initial value among the \code{N.random} starts that yielded the optimal result.
	}	
	\item{mse}{
		The mean squared error.
	}	
}
\references{
	Yamamoto, M. and Terada, Y. (2024). K-means clustering for sparsely observed longitudinal data. arXiv:2411.08256.
}
\author{
	Michio Yamamoto\cr
	\email{yamamoto.michio.hus@osaka-u.ac.jp}
}
% \note{
% %%  ~~further notes~~
% }

% %% ~Make other sections like Warning with \section{Warning }{....} ~

% \seealso{
% %% ~~objects to See Also as \code{\link{help}}, ~~~
% }
\examples{
	## Example data with dformat=1 (list format)
	dat <- SampleFunc(N.clust=2, reg=FALSE, dformat=1)
	fdata <- dat$data ## sample curves
	cls.t <- dat$cls ## true clusters

	## functional k-means clustering using B-spline basis functions with regularization
	res.fkm <- fkms(fdata=fdata, dformat=1, N.clust=2, N.bs=5, N.random=10, basis="B", lambda=0.1)

	## estimated assignments
	res.fkm$cluster

	## estimated cluster centers
	plot(res.fkm$x.centers, res.fkm$centers[, 1], type="l", xlim=c(0, 1), ylim=c(0, 1))
	lines(res.fkm$x.centers, res.fkm$centers[, 2], lty=2, lwd=4, col="red", xlim=c(0, 1), ylim=c(0, 1))



	## Example data with dformat=2 (data.frame format)
	dat <- SampleFunc(N.clust=3, reg=FALSE, dformat=2)
	fdata <- dat$data ## sample curves
	cls.t <- dat$cls ## true clusters
		
	## functional k-means clustering using B-spline basis functions with regularization
	res.fkm <- fkms(fdata=fdata, dformat=2, N.clust=3, N.bs=5, N.random=10, basis="B", lambda=0.1)

	## estimated assignments
	res.fkm$cluster

	## estimated cluster centers
	plot(res.fkm$x.centers, res.fkm$centers[, 1], type="l", xlim=c(0, 1), ylim=c(-1, 1))
	lines(res.fkm$x.centers, res.fkm$centers[, 2], lty=2, lwd=4, col="red", xlim=c(0, 1), ylim=c(-1, 1))
	lines(res.fkm$x.centers, res.fkm$centers[, 3], lty=3, lwd=4, col="blue", xlim=c(0, 1), ylim=c(-1, 1))
}

