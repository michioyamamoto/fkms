##☆★☆★☆★☆★☆★☆★☆★☆★☆★☆★☆★☆★☆★☆★☆★
##      Rパッケージでのサンプルデータ生成関数
##
##  ファイル名：SampleFunc.R
##  ファイル内容：
##  作成者：YAMAMOTO, Michio
##  作成日：2024年09月17日
##  最終更新日：2024年09月17日
##  コメント：funcyパッケージのsampleFuncy関数を参考にした
##☆★☆★☆★☆★☆★☆★☆★☆★☆★☆★☆★☆★☆★☆★☆★


##### 引数 #####
## *N.sub.clust: Sample size for each cluster.
## *N.clust: Number of clusters (1-5).
## *N.tp: Number of time points for regular datasets.
## *max.tp: Maximal number of time points for irregular datasets.
## *min.tp: Minimal number of time points for irregular datasets.
## *int.tp: Time interval where time points are drawn from. 長さ2のベクトルで指定する必要がある。
## *sd: Standard deviation from the cluster centers.
## *reg: TRUE=regular datasets, FALSE=irregular datasets.
## *dformat: data format; 1=list format, 2=data.frame format.

#### 詳細 ####
## クラスタ中心関数は順に、x^2, sqrt(x), sin(2*pi*x), x^3 としている

#### 戻り値 ####
## *data: 生成されたデータセット
## *cls: 各対象が所属するクラスタ番号


## 関数本体
## N.sub.clust <- 50; N.clust <- 3; N.tp <- 10; max.tp <- 10; min.tp <- 2; int.tp <- c(0, 1); sd <- 0.2; reg <- TRUE; format <- 1
SampleFunc <- function (N.sub.clust=50, N.clust=2, N.tp=10, max.tp=10, min.tp=2, int.tp=c(0, 1), sd=0.2, reg=TRUE, dformat=1) {
	if (N.clust < 2 | N.clust > 4)
		stop("Error: Check the argument 'N.clust'")

	N.sub <- N.sub.clust * N.clust ## 全体のサンプルサイズ
	cls.vec <- NULL
	for (n.clust in 1:N.clust)
		cls.vec <- c(cls.vec, rep(n.clust, len=N.sub.clust))
	
	if (reg) {
		## regular dataset
		grid.tp <- seq(int.tp[1], int.tp[2], len=N.tp)
		cls.centers <- matrix(0, N.tp, 5)

		cls.centers[, 1] <- grid.tp^2
		cls.centers[, 2] <- sqrt(grid.tp)
		cls.centers[, 3] <- sin(2*pi*grid.tp)
		cls.centers[, 4] <- grid.tp^3
		cls.centers[, 5] <- -grid.tp^2

		curves.mat <- NULL

		fdata <- rep(list(0), N.sub)
		for (n.sub in 1:N.sub) {
			tvec <- grid.tp
			fdata[[n.sub]] <- matrix(0, N.tp, 2)
			fdata[[n.sub]][, 1] <- grid.tp
			fdata[[n.sub]][, 2] <- cls.centers[, cls.vec[n.sub]] + rnorm(N.tp, 0, sd)
			colnames(fdata[[n.sub]]) <- c("time", "val")
		}
		res <- list()
		res$data <- fdata
		res$cls <- cls.vec

		if (dformat == 2) {
			## とりあえず順次連結して作成する
			nt.vec <- numeric(N.sub)
			x.vec <- numeric(0)
			y.vec <- numeric(0)
			id.vec <- numeric(0)
			cls.vec <- numeric(0)
			wt.vec <- numeric(0)
			for (n.sub in 1:N.sub) {
				nt.vec[n.sub] <- dim(fdata[[n.sub]])[1]
				x.vec <- c(x.vec, fdata[[n.sub]][, 1])
				y.vec <- c(y.vec, fdata[[n.sub]][, 2])
				id.vec <- c(id.vec, rep(n.sub, nt.vec[n.sub]))
				wt.vec <- c(wt.vec, rep(1 / nt.vec[n.sub], length=nt.vec[n.sub]))
			}
			fdata.df.tmp <- data.frame(id=id.vec, x=x.vec, y=y.vec, wt=wt.vec)

			res$data <- fdata.df.tmp
		}

	} else {
		## irregular dataset
		fdata <- rep(list(0), N.sub)

		for (n.sub in 1:N.sub) {
	    nt <- sample(min.tp:max.tp, 1) ## 測定回数
			tvec <- sort(runif(n=nt, min=0, max=1)) ##観測ポイントの発生 U[0, 1]
			fdata[[n.sub]] <- matrix(0, nt, 2)
			fdata[[n.sub]][, 1] <- tvec
			if (cls.vec[n.sub] == 1) {
				fdata[[n.sub]][, 2] <- tvec^2
			} else if (cls.vec[n.sub] == 2) {
				fdata[[n.sub]][, 2] <- sqrt(tvec)
			} else if (cls.vec[n.sub] == 3) {
				fdata[[n.sub]][, 2] <- sin(2*pi*tvec)
			} else if (cls.vec[n.sub] == 4) {
				fdata[[n.sub]][, 2] <- tvec^3
			} else if (cls.vec[n.sub] == 5) {
				fdata[[n.sub]][, 2] <- -tvec^2
			}
			colnames(fdata[[n.sub]]) <- c("x", "y")
		}
		res <- list()
		res$data <- fdata
		res$cls <- cls.vec

		if (dformat == 2) {
			## とりあえず順次連結して作成する
			nt.vec <- numeric(N.sub)
			x.vec <- numeric(0)
			y.vec <- numeric(0)
			id.vec <- numeric(0)
			cls.vec <- numeric(0)
			wt.vec <- numeric(0)
			for (n.sub in 1:N.sub) {
				nt.vec[n.sub] <- dim(fdata[[n.sub]])[1]
				x.vec <- c(x.vec, fdata[[n.sub]][, 1])
				y.vec <- c(y.vec, fdata[[n.sub]][, 2])
				id.vec <- c(id.vec, rep(n.sub, nt.vec[n.sub]))
				wt.vec <- c(wt.vec, rep(1 / nt.vec[n.sub], length=nt.vec[n.sub]))
			}
			fdata.df.tmp <- data.frame(id=id.vec, x=x.vec, y=y.vec, wt=wt.vec)

			res$data <- fdata.df.tmp
		}

	}

	return(res)
}

