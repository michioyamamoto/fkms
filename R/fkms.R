##### 引数 #####
## *fdata: データを指定する．list形式もしくは縦長のデータフレーム形式である必要がある．list形式の場合は、listの各要素がtime points\eqn{\times2} の行列で，1列目が測定時点を，2列目が測定値を表す．縦長のデータフレーム形式では，1列目に対象の番号を，2列目に測定時点を，3列目に測定値を格納しておく．さらに，4列目に各対象の各測定値に対して重みを指定することも可能である．4列目がない場合は重みは自動的に生成される．
## *dformat: data format; 1=list format, 2=data.frame format.
## *N.clust: クラスタ数
## *N.bs: the number of basis functions．Fourier基底の場合奇数で指定すること．偶数の場合、+1する値をN.bsとして用いる．
## *N.random: 初期値のランダムスタートの数
## *N.rep: アルゴリズムの反復の最大値
## *ini.cls: 初期クラスタベクトル．指定しない場合はランダムに割り当てる．
## *show.random.ite: TRUEならば実行中のランダムスタートの番号を出力する．
## *basis: 基底関数をしていする．"B"=B-spline basis, "F"=Fourier basis
## *lambda: 平滑化パラメータを指定する．現在は，basis="B"の場合にのみ有効である．


#### 戻り値 ####
## *cluster: 推定された各対象のクラスタを格納したベクトル
## *centers: 推定されたクラスタ中心関数からなる行列
## *x.centers: centersに格納されている関数の評価されている時点
## *res.fit: 各クラスタ中心関数の推定された回帰係数
## *n.opt: N.random個の初期値のうち，最適だった番号
## *mse: mean squared error



##### 依存関数 #####
# library(splines) ## bs関数のため
# library(mclust) ## classError関数のため
# library(Rcpp)
# library(RcppArmadillo)
# library(inline)
# library(fda) ## bsplinepen関数のため


# ## RcppArmadilloでの関数（srcフォルダ以下に別ファイルとして定義した 2024/09/26）
# src <- '
#  using namespace arma;

#  vec id_vec = Rcpp::as<vec>(id);
#  mat loss_mat = Rcpp::as<mat>(loss);
#  int N_sub = Rcpp::as<int>(N_sub_r);
#  int N_clust = Rcpp::as<int>(N_clust_r);

#  uvec cls_tmp(N_sub, fill::zeros);
 
#  double val=0.0;

#  for (uword n_sub = 0; n_sub < N_sub; n_sub++) { // If you want to access matrix or vector elements, define them in uword.
#   int n_sub_p1 = n_sub + 1;
#   uvec indices = find(id_vec == n_sub_p1);
#   vec loss_clust(N_clust, fill::zeros);
#   mat loss_clust_temp = loss_mat.rows(indices);

#   for (uword n_clust = 0; n_clust < N_clust; n_clust++) {
#    loss_clust(n_clust) = sum(loss_clust_temp.col(n_clust));
#   }

#   cls_tmp[n_sub] = index_min(loss_clust) + 1; // index_minは0が起点なので1を足しておく
#  }

#  return Rcpp::wrap(cls_tmp);
# '
# AssignClust <- cxxfunction(signature(id="int", loss="numeric", N_sub_r="int", N_clust_r="int"), src, plugin="RcppArmadillo")


## .cppファイルからの関数の読み込み（2024/09/26）
# Rcpp::sourceCpp("fkms_package/fkms/src/AssignClust.cpp")

# AssignClust_Cpp <- function(id_vec, loss_mat, N_sub, N_clust) {
#     .Call("AssignClust", 
#     as.double(id_vec), 
#     as.double(loss_mat), 
#     as.double(N_sub), 
#     as.double(N_clust), 
#     pacakge="fkms")
# }


## .cファイルからの関数の読み込み（2024/09/26）
# dyn.load("fkms_package/fkms/src/AssignClust.dll")
AssignClust_C <- function(id_vec, loss_mat, N_sub, N_clust) {
    .Call("AssignClust", 
    as.numeric(id_vec), 
    as.numeric(loss_mat), 
    as.integer(N_sub), 
    as.integer(N_clust),
    package="fkms")
}


##### 関数本体 #####
# fdata <- fdata; dformat <- 2; N.clust <- 3; N.bs <- 10; N.random <- 100; N.rep <- 100; ini.cls <- NA; show.random.ite <- FALSE; basis <- "B"; lambda <- 0.1
fkms <- function(fdata, dformat=1, N.clust=2, N.bs=NULL, N.random=1, N.rep=100, ini.cls=NA, show.random.ite=FALSE, basis="B", lambda=0){

  ##------------------------
  ##    データの整形
  ##------------------------
  if (dformat == 1) {
    ## データがlist形式の場合、data.frame形式に変更する
    N.sub <- length(fdata)
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
  } else if (dformat == 2) {
    if (dim(fdata)[2] < 4) {
      ## wt変数がない場合は作成しておく
      id.vec <- fdata[, 1]
      id.unique <- unique(id.vec)
      N.sub <- length(id.unique)

      wt.vec <- numeric(0)
      for (n.sub in 1:N.sub) {
        nt.sub <- length(which(id.vec == id.unique[n.sub]))
        wt.vec <- c(wt.vec, rep(1 / nt.sub, len=nt.sub))
      }
      fdata.df.tmp <- data.frame(id=fdata[, 1], x=fdata[, 2], y=fdata[, 3], wt=wt.vec)
    } else {
      ## wt変数がfdataに用意されている場合
      fdata.df.tmp <- data.frame(id=fdata[, 1], x=fdata[, 2], y=fdata[, 3], wt=fdata[, 4])
    }
  } else {
    stop("Error: check the argument 'dformat' in fkms.")
  }

  if (is.null(N.bs))
      stop("Error: check the argument 'N.bs' in fkms.")
  
  if (basis == "B") {
    N.knots <- N.bs - 3
  } else if (basis == "F"){
    ## 偶数の場合は奇数に変換しておく
    if (N.bs %% 2 == 0)
      N.bs <- N.bs + 1
  } else {
      stop("Error: check the argument 'basis' in fkms.")
  }


  ## 以下で必要なオブジェクトの生成
  N.df.row <- dim(fdata.df.tmp)[1] ## データフレームの行数
  id.vec <- fdata.df.tmp[, 1]
  id.unique <- unique(id.vec)
  x.vec <- fdata.df.tmp[, 2]
  N.sub <- length(id.unique)
  nt.vec <- numeric(N.sub)
  for (n.sub in 1:N.sub)
    nt.vec[n.sub] <- length(which(id.vec == id.unique[n.sub]))


  ##------------------------
  ##    基底関数の評価
  ##------------------------

  ## B-スプライン基底もしくはフーリエ基底の生成
  ## サンプル全体でのユニークな時点で基底を評価しておく必要あり
  eval.x <- sort(unique(fdata.df.tmp$x))

  if (basis == "B") {
    ## B-スプライン基底の生成
    S.all <- bs(eval.x, df=N.bs)
  } else if (basis == "F") {
    ## フーリエ基底の生成
    S.all <- fourier(x=eval.x, nbasis=N.bs)
  }
  S.list <- rep(list(0), N.sub)
  S.all.mat <- NULL

  for (n.sub in 1:N.sub) {
    # temp.x <- fdata[[n.sub]][, 1]
    temp.x <- x.vec[which(id.vec == id.unique[n.sub])]
    N.t.i <- length(temp.x)
    S.i <- matrix(0, N.t.i, N.bs)
    for (n.t.i in 1:N.t.i)
      S.i[n.t.i, ] <- S.all[which(eval.x == temp.x[n.t.i]),]
    S.list[[n.sub]] <- S.i
    S.all.mat <- rbind(S.all.mat, S.i)
  }

  colnames(S.all.mat) <- paste("s", 1:N.bs, sep="")
  fdata.df.tmp <- data.frame(fdata.df.tmp, S.all.mat)


  
  if (basis == "B" && lambda > 0) {
    basis.obj <- create.bspline.basis(rangeval=c(min(eval.x),max(eval.x)), nbasis=N.bs)
    Lfdobj <- int2Lfd(2)  # 2階導関数を指定
    P.mat <- bsplinepen(basis.obj, Lfdobj)
  } else {
    ## 一応，フーリエ基底の場合にも2階導関数の値を0行列で作成しておく（後々正則化回帰で実装するため）
    P.mat <- matrix(0, N.bs, N.bs)
  }



  ##-------------------------------------
  ## クラスタの初期値をランダムスタート
  ##-------------------------------------
  mse.res <- Inf ## クラスタ内MSEの値の初期値
  N.opt <- 1 ## 最適解を与えたn.randomの値
  fdata.df <- fdata.df.tmp ## データフレームの受け渡し
  

  ##------------------
  ##  random start
  ##------------------
  ## n.random <- 1
  for (n.random in 1:N.random) {

    if (show.random.ite)
      if (n.random %% (N.random / 10) == 0)
        cat(paste(n.random, "; ", sep=""))

    ##------------------------
    ##  クラスタの初期値設定
    ##------------------------
    ## もし ini.cls が指定されていればそれを初期値に含めて実行
    if (n.random == 1) {
      if (!all(is.na(ini.cls))) {
        tmp.cls <- matrix(0, N.sub, N.clust)
        for(k in 1:N.clust){
          tmp.cls[which(ini.cls==k), k] <- 1
        }
      } else {
        tmp.cls <- t(rmultinom(n=N.sub, size=1, prob=rep(1 / N.clust, N.clust))) ## 等確率で多項分布に従って生成
      }
    } else {
      tmp.cls <- t(rmultinom(n=N.sub, size=1, prob=rep(1 / N.clust, N.clust))) ## 等確率で多項分布に従って生成
    }

    cls.tmp <- which(tmp.cls==1, arr.ind=TRUE)
    cls.tmp.sort <- cls.tmp[order(cls.tmp[,"row"]),]
    cls.prev <- cls.tmp.sort[, 2] ## 終了判定用のオブジェクト

    ## 全ての対象が同じ時点数の場合は特別な処理が必要
    if (all(nt.vec == nt.vec[1])) {
      ## 全ての対象が同じ時点数の場合
      cls.vec <- c(apply(cls.tmp.sort, 1, FUN=function(x, nt.vec){rep(x[2], nt.vec[x[1]])}, nt.vec=nt.vec))
    } else {
      ## 異なる時点数が存在する場合
      cls.vec <- unlist(apply(cls.tmp.sort, 1, FUN=function(x, nt.vec){rep(x[2], nt.vec[x[1]])}, nt.vec=nt.vec))
    }

    fdata.df$cls <- cls.vec

    for(n.rep in 1:N.rep){
      ## print(paste(n.rep, "回目"))

      ##-----------------
      ##  平均関数の推定
      ##-----------------
      ## とりあえず，平均関数の推定は各回1度だけとする
      ## 重み付き最小二乗推定

      if (basis == "B") {
        res.fit.list <- list()
        pred.mat <- matrix(0, N.df.row, N.clust)
        loss.mat <- matrix(0, N.df.row, N.clust)

        for (k in 1:N.clust) {
          fdata.df.tmp <- subset(fdata.df, cls==k)

          ## 目的変数ベクトルと基底関数行列，ペナルティ行列構築から正則化回帰の実行
          y.tmp <- fdata.df.tmp$y
          S.tmp <- as.matrix(fdata.df.tmp[, paste0("s", 1:N.bs)])
          W.mat <- diag(fdata.df.tmp$wt)
          # S.2d.tmp <- as.matrix(fdata.df.tmp[, paste0("s2d", 1:N.bs)])

          ## 零列の存在の有無によって場合分け
          zero.cols <- which(colSums(abs(S.tmp) == 0) == nrow(S.tmp))
          if (length(zero.cols) == 0) {
            ## ペナルティ行列の構築
            # P.mat <- t(S.2d.tmp) %*% S.2d.tmp

            ## 正則化回帰の実行
            beta.hat <- solve(t(S.tmp) %*% W.mat %*% S.tmp + lambda * P.mat) %*% t(S.tmp) %*% W.mat %*% y.tmp

          } else {
            ## 基底関数行列の構築
            S.tmp2 <- S.tmp[, -zero.cols]
            # S.2d.tmp2 <- S.2d.tmp[, -zero.cols]

            ## ペナルティ行列の構築
            P.mat.tmp <- P.mat[-zero.cols, -zero.cols]
            # P.mat <- t(S.2d.tmp2) %*% S.2d.tmp2

            ## 正則化回帰の実行
            beta.hat.tmp <- solve(t(S.tmp2) %*% W.mat %*% S.tmp2 + lambda * P.mat.tmp) %*% t(S.tmp2) %*% W.mat %*% y.tmp
            beta.hat <- rep(0, N.bs)
            beta.hat[-zero.cols] <- beta.hat.tmp
          }

          res.fit.list <- append(res.fit.list, list(beta.hat))

          ## 対象全体での予測値の計算
          pred.mat[, k] <- as.matrix(fdata.df[, paste0("s", 1:N.bs)]) %*% beta.hat

          ## 対象全体でのlossの計算
          loss.mat[, k] <- (fdata.df$y - pred.mat[, k])^2
        }

        # fmla <- as.formula(paste("y ~ 0 +", paste(paste0("s", 1:N.bs), collapse="+")))  
        # res.fit.list <- list()

        # for (k in 1:N.clust) {
        #   fdata.df.tmp <- subset(fdata.df, cls==k)
        #   res.fit.tmp <- lm(fmla, data=fdata.df.tmp, weights=wt)
        #   res.fit.list <- append(res.fit.list, list(res.fit.tmp))
        # }

        # ##  対象全体での各クラスタでの予測値を求める
        # pred.mat <- matrix(0, N.df.row, N.clust)
        # loss.mat <- matrix(0, N.df.row, N.clust)
        # for (k in 1:N.clust) {
        #   pred.mat[, k] <- predict(res.fit.list[[k]], data.frame(fdata.df[, paste("s", 1:N.bs, sep="")]))
        #   loss.mat[, k] <- (fdata.df$y - pred.mat[, k])^2
        # }

      } else if (basis == "F") {
        fmla <- as.formula(paste("y ~ 0 +", paste(paste0("s", 1:N.bs), collapse="+")))
        res.fit.list <- list()

        for (k in 1:N.clust) {
          fdata.df.tmp <- subset(fdata.df, cls==k)
          res.fit.tmp <- lm(fmla, data=fdata.df.tmp, weights=wt)
          res.fit.list <- append(res.fit.list, list(res.fit.tmp))
        }

        ##  対象全体での各クラスタでの予測値を求める
        pred.mat <- matrix(0, N.df.row, N.clust)
        loss.mat <- matrix(0, N.df.row, N.clust)
        for (k in 1:N.clust) {
          pred.mat[, k] <- predict(res.fit.list[[k]], data.frame(fdata.df[, paste("s", 1:N.bs, sep="")]))
          loss.mat[, k] <- (fdata.df$y - pred.mat[, k])^2
        }
      }

      ##-----------------
      ## 各対象の割り当て
      ##-----------------
      cls.sub <- c(AssignClust_C(id_vec=fdata.df$id, loss_mat=loss.mat, N_sub=N.sub, N_clust=N.clust))
      # cls.sub <- c(AssignClust(id=fdata.df$id, loss=loss.mat, N_sub_r=N.sub, N_clust_r=N.clust))
      cls.tmp.sort <- cbind(1:N.sub, cls.sub)

      ## **全ての対象が同じ時点数の場合は特別な処理が必要**
      if (all(nt.vec == nt.vec[1])) {
        ## 全ての対象が同じ時点数の場合
        cls.vec <- c(apply(cls.tmp.sort, 1, FUN=function(x, nt.vec){rep(x[2], nt.vec[x[1]])}, nt.vec=nt.vec))
      } else {
        ## 異なる時点数が存在する場合
        cls.vec <- unlist(apply(cls.tmp.sort, 1, FUN=function(x, nt.vec){rep(x[2], nt.vec[x[1]])}, nt.vec=nt.vec))
      }

      ## 同じクラスタに割り当てられてしまった場合の対処
      ## ひとつ前のループの結果を最終結果として終了する
      ## そもそも基底の数よりもそのクラスタの観測数が少なかった場合も終了する
      if (length(table(cls.vec)) < N.clust) {
        cls.sub <- cls.prev
        break
      }

      nt.cls.vec <- table(cls.vec)
      if (sum(nt.cls.vec < N.bs) > 0) {
        cls.sub <- cls.prev
        break
      }


      ## 更新されたクラスタ情報を付加
      fdata.df$cls <- cls.vec

      ## 更新がなければ終了
      diff.prev <- classError(cls.sub, cls.prev)$errorRate
      if (diff.prev == 0) break
      cls.prev <- cls.sub

    } ## END of N.rep loop


    ##----------------------
    ## クラスタ内MSEの評価
    ##----------------------
    ## クラスタごとの平均関数と所属データの誤差二乗和を求める
    ## basis="B"と"F"で分けて計算しておく（2024/08/12）
    ## -- "B"では正則化回帰を利用しているため
    ## -- "F"ではpredict関数を利用したいので
    
    mse.vec <- numeric(N.clust)
    if (basis == "B") {
      for (k in 1:N.clust) {
        fdata.df.tmp <- subset(fdata.df, cls==k)
        y.tmp <- fdata.df.tmp$y
        x.tmp <- fdata.df.tmp$x
        nt.tmp <- dim(fdata.df.tmp)[1]

        y.pred <- as.matrix(fdata.df.tmp[, paste0("s", 1:N.bs)]) %*% res.fit.list[[k]]
        mse.vec[k] <- sum((y.tmp - y.pred)^2) / nt.tmp
      }
    } else if (basis == "F") {
      for (k in 1:N.clust) {
        fdata.df.tmp <- subset(fdata.df, cls==k)
        y.tmp <- fdata.df.tmp$y
        x.tmp <- fdata.df.tmp$x
        nt.tmp <- dim(fdata.df.tmp)[1]
        y.pred <- predict(res.fit.list[[k]], data.frame(fdata.df.tmp[, paste("s", 1:N.bs, sep="")]))
        mse.vec[k] <- sum((y.tmp - y.pred)^2) / nt.tmp
      }
    }

    mse.tmp <- sum(mse.vec)

    ##----------------------
    ## クラスタ内MSEが最小化か？
    ##----------------------
    if (!is.na(mse.tmp)) {
      if (mse.res > mse.tmp) {
        mse.res <- mse.tmp
        cls.res <- cls.sub
        res.fit.list.res <- res.fit.list
        n.opt <- n.random
      }
    }

  } ## "N.random" loop End


  if (show.random.ite)
    cat("\n")

  ## クラスタ中心関数を計算する
  mean.func.mat <- matrix(0, length(eval.x), N.clust)
  colnames(S.all) <- paste("s", 1:N.bs, sep="")
  for (n.clust in 1:N.clust) 
    mean.func.mat[, n.clust] <- as.matrix(S.all) %*% res.fit.list.res[[n.clust]]
  
  return(list(cluster=cls.res, centers=mean.func.mat, x.centers=eval.x, res.fit=res.fit.list.res, n.opt=n.opt, mse=mse.res))
}


