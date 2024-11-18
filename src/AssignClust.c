#include <R.h>
#include <Rinternals.h>

// Rのベクトル・行列の要素を取得するためのマクロ
#define INTEGER_ELT(x, i) (INTEGER(x)[i])
#define REAL_ELT(x, i) (REAL(x)[i])

// 関数のプロトタイプ宣言
SEXP AssignClust(SEXP id_vec, SEXP loss_mat, SEXP N_sub_r, SEXP N_clust_r);

SEXP AssignClust(SEXP id_vec, SEXP loss_mat, SEXP N_sub_r, SEXP N_clust_r) {
    // N_sub と N_clust を整数として取得
    unsigned int N_sub = INTEGER(N_sub_r)[0];  // サンプルサイズn
    unsigned int N_clust = INTEGER(N_clust_r)[0];  // クラスタの数

    // 結果を格納するベクトル
    SEXP cls_tmp = PROTECT(allocVector(INTSXP, N_sub));

    // Rのベクトル・行列データを直接参照
    double *id_vec_data = REAL(id_vec);
    double *loss_mat_data = REAL(loss_mat);

    // 各対象の測定時点数に基づいて処理を行う
    unsigned int total_measurements = LENGTH(id_vec);  // 全測定点数

    for (unsigned int n_sub = 0; n_sub < N_sub; n_sub++) {
        int n_sub_p1 = n_sub + 1;

        // id_vec内で対象IDが一致するインデックスを取得
        unsigned int start_idx = 0;
        unsigned int count = 0;

        for (unsigned int i = 0; i < total_measurements; i++) {
            if ((int)id_vec_data[i] == n_sub_p1) {
                if (count == 0) start_idx = i;  // 最初に出現する位置を記録
                count++;  // 対象n_subに対応する測定点数のカウント
            }
        }

        // loss_clust を初期化
        double *loss_clust = (double *) R_alloc(N_clust, sizeof(double));
        for (unsigned int i = 0; i < N_clust; i++) {
            loss_clust[i] = 0.0;
        }

        // 対象n_subに対応する測定点数 (count) についてクラスタごとの誤差を集計
        for (unsigned int i = 0; i < count; i++) {
            int index = start_idx + i;
            for (unsigned int n_clust = 0; n_clust < N_clust; n_clust++) {
                // Rの行列は列優先で保存されているため、インデックスを修正
                loss_clust[n_clust] += loss_mat_data[index + total_measurements * n_clust];
            }
        }

        // 最小値を持つクラスタを探す (index_min の代わりに手動で検索)
        unsigned int min_index = 0;
        for (unsigned int i = 1; i < N_clust; i++) {
            if (loss_clust[i] < loss_clust[min_index]) {
                min_index = i;
            }
        }

        // クラスタ番号を結果に格納
        INTEGER(cls_tmp)[n_sub] = min_index + 1;  // 1-based indexing
    }

    UNPROTECT(1);
    return cls_tmp;
}
