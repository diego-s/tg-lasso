#include <R.h>
#include <R_ext/BLAS.h>
#include <Rdefines.h>
#include <math.h>

// low level functions

int inline compute_backtracking_criterion(double g_W_c, double g_V_c, 
    double* g_gradient, double* W_c, double* V_c, double step, double* 
    delta_W_c, int M, int T) {
    for (int i = 0; i < M * T; i++) delta_W_c[i] = W_c[i] - V_c[i];
    int K = M * T;
    int one = 1;
    double right = g_V_c;
    right += F77_CALL(ddot)(&K, g_gradient, &one, delta_W_c, &one);
    //right += (1 / (2 * step)) * pow(F77_CALL(dnrm2)(&K, delta_W_c, &one), 2);
    right += 0.5 * pow(F77_CALL(dnrm2)(&K, delta_W_c, &one), 2);
    return (g_W_c <= right) ? 1 : 0;
}

void inline compute_dg_dW(double* dg_dW, double* R, double* X_c, double* W_c, 
    double lambda_1, double lambda_2, double* intervals, int N, int M, int T) {
    double one = 1;
    double zero = 0;
    F77_CALL(dgemm)("T", "N", &M, &T, &N, &one, X_c, &N, R, &N, &zero, dg_dW, 
        &M);
    int k;
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < T; j++) {
            k = i + (j * M);
            dg_dW[k] += lambda_1 * W_c[k];
            if (j > 0) dg_dW[k] -= lambda_2 * (1 / intervals[j - 1]) * 
                (W_c[k - M] - W_c[k]);
            if (j < (T - 1)) dg_dW[k] += lambda_2 * (1 / intervals[j]) * 
                (W_c[k] - W_c[k + M]);
        }
    }
}

double inline compute_difference(double* W_c, double* V_c, int M, int T) {
    double numerator_sum = 0;
    double denominator_sum = 0;
    for (int i = 0; i < M * T; i++) {
        numerator_sum += (W_c[i] - V_c[i]) * (W_c[i] - V_c[i]);
        denominator_sum += W_c[i] * W_c[i];
    }
    return (denominator_sum > 1) ? sqrt(numerator_sum) / sqrt(denominator_sum) :
        sqrt(numerator_sum);
}

double inline compute_g(double* R, double* W_c, double lambda_1, 
    double lambda_2, double* intervals, int N, int M, int T) {
    double g = 0;
    for (int i = 0; i < N * T; i++) g += R[i] * R[i];
    int k;
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < T; j++) {
            k = i + (j * M);
            // TODO: Check if using squared intervals causes big difference
            g += lambda_1 * W_c[k] * W_c[k];
            if (j < (T - 1)) g += lambda_2 * (1 / intervals[j]) * (W_c[k] - 
                W_c[k + M]) * (W_c[k] - W_c[k + M]);
        }
    }
    return g;
}

void inline compute_R(double* R, double* Y, double* Y_hat, int N, int T, 
    int allow_missing_Y) {
    for (int i = 0; i < N * T; i++) {
        R[i] = Y_hat[i] - Y[i];
        if (allow_missing_Y && ISNA(R[i])) R[i] = 0;
    }
}

void inline compute_V_c(double* V_c, double* W_c, double* W_c_old, int M, int T,
    int i) {
    for (int j = 0; j < M * T; j++) V_c[j] = W_c[j] + ((i - 2) / (i + 1)) * 
        (W_c[j] - W_c_old[j]);
}

void inline compute_Y_hat(double* Y_hat, double* X_c, double* W_c, int N, int M,
    int T) {
    double one = 1;
    double zero = 0;
    F77_CALL(dgemm)("N", "N", &N, &T, &M, &one, X_c, &N, W_c, &M, &zero, Y_hat, 
        &N);
}

void inline compute_Z_c(double* Z_c, double* V_c, double step, double* 
    g_gradient, int M, int T) {
    for (int i = 0; i < M * T; i++) Z_c[i] = V_c[i] - (step * g_gradient[i]);
}

void inline initialize_matrix(double* M, int k, int l, double v) {
    for (int i = 0; i < k * l; i++) M[i] = v;
}

void inline prox_l21(double* W_c, double lambda, double* Z_c, int M, int T) {
    double sum;
    double norm;
    double factor;
    int k;
    for (int i = 0; i < M; i++) {
        sum = 0;
        for (int j = 0; j < T; j++) {
            k = i + (j * M);
            sum += Z_c[k] * Z_c[k];
        }
        //norm = F77_CALL(dnrm2)(&K, Z_c + i, &M);
        norm = sqrt(sum);
        for (int j = 0; j < T; j++) {
            k = i + (j * M);
            factor = ((1 - (lambda / norm)) > 0) ? (1 - (lambda / norm)) : 0;
            //printf("f: %f.\n", factor);
            W_c[k] = Z_c[k] * factor;
        }
    }
}

// high level functions

double inline compute_g_W(double* X_c, double* W_c, double* Y_c, double* Y_hat, 
    double* R, double lambda_1, double lambda_2, double* intervals, int N, 
    int M, int T, int allow_missing_Y) {
    compute_Y_hat(Y_hat, X_c, W_c, N, M, T);
    compute_R(R, Y_c, Y_hat, N, T, allow_missing_Y);
    double g = compute_g(R, W_c, lambda_1, lambda_2, intervals, N, M, T);
    return g;
}

SEXP fit_tg_lasso(SEXP X_c_exp, SEXP Y_c_exp, SEXP W_c_exp, SEXP lambda_1_exp, 
    SEXP lambda_2_exp, SEXP lambda_3_exp, SEXP intervals_exp, 
    SEXP allow_missing_Y_exp, SEXP min_difference_exp) {
    // declarations
    int N = nrows(X_c_exp);
    int M = ncols(X_c_exp);
    int T = ncols(Y_c_exp);
    SEXP Y_hat_exp = PROTECT(allocMatrix(REALSXP, N, T));
    SEXP R_exp = PROTECT(allocMatrix(REALSXP, N, T));
    SEXP W_c_old_exp = PROTECT(allocMatrix(REALSXP, M, T));
    SEXP V_c_exp = PROTECT(allocMatrix(REALSXP, M, T));
    SEXP delta_W_c_exp = PROTECT(allocMatrix(REALSXP, M, T));
    SEXP g_gradient_exp = PROTECT(allocMatrix(REALSXP, M, T));
    SEXP Z_c_exp = PROTECT(allocMatrix(REALSXP, M, T));
    double* X_c = REAL(X_c_exp);
    double* Y_c = REAL(Y_c_exp);
    double lambda_1 = REAL(lambda_1_exp)[0];
    double lambda_2 = REAL(lambda_2_exp)[0];
    double lambda_3 = REAL(lambda_3_exp)[0];
    double* intervals = REAL(intervals_exp);
    int allow_missing_Y = LOGICAL(allow_missing_Y_exp)[0];
    double min_difference = REAL(min_difference_exp)[0];
    double* Y_hat = REAL(Y_hat_exp);
    double* R = REAL(R_exp);
    double* W_c = REAL(W_c_exp);
    double* W_c_old = REAL(W_c_old_exp);
    double* V_c = REAL(V_c_exp);
    double* delta_W_c = REAL(delta_W_c_exp);
    double g_W_c;
    double g_V_c;
    double* g_gradient = REAL(g_gradient_exp);
    double* Z_c = REAL(Z_c_exp);
    int max_iterations = 2000;
    int max_line_searches = 10000;
    int backtracking_criterion;
    double difference;
    // algorithm body
    copyMatrix(W_c_old_exp, W_c_exp, 0);
    g_W_c = compute_g_W(X_c, W_c, Y_c, Y_hat, R, lambda_1, lambda_2, intervals, 
        N, M, T, allow_missing_Y);
    double step = 1;
    for (int i = 1; i <= max_iterations; i++) {
        compute_V_c(V_c, W_c, W_c_old, M, T, i);
        compute_Y_hat(Y_hat, X_c, W_c, N, M, T);
        compute_R(R, Y_c, Y_hat, N, T, allow_missing_Y);
        g_V_c = compute_g(R, W_c, lambda_1, lambda_2, intervals, N, M, T);
        compute_dg_dW(g_gradient, R, X_c, V_c, lambda_1, lambda_2, intervals, N,
            M, T);
        copyMatrix(W_c_old_exp, W_c_exp, 0);
        for (int j = 0; j < max_line_searches; j++) {
            compute_Z_c(Z_c, V_c, step, g_gradient, M, T);
            prox_l21(W_c, step * lambda_3, Z_c, M, T);
            g_W_c = compute_g_W(X_c, W_c, Y_c, Y_hat, R, lambda_1, lambda_2, 
                intervals, N, M, T, allow_missing_Y);
            backtracking_criterion = compute_backtracking_criterion(g_W_c, 
                g_V_c, g_gradient, W_c, V_c, step, delta_W_c, M, T);
            if (backtracking_criterion) break;
            else step *= 0.1;
        }
        difference = compute_difference(W_c, V_c, M, T);
        if (difference < min_difference) break;
    }
    UNPROTECT(7);
    return(W_c_exp);
}
