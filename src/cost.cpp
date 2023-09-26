#include <Rcpp.h>
#include <iostream>
#include <vector>
#include <cmath>

using namespace Rcpp;


float fast_cor(const NumericVector& cell_1, const NumericVector& cell_2){
    int n = cell_1.size();
    double sum_cell_1 = 0.0, sum_cell_2 = 0.0, sum_cell_1_2 = 0.0;
    double sq_sum_cell_1 = 0, sq_sum_cell_2 = 0.0;
    for (int i = 0; i < n; i++){
        sum_cell_1 += cell_1[i];
        sum_cell_2 += cell_2[i];
        sum_cell_1_2 += cell_1[i] * cell_2[i];
        sq_sum_cell_1 += cell_1[i] * cell_1[i];
        sq_sum_cell_2 += cell_2[i] * cell_2[i];
    }
    float corr = (float)(n * sum_cell_1_2 - sum_cell_1 * sum_cell_2)
                  / sqrt((n * sq_sum_cell_1 - sum_cell_1 * sum_cell_1)
                      * (n * sq_sum_cell_2 - sum_cell_2 * sum_cell_2));
    return corr;
}

// [[Rcpp::export]]
NumericMatrix feature_dist_fast(const NumericMatrix seed, const NumericMatrix query) {
    int cell_seed = seed.ncol();
    int cell_query = query.ncol();
    NumericMatrix cost(cell_query, cell_seed);
    for (int i = 0 ; i < cell_seed; i++){
        for (int j = 0; j < cell_query; j++){
            NumericVector s = seed(_,i);
            NumericVector q = query(_,j);
            cost(j,i) = fast_cor(s,q);
        }
    }
    return cost;
}

