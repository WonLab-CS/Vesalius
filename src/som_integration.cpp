// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include <iostream>
#include <vector>
#include <cmath>


using namespace Rcpp;



Rcpp::NumericVector rescale(NumericVector& s, NumericVector q){
    double max_q = Rcpp::max(q);
    double min_q = Rcpp::min(q);
    double max_s = Rcpp::max(s);
    double min_s = Rcpp::min(s);
    for (int i = 0; i < q.size(); i++){
        q[i] = ((q[i] - min_q) / (max_q - min_q)) * (max_s - min_s) + min_s;
    }
    return q;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix rescale_to_seed(const NumericVector& seed_idx,
    const NumericVector& query_idx,
    const NumericMatrix& seed,
    NumericMatrix query) {
    for (int i = 0; i < seed_idx.size(); i++){
        NumericVector q = query(_,query_idx[i]);
        NumericVector s = seed(_, seed_idx[i]);
        query(_,i) = rescale(s,q);
    }
    return query;
}



