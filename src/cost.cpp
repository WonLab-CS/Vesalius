#include <Rcpp.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <unordered_set>

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
    float corr = static_cast<float>((n * sum_cell_1_2 - sum_cell_1 * sum_cell_2)
                  / sqrt((n * sq_sum_cell_1 - sum_cell_1 * sum_cell_1)
                      * (n * sq_sum_cell_2 - sum_cell_2 * sum_cell_2)));
    return corr;
}

float fast_jaccard(const CharacterVector& niche_1, const CharacterVector& niche_2){
    std::unordered_set<std::string> niche_set1(niche_1.begin(), niche_1.end());
    std::unordered_set<std::string> niche_set2(niche_2.begin(), niche_2.end());
    std::unordered_set<std::string> intersection;
    for (const std::string& element : niche_set1) {
        if (niche_set2.find(element) != niche_set2.end()) {
            intersection.insert(element);
        }
    }
    size_t union_size = niche_set1.size() + niche_2.size() - intersection.size();
    if (union_size == 0) {
        return 0.0;
    } else {
        return static_cast<float>(intersection.size()) / union_size;
    }
}

// [[Rcpp::export]]
Rcpp::NumericMatrix feature_cost(const NumericMatrix& seed,
    const NumericMatrix& query) {
    int cell_seed = seed.ncol();
    int cell_query = query.ncol();
    NumericMatrix score(cell_query, cell_seed);
    for (int i = 0 ; i < cell_seed; i++){
        for (int j = 0; j < cell_query; j++){
            NumericVector s = seed(_,i);
            NumericVector q = query(_,j);
            score(j,i) = fast_cor(s,q);
        }
    }
    return score;
}


// [[Rcpp::export]]
Rcpp::NumericMatrix compare_niche_fast(const List& seed,
    const List& query) {
    int cell_seed = seed.size();
    int cell_query = query.size();
    NumericMatrix composition(cell_query, cell_seed);
    for (int i = 0 ; i < cell_seed; i++){
        for (int j = 0; j < cell_query; j++){
            CharacterVector s = seed[i];
            CharacterVector q = query[j];
            composition(j,i) = fast_jaccard(s,q);
        }
    }
    return composition;
}




