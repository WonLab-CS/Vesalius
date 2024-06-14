#include <Rcpp.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <math.h>
#include <stdio.h>
#include <algorithm>
#include <string>
#include <unordered_map>


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




// Rcpp::CharacterVector intersection(CharacterVector& v1,
//     CharacterVector& v2){
//     Rcpp::CharacterVector v3;
//     std::sort(v1.begin(), v1.end());
//     std::sort(v2.begin(), v2.end());
//     std::set_intersection(v1.begin(),v1.end(),
//                           v2.begin(),v2.end(),
//                           std::back_inserter(v3));
//     return v3;
// }

// Moved to using std::vector instead of RcppCharacterVector. 
// Clang compiler on mac was throwing a hissy fit on this...
// which was strange since it had worked for months with no issue and 
// I don't remember making any updated to the compiler...
std::vector<std::string> intersection(std::vector<std::string>& v1, std::vector<std::string>& v2) {
    std::vector<std::string> v3;
    std::sort(v1.begin(), v1.end());
    std::sort(v2.begin(), v2.end());
    std::set_intersection(v1.begin(), v1.end(),
                          v2.begin(), v2.end(),
                          std::back_inserter(v3));
    return v3;
}


float jaccard_index(std::vector<std::string>& niche_1, std::vector<std::string>& niche_2){
    double size_niche_1 = niche_1.size();
    double size_niche_2 = niche_2.size();
    std::vector<std::string> intersect = intersection(niche_1, niche_2);
    double size_in = intersect.size();
    float jaccard_index = static_cast<float>(size_in
                        / (size_niche_1 + size_niche_2 - size_in));
    return jaccard_index;
}




// [[Rcpp::export]]
Rcpp::NumericMatrix pearson_cost(const List& seed,
    const List& query) {
    int cell_seed = seed.size();
    int cell_query = query.size();
    Rcpp::NumericMatrix score(cell_query, cell_seed);
    for (int i = 0 ; i < cell_seed; i++){
        for (int j = 0; j < cell_query; j++){
            Rcpp::NumericVector s = seed[i];
            Rcpp::NumericVector q = query[j];
            score(j,i) = fast_cor(s,q);
        }
    }
    return score;
}


// [[Rcpp::export]]
Rcpp::NumericMatrix jaccard_cost(const List& seed,
    const List& query) {
    int cell_seed = seed.size();
    int cell_query = query.size();
    NumericMatrix composition(cell_query, cell_seed);
    for (int i = 0 ; i < cell_seed; i++){
        for (int j = 0; j < cell_query; j++){
            std::vector<std::string> s = seed[i];
            std::vector<std::string> q = query[j];
            composition(j,i) = jaccard_index(s, q);
        }
    }
    return composition;
}




