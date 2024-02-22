// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include <iostream>
#include <vector>
#include <cmath>
<<<<<<< HEAD
#include <algorithm>
#include <string>
=======
#include <math.h>
#include <algorithm>
#include <string>
#include <unordered_map>
>>>>>>> develop


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


<<<<<<< HEAD
Rcpp::CharacterVector intersection(CharacterVector v1,
    CharacterVector v2){
    Rcpp::CharacterVector v3;
    std::sort(v1.begin(), v1.end());
    std::sort(v2.begin(), v2.end());
    std::set_intersection(v1.begin(),v1.end(),
                          v2.begin(),v2.end(),
                          std::back_inserter(v3));
    return v3;
}


=======

double RBO(CharacterVector& niche_1, CharacterVector& niche_2 , double p){
    Rcpp::CharacterVector short_list;
    Rcpp::CharacterVector long_list;
    int s;
    int l;
    double rbo = 0;
    bool even;
    // setting list length
    int lenN1 = niche_1.size();
    int lenN2 = niche_2.size();
    if(lenN1 == lenN2){
        even = true ;
        short_list = niche_1;
        long_list = niche_2;
    } else if (lenN1 < lenN2){
        even = false;
        short_list = niche_1;
        long_list = niche_2;
    } else {
        even = false;
        short_list = niche_2;
        long_list = niche_1;
    }
    
    s = short_list.size();
    l = long_list.size();
    // overlap
    // initalize overlap vector
    Rcpp::NumericVector overlap(l);
    for(int i = 0; i < l;i++){
        Rcpp::CharacterVector g1(i + 1);
        Rcpp::CharacterVector g2(i + 1);
        Rcpp::CharacterVector tmp;
        for(int j=0; j < g1.size() ;j++){
            g1[j] = long_list[j];
            int at = std::min((s -1),j);
            if ( (s - 1) < j) {
                g2[j] = "";
            } else {
                g2[j] = short_list[at];
            }
            
        }
        std::sort(g1.begin(), g1.end());
        std::sort(g2.begin(), g2.end());
        std::set_intersection(
            g1.begin(),g1.end(),
            g2.begin(),g2.end(),
            std::back_inserter(tmp));
        double size = tmp.size();
        overlap[i] = size;
    }
    // comp RBO
    
    if(!even){
        double odd_sum_long = 0;
        double odd_sum_diff = 0;
        for(int i = 0 ; i < overlap.size();i++){
            odd_sum_long += overlap[i] / (i + 1) * pow(p,(i + 1));
        }

        for(int i = s; i < overlap.size();i++){
            odd_sum_diff += (overlap[(s -1)] * ((i + 1 ) - s) / (s * (i + 1)) * pow(p,(i+1)));
        }

        rbo = ((1-p)/p) * 
            ((odd_sum_long) + 
            (odd_sum_diff)) +
            ((overlap.at((l - 1)) - overlap.at((s - 1))) / l + (overlap.at((s - 1)) / s)) *
            pow(p,l);

    } else {
        double es =0;
        int last = overlap.at(l - 1);
        for(int i=0; i < overlap.size(); i++){
            es += (overlap[i] / (i+1)) * pow(p,(i+1));
        }
        rbo = (last/l) * pow(p,l) + (((1-p)/p) * es);
    }
    return rbo;
}



Rcpp::CharacterVector intersection(CharacterVector v1,
    CharacterVector v2){
    Rcpp::CharacterVector v3;
    std::sort(v1.begin(), v1.end());
    std::sort(v2.begin(), v2.end());
    std::set_intersection(v1.begin(),v1.end(),
                          v2.begin(),v2.end(),
                          std::back_inserter(v3));
    return v3;
}


>>>>>>> develop
float jaccard_index(const CharacterVector& niche_1, const CharacterVector& niche_2){
    double size_niche_1 = niche_1.size();
    double size_niche_2 = niche_2.size();
    Rcpp::CharacterVector intersect = intersection(niche_1, niche_2);
    double size_in = intersect.size();
    float jaccard_index = static_cast<float>(size_in
                        / (size_niche_1 + size_niche_2 - size_in));
    return jaccard_index;
}

<<<<<<< HEAD
// float fast_jaccard(){
//     std::unordered_set<String> niche_set1(niche_1.begin(), niche_1.end());
//     std::unordered_set<String> niche_set2(niche_2.begin(), niche_2.end());
//     std::unordered_set<String> intersection;
    
//     size_t union_size = niche_set1.size() + niche_2.size() - intersection.size();
//     if (union_size == 0) {
//         return 0.0;
//     } else {
//         return static_cast<float>(intersection.size()) / union_size;
//     }
// }

// [[Rcpp::export]]
Rcpp::NumericMatrix feature_cost(const NumericMatrix& seed,
    const NumericMatrix& query) {
    int cell_seed = seed.ncol();
    int cell_query = query.ncol();
    Rcpp::NumericMatrix score(cell_query, cell_seed);
    for (int i = 0 ; i < cell_seed; i++){
        for (int j = 0; j < cell_query; j++){
            Rcpp::NumericVector s = seed(_,i);
            Rcpp::NumericVector q = query(_,j);
=======



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
>>>>>>> develop
            score(j,i) = fast_cor(s,q);
        }
    }
    return score;
}


// [[Rcpp::export]]
<<<<<<< HEAD
Rcpp::NumericMatrix compare_niche_fast(const List& seed,
=======
Rcpp::NumericMatrix jaccard_cost(const List& seed,
>>>>>>> develop
    const List& query) {
    int cell_seed = seed.size();
    int cell_query = query.size();
    NumericMatrix composition(cell_query, cell_seed);
    for (int i = 0 ; i < cell_seed; i++){
        for (int j = 0; j < cell_query; j++){
<<<<<<< HEAD
            CharacterVector s = seed[i];
            CharacterVector q = query[j];
            composition(j,i) = jaccard_index(s,q);
=======
            Rcpp::CharacterVector s = seed[i];
            Rcpp::CharacterVector q = query[j];
            composition(j,i) = jaccard_index(s, q);
>>>>>>> develop
        }
    }
    return composition;
}




