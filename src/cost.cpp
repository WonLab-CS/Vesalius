#include <Rcpp.h>
#include <RcppEigen.h>
#include <vector>
#include <cmath>
#include <math.h>
#include <stdio.h>
#include <algorithm>
#include <string>
#include <unordered_map>


using namespace Rcpp;



Eigen::VectorXd col_means(const Eigen::MatrixXd& mat) {
    return mat.colwise().mean();
}


Eigen::VectorXd col_stds(const Eigen::MatrixXd& mat, const Eigen::VectorXd& means) {
    Eigen::VectorXd stds(mat.cols());
    for (int i = 0; i < mat.cols(); i++) {
        stds[i] = sqrt((mat.col(i).array() - means[i]).square().sum() / (mat.rows() - 1));
    }
    return stds;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix pearson_cost(const Eigen::Map<Eigen::MatrixXd>& seed,
                               const Eigen::Map<Eigen::MatrixXd>& query) {
    int seed_cols = seed.cols(); // Number of columns in matA
    int query_cols = query.cols(); // Number of columns in matB
    NumericMatrix result(query_cols, seed_cols);

    // Compute means and standard deviations for both matrices
    Eigen::VectorXd mean_seed = col_means(seed);
    Eigen::VectorXd mean_query = col_means(query);
    Eigen::VectorXd std_seed = col_stds(seed, mean_seed);
    Eigen::VectorXd std_query = col_stds(query, mean_query);

    // Center the matrices
    Eigen::MatrixXd centered_seed = seed.rowwise() - mean_seed.transpose();
    Eigen::MatrixXd centered_query = query.rowwise() - mean_query.transpose();

    // Compute correlations sequentially
    for (int i = 0; i < seed_cols; i++) {
        for (int j = 0; j < query_cols; j++) {
            double numerator = (centered_seed.col(i).array() * centered_query.col(j).array()).sum();
            double denominator = (std_seed[i] * std_query[j]) * (seed.rows() - 1);
            result(j, i) = numerator / denominator;
        }
    }

    return result;
}



;


inline double euclidean_distance(const Eigen::VectorXd& vec1, const Eigen::VectorXd& vec2) {
    return (vec1 - vec2).norm();
}

// [[Rcpp::export]]
Eigen::MatrixXd distance_cost(const Eigen::MatrixXd& seed, const Eigen::MatrixXd& query) {
    size_t cols1 = seed.cols(); 
    size_t cols2 = query.cols(); 
    Eigen::MatrixXd distance_matrix(cols2, cols1);
    for (size_t i = 0; i < cols2; ++i) {  
        for (size_t j = 0; j < cols1; ++j) {  
            distance_matrix(i, j) = euclidean_distance(seed.col(j), query.col(i));
        }
    }
    return distance_matrix;
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


double jaccard_index(std::vector<std::string>& niche_1, std::vector<std::string>& niche_2){
    double size_niche_1 = niche_1.size();
    double size_niche_2 = niche_2.size();
    std::vector<std::string> intersect = intersection(niche_1, niche_2);
    double size_in = intersect.size();
    double jaccard_index = static_cast<float>(size_in
                        / (size_niche_1 + size_niche_2 - size_in));
    return jaccard_index;
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




