#include <Rcpp.h>
#include <RcppEigen.h>
#include <vector>
#include <cmath>
#include <math.h>
#include <stdio.h>
#include <algorithm>
#include <string>
#include <unordered_map>
#include <unordered_set>



using namespace Rcpp;



// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
Eigen::MatrixXd pearson_fast(const Eigen::MatrixXd& seed, const Eigen::MatrixXd& query) {
  int n = seed.rows();

  Eigen::VectorXd mean_seed = seed.colwise().mean();
  Eigen::VectorXd mean_query = query.colwise().mean();

  Eigen::MatrixXd centered_seed = seed.colwise() - mean_seed;
  Eigen::MatrixXd centered_query = query.colwise() - mean_query;

  Eigen::VectorXd std_seed = centered_seed.colwise().norm() / std::sqrt(n - 1);
  Eigen::VectorXd std_query = centered_query.colwise().norm() / std::sqrt(n - 1);

  Eigen::MatrixXd correlation = (centered_query.transpose() * centered_seed) 
                             .cwiseQuotient(std_query * std_seed.transpose());

  return correlation;
}



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
Rcpp::NumericMatrix pearson_exact(const Eigen::Map<Eigen::MatrixXd>& seed,
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




std::vector<std::string> intersection(const std::vector<std::string>& v1, const std::vector<std::string>& v2) {
    std::unordered_set<std::string> set_v1(v1.begin(), v1.end());
    std::vector<std::string> v3;
    for (const auto& item : v2) {
        if (set_v1.count(item)) {
            v3.push_back(item);
        }
    }
    return v3;
}


double jaccard_index(const std::vector<std::string>& niche_1, const std::vector<std::string>& niche_2) {
    std::unordered_set<std::string> set_niche_1, set_niche_2;

    // Add non-NA elements to hash sets
    for (const auto& item : niche_1) {
        if (item != "NA") set_niche_1.insert(item);
    }
    for (const auto& item : niche_2) {
        if (item != "NA") set_niche_2.insert(item);
    }

    // Compute intersection size
    size_t intersection_size = 0;
    for (const auto& item : set_niche_2) {
        if (set_niche_1.count(item)) {
            ++intersection_size;
        }
    }

    // Compute union size
    size_t union_size = set_niche_1.size() + set_niche_2.size() - intersection_size;

    // Return Jaccard index
    return (union_size > 0) ? static_cast<double>(intersection_size) / union_size : 0.0;
}



// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
Eigen::MatrixXd jaccard_cost(const Rcpp::CharacterMatrix& seed,
                             const Rcpp::CharacterMatrix& query) {
    int cell_seed = seed.ncol();
    int cell_query = query.ncol();
    Eigen::MatrixXd composition(cell_query, cell_seed);

    // Process each pair of columns
    for (int i = 0; i < cell_seed; ++i) {
        for (int j = 0; j < cell_query; ++j) {
            std::vector<std::string> s(seed(_, i).begin(), seed(_, i).end());
            std::vector<std::string> q(query(_, j).begin(), query(_, j).end());
            composition(j, i) = jaccard_index(s, q);
        }
    }

    return composition;
}



