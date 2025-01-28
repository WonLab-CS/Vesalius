#include <RcppEigen.h>
#include <queue>
#include <unordered_set>
#include <random>
#include <cmath>

using namespace Rcpp;

// Function to compute squared Euclidean distance
inline double squared_distance(const Eigen::VectorXd& a, const Eigen::VectorXd& b) {
    return (a - b).squaredNorm();
}

// Random sampling helper
int random_sample(const std::vector<int>& vec) {
    static std::random_device rd;
    static std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis(0, vec.size() - 1);
    return vec[dis(gen)];
}

// [[Rcpp::export]]
IntegerVector bubble_stack_optimized(const Eigen::MatrixXd& coordinates, int n_centers = 500, int max_iter = 500) {
    int n = coordinates.rows();
    if (n < 4) stop("At least 4 points are required.");

    // Initial radius estimate
    double radius = 0.0;
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            radius = std::max(radius, (coordinates.row(i) - coordinates.row(j)).norm());
        }
    }
    radius /= 2.0;

    // Main loop variables
    bool convergence = false;
    int iter = 1;
    std::vector<int> background_grid;

    while (!convergence) {
        background_grid.clear();
        std::vector<bool> active(n, true);

        // Inner active point selection loop
        for (int i = 0; i < n; ++i) {
            if (!active[i]) continue;
            background_grid.push_back(i);

            // Deactivate points within the radius
            for (int j = 0; j < n; ++j) {
                if (active[j] && squared_distance(coordinates.row(i), coordinates.row(j)) <= radius * radius) {
                    active[j] = false;
                }
            }
        }

        // Check convergence
        if (background_grid.size() == n_centers) {
            convergence = true;
        } else if (iter >= max_iter) {
            convergence = true;
            warning("Max iterations reached with no convergence. Returning approximation.");
        } else {
            // Adjust radius based on the number of points
            if (background_grid.size() < n_centers) {
                radius *= 0.75;  // Reduce radius to select more points
            } else {
                radius *= 1.25;  // Increase radius to select fewer points
            }
        }

        ++iter;
    }

    return wrap(background_grid);
}
