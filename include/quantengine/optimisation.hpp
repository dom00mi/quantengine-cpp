#pragma once

#include <cstddef>
#include <functional>
#include <vector>

namespace quantengine::optimisation {

struct Minimum {
  std::vector<double> parameters;
  double objective;
  std::size_t iterations;
  bool converged;
};

// Calibrators transform
// their parameters to satisfy positivity and correlation constraints.
[[nodiscard]] Minimum nelder_mead(
    const std::function<double(const std::vector<double>&)>& objective,
    const std::vector<double>& initial, const std::vector<double>& steps,
    std::size_t max_iterations = 2000, double tolerance = 1e-10);

}  
