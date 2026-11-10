#include "quantengine/numerics.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <numeric>
#include <stdexcept>

namespace quantengine::numerics {

Quadratic fit_quadratic(const std::vector<double>& x,
                        const std::vector<double>& y) {
  if (x.empty() || x.size() != y.size())
    throw std::invalid_argument("Invalid regression samples");
  for (std::size_t i = 0; i < x.size(); ++i)
    if (!std::isfinite(x[i]) || !std::isfinite(y[i]))
      throw std::invalid_argument("Non-finite regression sample");
  const double mean = std::accumulate(y.begin(), y.end(), 0.0) /
                      static_cast<double>(y.size());
  if (x.size() < 3) return {mean, 0, 0};

  // Normal equations for a small basis; pivoting and a constant-model fallback
  // prevent unstable exercise decisions on near-singular simulated samples.
  std::array<std::array<double, 4>, 3> system{};
  for (std::size_t j = 0; j < x.size(); ++j) {
    const std::array<double, 3> basis{1.0, x[j], x[j] * x[j]};
    for (std::size_t row = 0; row < 3; ++row) {
      for (std::size_t col = 0; col < 3; ++col)
        system[row][col] += basis[row] * basis[col];
      system[row][3] += basis[row] * y[j];
    }
  }
  for (std::size_t pivot = 0; pivot < 3; ++pivot) {
    std::size_t best = pivot;
    for (std::size_t row = pivot + 1; row < 3; ++row)
      if (std::abs(system[row][pivot]) > std::abs(system[best][pivot])) best = row;
    if (std::abs(system[best][pivot]) < 1e-12 * x.size())
      return {mean, 0, 0};
    std::swap(system[pivot], system[best]);
    const double scale = system[pivot][pivot];
    for (std::size_t col = pivot; col < 4; ++col) system[pivot][col] /= scale;
    for (std::size_t row = 0; row < 3; ++row) {
      if (row == pivot) continue;
      const double multiple = system[row][pivot];
      for (std::size_t col = pivot; col < 4; ++col)
        system[row][col] -= multiple * system[pivot][col];
    }
  }
  return {system[0][3], system[1][3], system[2][3]};
}

}
