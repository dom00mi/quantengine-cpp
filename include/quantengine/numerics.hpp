#pragma once

#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>

namespace quantengine::numerics {

template <class F>
[[nodiscard]] double bisection(F&& function, double lo, double hi,
                               double tolerance = 1e-12,
                               std::size_t max_iterations = 250) {
  if (!std::isfinite(lo) || !std::isfinite(hi) || lo >= hi ||
      !std::isfinite(tolerance) || tolerance <= 0 || max_iterations == 0)
    throw std::invalid_argument("Invalid root-search interval");
  double left_value = function(lo);
  double right_value = function(hi);
  if (!std::isfinite(left_value) || !std::isfinite(right_value))
    throw std::domain_error("Non-finite endpoint value");
  if (left_value == 0) return lo;
  if (right_value == 0) return hi;
  if (std::signbit(left_value) == std::signbit(right_value))
    throw std::domain_error("Root must be bracketed");
  for (std::size_t i = 0; i < max_iterations; ++i) {
    const double mid = lo + 0.5 * (hi - lo);
    const double val = function(mid);
    if (!std::isfinite(val)) throw std::domain_error("Non-finite root residual");
    if (std::abs(val) <= tolerance || hi - lo <= tolerance ||
        mid == lo || mid == hi)
      return mid;
    if (std::signbit(val) == std::signbit(left_value)) {
      lo = mid;
      left_value = val;
    } else {
      hi = mid;
    }
  }
  throw std::runtime_error("Root solver did not converge");
}

// Least-squares continuation: y ~= a + b*x + c*x*x.
struct Quadratic {
  double a = 0.0;
  double b = 0.0;
  double c = 0.0;
  [[nodiscard]] double operator()(double x) const {
    return a + b * x + c * x * x;
  }
};

[[nodiscard]] Quadratic fit_quadratic(const std::vector<double>& x,
                                       const std::vector<double>& y);

}  
