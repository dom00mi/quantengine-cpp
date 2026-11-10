#pragma once

#include <cstddef>
#include <vector>

namespace quantengine::market {

// Positive discount factors interpolated linearly in log space.
// The origin (0, 1) is implicit; last observed forward slope extrapolates.
class DiscountCurve {
 public:
  DiscountCurve(std::vector<double> maturities,
                std::vector<double> discount_factors);
  [[nodiscard]] static DiscountCurve flat(double continuously_compounded_rate);
  [[nodiscard]] double discount(double time) const;
  [[nodiscard]] double zero_rate(double maturity) const;
  [[nodiscard]] double forward_rate(double start, double end) const;

 private:
  std::vector<double> times_;
  std::vector<double> log_discount_factors_;
};

// Straight-line interpolation with constant endpoint extrapolation.
class LinearCurve {
 public:
  LinearCurve(std::vector<double> times, std::vector<double> values);
  [[nodiscard]] double value(double time) const;

 private:
  std::vector<double> times_;
  std::vector<double> values_;
};

}
