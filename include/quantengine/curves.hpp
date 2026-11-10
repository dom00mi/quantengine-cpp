#pragma once

#include <cstddef>
#include <vector>

#include "quantengine/market.hpp"

namespace quantengine::curves {

struct DepositQuote {
  double maturity;
  double simple_rate;
};

struct ParSwapQuote {
  double maturity;
  double fixed_rate;
  std::size_t coupons_per_year = 1;
};

// Sequential single-curve bootstrap from simple deposits and par fixed swaps.
// Floating leg is represented by 1 - P(0,T); fixed coupons have uniform periods.
[[nodiscard]] market::DiscountCurve bootstrap_discount_curve(
    std::vector<DepositQuote> deposits, std::vector<ParSwapQuote> swaps);

// Shape-preserving cubic Hermite interpolation for market tenors and smiles.
class MonotoneCubicCurve {
 public:
  MonotoneCubicCurve(std::vector<double> x, std::vector<double> y);
  [[nodiscard]] double value(double x) const;

 private:
  std::vector<double> x_;
  std::vector<double> y_;
  std::vector<double> slopes_;
};

}
