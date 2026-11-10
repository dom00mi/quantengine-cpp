#pragma once

#include <vector>

#include "quantengine/analytics.hpp"

namespace quantengine::calibration {

struct OptionQuote {
  double strike;
  double maturity;
  double price;
  double weight = 1.0;
  payoffs::Right right = payoffs::Right::Call;
};

struct FlatVolFit {
  double volatility;
  double weighted_rmse;
};

// Least-squares constant-volatility fit across a cross-section of quotes.
// This is intentionally simple; not an arbitrage-free volatility-surface fit.
[[nodiscard]] FlatVolFit calibrate_flat_black_vol(
    double spot, double rate, double dividend_yield,
    const std::vector<OptionQuote>& quotes,
    double maximum_volatility = 5.0);

}
