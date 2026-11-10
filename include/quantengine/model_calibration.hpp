#pragma once

#include <vector>

#include "quantengine/market.hpp"

namespace quantengine::model_calibration {

struct CapletMarketQuote {
  double start;
  double end;
  double strike;
  double market_price;
  double weight = 1;
};

struct HullWhiteFit {
  double mean_reversion;
  double volatility;
  double weighted_rmse;
  bool converged;
};

// Fit constant mean reversion and short-rate volatility to a collection of
// caplet prices with different expiries and tenors, holding the initial curve
// fixed. Uses an unconstrained log-parameter Nelder–Mead search.
[[nodiscard]] HullWhiteFit fit_hull_white_caplets(
    const market::DiscountCurve& curve,
    const std::vector<CapletMarketQuote>& quotes,
    double initial_mean_reversion = 0.5,
    double initial_volatility = 0.01);

}
