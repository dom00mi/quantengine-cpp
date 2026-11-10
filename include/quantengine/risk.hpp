#pragma once

#include <functional>
#include <vector>

#include "quantengine/analytics.hpp"

namespace quantengine::risk {

struct BumpSettings {
  double relative_spot = 1e-3;
  double absolute_volatility = 1e-3;
  double absolute_rate = 1e-4;
  double absolute_maturity = 1e-3;
};

// Bump-and-revalue Greeks with Richardson extrapolation. The callback may
// price analytically by PDE or by simulation
[[nodiscard]] analytics::Sensitivities numerical_greeks(
    const std::function<double(const analytics::EuropeanMarket&)>& pricer,
    analytics::EuropeanMarket market, BumpSettings bumps = {});

struct VanillaPosition {
  analytics::EuropeanMarket market;
  payoffs::Right right;
  double quantity;
};

struct PortfolioRisk {
  double value;
  analytics::Sensitivities sensitivities;
};

[[nodiscard]] PortfolioRisk aggregate_black_scholes(
    const std::vector<VanillaPosition>& positions);

struct TailRisk {
  double value_at_risk;
  double expected_shortfall;
  double confidence;
};

// P&L input is scenario profit (positive) or loss (negative). Both returned
// metrics are positive loss amounts, clipped to zero when there is no loss.
[[nodiscard]] TailRisk historical_tail_risk(
    std::vector<double> scenario_pnl, double confidence = 0.99);

}
