#include "quantengine/risk.hpp"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <stdexcept>
#include <utility>

namespace quantengine::risk {
namespace {

double derivative(const std::function<double(double)>& p,
                  double value, double step) {
  const auto central = [&](double h) {
    return (p(value + h) - p(value - h)) / (2 * h);
  };
  return (4 * central(step / 2) - central(step)) / 3;
}

double second_derivative(const std::function<double(double)>& p,
                         double value, double step) {
  const double center = p(value);
  const auto central = [&](double h) {
    return (p(value + h) - 2 * center + p(value - h)) / (h * h);
  };
  return (4 * central(step / 2) - central(step)) / 3;
}

}

analytics::Sensitivities numerical_greeks(
    const std::function<double(const analytics::EuropeanMarket&)>& pricer,
    analytics::EuropeanMarket m, BumpSettings b) {
  if (!std::isfinite(m.spot) || !std::isfinite(m.volatility) ||
      !std::isfinite(m.maturity) ||
      m.spot <= 0 || m.volatility <= 0 || m.maturity <= 0 ||
      !std::isfinite(b.relative_spot) || b.relative_spot <= 0 ||
      !std::isfinite(b.absolute_rate) || b.absolute_rate <= 0 ||
      !std::isfinite(b.absolute_volatility) || b.absolute_volatility <= 0 ||
      !std::isfinite(b.absolute_maturity) || b.absolute_maturity <= 0)
    throw std::invalid_argument("Invalid risk bump configuration");
  const double spot_step = m.spot * std::min(b.relative_spot, 0.5);
  const double vol_step = std::min(b.absolute_volatility, 0.5 * m.volatility);
  const double time_step = std::min(b.absolute_maturity, 0.5 * m.maturity);
  const auto spot_pricer = [&](double x) { auto p = m; p.spot = x; return pricer(p); };
  const auto vol_pricer = [&](double x) { auto p = m; p.volatility = x; return pricer(p); };
  const auto rate_pricer = [&](double x) { auto p = m; p.rate = x; return pricer(p); };
  const auto time_pricer = [&](double x) { auto p = m; p.maturity = x; return pricer(p); };
  return {
      derivative(spot_pricer, m.spot, spot_step),
      second_derivative(spot_pricer, m.spot, spot_step),
      derivative(vol_pricer, m.volatility, vol_step),
      -derivative(time_pricer, m.maturity, time_step),
      derivative(rate_pricer, m.rate, b.absolute_rate),
  };
}

PortfolioRisk aggregate_black_scholes(
    const std::vector<VanillaPosition>& positions) {
  PortfolioRisk total{0, {0, 0, 0, 0, 0}};
  for (const auto& position : positions) {
    if (!std::isfinite(position.quantity))
      throw std::invalid_argument("Invalid portfolio quantity");
    total.value += position.quantity *
                   analytics::black_scholes(position.market, position.right);
    const auto g = analytics::black_scholes_risk(position.market, position.right);
    total.sensitivities.delta += position.quantity * g.delta;
    total.sensitivities.gamma += position.quantity * g.gamma;
    total.sensitivities.vega += position.quantity * g.vega;
    total.sensitivities.theta += position.quantity * g.theta;
    total.sensitivities.rho += position.quantity * g.rho;
  }
  return total;
}

TailRisk historical_tail_risk(std::vector<double> pnl, double confidence) {
  if (pnl.size() < 2 || !std::isfinite(confidence) || confidence <= 0 ||
      confidence >= 1)
    throw std::invalid_argument("Invalid tail-risk sample or confidence level");
  for (double p : pnl)
    if (!std::isfinite(p)) throw std::invalid_argument("Non-finite scenario P&L");
  std::sort(pnl.begin(), pnl.end());
  const std::size_t n = pnl.size();
  const double tail_size = (1 - confidence) * n;
  const std::size_t tail_count = std::max<std::size_t>(
      1, static_cast<std::size_t>(std::ceil(tail_size)));
  const double var = std::max(-pnl[tail_count - 1], 0.0);
  const double es = -std::accumulate(pnl.begin(),
                                    pnl.begin() + tail_count, 0.0) /
                    static_cast<double>(tail_count);
  return {var, std::max(es, 0.0), confidence};
}

}
