#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

#include "quantengine/quantengine.hpp"

int main() {
  namespace q = quantengine;
  using q::payoffs::Right;
  std::cout << std::fixed << std::setprecision(6);

  const q::analytics::EuropeanMarket equity{100, 100, 0.05, 0.02, 0.2, 1};
  const auto grid_result = q::pde::price(equity, Right::Call);
  std::cout << "European call, PDE:      " << grid_result.price << '\n'
            << "European call, analytic: "
            << q::analytics::black_scholes(equity, Right::Call) << '\n';

  const q::multifactor::CorrelatedGBM assets(
      {100, 80}, {0.01, 0.02}, {0.20, 0.30}, 0.04,
      {{1, 0.6}, {0.6, 1}});
  const auto basket = q::multifactor::price(
      assets, q::models::TimeGrid::uniform(1, 24),
      q::multifactor::BasketCall{{0.6, 0.4}, 95}, 30'000);
  std::cout << "Correlated basket call:  " << basket.price
            << " +/- " << basket.standard_error << '\n';

  const auto discount = q::curves::bootstrap_discount_curve(
      {{1, 0.03}}, {{2, 0.04, 1}, {3, 0.045, 1}});
  std::cout << "Bootstrapped 3y DF:     " << discount.discount(3) << '\n';

  const q::rates::SwapSchedule swap{1, 3, 2};
  std::cout << "Payer swaption (Black): "
            << q::rates::swaption(discount, swap, 0.04, 0.25,
                                 q::rates::VolatilityQuote::Lognormal) << '\n';
  std::cout << "Hull-White caplet:      "
            << q::rates::hull_white_caplet(discount, 1, 1.5, 0.04, 0.5, 0.01)
            << '\n';

  const q::volatility::SviParameters parameters{0.01, 0.12, -0.3, 0, 0.2};
  std::vector<q::volatility::SmileQuote> quotes;
  for (int i = -5; i <= 5; ++i) {
    const double k = i / 10.0;
    quotes.push_back({k, std::sqrt(
                             q::volatility::svi_total_variance(k, parameters))});
  }
  const auto svi = q::volatility::calibrate_svi(1, quotes);
  std::cout << "SVI fit variance RMSE:  " << svi.total_variance_rmse << '\n';

  const auto greeks = q::risk::numerical_greeks(
      [](const auto& inputs) {
        return q::analytics::black_scholes(inputs, Right::Call);
      }, equity);
  std::cout << "Richardson delta / vega: " << greeks.delta << " / "
            << greeks.vega << '\n';
}
