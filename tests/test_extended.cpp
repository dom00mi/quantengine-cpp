#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

#include "quantengine/quantengine.hpp"

namespace {
void near(double actual, double expected, double tolerance = 1e-9) {
  if (!std::isfinite(actual) || std::abs(actual - expected) > tolerance) {
    std::cerr << "Expected " << expected << " but received " << actual << '\n';
    throw std::runtime_error("Extended numerical benchmark failed");
  }
}

template <class Fun>
void must_throw(Fun&& call) {
  try {
    call();
  } catch (const std::exception&) {
    return;
  }
  throw std::runtime_error("Missing validation exception");
}
}

int main() {
  using namespace quantengine;
  using payoffs::Right;

  // 1. Two-factor correlation and independently priced basket.
  const multifactor::CorrelatedGBM basket_model(
      {100, 80}, {0.01, 0.02}, {0.2, 0.3}, 0.04,
      {{1, 0.6}, {0.6, 1}});
  const std::array<double, 2> input{1, 0};
  const auto correlated = basket_model.correlate(input);
  near(correlated[0], 1);
  near(correlated[1], 0.6);
  const auto basket_result = multifactor::price(
      basket_model, models::TimeGrid::uniform(1, 24),
      multifactor::BasketCall{{0.5, 0.5}, 0}, 100'000, 47);
  const double exact_basket = 50 * std::exp(-0.01) + 40 * std::exp(-0.02);
  if (std::abs(basket_result.price - exact_basket) >
      5 * basket_result.standard_error + 0.03)
    throw std::runtime_error("Correlated basket forward benchmark failed");
  must_throw([] {
    (void)multifactor::CorrelatedGBM({100, 100}, {0, 0}, {0.2, 0.2}, 0,
                                     {{1, 1.2}, {1.2, 1}});
  });

  // 2. European and American PDE vs exact Black-Scholes.
  const analytics::EuropeanMarket m{100, 100, 0.05, 0.02, 0.2, 1};
  const auto pde_call = pde::price(m, Right::Call);
  const auto bs_call = analytics::black_scholes(m, Right::Call);
  const auto bs_risk = analytics::black_scholes_risk(m, Right::Call);
  near(pde_call.price, bs_call, 0.06);
  near(pde_call.delta, bs_risk.delta, 0.02);
  near(pde_call.gamma, bs_risk.gamma, 0.004);
  const auto eu_put = pde::price(m, Right::Put);
  const auto am_put = pde::price(m, Right::Put, true);
  if (am_put.price + 1e-6 < eu_put.price ||
      am_put.price < analytics::black_scholes(m, Right::Put))
    throw std::runtime_error("American PDE price is below European price");

  // 3. Deposits, par swaps and shape-preserving interpolation.
  const auto bootstrapped = curves::bootstrap_discount_curve(
      {{1, 0.03}}, {{2, 0.04, 1}});
  const double p1 = 1 / 1.03;
  const double p2 = (1 - 0.04 * p1) / 1.04;
  near(bootstrapped.discount(1), p1);
  near(bootstrapped.discount(2), p2);
  near(rates::forward_swap_rate(bootstrapped, {0, 2, 1}), 0.04, 1e-10);
  const curves::MonotoneCubicCurve monotone({0, 1, 2, 3}, {0, 0.3, 0.6, 1});
  for (int i = 0; i <= 300; ++i) {
    const double value = monotone.value(i / 100.0);
    if (value < -1e-12 || value > 1 + 1e-12)
      throw std::runtime_error("Monotone interpolant overshot");
  }

  // 4. Caplet / floorlet and payer / receiver swaption parity.
  const auto curve = market::DiscountCurve::flat(0.035);
  const double start = 1.0, end = 1.5, strike = 0.034;
  const double cap = rates::caplet(curve, start, end, strike, 0.2,
                                  rates::VolatilityQuote::Lognormal);
  const double floorlet = rates::caplet(curve, start, end, strike, 0.2,
                                       rates::VolatilityQuote::Lognormal,
                                       Right::Put);
  near(cap - floorlet, (end - start) * curve.discount(end) *
                           (rates::forward_libor(curve, start, end) - strike));
  const rates::SwapSchedule schedule{1, 4, 2};
  const double annuity = rates::swap_annuity(curve, schedule);
  const double swap_fwd = rates::forward_swap_rate(curve, schedule);
  const double payer = rates::swaption(curve, schedule, strike, 0.25,
                                      rates::VolatilityQuote::Lognormal);
  const double receiver = rates::swaption(curve, schedule, strike, 0.25,
                                         rates::VolatilityQuote::Lognormal,
                                         Right::Put);
  near(payer - receiver, annuity * (swap_fwd - strike));
  near(rates::payer_swap_pv(curve, schedule, strike),
       annuity * (swap_fwd - strike));
  const double intrinsic_cap = rates::caplet(curve, start, end, strike, 0,
                                             rates::VolatilityQuote::Normal);
  near(rates::hull_white_caplet(curve, start, end, strike, 0.6, 0),
       intrinsic_cap);

  // 5a. Calibrate a Hull-White model to synthetic caplet market prices.
  std::vector<model_calibration::CapletMarketQuote> caplets;
  for (double t : {0.5, 1.0, 1.5, 2.5, 4.0}) {
    const double k = rates::forward_libor(curve, t, t + 0.5) + 0.001;
    caplets.push_back({t, t + 0.5, k,
                       rates::hull_white_caplet(curve, t, t + 0.5, k,
                                                0.65, 0.015)});
  }
  const auto hw_fit = model_calibration::fit_hull_white_caplets(
      curve, caplets, 0.45, 0.01);
  near(hw_fit.weighted_rmse, 0, 1e-6);
  near(hw_fit.mean_reversion, 0.65, 0.12);
  near(hw_fit.volatility, 0.015, 0.002);

  // 5b. Raw-SVI smile fit and interpolated total variance surface.
  const volatility::SviParameters raw{0.01, 0.12, -0.35, 0.02, 0.20};
  std::vector<volatility::SmileQuote> smile;
  for (int i = -6; i <= 6; ++i) {
    const double k = i / 15.0;
    smile.push_back({k, std::sqrt(volatility::svi_total_variance(k, raw))});
  }
  const auto smile_fit = volatility::calibrate_svi(
      1, smile, {0.015, 0.1, -0.2, 0, 0.25});
  near(smile_fit.total_variance_rmse, 0, 2e-4);
  const volatility::SviParameters second{0.02, 0.24, -0.35, 0.02, 0.2};
  const volatility::SviSurface surface({{1, raw}, {2, second}});
  if (!surface.calendar_monotone())
    throw std::runtime_error("Synthetic SVI calendar monotonicity failed");
  const double k = 0.05;
  near(surface.total_variance(k, 1.5),
       0.5 * (volatility::svi_total_variance(k, raw) +
              volatility::svi_total_variance(k, second)));

  // 6. Richardson Greeks, portfolio aggregation and loss-tail risk.
  const auto numerical = risk::numerical_greeks(
      [](const analytics::EuropeanMarket& market) {
        return analytics::black_scholes(market, Right::Call);
      }, m);
  near(numerical.delta, bs_risk.delta, 1e-6);
  near(numerical.gamma, bs_risk.gamma, 1e-5);
  near(numerical.vega, bs_risk.vega, 1e-5);
  near(numerical.rho, bs_risk.rho, 1e-5);
  near(numerical.theta, bs_risk.theta, 1e-5);
  const auto portfolio = risk::aggregate_black_scholes(
      {{m, Right::Call, 2}, {m, Right::Put, -1}});
  near(portfolio.value, 2 * bs_call - analytics::black_scholes(m, Right::Put));
  const auto tail = risk::historical_tail_risk({-5, -3, -1, 0, 2, 4}, 0.5);
  near(tail.value_at_risk, 1);
  near(tail.expected_shortfall, 3);

  std::cout << "All extended quantitative-finance tests passed\n";
}
