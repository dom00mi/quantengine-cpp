#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

#include "quantengine/quantengine.hpp"

namespace {

void near(double a, double b, double tolerance = 1e-9) {
  if (!std::isfinite(a) || std::abs(a - b) > tolerance) {
    std::cerr << "Expected " << b << ", observed " << a << '\n';
    throw std::runtime_error("Numerical regression failed");
  }
}

template <class F>
void throws(F&& fun) {
  try {
    fun();
  } catch (const std::exception&) {
    return;
  }
  throw std::runtime_error("Expected exception was not thrown");
}

}

int main() {
  using namespace quantengine;
  using payoffs::Right;

  const market::DiscountCurve curve({1, 2}, {std::exp(-0.03), std::exp(-0.08)});
  near(curve.discount(0), 1);
  near(curve.discount(1.5), std::exp(-0.055));
  near(curve.discount(3), std::exp(-0.13));
  near(curve.forward_rate(1, 2), 0.05);
  near(market::DiscountCurve::flat(0.04).zero_rate(2), 0.04);
  const market::LinearCurve vol_curve({0, 1, 2}, {0.2, 0.3, 0.5});
  near(vol_curve.value(1.5), 0.4);
  near(vol_curve.value(4), 0.5);
  throws([] { (void)market::DiscountCurve({1, 1}, {0.9, 0.8}); });

  const analytics::EuropeanMarket m{100, 100, 0.05, 0.02, 0.2, 1};
  const double call = analytics::black_scholes(m, Right::Call);
  const double put = analytics::black_scholes(m, Right::Put);
  near(call, 9.227005508154036);
  near(put, 6.330080627549918);
  near(call - put, 100 * std::exp(-0.02) - 100 * std::exp(-0.05));
  near(analytics::implied_black_vol(call, m, Right::Call), 0.2, 1e-8);
  const auto risk = analytics::black_scholes_risk(m, Right::Call);
  auto up = m;
  auto down = m;
  up.spot += 1e-4;
  down.spot -= 1e-4;
  near((analytics::black_scholes(up, Right::Call) -
        analytics::black_scholes(down, Right::Call)) / 2e-4,
       risk.delta, 1e-8);
  near(analytics::bachelier(100, 100, 1, 10, 1, Right::Call),
       10 / std::sqrt(2 * std::acos(-1.0)));
  near(analytics::bachelier(-2, -1, 0.95, 0, 1, Right::Put), 0.95);

  near(numerics::bisection([](double x) { return x * x - 2; }, 0, 2),
       std::sqrt(2.0));
  const std::vector<double> x{-2, -1, 0, 1, 2};
  const std::vector<double> y{9, 4, 1, 0, 1};
  const auto quadratic = numerics::fit_quadratic(x, y);
  near(quadratic.a, 1);
  near(quadratic.b, -2);
  near(quadratic.c, 1);

  std::vector<calibration::OptionQuote> quotes;
  for (double strike : {80.0, 100.0, 120.0}) {
    for (double maturity : {0.5, 1.0, 2.0}) {
      quotes.push_back({strike, maturity,
                        analytics::black_scholes({100, strike, 0.05, 0.02,
                                                  0.24, maturity}, Right::Call)});
    }
  }
  const auto fitted = calibration::calibrate_flat_black_vol(100, 0.05, 0.02, quotes);
  near(fitted.volatility, 0.24, 1e-7);
  near(fitted.weighted_rmse, 0, 1e-7);

  const auto grid = models::TimeGrid::uniform(1.0, 24);
  const models::GeometricBrownian gbm{100, 0.05, 0.02, 0.2};
  const payoffs::European european{100, Right::Call};
  const auto mc = simulation::price(gbm, grid, european, 100'000, 123);
  if (std::abs(mc.price - call) > 5 * mc.standard_error + 0.02)
    throw std::runtime_error("GBM Monte Carlo does not match analytical price");
  near(mc.price, simulation::price(gbm, grid, european, 100'000, 123).price, 0);
  const auto asian = simulation::price(gbm, grid,
                                     payoffs::Asian{100, Right::Call}, 40'000, 44);
  const auto barrier = simulation::price(gbm, grid,
                                        payoffs::UpAndOut{100, 130, Right::Call},
                                        40'000, 44);
  if (asian.price <= 0 || barrier.price < 0 || barrier.price >= call)
    throw std::runtime_error("Unexpected path-dependent option price");

  const models::ArithmeticBrownian normal{100, 0, 10, 0.03};
  const auto normal_mc = simulation::price(normal, grid,
                                          payoffs::European{100, Right::Call},
                                          80'000, 77);
  const double normal_closed = analytics::bachelier(
      100, 100, std::exp(-0.03), 10, 1, Right::Call);
  if (std::abs(normal_mc.price - normal_closed) > 5 * normal_mc.standard_error + 0.01)
    throw std::runtime_error("Arithmetic Brownian price mismatch");

  const models::GeometricBrownian put_model{90, 0.05, 0, 0.2};
  const auto amer = simulation::american_option(
      put_model, models::TimeGrid::uniform(1, 32), 100, Right::Put, 24'000, 8);
  const double euro_put = analytics::black_scholes(
      {90, 100, 0.05, 0, 0.2, 1}, Right::Put);
  if (amer.price < euro_put - 0.15)
    throw std::runtime_error("American put violates European lower benchmark");

  const models::MeanRevertingShortRate short_rate{0.03, 0.7, 0.04, 0.01};
  const double maturity = 2.0;
  const double a = short_rate.mean_reversion;
  const double b = (1 - std::exp(-a * maturity)) / a;
  const double expected_integral = short_rate.long_run_level * maturity +
      (short_rate.initial_rate - short_rate.long_run_level) * b;
  const double integrated_variance = std::pow(short_rate.volatility / a, 2) *
      (maturity - 2 * b + (1 - std::exp(-2 * a * maturity)) / (2 * a));
  const double bond_exact = std::exp(-expected_integral + 0.5 * integrated_variance);
  const auto bond_mc = simulation::price(short_rate,
                                        models::TimeGrid::uniform(maturity, 80),
                                        payoffs::ZeroCouponBond{}, 35'000, 5);
  if (std::abs(bond_mc.price - bond_exact) > 5 * bond_mc.standard_error + 0.0003)
    throw std::runtime_error("Short-rate zero-coupon bond benchmark failed");

  const auto initial_term_structure = market::DiscountCurve::flat(0.035);
  const models::HullWhiteOneFactor hw(initial_term_structure, 0.8, 0.012);
  const auto hw_bond = simulation::price(
      hw, models::TimeGrid::uniform(2, 100), payoffs::ZeroCouponBond{},
      36'000, 888);
  if (std::abs(hw_bond.price - initial_term_structure.discount(2)) >
      5 * hw_bond.standard_error + 0.0003)
    throw std::runtime_error("Hull-White initial curve reconstruction failed");

  throws([] { (void)models::TimeGrid::uniform(0, 4); });
  throws([&] { (void)simulation::price(gbm, grid, european, 3); });
  throws([&] { (void)analytics::implied_black_vol(101, m, Right::Call); });

  std::cout << "All QuantEngine regression tests passed\n";
}
