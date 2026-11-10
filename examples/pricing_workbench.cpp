#include <iomanip>
#include <iostream>

#include "quantengine/quantengine.hpp"

int main() {
  using namespace quantengine;
  const models::GeometricBrownian equity{100, 0.04, 0.01, 0.2};
  const auto timeline = models::TimeGrid::uniform(1.0, 50);

  const auto european = simulation::price(
      equity, timeline, payoffs::European{105, payoffs::Right::Call}, 100'000);
  const auto asian = simulation::price(
      equity, timeline, payoffs::Asian{105, payoffs::Right::Call}, 100'000);
  const auto american_put = simulation::american_option(
      equity, timeline, 105, payoffs::Right::Put, 25'000);

  std::cout << std::fixed << std::setprecision(5)
            << "European call (Monte Carlo): " << european.price
            << " +/- " << european.standard_error << "\n"
            << "Asian call (Monte Carlo):    " << asian.price
            << " +/- " << asian.standard_error << "\n"
            << "American put (LSM):         " << american_put.price
            << " +/- " << american_put.standard_error << "\n";

  const auto market = analytics::EuropeanMarket{100, 105, 0.04, 0.01, 0.2, 1};
  std::cout << "European call (closed form): "
            << analytics::black_scholes(market, payoffs::Right::Call) << '\n';
}
