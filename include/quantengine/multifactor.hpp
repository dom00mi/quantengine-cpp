#pragma once

#include <cstddef>
#include <cstdint>
#include <cmath>
#include <random>
#include <span>
#include <stdexcept>
#include <vector>

#include "quantengine/models.hpp"
#include "quantengine/monte_carlo.hpp"

namespace quantengine::multifactor {

using Matrix = std::vector<std::vector<double>>;

struct Observation {
  double time;
  std::vector<double> assets;
  double discount_factor;
};

class CorrelatedGBM {
 public:
  CorrelatedGBM(std::vector<double> initial_spots,
                std::vector<double> dividend_yields,
                std::vector<double> volatilities, double risk_free_rate,
                Matrix correlation);

  [[nodiscard]] std::size_t dimension() const noexcept { return spots_.size(); }
  [[nodiscard]] const std::vector<double>& initial_spots() const noexcept {
    return spots_;
  }
  [[nodiscard]] double rate() const noexcept { return rate_; }
  // Transforms independent standard normals into correlated innovations
  [[nodiscard]] std::vector<double> correlate(
      std::span<const double> independent_normals) const;
  [[nodiscard]] std::vector<double> evolve(
      std::span<const double> current, double dt,
      std::span<const double> independent_normals) const;

 private:
  std::vector<double> spots_;
  std::vector<double> yields_;
  std::vector<double> vols_;
  double rate_;
  Matrix cholesky_;
};

struct BasketCall {
  std::vector<double> weights;
  double strike;

  [[nodiscard]] double operator()(std::span<const Observation> path) const;
};

// Antithetic multi-factor paths: a payoff receives the full joint history
template <class Payoff>
[[nodiscard]] simulation::Estimate price(const CorrelatedGBM& model,
                                           const models::TimeGrid& times,
                                           const Payoff& payoff,
                                           std::size_t paths = 100'000,
                                           std::uint64_t seed = 42) {
  if (paths < 4 || paths % 2)
    throw std::invalid_argument("Multi-factor path count must be even and >= 4");

  const auto& grid = times.times();
  const std::size_t dimensions = model.dimension();
  std::mt19937_64 engine(seed);
  std::normal_distribution<double> normal(0, 1);
  std::vector<std::vector<double>> shocks(grid.size() - 1,
                                          std::vector<double>(dimensions));
  simulation::RunningStatistics estimator;
  for (std::size_t pair = 0; pair < paths / 2; ++pair) {
    for (auto& step : shocks)
      for (double& shock : step) shock = normal(engine);

    const auto discounted_payoff = [&](double sign) {
      std::vector<double> current = model.initial_spots();
      std::vector<Observation> history;
      history.reserve(grid.size());
      history.push_back({0, current, 1});
      for (std::size_t k = 1; k < grid.size(); ++k) {
        std::vector<double> signed_shocks = shocks[k - 1];
        if (sign < 0)
          for (double& shock : signed_shocks) shock = -shock;
        current = model.evolve(current, grid[k] - grid[k - 1], signed_shocks);
        history.push_back({grid[k], current,
                           std::exp(-model.rate() * grid[k])});
      }
      return payoff(history) * history.back().discount_factor;
    };
    const double a = discounted_payoff(1);
    const double b = discounted_payoff(-1);
    if (!std::isfinite(a) || !std::isfinite(b))
      throw std::overflow_error("Non-finite basket payoff");
    estimator.add(0.5 * (a + b));
  }
  return estimator.estimate(paths);
}

}  
