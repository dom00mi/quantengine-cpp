#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <random>
#include <span>
#include <stdexcept>
#include <vector>

#include "quantengine/models.hpp"
#include "quantengine/numerics.hpp"
#include "quantengine/payoffs.hpp"

namespace quantengine::simulation {

struct Estimate {
  double price;
  double standard_error;
  std::size_t paths;
};

class RunningStatistics {
 public:
  void add(double observation) {
    ++count_;
    const double delta = observation - mean_;
    mean_ += delta / static_cast<double>(count_);
    sum_squares_ += delta * (observation - mean_);
  }

  [[nodiscard]] Estimate estimate(std::size_t actual_paths) const {
    if (count_ == 0) throw std::logic_error("No Monte Carlo observations");
    const double variance = count_ > 1
                                ? sum_squares_ / static_cast<double>(count_ - 1)
                                : 0;
    return {mean_, std::sqrt(std::max(variance, 0.0) / count_), actual_paths};
  }

 private:
  std::size_t count_ = 0;
  double mean_ = 0;
  double sum_squares_ = 0;
};

// Path generation stores observations, not state-indexed value functions.
template <models::PathModel Model>
[[nodiscard]] std::vector<models::Observation> generate_path(
    const Model& model, const models::TimeGrid& grid,
    std::span<const double> standard_normals) {
  const auto& times = grid.times();
  if (standard_normals.size() != times.size() - 1)
    throw std::invalid_argument("Incorrect number of Brownian shocks");
  model.validate();
  auto state = model.initial();
  std::vector<models::Observation> observations;
  observations.reserve(times.size());
  observations.push_back(model.observe(state, times.front()));
  for (std::size_t j = 1; j < times.size(); ++j) {
    state = model.evolve(state, times[j - 1], times[j] - times[j - 1],
                         standard_normals[j - 1]);
    observations.push_back(model.observe(state, times[j]));
  }
  return observations;
}

// A generic functional-pricing engine: products consume the full simulated
// history and return an undiscounted maturity payoff.
template <models::PathModel Model, class Payoff>
[[nodiscard]] Estimate price(const Model& model, const models::TimeGrid& grid,
                             const Payoff& payoff, std::size_t paths = 100'000,
                             std::uint64_t seed = 42) {
  if (paths < 4 || paths % 2)
    throw std::invalid_argument("Antithetic simulation needs an even path count >= 4");
  model.validate();
  std::mt19937_64 engine(seed);
  std::normal_distribution<double> normal(0, 1);
  std::vector<double> shocks(grid.times().size() - 1);
  RunningStatistics statistics;
  for (std::size_t pair = 0; pair < paths / 2; ++pair) {
    for (double& shock : shocks) shock = normal(engine);
    const auto first = generate_path(model, grid, shocks);
    for (double& shock : shocks) shock = -shock;
    const auto second = generate_path(model, grid, shocks);
    const double pv_a = first.back().discount_factor * payoff(first);
    const double pv_b = second.back().discount_factor * payoff(second);
    if (!std::isfinite(pv_a) || !std::isfinite(pv_b))
      throw std::overflow_error("Non-finite discounted payoff");
    // Antithetic pair means are independent; individual payoffs are not.
    statistics.add(0.5 * (pv_a + pv_b));
  }
  return statistics.estimate(paths);
}

// Least-squares Monte Carlo (Longstaff–Schwartz), with a quadratic continuation
// regression. Exercise is permitted only on the dates in the supplied grid.
template <models::PathModel Model>
[[nodiscard]] Estimate american_option(const Model& model,
                                       const models::TimeGrid& grid,
                                       double strike, payoffs::Right right,
                                       std::size_t paths = 50'000,
                                       std::uint64_t seed = 42) {
  if (paths < 100 || !std::isfinite(strike))
    throw std::invalid_argument("Invalid American option configuration");
  model.validate();
  std::mt19937_64 engine(seed);
  std::normal_distribution<double> normal(0, 1);
  std::vector<double> shocks(grid.times().size() - 1);
  std::vector<std::vector<models::Observation>> histories;
  histories.reserve(paths);
  std::vector<double> values(paths);
  for (std::size_t i = 0; i < paths; ++i) {
    for (double& shock : shocks) shock = normal(engine);
    histories.push_back(generate_path(model, grid, shocks));
    const auto& last = histories.back().back();
    values[i] = last.discount_factor *
                payoffs::intrinsic(last.value, strike, right);
  }

  const double scale = std::max(std::abs(strike), 1.0);
  for (std::size_t step = grid.times().size() - 2; step > 0; --step) {
    std::vector<double> x;
    std::vector<double> y;
    x.reserve(paths);
    y.reserve(paths);
    for (std::size_t i = 0; i < paths; ++i) {
      const auto& observation = histories[i][step];
      if (payoffs::intrinsic(observation.value, strike, right) > 0) {
        x.push_back(observation.value / scale);
        y.push_back(values[i] / observation.discount_factor);
      }
    }
    if (x.empty()) continue;
    const auto continuation = numerics::fit_quadratic(x, y);
    for (std::size_t i = 0; i < paths; ++i) {
      const auto& observation = histories[i][step];
      const double exercise = payoffs::intrinsic(observation.value, strike, right);
      if (exercise > 0 &&
          exercise > std::max(0.0, continuation(observation.value / scale)))
        values[i] = exercise * observation.discount_factor;
    }
  }
  RunningStatistics statistics;
  for (double value : values) statistics.add(value);
  auto result = statistics.estimate(paths);
  const double immediate = payoffs::intrinsic(model.initial().value, strike, right);
  if (immediate > result.price) {
    result.price = immediate;
    result.standard_error = 0;
  }
  return result;
}

} 
