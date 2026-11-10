#pragma once

#include <cmath>
#include <concepts>
#include <cstddef>
#include <vector>

#include "quantengine/market.hpp"

namespace quantengine::models {

struct Observation {
  double time;
  double value;
  double short_rate;
  double discount_factor;  // Numeraire inverse accumulated along the path.
};

struct State {
  double value;
  double discount_factor;
};

class TimeGrid {
 public:
  explicit TimeGrid(std::vector<double> times);
  [[nodiscard]] static TimeGrid uniform(double maturity, std::size_t steps);
  [[nodiscard]] const std::vector<double>& times() const noexcept {
    return times_;
  }

 private:
  std::vector<double> times_;
};

// Risk-neutral GBM under continuous dividend yield.
struct GeometricBrownian {
  double spot;
  double risk_free_rate;
  double dividend_yield;
  double volatility;

  void validate() const;
  [[nodiscard]] State initial() const;
  [[nodiscard]] State evolve(State current, double time, double dt,
                             double normal_shock) const;
  [[nodiscard]] Observation observe(State current, double time) const;
};

// Arithmetic Brownian motion with a constant risk-neutral drift
struct ArithmeticBrownian {
  double initial_forward;
  double drift;
  double normal_volatility;
  double risk_free_rate;

  void validate() const;
  [[nodiscard]] State initial() const;
  [[nodiscard]] State evolve(State current, double time, double dt,
                             double normal_shock) const;
  [[nodiscard]] Observation observe(State current, double time) const;
};

// One-factor Ornstein–Uhlenbeck short-rate model (Vasicek)
// The Hull–White family is represented by this constant-parameter case.
struct MeanRevertingShortRate {
  double initial_rate;
  double mean_reversion;
  double long_run_level;
  double volatility;

  void validate() const;
  [[nodiscard]] State initial() const;
  [[nodiscard]] State evolve(State current, double time, double dt,
                             double normal_shock) const;
  [[nodiscard]] Observation observe(State current, double time) const;
};

// One-factor Gaussian Hull–White model fitted to an initial discount curve.
// Expressed as r(t) = deterministic_shift(t) + mean-zero OU factor(t).
class HullWhiteOneFactor {
 public:
  HullWhiteOneFactor(market::DiscountCurve initial_discount_curve,
                     double mean_reversion, double volatility);
  void validate() const;
  [[nodiscard]] State initial() const;
  [[nodiscard]] State evolve(State current, double time, double dt,
                             double normal_shock) const;
  [[nodiscard]] Observation observe(State current, double time) const;

 private:
  [[nodiscard]] double deterministic_shift(double time) const;
  market::DiscountCurve curve_;
  double mean_reversion_;
  double volatility_;
};

template <class Model>
concept PathModel = requires(const Model& model, State state, double time,
                             double dt, double normal) {
  { model.validate() } -> std::same_as<void>;
  { model.initial() } -> std::same_as<State>;
  { model.evolve(state, time, dt, normal) } -> std::same_as<State>;
  { model.observe(state, time) } -> std::same_as<Observation>;
};

}  
