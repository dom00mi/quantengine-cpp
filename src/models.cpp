#include "quantengine/models.hpp"

#include <cmath>
#include <stdexcept>
#include <utility>

namespace quantengine::models {
namespace {
void validate_step(double dt, double z) {
  if (!std::isfinite(dt) || dt <= 0 || !std::isfinite(z))
    throw std::invalid_argument("Invalid process step");
}
}

TimeGrid::TimeGrid(std::vector<double> times) : times_(std::move(times)) {
  if (times_.size() < 2 || times_.front() != 0)
    throw std::invalid_argument("Time grid must begin at 0 and contain steps");
  for (std::size_t i = 1; i < times_.size(); ++i) {
    if (!std::isfinite(times_[i]) || times_[i] <= times_[i - 1])
      throw std::invalid_argument("Time grid must increase strictly");
  }
}

TimeGrid TimeGrid::uniform(double maturity, std::size_t steps) {
  if (!std::isfinite(maturity) || maturity <= 0 || steps == 0)
    throw std::invalid_argument("Invalid simulation horizon");
  std::vector<double> times(steps + 1);
  for (std::size_t i = 0; i <= steps; ++i)
    times[i] = maturity * static_cast<double>(i) / static_cast<double>(steps);
  return TimeGrid(std::move(times));
}

void GeometricBrownian::validate() const {
  if (!std::isfinite(spot) || spot <= 0 || !std::isfinite(risk_free_rate) ||
      !std::isfinite(dividend_yield) || !std::isfinite(volatility) || volatility < 0)
    throw std::invalid_argument("Invalid geometric Brownian model");
}

State GeometricBrownian::initial() const { return {spot, 1.0}; }

State GeometricBrownian::evolve(State s, double /*time*/, double dt,
                                double z) const {
  validate_step(dt, z);
  return {s.value * std::exp((risk_free_rate - dividend_yield -
                              0.5 * volatility * volatility) * dt +
                             volatility * std::sqrt(dt) * z),
          s.discount_factor * std::exp(-risk_free_rate * dt)};
}

Observation GeometricBrownian::observe(State s, double time) const {
  return {time, s.value, risk_free_rate, s.discount_factor};
}

void ArithmeticBrownian::validate() const {
  if (!std::isfinite(initial_forward) || !std::isfinite(drift) ||
      !std::isfinite(normal_volatility) || normal_volatility < 0 ||
      !std::isfinite(risk_free_rate))
    throw std::invalid_argument("Invalid arithmetic Brownian model");
}

State ArithmeticBrownian::initial() const { return {initial_forward, 1.0}; }

State ArithmeticBrownian::evolve(State s, double /*time*/, double dt,
                                 double z) const {
  validate_step(dt, z);
  return {s.value + drift * dt + normal_volatility * std::sqrt(dt) * z,
          s.discount_factor * std::exp(-risk_free_rate * dt)};
}

Observation ArithmeticBrownian::observe(State s, double time) const {
  return {time, s.value, risk_free_rate, s.discount_factor};
}

void MeanRevertingShortRate::validate() const {
  if (!std::isfinite(initial_rate) || !std::isfinite(mean_reversion) ||
      mean_reversion <= 0 || !std::isfinite(long_run_level) ||
      !std::isfinite(volatility) || volatility < 0)
    throw std::invalid_argument("Invalid mean-reverting rate model");
}

State MeanRevertingShortRate::initial() const { return {initial_rate, 1.0}; }

State MeanRevertingShortRate::evolve(State s, double /*time*/, double dt,
                                    double z) const {
  validate_step(dt, z);
  const double decay = std::exp(-mean_reversion * dt);
  const double next_rate = long_run_level + (s.value - long_run_level) * decay +
                           volatility * std::sqrt(-std::expm1(-2 * mean_reversion * dt) /
                                                  (2 * mean_reversion)) * z;
  // Trapezoidal pathwise approximation to exp(-integral r dt). This does not
  // simulate the exact joint law of the OU endpoint and integrated short rate.
  return {next_rate, s.discount_factor *
                         std::exp(-0.5 * (s.value + next_rate) * dt)};
}

Observation MeanRevertingShortRate::observe(State s, double time) const {
  return {time, s.value, s.value, s.discount_factor};
}

HullWhiteOneFactor::HullWhiteOneFactor(market::DiscountCurve initial_discount_curve,
                                       double mean_reversion, double volatility)
    : curve_(std::move(initial_discount_curve)),
      mean_reversion_(mean_reversion), volatility_(volatility) {
  validate();
}

void HullWhiteOneFactor::validate() const {
  if (!std::isfinite(mean_reversion_) || mean_reversion_ <= 0 ||
      !std::isfinite(volatility_) || volatility_ < 0)
    throw std::invalid_argument("Invalid Hull-White model parameters");
}

double HullWhiteOneFactor::deterministic_shift(double time) const {
  constexpr double forward_width = 1e-5;
  const double forward = curve_.forward_rate(time, time + forward_width);
  const double adjustment =
      0.5 * std::pow(volatility_ / mean_reversion_, 2) *
      std::pow(-std::expm1(-mean_reversion_ * time), 2);
  return forward + adjustment;
}

State HullWhiteOneFactor::initial() const {
  return {deterministic_shift(0), 1.0};
}

State HullWhiteOneFactor::evolve(State s, double time, double dt,
                                 double z) const {
  validate_step(dt, z);
  const double factor = s.value - deterministic_shift(time);
  const double decay = std::exp(-mean_reversion_ * dt);
  const double next_factor = factor * decay +
      volatility_ * std::sqrt(-std::expm1(-2 * mean_reversion_ * dt) /
                              (2 * mean_reversion_)) * z;
  const double next_rate = next_factor + deterministic_shift(time + dt);
  return {next_rate,
          s.discount_factor * std::exp(-0.5 * (s.value + next_rate) * dt)};
}

Observation HullWhiteOneFactor::observe(State s, double time) const {
  return {time, s.value, s.value, s.discount_factor};
}

}
