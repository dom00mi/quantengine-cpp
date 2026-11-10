#pragma once

#include <algorithm>
#include <cmath>
#include <span>
#include <stdexcept>

#include "quantengine/models.hpp"

namespace quantengine::payoffs {

enum class Right { Call, Put };

inline double intrinsic(double underlying, double strike, Right right) {
  const double direction = right == Right::Call ? 1.0 : -1.0;
  return std::max(direction * (underlying - strike), 0.0);
}

struct European {
  double strike;
  Right right;
  [[nodiscard]] double operator()(std::span<const models::Observation> path) const {
    return intrinsic(path.back().value, strike, right);
  }
};

// Arithmetic average of observed values after initial time
struct Asian {
  double strike;
  Right right;
  [[nodiscard]] double operator()(std::span<const models::Observation> path) const {
    if (path.size() < 2) throw std::invalid_argument("Asian payoff needs time steps");
    double sum = 0;
    for (std::size_t i = 1; i < path.size(); ++i) sum += path[i].value;
    return intrinsic(sum / static_cast<double>(path.size() - 1), strike, right);
  }
};

// Discretely monitored up-and-out European call/put.
struct UpAndOut {
  double strike;
  double barrier;
  Right right;
  [[nodiscard]] double operator()(std::span<const models::Observation> path) const {
    if (std::any_of(path.begin(), path.end(), [this](const auto& obs) {
          return obs.value >= barrier;
        })) return 0;
    return intrinsic(path.back().value, strike, right);
  }
};

// Pays one unit at maturity, discounted through the simulated short rate.
struct ZeroCouponBond {
  [[nodiscard]] double operator()(std::span<const models::Observation>) const {
    return 1.0;
  }
};

} 
