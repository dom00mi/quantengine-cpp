#include "quantengine/multifactor.hpp"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <stdexcept>
#include <utility>

namespace quantengine::multifactor {

CorrelatedGBM::CorrelatedGBM(std::vector<double> initial_spots,
                             std::vector<double> dividend_yields,
                             std::vector<double> volatilities,
                             double risk_free_rate, Matrix correlation)
    : spots_(std::move(initial_spots)),
      yields_(std::move(dividend_yields)),
      vols_(std::move(volatilities)), rate_(risk_free_rate) {
  const auto dim = spots_.size();
  if (dim == 0 || yields_.size() != dim || vols_.size() != dim ||
      correlation.size() != dim || !std::isfinite(rate_))
    throw std::invalid_argument("Incompatible multi-factor model dimensions");
  for (const auto& row : correlation)
    if (row.size() != dim)
      throw std::invalid_argument("Correlation matrix must be square");
  for (std::size_t i = 0; i < dim; ++i) {
    if (!std::isfinite(spots_[i]) || spots_[i] <= 0 ||
        !std::isfinite(yields_[i]) || !std::isfinite(vols_[i]) || vols_[i] < 0 ||
        correlation[i].size() != dim)
      throw std::invalid_argument("Invalid multi-factor market parameters");
    for (std::size_t j = 0; j < dim; ++j) {
      if (!std::isfinite(correlation[i][j]) ||
          std::abs(correlation[i][j]) > 1 + 1e-12 ||
          std::abs(correlation[i][j] - correlation[j][i]) > 1e-10)
        throw std::invalid_argument("Correlation matrix must be symmetric");
    }
    if (std::abs(correlation[i][i] - 1) > 1e-10)
      throw std::invalid_argument("Correlation diagonal must equal one");
  }
  cholesky_ = Matrix(dim, std::vector<double>(dim, 0));
  for (std::size_t i = 0; i < dim; ++i) {
    for (std::size_t j = 0; j <= i; ++j) {
      double value = correlation[i][j];
      for (std::size_t k = 0; k < j; ++k)
        value -= cholesky_[i][k] * cholesky_[j][k];
      if (i == j) {
        if (value <= 1e-12)
          throw std::invalid_argument("Correlation matrix must be positive definite");
        cholesky_[i][j] = std::sqrt(value);
      } else {
        cholesky_[i][j] = value / cholesky_[j][j];
      }
    }
  }
}

std::vector<double> CorrelatedGBM::correlate(
    std::span<const double> shocks) const {
  if (shocks.size() != dimension())
    throw std::invalid_argument("Incorrect number of stochastic factors");
  std::vector<double> result(dimension());
  for (std::size_t i = 0; i < dimension(); ++i) {
    for (std::size_t j = 0; j <= i; ++j)
      result[i] += cholesky_[i][j] * shocks[j];
  }
  return result;
}

std::vector<double> CorrelatedGBM::evolve(
    std::span<const double> current, double dt,
    std::span<const double> shocks) const {
  if (current.size() != dimension() || !std::isfinite(dt) || dt <= 0)
    throw std::invalid_argument("Invalid multi-factor evolution step");
  const auto correlated = correlate(shocks);
  const double sqrt_t = std::sqrt(dt);
  std::vector<double> next(dimension());
  for (std::size_t i = 0; i < dimension(); ++i) {
    next[i] = current[i] *
        std::exp((rate_ - yields_[i] - 0.5 * vols_[i] * vols_[i]) * dt +
                 vols_[i] * sqrt_t * correlated[i]);
  }
  return next;
}

double BasketCall::operator()(std::span<const Observation> path) const {
  if (path.empty() || weights.size() != path.back().assets.size() ||
      !std::isfinite(strike))
    throw std::invalid_argument("Invalid basket payoff inputs");
  double basket = 0;
  for (std::size_t i = 0; i < weights.size(); ++i)
    basket += weights[i] * path.back().assets[i];
  return std::max(basket - strike, 0.0);
}

}
