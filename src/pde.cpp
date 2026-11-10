#include "quantengine/pde.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace quantengine::pde {
namespace {

double boundary(const analytics::EuropeanMarket& m, payoffs::Right right,
                bool american, double tau, double spot) {
  if (right == payoffs::Right::Call) {
    if (spot == 0) return 0;
    const double european = std::max(spot * std::exp(-m.dividend_yield * tau) -
                                      m.strike * std::exp(-m.rate * tau), 0.0);
    return american ? std::max(european, spot - m.strike) : european;
  }
  if (spot > 0) return 0;
  const double euro = m.strike * std::exp(-m.rate * tau);
  return american ? std::max(euro, m.strike) : euro;
}

std::vector<double> tridiagonal(const std::vector<double>& lower,
                                 const std::vector<double>& diagonal,
                                 const std::vector<double>& upper,
                                 const std::vector<double>& rhs) {
  const std::size_t n = rhs.size();
  std::vector<double> c(n), d(n), x(n);
  double pivot = diagonal[0];
  if (std::abs(pivot) < 1e-14) throw std::runtime_error("Singular PDE matrix");
  c[0] = upper[0] / pivot;
  d[0] = rhs[0] / pivot;
  for (std::size_t i = 1; i < n; ++i) {
    pivot = diagonal[i] - lower[i] * c[i - 1];
    if (std::abs(pivot) < 1e-14) throw std::runtime_error("Singular PDE matrix");
    c[i] = upper[i] / pivot;
    d[i] = (rhs[i] - lower[i] * d[i - 1]) / pivot;
  }
  x[n - 1] = d[n - 1];
  for (std::size_t i = n - 1; i-- > 0;)
    x[i] = d[i] - c[i] * x[i + 1];
  return x;
}

}

Result price(const analytics::EuropeanMarket& m, payoffs::Right right,
             bool american, GridSettings settings) {
  if (!std::isfinite(m.spot) || !std::isfinite(m.strike) ||
      !std::isfinite(m.volatility) || !std::isfinite(m.maturity) ||
      !std::isfinite(m.rate) || !std::isfinite(m.dividend_yield) ||
      m.spot <= 0 || m.strike <= 0 || m.volatility <= 0 || m.maturity < 0 ||
      settings.price_steps < 4 || settings.time_steps < 1 ||
      !std::isfinite(settings.upper_spot) || settings.upper_spot < 0 ||
      !std::isfinite(settings.theta) ||
      settings.theta < 0.5 || settings.theta > 1 ||
      !std::isfinite(settings.tolerance) ||
      !std::isfinite(settings.sor_relaxation) ||
      settings.tolerance <= 0 || settings.sor_relaxation <= 0 ||
      settings.sor_relaxation >= 2 || settings.sor_max_iterations == 0)
    throw std::invalid_argument("Invalid PDE market or grid parameters");

  if (m.maturity == 0) {
    return {payoffs::intrinsic(m.spot, m.strike, right), 0, 0, 0};
  }
  const std::size_t N = settings.price_steps;
  const double spot_max = settings.upper_spot > 0
      ? settings.upper_spot : 4.0 * std::max(m.spot, m.strike);
  if (!std::isfinite(spot_max) || spot_max <= m.spot ||
      spot_max <= m.strike)
    throw std::invalid_argument("Truncated price grid is too narrow");

  const double dS = spot_max / static_cast<double>(N);
  const double dt = m.maturity / static_cast<double>(settings.time_steps);
  const double theta = settings.theta;
  const double drift = m.rate - m.dividend_yield;
  std::vector<double> old(N + 1), next(N + 1);
  for (std::size_t i = 0; i <= N; ++i)
    old[i] = payoffs::intrinsic(i * dS, m.strike, right);

  std::vector<double> lower(N - 1), diag(N - 1), upper(N - 1), rhs(N - 1);
  for (std::size_t step = 1; step <= settings.time_steps; ++step) {
    const double tau = step * dt;
    next.front() = boundary(m, right, american, tau, 0);
    next.back() = boundary(m, right, american, tau, spot_max);

    for (std::size_t i = 1; i < N; ++i) {
      const double index = static_cast<double>(i);
      const double variance = m.volatility * m.volatility * index * index;
      const double a = 0.5 * (variance - drift * index);
      const double b = -(variance + m.rate);
      const double c = 0.5 * (variance + drift * index);
      const std::size_t j = i - 1;
      lower[j] = -theta * dt * a;
      diag[j] = 1 - theta * dt * b;
      upper[j] = -theta * dt * c;
      rhs[j] = (1 + (1 - theta) * dt * b) * old[i] +
               (1 - theta) * dt * (a * old[i - 1] + c * old[i + 1]);
      if (i == 1) rhs[j] -= lower[j] * next.front();
      if (i == N - 1) rhs[j] -= upper[j] * next.back();
    }

    if (!american) {
      const auto interior = tridiagonal(lower, diag, upper, rhs);
      for (std::size_t j = 0; j < interior.size(); ++j)
        next[j + 1] = interior[j];
    } else {
      for (std::size_t i = 1; i < N; ++i) next[i] = old[i];
      bool converged = false;
      for (std::size_t iter = 0; iter < settings.sor_max_iterations; ++iter) {
        double max_change = 0;
        for (std::size_t i = 1; i < N; ++i) {
          const std::size_t j = i - 1;
          // Boundary effects have already been added to rhs; only unknown
          // interior values enter the matrix-vector iteration.
          const double left = i > 1 ? lower[j] * next[i - 1] : 0;
          const double right_neighbor = i < N - 1 ? upper[j] * next[i + 1] : 0;
          const double raw = (rhs[j] - left - right_neighbor) / diag[j];
          const double updated = std::max(
              payoffs::intrinsic(i * dS, m.strike, right),
              next[i] + settings.sor_relaxation * (raw - next[i]));
          max_change = std::max(max_change, std::abs(updated - next[i]));
          next[i] = updated;
        }
        if (max_change < settings.tolerance) {
          converged = true;
          break;
        }
      }
      if (!converged)
        throw std::runtime_error("Projected SOR failed to converge");
    }
    old.swap(next);
  }

  const double position = m.spot / dS;
  const std::size_t left_index = std::min(
      static_cast<std::size_t>(position), N - 1);
  const double weight = position - left_index;
  const double option_value = (1 - weight) * old[left_index] +
                              weight * old[left_index + 1];
  const std::size_t center = std::clamp(left_index, std::size_t{1}, N - 1);
  return {option_value,
          (old[center + 1] - old[center - 1]) / (2 * dS),
          (old[center + 1] - 2 * old[center] + old[center - 1]) / (dS * dS),
          dS};
}

}
