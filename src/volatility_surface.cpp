#include "quantengine/volatility_surface.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

#include "quantengine/optimisation.hpp"

namespace quantengine::volatility {
namespace {

void validate(SviParameters p) {
  if (!std::isfinite(p.a) || p.a < 0 || !std::isfinite(p.b) || p.b <= 0 ||
      !std::isfinite(p.rho) || std::abs(p.rho) >= 1 ||
      !std::isfinite(p.m) || !std::isfinite(p.sigma) || p.sigma <= 0)
    throw std::invalid_argument("SVI requires a >= 0, b/sigma > 0 and |rho| < 1");
}

SviParameters decode(const std::vector<double>& x) {
  return {std::exp(x[0]), std::exp(x[1]), std::tanh(x[2]),
          x[3], std::exp(x[4])};
}

std::vector<double> encode(SviParameters p) {
  validate(p);
  return {std::log(std::max(p.a, 1e-12)), std::log(p.b),
          std::atanh(p.rho), p.m, std::log(p.sigma)};
}

}

double svi_total_variance(double k, const SviParameters& p) {
  validate(p);
  if (!std::isfinite(k)) throw std::invalid_argument("Non-finite log moneyness");
  const double shifted = k - p.m;
  return p.a + p.b * (p.rho * shifted +
                      std::hypot(shifted, p.sigma));
}

SviFit calibrate_svi(double maturity, const std::vector<SmileQuote>& quotes,
                     SviParameters initial) {
  if (!std::isfinite(maturity) || maturity <= 0 || quotes.size() < 5)
    throw std::invalid_argument("SVI calibration requires >=5 smile quotes");
  validate(initial);
  double total_weight = 0;
  for (const auto& quote : quotes) {
    if (!std::isfinite(quote.log_moneyness) ||
        !std::isfinite(quote.implied_volatility) || quote.implied_volatility <= 0 ||
        !std::isfinite(quote.weight) || quote.weight <= 0)
      throw std::invalid_argument("Invalid SVI quote");
    total_weight += quote.weight;
  }

  const auto loss = [&](const std::vector<double>& x) {
    const auto params = decode(x);
    if (!std::isfinite(params.a) || !std::isfinite(params.b) ||
        !std::isfinite(params.sigma) || params.b <= 0 || params.sigma <= 0)
      return std::numeric_limits<double>::infinity();
    double error = 0;
    for (const auto& quote : quotes) {
      const double fitted = svi_total_variance(quote.log_moneyness, params);
      const double observed = maturity * quote.implied_volatility *
                              quote.implied_volatility;
      error += quote.weight * std::pow(fitted - observed, 2);
    }
    return error / total_weight;
  };
  const std::vector<double> steps{0.15, 0.15, 0.15, 0.05, 0.15};
  const auto first = optimisation::nelder_mead(
      loss, encode(initial), steps, 5000, 1e-15);
  const auto second = optimisation::nelder_mead(
      loss, encode({0.01, 0.2, 0, 0, 0.3}), steps, 5000, 1e-15);
  const auto& best = first.objective <= second.objective ? first : second;
  return {decode(best.parameters), std::sqrt(best.objective), best.converged};
}

SviSurface::SviSurface(std::vector<SviSlice> smiles)
    : smiles_(std::move(smiles)) {
  if (smiles_.empty()) throw std::invalid_argument("Empty volatility surface");
  for (std::size_t i = 0; i < smiles_.size(); ++i) {
    validate(smiles_[i].parameters);
    if (!std::isfinite(smiles_[i].maturity) || smiles_[i].maturity <= 0 ||
        (i > 0 && smiles_[i].maturity <= smiles_[i - 1].maturity))
      throw std::invalid_argument("SVI smile maturities must increase");
  }
}

double SviSurface::total_variance(double k, double maturity) const {
  if (!std::isfinite(maturity) || maturity <= 0)
    throw std::invalid_argument("Invalid volatility maturity");
  if (maturity <= smiles_.front().maturity)
    return maturity / smiles_.front().maturity *
           svi_total_variance(k, smiles_.front().parameters);
  if (maturity >= smiles_.back().maturity)
    return maturity / smiles_.back().maturity *
           svi_total_variance(k, smiles_.back().parameters);
  const auto it = std::lower_bound(
      smiles_.begin(), smiles_.end(), maturity,
      [](const SviSlice& smile, double t) { return smile.maturity < t; });
  const std::size_t high = static_cast<std::size_t>(it - smiles_.begin());
  const std::size_t low = high - 1;
  const double weight = (maturity - smiles_[low].maturity) /
                        (smiles_[high].maturity - smiles_[low].maturity);
  return (1 - weight) * svi_total_variance(k, smiles_[low].parameters) +
         weight * svi_total_variance(k, smiles_[high].parameters);
}

double SviSurface::implied_volatility(double k, double maturity) const {
  return std::sqrt(total_variance(k, maturity) / maturity);
}

bool SviSurface::calendar_monotone(double min_k, double max_k,
                                   std::size_t samples) const {
  if (!std::isfinite(min_k) || !std::isfinite(max_k) || min_k >= max_k ||
      samples < 2)
    throw std::invalid_argument("Invalid SVI calendar check range");
  for (std::size_t j = 0; j < samples; ++j) {
    const double k = min_k + (max_k - min_k) * j / (samples - 1);
    for (std::size_t i = 1; i < smiles_.size(); ++i)
      if (svi_total_variance(k, smiles_[i].parameters) + 1e-12 <
          svi_total_variance(k, smiles_[i - 1].parameters))
        return false;
  }
  return true;
}

}
