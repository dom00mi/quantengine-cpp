#include "quantengine/analytics.hpp"

#include <algorithm>
#include <cmath>
#include <numbers>
#include <stdexcept>

#include "quantengine/numerics.hpp"

namespace quantengine::analytics {
namespace {

double direction(payoffs::Right right) {
  return right == payoffs::Right::Call ? 1.0 : -1.0;
}

void check(const EuropeanMarket& p) {
  if (!std::isfinite(p.spot) || p.spot <= 0 || !std::isfinite(p.strike) ||
      p.strike <= 0 || !std::isfinite(p.rate) ||
      !std::isfinite(p.dividend_yield) || !std::isfinite(p.volatility) ||
      p.volatility < 0 || !std::isfinite(p.maturity) || p.maturity < 0)
    throw std::invalid_argument("Invalid lognormal option inputs");
}

struct Deviates { double d1; double d2; };

Deviates deviates(const EuropeanMarket& p) {
  const double sig_sqrt_t = p.volatility * std::sqrt(p.maturity);
  const double d1 = (std::log(p.spot / p.strike) +
                     (p.rate - p.dividend_yield +
                      0.5 * p.volatility * p.volatility) * p.maturity) /
                    sig_sqrt_t;
  return {d1, d1 - sig_sqrt_t};
}

}

double normal_pdf(double x) noexcept {
  return std::exp(-0.5 * x * x) / std::sqrt(2 * std::numbers::pi);
}

double normal_cdf(double x) noexcept {
  return 0.5 * std::erfc(-x / std::numbers::sqrt2);
}

double black_scholes(const EuropeanMarket& p, payoffs::Right right) {
  check(p);
  const double discounted_spot = p.spot * std::exp(-p.dividend_yield * p.maturity);
  const double discounted_strike = p.strike * std::exp(-p.rate * p.maturity);
  if (p.volatility == 0 || p.maturity == 0)
    return std::max(direction(right) * (discounted_spot - discounted_strike), 0.0);
  const auto [d1, d2] = deviates(p);
  const double sign = direction(right);
  return sign * (discounted_spot * normal_cdf(sign * d1) -
                 discounted_strike * normal_cdf(sign * d2));
}

Sensitivities black_scholes_risk(const EuropeanMarket& p, payoffs::Right right) {
  check(p);
  if (p.volatility == 0 || p.maturity == 0)
    throw std::domain_error("Zero-volatility/expiry Greeks are singular");
  const auto [d1, d2] = deviates(p);
  const double sign = direction(right);
  const double sqrt_t = std::sqrt(p.maturity);
  const double df_dividend = std::exp(-p.dividend_yield * p.maturity);
  const double df_rate = std::exp(-p.rate * p.maturity);
  const double cdf1 = normal_cdf(sign * d1);
  const double cdf2 = normal_cdf(sign * d2);
  return {
      sign * df_dividend * cdf1,
      df_dividend * normal_pdf(d1) / (p.spot * p.volatility * sqrt_t),
      p.spot * df_dividend * normal_pdf(d1) * sqrt_t,
      -p.spot * df_dividend * normal_pdf(d1) * p.volatility / (2 * sqrt_t) -
          sign * p.rate * p.strike * df_rate * cdf2 +
          sign * p.dividend_yield * p.spot * df_dividend * cdf1,
      sign * p.strike * p.maturity * df_rate * cdf2,
  };
}

double bachelier(double forward, double strike, double discount_factor,
                 double normal_volatility, double maturity,
                 payoffs::Right right) {
  if (!std::isfinite(forward) || !std::isfinite(strike) ||
      !std::isfinite(discount_factor) || discount_factor <= 0 ||
      !std::isfinite(normal_volatility) || normal_volatility < 0 ||
      !std::isfinite(maturity) || maturity < 0)
    throw std::invalid_argument("Invalid normal-model inputs");
  const double stddev = normal_volatility * std::sqrt(maturity);
  if (stddev == 0)
    return discount_factor * payoffs::intrinsic(forward, strike, right);
  const double sign = direction(right);
  const double d = (forward - strike) / stddev;
  return discount_factor * (sign * (forward - strike) * normal_cdf(sign * d) +
                            stddev * normal_pdf(d));
}

double implied_black_vol(double price, EuropeanMarket p, payoffs::Right right) {
  check(p);
  if (p.maturity == 0 || !std::isfinite(price))
    throw std::invalid_argument("Invalid implied volatility inversion");
  p.volatility = 0;
  const double intrinsic = black_scholes(p, right);
  const double upper = right == payoffs::Right::Call
                           ? p.spot * std::exp(-p.dividend_yield * p.maturity)
                           : p.strike * std::exp(-p.rate * p.maturity);
  if (price < intrinsic - 1e-10 || price > upper + 1e-10)
    throw std::domain_error("Option quote violates no-arbitrage bounds");
  if (price <= intrinsic + 1e-10) return 0;
  const auto residual = [&](double vol) {
    p.volatility = vol;
    return black_scholes(p, right) - price;
  };
  double high = 1;
  while (residual(high) < 0 && high < 16) high *= 2;
  if (residual(high) < 0)
    throw std::domain_error("Implied volatility out of search range");
  return numerics::bisection(residual, 0, high, 1e-11);
}

}
