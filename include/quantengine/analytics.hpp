#pragma once

#include "quantengine/payoffs.hpp"

namespace quantengine::analytics {

struct EuropeanMarket {
  double spot;
  double strike;
  double rate;
  double dividend_yield;
  double volatility;
  double maturity;
};

struct Sensitivities {
  double delta;
  double gamma;
  double vega;
  double theta;
  double rho;
};

[[nodiscard]] double normal_pdf(double x) noexcept;
[[nodiscard]] double normal_cdf(double x) noexcept;
[[nodiscard]] double black_scholes(const EuropeanMarket& market,
                                   payoffs::Right right);
[[nodiscard]] Sensitivities black_scholes_risk(const EuropeanMarket& market,
                                               payoffs::Right right);
[[nodiscard]] double bachelier(double forward, double strike,
                               double discount_factor, double normal_volatility,
                               double maturity, payoffs::Right right);
[[nodiscard]] double implied_black_vol(double option_price,
                                       EuropeanMarket market,
                                       payoffs::Right right);

}
