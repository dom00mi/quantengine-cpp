#include "quantengine/calibration.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace quantengine::calibration {

FlatVolFit calibrate_flat_black_vol(double spot, double rate,
                                    double dividend_yield,
                                    const std::vector<OptionQuote>& quotes,
                                    double maximum_volatility) {
  if (quotes.empty() || !std::isfinite(maximum_volatility) ||
      maximum_volatility <= 0)
    throw std::invalid_argument("Invalid calibration configuration");
  double total_weight = 0;
  for (const auto& quote : quotes) {
    if (!std::isfinite(quote.price) || quote.price < 0 ||
        !std::isfinite(quote.weight) || quote.weight <= 0 ||
        !std::isfinite(quote.strike) || quote.strike <= 0 ||
        !std::isfinite(quote.maturity) || quote.maturity <= 0)
      throw std::invalid_argument("Invalid calibration quote");
    total_weight += quote.weight;
  }

  const auto squared_error = [&](double vol) {
    double error = 0;
    for (const auto& quote : quotes) {
      const analytics::EuropeanMarket inputs{
          spot, quote.strike, rate, dividend_yield, vol, quote.maturity};
      const double residual = analytics::black_scholes(inputs, quote.right) -
                              quote.price;
      error += quote.weight * residual * residual;
    }
    return error;
  };

  // Bounded, derivative-free minimisation over one common volatility.
  constexpr double phi = 0.6180339887498948482;
  double left = 0;
  double right = maximum_volatility;
  double x1 = right - phi * (right - left);
  double x2 = left + phi * (right - left);
  double y1 = squared_error(x1);
  double y2 = squared_error(x2);
  for (std::size_t step = 0; step < 250 && right - left > 1e-11; ++step) {
    if (y1 < y2) {
      right = x2;
      x2 = x1;
      y2 = y1;
      x1 = right - phi * (right - left);
      y1 = squared_error(x1);
    } else {
      left = x1;
      x1 = x2;
      y1 = y2;
      x2 = left + phi * (right - left);
      y2 = squared_error(x2);
    }
  }
  const double vol = 0.5 * (left + right);
  return {vol, std::sqrt(squared_error(vol) / total_weight)};
}

}
