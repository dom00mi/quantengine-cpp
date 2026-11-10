#include "quantengine/model_calibration.hpp"

#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

#include "quantengine/optimisation.hpp"
#include "quantengine/rates.hpp"

namespace quantengine::model_calibration {

HullWhiteFit fit_hull_white_caplets(const market::DiscountCurve& curve,
                                    const std::vector<CapletMarketQuote>& quotes,
                                    double initial_mean_reversion,
                                    double initial_volatility) {
  if (quotes.size() < 2 || !std::isfinite(initial_mean_reversion) ||
      !std::isfinite(initial_volatility) || initial_mean_reversion <= 0 ||
      initial_volatility <= 0)
    throw std::invalid_argument("Insufficient Hull-White calibration inputs");
  double weight_sum = 0;
  for (const auto& q : quotes) {
    if (!std::isfinite(q.start) || !std::isfinite(q.end) ||
        !std::isfinite(q.strike) || !std::isfinite(q.market_price) ||
        !std::isfinite(q.weight) || q.start <= 0 || q.end <= q.start ||
        q.market_price < 0 || q.weight <= 0)
      throw std::invalid_argument("Invalid caplet calibration quote");
    weight_sum += q.weight;
  }

  const auto loss = [&](const std::vector<double>& x) {
    const double a = std::exp(x[0]);
    const double sigma = std::exp(x[1]);
    if (!std::isfinite(a) || !std::isfinite(sigma) ||
        a < 1e-4 || a > 10 || sigma < 1e-6 || sigma > 1)
      return std::numeric_limits<double>::infinity();
    double error = 0;
    for (const auto& q : quotes) {
      const double fitted = rates::hull_white_caplet(
          curve, q.start, q.end, q.strike, a, sigma);
      error += q.weight * (fitted - q.market_price) * (fitted - q.market_price);
    }
    return error / weight_sum;
  };

  const auto fit = optimisation::nelder_mead(
      loss, {std::log(initial_mean_reversion), std::log(initial_volatility)},
      {0.1, 0.1}, 2500, 1e-16);
  return {std::exp(fit.parameters[0]), std::exp(fit.parameters[1]),
          std::sqrt(fit.objective), fit.converged};
}

}
