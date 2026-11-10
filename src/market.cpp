#include "quantengine/market.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>

namespace quantengine::market {
namespace {
void validate_times(const std::vector<double>& times, bool require_positive) {
  if (times.empty()) throw std::invalid_argument("Empty market curve");
  for (std::size_t i = 0; i < times.size(); ++i) {
    if (!std::isfinite(times[i]) ||
        (require_positive ? times[i] <= 0 : times[i] < 0) ||
        (i > 0 && times[i] <= times[i - 1]))
      throw std::invalid_argument("Curve times must be strictly increasing");
  }
}
}

DiscountCurve::DiscountCurve(std::vector<double> maturities,
                             std::vector<double> discount_factors)
    : times_(std::move(maturities)) {
  validate_times(times_, true);
  if (discount_factors.size() != times_.size())
    throw std::invalid_argument("Curve maturities/factors have different sizes");
  log_discount_factors_.reserve(discount_factors.size());
  for (double factor : discount_factors) {
    if (!std::isfinite(factor) || factor <= 0)
      throw std::invalid_argument("Discount factors must be positive and finite");
    log_discount_factors_.push_back(std::log(factor));
  }
}

DiscountCurve DiscountCurve::flat(double rate) {
  if (!std::isfinite(rate)) throw std::invalid_argument("Non-finite interest rate");
  return DiscountCurve({1.0}, {std::exp(-rate)});
}

double DiscountCurve::discount(double time) const {
  if (!std::isfinite(time) || time < 0)
    throw std::invalid_argument("Invalid curve evaluation time");
  if (time == 0) return 1;
  // With one pillar, extrapolate the corresponding continuously compounded rate.
  if (times_.size() == 1 || time <= times_.front())
    return std::exp(log_discount_factors_.front() * time / times_.front());

  const auto bound = std::lower_bound(times_.begin(), times_.end(), time);
  const std::size_t hi = bound == times_.end()
                             ? times_.size() - 1
                             : static_cast<std::size_t>(bound - times_.begin());
  const std::size_t lo = hi - 1;
  const double weight = (time - times_[lo]) / (times_[hi] - times_[lo]);
  return std::exp(log_discount_factors_[lo] +
                  weight * (log_discount_factors_[hi] - log_discount_factors_[lo]));
}

double DiscountCurve::zero_rate(double maturity) const {
  if (!std::isfinite(maturity) || maturity < 0)
    throw std::invalid_argument("Invalid maturity");
  if (maturity == 0) return -log_discount_factors_.front() / times_.front();
  return -std::log(discount(maturity)) / maturity;
}

double DiscountCurve::forward_rate(double start, double end) const {
  if (!std::isfinite(start) || !std::isfinite(end) || start < 0 || end <= start)
    throw std::invalid_argument("Invalid forward interval");
  return std::log(discount(start) / discount(end)) / (end - start);
}

LinearCurve::LinearCurve(std::vector<double> times, std::vector<double> values)
    : times_(std::move(times)), values_(std::move(values)) {
  validate_times(times_, false);
  if (times_.size() != values_.size())
    throw std::invalid_argument("Curve times/values have different sizes");
  for (double value : values_)
    if (!std::isfinite(value)) throw std::invalid_argument("Non-finite curve value");
}

double LinearCurve::value(double time) const {
  if (!std::isfinite(time) || time < 0)
    throw std::invalid_argument("Invalid curve evaluation time");
  if (time <= times_.front()) return values_.front();
  if (time >= times_.back()) return values_.back();
  const auto bound = std::lower_bound(times_.begin(), times_.end(), time);
  const std::size_t hi = static_cast<std::size_t>(bound - times_.begin());
  const double weight = (time - times_[hi - 1]) / (times_[hi] - times_[hi - 1]);
  return values_[hi - 1] + weight * (values_[hi] - values_[hi - 1]);
}

}
