#include "quantengine/curves.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>

#include "quantengine/numerics.hpp"

namespace quantengine::curves {
namespace {

double endpoint_slope(double h0, double h1, double d0, double d1) {
  double m = ((2 * h0 + h1) * d0 - h0 * d1) / (h0 + h1);
  if (std::signbit(m) != std::signbit(d0)) return 0;
  if (std::signbit(d0) != std::signbit(d1) && std::abs(m) > 3 * std::abs(d0))
    return 3 * d0;
  return m;
}

}

MonotoneCubicCurve::MonotoneCubicCurve(std::vector<double> x,
                                       std::vector<double> y)
    : x_(std::move(x)), y_(std::move(y)) {
  if (x_.size() < 2 || x_.size() != y_.size())
    throw std::invalid_argument("Cubic curve needs at least two matching nodes");
  const std::size_t n = x_.size();
  std::vector<double> h(n - 1), d(n - 1);
  for (std::size_t i = 0; i < n; ++i) {
    if (!std::isfinite(x_[i]) || !std::isfinite(y_[i]) ||
        (i > 0 && x_[i] <= x_[i - 1]))
      throw std::invalid_argument("Invalid cubic interpolation nodes");
    if (i > 0) {
      h[i - 1] = x_[i] - x_[i - 1];
      d[i - 1] = (y_[i] - y_[i - 1]) / h[i - 1];
    }
  }
  slopes_.resize(n);
  if (n == 2) {
    slopes_[0] = slopes_[1] = d[0];
    return;
  }
  slopes_.front() = endpoint_slope(h[0], h[1], d[0], d[1]);
  slopes_.back() = endpoint_slope(h[n - 2], h[n - 3], d[n - 2], d[n - 3]);
  for (std::size_t i = 1; i + 1 < n; ++i) {
    if (d[i - 1] * d[i] <= 0) {
      slopes_[i] = 0;
    } else {
      const double w1 = 2 * h[i] + h[i - 1];
      const double w2 = h[i] + 2 * h[i - 1];
      slopes_[i] = (w1 + w2) / (w1 / d[i - 1] + w2 / d[i]);
    }
  }
}

double MonotoneCubicCurve::value(double t) const {
  if (!std::isfinite(t)) throw std::invalid_argument("Invalid interpolation point");
  if (t <= x_.front()) return y_.front();
  if (t >= x_.back()) return y_.back();
  const auto it = std::upper_bound(x_.begin(), x_.end(), t);
  const std::size_t right = static_cast<std::size_t>(it - x_.begin());
  const std::size_t left = right - 1;
  const double h = x_[right] - x_[left];
  const double u = (t - x_[left]) / h;
  return (2 * u * u * u - 3 * u * u + 1) * y_[left] +
         (u * u * u - 2 * u * u + u) * h * slopes_[left] +
         (-2 * u * u * u + 3 * u * u) * y_[right] +
         (u * u * u - u * u) * h * slopes_[right];
}

market::DiscountCurve bootstrap_discount_curve(
    std::vector<DepositQuote> deposits, std::vector<ParSwapQuote> swaps) {
  if (deposits.empty())
    throw std::invalid_argument("At least one deposit anchor is required");
  for (const auto& q : deposits)
    if (!std::isfinite(q.maturity) || !std::isfinite(q.simple_rate) ||
        q.maturity <= 0)
      throw std::invalid_argument("Invalid deposit tenor or rate");
  for (const auto& q : swaps)
    if (!std::isfinite(q.maturity) || !std::isfinite(q.fixed_rate) ||
        q.maturity <= 0 || q.coupons_per_year == 0)
      throw std::invalid_argument("Invalid par-swap tenor or rate");
  std::sort(deposits.begin(), deposits.end(),
            [](auto a, auto b) { return a.maturity < b.maturity; });
  std::sort(swaps.begin(), swaps.end(),
            [](auto a, auto b) { return a.maturity < b.maturity; });
  std::vector<double> times;
  std::vector<double> discounts;

  for (const auto& quote : deposits) {
    const double denominator = 1 + quote.simple_rate * quote.maturity;
    if (!std::isfinite(quote.maturity) || !std::isfinite(quote.simple_rate) ||
        quote.maturity <= 0 || denominator <= 0 ||
        (!times.empty() && quote.maturity <= times.back()))
      throw std::invalid_argument("Invalid deposit bootstrapping quote");
    times.push_back(quote.maturity);
    discounts.push_back(1 / denominator);
  }

  for (const auto& quote : swaps) {
    if (!std::isfinite(quote.maturity) || !std::isfinite(quote.fixed_rate) ||
        quote.maturity <= times.back() || quote.coupons_per_year == 0)
      throw std::invalid_argument("Invalid swap bootstrapping quote");
    const double exact_payments = quote.maturity * quote.coupons_per_year;
    const auto payments = static_cast<std::size_t>(std::round(exact_payments));
    if (payments == 0 || std::abs(exact_payments - payments) > 1e-8)
      throw std::invalid_argument("Swap maturity must lie on payment schedule");
    const double accrual = 1.0 / quote.coupons_per_year;
    const market::DiscountCurve prior(times, discounts);
    const double last_time = times.back();
    const double last_discount = discounts.back();

    const auto equation = [&](double candidate) {
      double fixed_pv = 0;
      for (std::size_t i = 1; i <= payments; ++i) {
        const double payment_time = i * accrual;
        double payment_discount = candidate;
        if (payment_time <= last_time) {
          payment_discount = prior.discount(payment_time);
        } else if (payment_time < quote.maturity) {
          const double fraction = (payment_time - last_time) /
                                  (quote.maturity - last_time);
          payment_discount = std::exp(std::log(last_discount) * (1 - fraction) +
                                      std::log(candidate) * fraction);
        }
        fixed_pv += accrual * payment_discount;
      }
      return quote.fixed_rate * fixed_pv + candidate - 1;
    };
    const double low = 1e-9;
    const double high = 100.0;
    if (equation(low) * equation(high) > 0)
      throw std::domain_error("Swap quote cannot be bootstrapped on this interval");
    const double final_discount = numerics::bisection(equation, low, high);
    times.push_back(quote.maturity);
    discounts.push_back(final_discount);
  }
  return market::DiscountCurve(std::move(times), std::move(discounts));
}

}
