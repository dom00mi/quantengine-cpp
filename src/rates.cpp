#include "quantengine/rates.hpp"

#include <cmath>
#include <stdexcept>

#include "quantengine/analytics.hpp"

namespace quantengine::rates {
namespace {

double option_on_forward(double forward, double strike, double volatility,
                         double expiry, VolatilityQuote convention,
                         payoffs::Right right) {
  if (!std::isfinite(forward) || !std::isfinite(strike) ||
      !std::isfinite(volatility) || volatility < 0 ||
      !std::isfinite(expiry) || expiry < 0)
    throw std::invalid_argument("Invalid forward option parameters");
  if (convention == VolatilityQuote::Normal)
    return analytics::bachelier(forward, strike, 1, volatility, expiry, right);
  if (forward <= 0 || strike <= 0)
    throw std::invalid_argument("Lognormal quotes require positive forward/strike");
  return analytics::black_scholes(
      {forward, strike, 0, 0, volatility, expiry}, right);
}

void validate_schedule(SwapSchedule schedule) {
  if (!std::isfinite(schedule.start) || !std::isfinite(schedule.end) ||
      schedule.start < 0 || schedule.end <= schedule.start ||
      schedule.coupons_per_year == 0)
    throw std::invalid_argument("Invalid fixed swap schedule");
  const double exact = (schedule.end - schedule.start) *
                       schedule.coupons_per_year;
  if (std::abs(exact - std::round(exact)) > 1e-8 || exact < 1)
    throw std::invalid_argument("Swap tenor must divide into uniform coupons");
}

}

double forward_libor(const market::DiscountCurve& curve,
                     double start, double end) {
  if (!std::isfinite(start) || !std::isfinite(end) || start < 0 || end <= start)
    throw std::invalid_argument("Invalid forward rate dates");
  return (curve.discount(start) / curve.discount(end) - 1) / (end - start);
}

double swap_annuity(const market::DiscountCurve& curve,
                    SwapSchedule schedule) {
  validate_schedule(schedule);
  const auto count = static_cast<std::size_t>(std::round(
      (schedule.end - schedule.start) * schedule.coupons_per_year));
  const double interval = 1.0 / schedule.coupons_per_year;
  double annuity = 0;
  for (std::size_t i = 1; i <= count; ++i)
    annuity += interval * curve.discount(schedule.start + i * interval);
  return annuity;
}

double forward_swap_rate(const market::DiscountCurve& curve,
                         SwapSchedule schedule) {
  return (curve.discount(schedule.start) - curve.discount(schedule.end)) /
         swap_annuity(curve, schedule);
}

double payer_swap_pv(const market::DiscountCurve& curve,
                     SwapSchedule schedule, double fixed_rate, double notional) {
  if (!std::isfinite(fixed_rate) || !std::isfinite(notional) || notional <= 0)
    throw std::invalid_argument("Invalid swap notional or fixed rate");
  return notional * ((curve.discount(schedule.start) -
                      curve.discount(schedule.end)) -
                     fixed_rate * swap_annuity(curve, schedule));
}

double caplet(const market::DiscountCurve& curve, double start, double end,
              double strike, double volatility, VolatilityQuote convention,
              payoffs::Right right, double notional) {
  if (!std::isfinite(notional) || notional <= 0)
    throw std::invalid_argument("Invalid caplet notional");
  const double accrual = end - start;
  const double forward = forward_libor(curve, start, end);
  return notional * accrual * curve.discount(end) *
         option_on_forward(forward, strike, volatility, start, convention, right);
}

double swaption(const market::DiscountCurve& curve, SwapSchedule schedule,
                double strike, double volatility, VolatilityQuote convention,
                payoffs::Right right, double notional) {
  if (!std::isfinite(notional) || notional <= 0)
    throw std::invalid_argument("Invalid swaption notional");
  const double annuity = swap_annuity(curve, schedule);
  const double fwd = forward_swap_rate(curve, schedule);
  return notional * annuity * option_on_forward(
      fwd, strike, volatility, schedule.start, convention, right);
}

double hull_white_caplet(const market::DiscountCurve& curve,
                         double start, double end, double strike,
                         double mean_reversion, double volatility,
                         double notional) {
  if (!std::isfinite(mean_reversion) || mean_reversion <= 0 ||
      !std::isfinite(volatility) || volatility < 0 ||
      !std::isfinite(notional) || notional <= 0 ||
      !std::isfinite(strike) || !(end > start) || start < 0)
    throw std::invalid_argument("Invalid Hull-White caplet parameters");

  const double period = end - start;
  const double repayment = 1 + period * strike;
  if (repayment <= 0)
    throw std::domain_error("Caplet bond-put mapping requires positive repayment");
  const double p_start = curve.discount(start);
  const double p_end = curve.discount(end);
  const double bond_strike = 1 / repayment;
  const double bond_vol = volatility / mean_reversion *
      (-std::expm1(-mean_reversion * period)) *
      std::sqrt((-std::expm1(-2 * mean_reversion * start)) /
                (2 * mean_reversion));
  const double bond_put = option_on_forward(
      p_end / p_start, bond_strike, bond_vol, 1,
      VolatilityQuote::Lognormal, payoffs::Right::Put);
  return notional * repayment * p_start * bond_put;
}

}
