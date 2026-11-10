#pragma once

#include <cstddef>

#include "quantengine/market.hpp"
#include "quantengine/payoffs.hpp"

namespace quantengine::rates {

enum class VolatilityQuote { Lognormal, Normal };

struct SwapSchedule {
  double start;
  double end;
  std::size_t coupons_per_year;
};

[[nodiscard]] double forward_libor(const market::DiscountCurve& curve,
                                   double start, double end);
[[nodiscard]] double swap_annuity(const market::DiscountCurve& curve,
                                  SwapSchedule schedule);
[[nodiscard]] double forward_swap_rate(const market::DiscountCurve& curve,
                                       SwapSchedule schedule);
[[nodiscard]] double payer_swap_pv(const market::DiscountCurve& curve,
                                   SwapSchedule schedule, double fixed_rate,
                                   double notional = 1.0);

// Caplets pay at end, with the forward fixed at start
[[nodiscard]] double caplet(const market::DiscountCurve& curve,
                            double start, double end, double strike,
                            double volatility, VolatilityQuote convention,
                            payoffs::Right right = payoffs::Right::Call,
                            double notional = 1.0);

// Payer call / receiver put on the forward swap rate
[[nodiscard]] double swaption(const market::DiscountCurve& curve,
                              SwapSchedule schedule, double strike,
                              double volatility, VolatilityQuote convention,
                              payoffs::Right right = payoffs::Right::Call,
                              double notional = 1.0);

// One-factor Hull–White caplet from the closed-form bond-put equivalent
[[nodiscard]] double hull_white_caplet(const market::DiscountCurve& curve,
                                       double start, double end, double strike,
                                       double mean_reversion, double volatility,
                                       double notional = 1.0);

}  
