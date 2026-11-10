#pragma once

#include <cstddef>

#include "quantengine/analytics.hpp"

namespace quantengine::pde {

struct GridSettings {
  std::size_t price_steps = 240;
  std::size_t time_steps = 240;
  double upper_spot = 0;       // Zero selects an automatic truncation.
  double theta = 0.5;         // 0.5 = Crank-Nicolson; 1 = implicit Euler.
  double sor_relaxation = 1.2;
  double tolerance = 1e-9;
  std::size_t sor_max_iterations = 20'000;
};

struct Result {
  double price;
  double delta;
  double gamma;
  double spatial_step;
};

// Black-Scholes 1D PDE on a truncated price grid. American exercise uses
// projected SOR; European options use a tridiagonal direct solve.
[[nodiscard]] Result price(const analytics::EuropeanMarket& market,
                            payoffs::Right right, bool american = false,
                            GridSettings settings = {});

}  
