#pragma once

#include <cstddef>
#include <vector>

namespace quantengine::volatility {

// Raw SVI total implied variance: a + b[rho(k-m)+sqrt((k-m)^2+sigma^2)]
struct SviParameters {
  double a;
  double b;
  double rho;
  double m;
  double sigma;
};

struct SmileQuote {
  double log_moneyness;
  double implied_volatility;
  double weight = 1.0;
};

struct SviFit {
  SviParameters parameters;
  double total_variance_rmse;
  bool converged;
};

[[nodiscard]] double svi_total_variance(double log_moneyness,
                                         const SviParameters& parameters);

[[nodiscard]] SviFit calibrate_svi(
    double maturity, const std::vector<SmileQuote>& quotes,
    SviParameters initial = {0.02, 0.1, -0.3, 0, 0.2});

struct SviSlice {
  double maturity;
  SviParameters parameters;
};

class SviSurface {
 public:
  explicit SviSurface(std::vector<SviSlice> smiles);
  [[nodiscard]] double total_variance(double k, double maturity) const;
  [[nodiscard]] double implied_volatility(double k, double maturity) const;
  // Checks calendar monotonicity numerically on the specified strike range.
  [[nodiscard]] bool calendar_monotone(double min_k = -1, double max_k = 1,
                                       std::size_t samples = 101) const;

 private:
  std::vector<SviSlice> smiles_;
};

} 
