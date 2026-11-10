# Numerical and financial limitations

This library is for **education, prototyping and research**.

1. **Monte Carlo:** Antithetic pairs can reduce variance but do not remove
   simulation error; pricing returns the sampling standard error only. The
   pseudo-random generator is deterministic for a chosen toolchain/seed, but
   standard-library distribution implementations may differ across platforms.
2. **Barrier monitoring:** Barriers are checked on discrete dates, so continuous
   barrier pricing has monitoring bias. No Brownian-bridge correction yet.
3. **Asian averaging:** The average uses simulated observations after time 0,
   equally weighted. It does not automatically weight irregular intervals.
4. **American exercise:** Longstaff–Schwartz uses an in-sample quadratic
   regression. It can have exercise-policy bias; reported SE ignores the
   regression uncertainty. Validate with trees/PDEs before practical use.
5. **Hull–White:** The model shifts a Gaussian OU factor to the initial curve,
   but approximates the discount integral by the trapezoidal rule. Curve
   forward rates derive from finite-width adjacent discount factors; at pillar
   kinks or under a coarse grid, fit to the initial curve is approximate.
6. **Multi-factor simulation:** Correlated GBM is implemented for a positive-
   definite correlation matrix. Exactly singular matrices are rejected. There
   is no stochastic covariance structure, stochastic volatility, local
   volatility or multi-factor stochastic discounting in this module yet.
7. **PDE:** The Black–Scholes solver is one-dimensional and uses uniform
   stock/time grids with a truncated upper spot boundary. Crank–Nicolson can
   oscillate around nonsmooth payoffs on coarse grids; no Rannacher smoothing,
   nonuniform meshes or PDE convergence certification is included. American
   exercise uses projected SOR and may fail to converge for unusual settings.
8. **Interest rates:** Caplet and swaption analytics use a single discount /
   projection curve, constant quoted volatility, regular payment periods,
   continuously annualised input volatilities and year fractions. They do
   not implement market day counts, curves under collateralisation or smile
   conventions. Hull–White caplets use the Gaussian bond-option equivalence.
9. **Market conventions:** Deposits use simple compounding; par swaps use
   uniform fixed coupons and a basic one-curve floating-leg identity.
   There are no date calendars, OIS/IBOR multi-curve conventions, day-count
   schedules or complete arbitrage repairs. The monotone cubic interpolant is
   a scalar interpolant, not a positive-discount-curve guarantee.
10. **Calibration:** SVI fits a restricted raw parameterisation with positive
    `a`, `b` and `sigma` and `|rho| < 1`. Calendar monotonicity is checked
    only over a user-specified log-moneyness mesh. No guarantee of butterfly
    arbitrage-freeness or robust global identification is provided.
    Hull–White parameters are fitted to a small caplet dataset, which may
    be weakly identifiable or sensitive to quote selection.
11. **Risk:** Bump-and-revalue estimates depend on bump size and pricing
    noise. Historical VaR/ES depend on the scenario distribution and simple
    empirical quantiles.

Use unit tests, convergence studies and independent financial benchmarks when
introducing each new model or contract.
