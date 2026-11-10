# QuantEngine C++

**A foundational computational-finance research library in C++20.**

I thought to share this C++ library for quantitative finance, which I named *QuantEngine*. It takes a different approach from classic state-grid / classic discounting frameworks under the risk-neutral measure $\mathbb{Q}$: stochastic models produce paths, products evaluate those paths,
and generic engines aggregate discounted payoffs. For cross-validation and
numerical breadth, it also contains an independent one-dimensional PDE solver,
rate-option analytics, market calibration and scenario-risk tools.

The scope of the financial problems to address is (for now) limited, as the library grows!

## Components

| Module | Included functionality |
| --- | --- |
| **Market** | Discount factors, continuous zero/forward rates, log-linear discount interpolation, scalar curves |
| **Models** | Geometric Brownian motion, arithmetic Brownian motion, Vasicek/OU short rates, term-structure-aware one-factor Hull–White |
| **Products** | European call/put, arithmetic Asian, up-and-out barrier, zero-coupon bond |
| **Pricing** | Antithetic Monte Carlo, pathwise stochastic discounting, Longstaff–Schwartz American option approximation |
| **Analytics** | Black–Scholes and Bachelier formulas, Greeks, implied Black volatility |
| **Calibration** | Least-squares flat volatility across European option quotes |
| **Numerics** | Bracketed root solver, quadratic least-squares regression |
| **Multi-factor** | Correlated multi-asset GBM, Cholesky factorisation and basket payoffs |
| **PDE** | Theta-method Black–Scholes pricing; Crank–Nicolson European and projected-SOR American exercise |
| **Curves** | Simple deposits/par-swap discount bootstrapping and monotone cubic interpolation |
| **Rates** | Forward rates, swap valuation, caplet/floorlet and payer/receiver swaption pricing in Black or normal conventions; Hull–White caplets |
| **Volatility** | Raw-SVI smile fitting, total-variance term interpolation and calendar checks |
| **Model calibration** | Hull–White mean reversion/volatility fitting from caplet prices; Nelder–Mead optimiser |
| **Risk** | Richardson finite-difference Greeks, portfolio sensitivities, historical VaR/expected shortfall |

---

## Development Environment

The project follows the same toolchain and workflow as
[**Numerical Linear Algebra in C++20**](https://github.com/dom00mi/numerical-linear-algebra-cpp/tree/main).

| Component | Technology |
| --- | --- |
| Programming language | C++20 |
| Compiler | GCC / G++ |
| Build system | CMake + Ninja |
| Testing framework | CTest |
| Debugger | GDB |
| Editor | Visual Studio Code |
| Operating system | Ubuntu (native or WSL2) |

The quantitative-finance algorithms use the C++ Standard Library without an
external numerical framework such as GSL. CMake 3.21 or newer is required for
the version-controlled build presets.

## Setup and Installation

### 1. Install the development tools

**If you already set up Ubuntu/WSL for the numerical-linear-algebra project,
skip this installation step.** Otherwise, in an Ubuntu terminal:

```bash
sudo apt update
sudo apt install -y build-essential gdb cmake ninja-build git
```

Verify the toolchain when needed:

```bash
g++ --version
cmake --version
ninja --version
gdb --version
```

### 2. Configure Visual Studio Code

For Windows users, develop inside Ubuntu on **WSL2** using the following
VS Code extensions:

- WSL (Microsoft)
- C/C++ Extension Pack (Microsoft)
- CMake Tools (Microsoft)

This repository contains `.vscode/extensions.json` and
`.vscode/settings.json` to recommend those extensions and enable CMake presets.

### 3. Open the project in Ubuntu

Copy/extract the `quantengine-cpp` directory into your Ubuntu project folder.
For example, from your Ubuntu terminal:

```bash
mkdir -p ~/projects/cpp
cd ~/projects/cpp/quantengine-cpp
code .
```

The folder must already contain the extracted project files before running
`cd`. Once you create and push a GitHub repository, you can use `git clone`
instead of copying the files. Check the bottom-left VS Code status area for
the WSL connection.

### 4. Configure CMake

The checked-in `CMakePresets.json` provides identical development preset names
and output layout to your other C++ project.

For the Debug configuration:

```bash
cmake --preset wsl-debug
```

Build files are generated in:

```text
build/wsl-debug/
```

Alternatively, in VS Code:

1. Open the Command Palette (`Ctrl + Shift + P`).
2. Select `CMake: Select Configure Preset`.
3. Choose `wsl-debug`.
4. Run `CMake: Configure`.

### 5. Build the library

```bash
cmake --build --preset wsl-debug
```

The main outputs are:

```text
build/wsl-debug/
├── libquantengine.a
├── pricing_workbench
├── expanded_workbench
├── quantengine_tests
└── quantengine_extended_tests
```

### 6. Run the tests

```bash
ctest --preset wsl-debug --output-on-failure
```

CTest executes the foundation and expanded numerical regression suites.

### 7. Run the examples

```bash
./build/wsl-debug/pricing_workbench
./build/wsl-debug/expanded_workbench
```

The first example covers European, Asian and American options. The second
demonstrates correlated baskets, PDE pricing, curve bootstrapping, rate options,
SVI calibration and risk sensitivities.

### 8. Debugging in Visual Studio Code

1. Open `examples/pricing_workbench.cpp` or `examples/expanded_workbench.cpp`.
2. Set a breakpoint next to a line of C++ code.
3. Open the Command Palette (`Ctrl + Shift + P`).
4. Select `CMake: Set Launch/Debug Target`.
5. Choose `pricing_workbench` or `expanded_workbench`.
6. Select `CMake: Debug` to start GDB under Ubuntu.

### 9. Release configuration

For optimised builds:

```bash
cmake --preset wsl-release
cmake --build --preset wsl-release
ctest --preset wsl-release --output-on-failure
```

Use the Debug preset during development and Release when measuring execution
time or running larger quantitative experiments.

---

## Minimal example

```cpp
#include <iostream>
#include "quantengine/quantengine.hpp"

int main() {
    namespace qe = quantengine;

    const qe::models::GeometricBrownian model{100.0, 0.04, 0.01, 0.20};
    const auto dates = qe::models::TimeGrid::uniform(1.0, 50);
    const auto call = qe::payoffs::European{105.0, qe::payoffs::Right::Call};

    const auto result = qe::simulation::price(model, dates, call, 100'000);
    std::cout << "Price: " << result.price
              << " / Monte Carlo SE: " << result.standard_error << '\n';
}
```

The same `simulation::price` engine accepts an Asian, barrier or bond payoff.
The same model can be passed to `simulation::american_option` for discrete
early-exercise valuation by least-squares regression.

## Six expanded research modules

1. **Correlated simulation:** `multifactor::CorrelatedGBM` accepts asset-wise
   volatilities and a positive-definite correlation matrix; its basket pricer
   uses joint shocks and antithetic pairs.
2. **PDEs:** `pde::price` solves a one-dimensional Black–Scholes PDE for a
   European or discretised American exercise obstacle, returning price, delta
   and gamma from the numerical grid.
3. **Curves:** `curves::bootstrap_discount_curve` combines simple deposits and
   single-curve par swaps; `MonotoneCubicCurve` supplies shape-preserving scalar
   interpolation.
4. **Interest rates:** `rates` includes caplets, floorlets, swaptions, swaps and
   Hull–White Gaussian bond-option-equivalent caplets.
5. **Calibrations:** `volatility::calibrate_svi` fits a five-parameter SVI smile;
   `model_calibration::fit_hull_white_caplets` fits mean reversion and volatility.
6. **Risk:** `risk::numerical_greeks` uses Richardson extrapolation; the module
   also aggregates option positions and calculates scenario-tail measures.

The extended numerical tests benchmark all six areas. See
[Numerical limitations](docs/LIMITATIONS.md) before interpreting the outputs.

## Project structure

```text
quantengine-cpp/
├── include/quantengine/
│   ├── analytics.hpp
│   ├── calibration.hpp
│   ├── curves.hpp
│   ├── market.hpp
│   ├── model_calibration.hpp
│   ├── models.hpp
│   ├── monte_carlo.hpp
│   ├── multifactor.hpp
│   ├── numerics.hpp
│   ├── optimisation.hpp
│   ├── payoffs.hpp
│   ├── pde.hpp
│   ├── rates.hpp
│   ├── risk.hpp
│   ├── volatility_surface.hpp
│   └── quantengine.hpp
├── src/                     # Implementations
├── tests/                   # Deterministic regression benchmarks
├── examples/                # Pricing workbench
├── docs/                    # Architecture, model limitations, roadmap
├── cmake/                   # Consumer package config
├── .vscode/                 # WSL/CMake Tools recommendations and settings
├── CMakePresets.json        # wsl-debug and wsl-release
├── LICENSE                  # MIT license
└── CMakeLists.txt
```

## Using as a library

After installation with `cmake --install build`, a CMake consumer can use:

```cmake
find_package(QuantEngine CONFIG REQUIRED)
target_link_libraries(my_pricer PRIVATE QuantEngine::QuantEngine)
```

## Financial conventions

- Continuous annualised rates and volatilities (0.20 = 20%).
- Year-fraction times, strictly increasing from 0 for simulation grids.
- `payoffs` return undiscounted amounts, while models accumulate stochastic
  discount factors and pricing engines return present values.
- The normal model allows negative forwards/strikes.
- Monte Carlo standard errors account for antithetic pairing; American option
  standard errors do **not** capture LSM regression or exercise-policy bias.
- Hull–White paths use exact OU factor transitions and trapezoidal integration
  for stochastic discounting. Accuracy of bond reconstruction depends on time
  discretisation and curve smoothness.
- PDE values and sensitivities depend on the spatial/time grids; American PDE
  exercise is handled by a projected iterative solver.
- SVI total-variance fitting and calendar checks **do not** guarantee absence
  of butterfly or calendar arbitrage between sampled strikes and maturities.
- The rates module assumes single-curve discounting and regular fixed coupons.

## Project maturity and scope

This is a **tested educational/research foundation**, and a starting point to obtain some more in-depth models and possibly satisfactory results, spanning across a variety of areas, from stochastic interest rate modelling to multi-asset approaches. The
current PDE engine is one-dimensional; correlated paths cover GBM rather than
every stochastic process. Calibration routines are constrained illustrative
models and need further robustness, arbitrage and benchmark validation.

See [Architecture](docs/ARCHITECTURE.md), [Numerical limitations](docs/LIMITATIONS.md)
and [Roadmap](docs/ROADMAP.md).

---

## Future Development

Planned improvements include multi-factor interest-rate dynamics, more flexible
curve conventions, PDE convergence studies, volatility-surface arbitrage
checks, robust calibration diagnostics and risk benchmarking.

## Applications

The library can support future mathematical-modelling and quantitative-finance
projects. Its stand-alone CMake package can also be linked from separate C++
applications, including the existing numerical-linear-algebra project if needed.

## Project Status

**Work in progress.** The library serves only the purpose only in the context of education, experimentation,
and reproducible quantitative-finance research.

## License

Copyright (c) 2026 Domenico Nava.

This project is licensed under the MIT License.
See [LICENSE](LICENSE) for details.
