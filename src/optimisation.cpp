#include "quantengine/optimisation.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

namespace quantengine::optimisation {
namespace {

struct Vertex {
  std::vector<double> x;
  double loss;
};

}

Minimum nelder_mead(
    const std::function<double(const std::vector<double>&)>& objective,
    const std::vector<double>& initial, const std::vector<double>& steps,
    std::size_t max_iterations, double tolerance) {
  const std::size_t dim = initial.size();
  if (dim == 0 || steps.size() != dim || max_iterations == 0 ||
      !std::isfinite(tolerance) || tolerance <= 0)
    throw std::invalid_argument("Invalid Nelder-Mead configuration");
  for (std::size_t i = 0; i < dim; ++i)
    if (!std::isfinite(initial[i]) || !std::isfinite(steps[i]) || steps[i] <= 0)
      throw std::invalid_argument("Invalid initial simplex");

  const auto evaluate = [&](std::vector<double> point) {
    double loss = objective(point);
    if (!std::isfinite(loss)) loss = std::numeric_limits<double>::infinity();
    return Vertex{std::move(point), loss};
  };
  std::vector<Vertex> simplex;
  simplex.reserve(dim + 1);
  simplex.push_back(evaluate(initial));
  for (std::size_t i = 0; i < dim; ++i) {
    auto point = initial;
    point[i] += steps[i];
    simplex.push_back(evaluate(std::move(point)));
  }

  for (std::size_t iteration = 0; iteration < max_iterations; ++iteration) {
    std::sort(simplex.begin(), simplex.end(),
              [](const auto& a, const auto& b) { return a.loss < b.loss; });
    double spread = 0;
    double distance = 0;
    for (std::size_t k = 1; k <= dim; ++k) {
      spread = std::max(spread, std::abs(simplex[k].loss - simplex[0].loss));
      for (std::size_t j = 0; j < dim; ++j)
        distance = std::max(distance,
                            std::abs(simplex[k].x[j] - simplex[0].x[j]));
    }
    if (spread < tolerance && distance < std::sqrt(tolerance))
      return {simplex[0].x, simplex[0].loss, iteration, true};

    std::vector<double> centroid(dim, 0);
    for (std::size_t k = 0; k < dim; ++k)
      for (std::size_t j = 0; j < dim; ++j)
        centroid[j] += simplex[k].x[j] / dim;
    const auto project = [&](double multiple, const std::vector<double>& x) {
      std::vector<double> point(dim);
      for (std::size_t j = 0; j < dim; ++j)
        point[j] = centroid[j] + multiple * (x[j] - centroid[j]);
      return evaluate(std::move(point));
    };
    auto reflected = project(-1, simplex.back().x);
    if (reflected.loss < simplex.front().loss) {
      auto expanded = project(2, reflected.x);
      simplex.back() = expanded.loss < reflected.loss
                           ? std::move(expanded) : std::move(reflected);
    } else if (reflected.loss < simplex[dim - 1].loss) {
      simplex.back() = std::move(reflected);
    } else {
      // Outside contraction if reflection improves the worst vertex, otherwise
      // inside contraction towards the centroid.
      auto contracted = reflected.loss < simplex.back().loss
                            ? project(0.5, reflected.x)
                            : project(0.5, simplex.back().x);
      const double limit = std::min(reflected.loss, simplex.back().loss);
      if (contracted.loss < limit) {
        simplex.back() = std::move(contracted);
      } else {
        for (std::size_t k = 1; k <= dim; ++k) {
          std::vector<double> shrunk(dim);
          for (std::size_t j = 0; j < dim; ++j)
            shrunk[j] = simplex[0].x[j] +
                        0.5 * (simplex[k].x[j] - simplex[0].x[j]);
          simplex[k] = evaluate(std::move(shrunk));
        }
      }
    }
  }
  std::sort(simplex.begin(), simplex.end(),
            [](const auto& a, const auto& b) { return a.loss < b.loss; });
  return {simplex.front().x, simplex.front().loss, max_iterations, false};
}

}
