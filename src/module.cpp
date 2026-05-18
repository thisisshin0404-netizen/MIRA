#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <cmath>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

namespace py = pybind11;

struct Vec2 {
    double x;
    double y;
};

struct ForceResult {
    std::vector<Vec2> forces;
    double potential_energy;
};

static std::vector<Vec2> py_to_vec2(const std::vector<std::vector<double>>& input,
                                    const std::string& name) {
    std::vector<Vec2> output;
    output.reserve(input.size());

    for (std::size_t i = 0; i < input.size(); ++i) {
        if (input[i].size() != 2) {
            throw std::invalid_argument(name + " must be a list of [x, y] pairs.");
        }
        output.push_back(Vec2{input[i][0], input[i][1]});
    }

    return output;
}

static std::vector<std::vector<double>> vec2_to_py(const std::vector<Vec2>& input) {
    std::vector<std::vector<double>> output;
    output.reserve(input.size());

    for (const auto& v : input) {
        output.push_back({v.x, v.y});
    }

    return output;
}

static void validate_positive(double value, const std::string& name) {
    if (value <= 0.0) {
        throw std::invalid_argument(name + " must be positive.");
    }
}

static void apply_pbc(Vec2& r, double box_length) {
    r.x -= box_length * std::floor(r.x / box_length);
    r.y -= box_length * std::floor(r.y / box_length);

    // Protect against rare floating-point edge cases.
    if (r.x >= box_length) r.x -= box_length;
    if (r.y >= box_length) r.y -= box_length;
    if (r.x < 0.0) r.x += box_length;
    if (r.y < 0.0) r.y += box_length;
}

static Vec2 minimum_image(Vec2 dr, double box_length) {
    dr.x -= box_length * std::round(dr.x / box_length);
    dr.y -= box_length * std::round(dr.y / box_length);
    return dr;
}

static bool is_perfect_square(std::size_t n) {
    const auto root = static_cast<std::size_t>(std::llround(std::sqrt(static_cast<double>(n))));
    return root * root == n;
}

static std::vector<std::vector<double>> initialize_square_positions(std::size_t n_particles,
                                                                    double box_length) {
    if (n_particles == 0) {
        throw std::invalid_argument("n_particles must be greater than zero.");
    }
    if (!is_perfect_square(n_particles)) {
        throw std::invalid_argument(
            "For position_rule='square', n_particles must be a perfect square: 4, 9, 16, 25, ...");
    }

    const auto n_side = static_cast<std::size_t>(std::llround(std::sqrt(static_cast<double>(n_particles))));
    const double spacing = box_length / static_cast<double>(n_side);

    std::vector<std::vector<double>> positions;
    positions.reserve(n_particles);

    for (std::size_t iy = 0; iy < n_side; ++iy) {
        for (std::size_t ix = 0; ix < n_side; ++ix) {
            const double x = (static_cast<double>(ix) + 0.5) * spacing;
            const double y = (static_cast<double>(iy) + 0.5) * spacing;
            positions.push_back({x, y});
        }
    }

    return positions;
}

static std::vector<std::vector<double>> initialize_random_positions(std::size_t n_particles,
                                                                    double box_length,
                                                                    unsigned int seed) {
    if (n_particles == 0) {
        throw std::invalid_argument("n_particles must be greater than zero.");
    }

    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> dist(0.0, box_length);

    std::vector<std::vector<double>> positions;
    positions.reserve(n_particles);

    for (std::size_t i = 0; i < n_particles; ++i) {
        positions.push_back({dist(rng), dist(rng)});
    }

    return positions;
}

static std::vector<std::vector<double>> initialize_positions(std::size_t n_particles,
                                                            double box_length,
                                                            const std::string& position_rule,
                                                            unsigned int seed) {
    validate_positive(box_length, "box_length");

    if (position_rule == "square") {
        return initialize_square_positions(n_particles, box_length);
    }

    if (position_rule == "random") {
        return initialize_random_positions(n_particles, box_length, seed);
    }

    throw std::invalid_argument("Unknown position_rule. Use 'square' or 'random'.");
}

static std::vector<std::vector<double>> initialize_velocities(std::size_t n_particles) {
    if (n_particles == 0) {
        throw std::invalid_argument("n_particles must be greater than zero.");
    }
    return std::vector<std::vector<double>>(n_particles, std::vector<double>{0.0, 0.0});
}

static py::dict initialize_system(std::size_t n_particles,
                                  double box_length,
                                  const std::string& position_rule = "square",
                                  unsigned int seed = 42) {
    auto positions = initialize_positions(n_particles, box_length, position_rule, seed);
    auto velocities = initialize_velocities(n_particles);

    py::dict result;
    result["positions"] = positions;
    result["velocities"] = velocities;
    result["box_length"] = box_length;
    result["n_particles"] = n_particles;
    result["position_rule"] = position_rule;
    result["seed"] = seed;
    return result;
}

static ForceResult compute_lj_forces(const std::vector<Vec2>& positions,
                                     double box_length,
                                     double epsilon,
                                     double sigma) {
    validate_positive(box_length, "box_length");
    validate_positive(epsilon, "epsilon");
    validate_positive(sigma, "sigma");

    const std::size_t n = positions.size();
    std::vector<Vec2> forces(n, Vec2{0.0, 0.0});
    double potential_energy = 0.0;

    const double sigma2 = sigma * sigma;

    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = i + 1; j < n; ++j) {
            Vec2 dr{positions[j].x - positions[i].x,
                    positions[j].y - positions[i].y};
            dr = minimum_image(dr, box_length);

            const double r2 = dr.x * dr.x + dr.y * dr.y;
            if (r2 <= 0.0) {
                throw std::runtime_error("Two particles occupy the same position or are numerically overlapping.");
            }

            const double inv_r2_scaled = sigma2 / r2;
            const double inv_r6 = inv_r2_scaled * inv_r2_scaled * inv_r2_scaled;
            const double inv_r12 = inv_r6 * inv_r6;

            potential_energy += 4.0 * epsilon * (inv_r12 - inv_r6);

            const double force_over_r = 24.0 * epsilon * (2.0 * inv_r12 - inv_r6) / r2;
            const double fx = force_over_r * dr.x;
            const double fy = force_over_r * dr.y;

            forces[i].x += fx;
            forces[i].y += fy;
            forces[j].x -= fx;
            forces[j].y -= fy;
        }
    }

    return ForceResult{forces, potential_energy};
}

static double compute_kinetic_energy(const std::vector<Vec2>& velocities) {
    double kinetic_energy = 0.0;
    for (const auto& v : velocities) {
        kinetic_energy += 0.5 * (v.x * v.x + v.y * v.y); // mass = 1
    }
    return kinetic_energy;
}

static py::dict step(const std::vector<std::vector<double>>& positions_in,
                     const std::vector<std::vector<double>>& velocities_in,
                     double dt,
                     double box_length,
                     double epsilon,
                     double sigma) {
    validate_positive(dt, "dt");
    validate_positive(box_length, "box_length");

    auto positions = py_to_vec2(positions_in, "positions");
    auto velocities = py_to_vec2(velocities_in, "velocities");

    if (positions.size() != velocities.size()) {
        throw std::invalid_argument("positions and velocities must contain the same number of particles.");
    }
    if (positions.empty()) {
        throw std::invalid_argument("At least one particle is required.");
    }

    auto force_result = compute_lj_forces(positions, box_length, epsilon, sigma);

    for (std::size_t i = 0; i < positions.size(); ++i) {
        // mass = 1, so acceleration = force
        velocities[i].x += force_result.forces[i].x * dt;
        velocities[i].y += force_result.forces[i].y * dt;

        positions[i].x += velocities[i].x * dt;
        positions[i].y += velocities[i].y * dt;

        apply_pbc(positions[i], box_length);
    }

    const double kinetic_energy = compute_kinetic_energy(velocities);

    py::dict result;
    result["positions"] = vec2_to_py(positions);
    result["velocities"] = vec2_to_py(velocities);
    result["potential_energy"] = force_result.potential_energy;
    result["kinetic_energy"] = kinetic_energy;
    result["total_energy"] = force_result.potential_energy + kinetic_energy;
    return result;
}

PYBIND11_MODULE(MIRA, m) {
    m.doc() = "Minimal 2D Lennard-Jones molecular dynamics module";

    m.def("initialize_positions", &initialize_positions,
          py::arg("n_particles"),
          py::arg("box_length"),
          py::arg("position_rule") = "square",
          py::arg("seed") = 42,
          "Initialize 2D particle positions using 'square' or 'random'.");

    m.def("initialize_system", &initialize_system,
          py::arg("n_particles"),
          py::arg("box_length"),
          py::arg("position_rule") = "square",
          py::arg("seed") = 42,
          "Initialize positions and zero velocities.");

    m.def("step", &step,
          py::arg("positions"),
          py::arg("velocities"),
          py::arg("dt"),
          py::arg("box_length"),
          py::arg("epsilon"),
          py::arg("sigma"),
          "Advance the system by one semi-implicit Euler step.");
}
