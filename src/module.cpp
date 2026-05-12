#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <array>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace py = pybind11;

struct Vec2 {
    double x;
    double y;
};

static inline Vec2 add(const Vec2& a, const Vec2& b) {
    return {a.x + b.x, a.y + b.y};
}

static inline Vec2 sub(const Vec2& a, const Vec2& b) {
    return {a.x - b.x, a.y - b.y};
}

static inline Vec2 mul(const Vec2& a, double s) {
    return {a.x * s, a.y * s};
}

static inline double norm2(const Vec2& a) {
    return a.x * a.x + a.y * a.y;
}

static inline void apply_pbc(Vec2& r, double L) {
    r.x -= L * std::floor(r.x / L);
    r.y -= L * std::floor(r.y / L);
}

static inline Vec2 minimum_image(Vec2 dr, double L) {
    dr.x -= L * std::round(dr.x / L);
    dr.y -= L * std::round(dr.y / L);
    return dr;
}

py::dict step(
    std::vector<std::array<double, 2>> positions_in,
    std::vector<std::array<double, 2>> velocities_in,
    double dt,
    double box_length,
    double epsilon = 1.0,
    double sigma = 1.0
) {
    const std::size_t n = positions_in.size();

    if (velocities_in.size() != n) {
        throw std::runtime_error("positions와 velocities의 길이가 다릅니다.");
    }
    if (dt <= 0.0) {
        throw std::runtime_error("dt는 양수여야 합니다.");
    }
    if (box_length <= 0.0) {
        throw std::runtime_error("box_length는 양수여야 합니다.");
    }

    std::vector<Vec2> r(n), v(n), f(n, {0.0, 0.0});

    for (std::size_t i = 0; i < n; ++i) {
        r[i] = {positions_in[i][0], positions_in[i][1]};
        v[i] = {velocities_in[i][0], velocities_in[i][1]};
        apply_pbc(r[i], box_length);
    }

    double potential_energy = 0.0;

    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = i + 1; j < n; ++j) {
            Vec2 dr = sub(r[j], r[i]);
            dr = minimum_image(dr, box_length);

            const double r2 = norm2(dr);
            if (r2 < 1e-12) {
                continue;
            }

            const double inv_r2 = 1.0 / r2;
            const double sig2_over_r2 = (sigma * sigma) * inv_r2;
            const double sr6 = sig2_over_r2 * sig2_over_r2 * sig2_over_r2;
            const double sr12 = sr6 * sr6;

            const double uij = 4.0 * epsilon * (sr12 - sr6);
            potential_energy += uij;

            const double coeff = 24.0 * epsilon * inv_r2 * (2.0 * sr12 - sr6);
            Vec2 fij = mul(dr, coeff);

            f[i] = sub(f[i], fij);
            f[j] = add(f[j], fij);
        }
    }

    for (std::size_t i = 0; i < n; ++i) {
        v[i] = add(v[i], mul(f[i], dt));
        r[i] = add(r[i], mul(v[i], dt));
        apply_pbc(r[i], box_length);
    }

    double kinetic_energy = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        kinetic_energy += 0.5 * norm2(v[i]);
    }

    std::vector<std::array<double, 2>> positions_out(n), velocities_out(n);
    for (std::size_t i = 0; i < n; ++i) {
        positions_out[i] = {r[i].x, r[i].y};
        velocities_out[i] = {v[i].x, v[i].y};
    }

    py::dict out;
    out["positions"] = positions_out;
    out["velocities"] = velocities_out;
    out["potential_energy"] = potential_energy;
    out["kinetic_energy"] = kinetic_energy;
    return out;
}

PYBIND11_MODULE(MIRA, m) {
    m.def("step",
          &step,
          py::arg("positions"),
          py::arg("velocities"),
          py::arg("dt"),
          py::arg("box_length"),
          py::arg("epsilon") = 1.0,
          py::arg("sigma") = 1.0);
}