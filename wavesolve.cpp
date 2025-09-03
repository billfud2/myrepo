#include <iostream>
#include <vector>

int main() {
    // Create WaveOrthotope
    auto rows = 25, cols = 50;
    auto c = 0.01;
    auto u = std::vector<double>(rows * cols, 0.0);
    auto v = std::vector<double>(rows * cols, 0.0);
    auto w = WaveOrthotope(rows, cols, c);
    for (size_t i = 1; i < rows - 1; ++i) {
        for (size_t j = 1; j < cols - 1; ++j) {
            v[i*cols + j] = 0.1;
        }
    }
    w.solve();

    std::cout << w.sim_time() << std::endl;
    return 0;
}

class WaveOrthotope {
protected:
    const size_t rows, cols;  // size
    const double c;           // damping coefficient
    double t;                 // simulation time
    std::vector<double> u, v; // displacement and velocity; size is rows*cols

public:
    WaveOrthotope(auto rows, auto cols, auto damping_coefficient);

    auto &displacement(auto i, auto j) { return u[i*cols+j]; }
    auto &velocity(    auto i, auto j) { return v[i*cols+j]; }

    auto sim_time() const { return t; }

    double energy(){
        double E = 0.0;
        // Dynamic energy
        for (size_t i = 1; i < rows - 1; ++i) {
            for (size_t j = 1; j < cols - 1; ++j) {
                E += v[i * cols + j] * v[i * cols + j] / 2.0;
            }
        }
        // Potential energy along x axis
        for (size_t i = 0; i < rows - 1; ++i) {
            for (size_t j = 1; j < cols - 1; ++j) {
                double diff = u[i * cols + j] - u[(i + 1) * cols + j];
                E += (diff * diff) / 4.0;
            }
        }
        // Potential energy along y axis
        for (size_t i = 1; i < rows - 1; ++i) {
            for (size_t j = 0; j < cols - 1; ++j) {
                double diff = u[i * cols + j] - u[i * cols + (j + 1)];
                E += (diff * diff) / 4.0;
            }
        }
        return E;
    }

    double step(double dt) {
        // Update velocity v
        for (size_t i = 1; i < rows - 1; ++i) {
            for (size_t j = 1; j < cols - 1; ++j) {
                double L = (u[(i-1)*cols + j] + u[(i+1)*cols + j] + u[i*cols + (j-1)] + u[i*cols + (j+1)]) / 2.0
                           - 2.0 * u[i*cols + j];
                v[i*cols + j] = (1.0 - dt * c) * v[i*cols + j] + dt * L;
            }
        }
        // Update displacement u
        for (size_t i = 1; i < rows - 1; ++i) {
            for (size_t j = 1; j < cols - 1; ++j) {
                u[i*cols + j] += v[i*cols + j] * dt;
            }
        }
        // Enforce boundary conditions (u and v are zero at boundaries)
        for (size_t i = 0; i < rows; ++i) {
            // Set first and last column to 0.0
            u[i*cols + 0] = 0.0;
            u[i*cols + (cols - 1)] = 0.0;
            v[i*cols + 0] = 0.0;
            v[i*cols + (cols - 1)] = 0.0;
        }
        // Set first and last row to 0.0
        for (size_t j = 0; j < cols; ++j) {
            u[0 * cols + j] = 0.0;
            u[(rows - 1) * cols + j] = 0.0;
            v[0 * cols + j] = 0.0;
            v[(rows - 1) * cols + j] = 0.0;
        }
        t += dt;
        return t;
    }

    double solve(){
        double dt = 0.01;
        double stopping_energy = (rows - 2) * (cols - 2) * 0.001;
        while (energy() > stopping_energy) {
            step(dt);
        }
        return t;
    }
}
