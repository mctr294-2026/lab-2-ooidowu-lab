#include "roots.hpp"

#include <cmath>
#include <limits>

static constexpr double TOL = 1e-6;
static constexpr int MAX_ITER = 1'000'000;

// BISECTION METHOD

bool bisection(std::function<double(double)> f,
               double a, double b,
               double *root)
{
    if (!root) return false;

    double fa = f(a);
    double fb = f(b);

    if (std::abs(fa) < TOL) { *root = a; return true; }
    if (std::abs(fb) < TOL) { *root = b; return true; }

    if (fa * fb > 0.0) return false;

    for (int i = 0; i < MAX_ITER; ++i) {
        double m  = 0.5 * (a + b);
        double fm = f(m);

        if (std::abs(fm) < TOL || std::abs(b - a) < TOL) {
            *root = m;
            return true;
        }

        if (fa * fm < 0.0) {
            b = m; fb = fm;
        } else {
            a = m; fa = fm;
        }
    }
    return false;
}


// REGULA FALSI

bool regula_falsi(std::function<double(double)> f,
                  double a, double b,
                  double *root)
{
    if (!root) return false;

    double fa = f(a);
    double fb = f(b);

    if (std::abs(fa) < TOL) { *root = a; return true; }
    if (std::abs(fb) < TOL) { *root = b; return true; }

    if (fa * fb > 0.0) return false;

    for (int i = 0; i < MAX_ITER; ++i) {
        double denom = fb - fa;
        if (std::abs(denom) < std::numeric_limits<double>::epsilon())
            return false;

        double c  = (a * fb - b * fa) / denom;
        double fc = f(c);

        if (std::abs(fc) < TOL) {
            *root = c;
            return true;
        }

        if (fa * fc < 0.0) {
            b = c; fb = fc;
        } else {
            a = c; fa = fc;
        }

        if (std::abs(b - a) < TOL) {
            *root = 0.5 * (a + b);
            return true;
        }
    }
    return false;
}

// NEWTON–RAPHSON

bool newton_raphson(std::function<double(double)> f,
                    std::function<double(double)> g,
                    double a, double b, double c,
                    double *root)
{
    if (!root) return false;

    double x = c;
    if (x < a) x = a;
    if (x > b) x = b;

    for (int i = 0; i < MAX_ITER; ++i) {
        double fx = f(x);
        if (std::abs(fx) < TOL) {
            *root = x;
            return true;
        }

        double gx = g(x);
        if (std::abs(gx) < 1e-12) return false;

        double x_next = x - fx / gx;

        if (x_next < a) x_next = a;
        if (x_next > b) x_next = b;

        if (std::abs(x_next - x) < TOL) {
            *root = x_next;
            return true;
        }

        x = x_next;
    }
    return false;
}


// SECANT METHOD

bool secant(std::function<double(double)> f,
            double a, double b, double c,
            double *root)
{
    if (!root) return false;

    double x0 = c;
    if (x0 < a) x0 = a;
    if (x0 > b) x0 = b;

    double x1 = x0 + 1e-3;
    if (x1 > b) x1 = x0 - 1e-3;
    if (x1 < a) x1 = a;

    double f0 = f(x0);
    double f1 = f(x1);

    for (int i = 0; i < MAX_ITER; ++i) {
        if (std::abs(f1) < TOL) {
            *root = x1;
            return true;
        }

        double denom = f1 - f0;
        if (std::abs(denom) < 1e-12) return false;

        double x2 = x1 - f1 * (x1 - x0) / denom;

        if (x2 < a) x2 = a;
        if (x2 > b) x2 = b;

        if (std::abs(x2 - x1) < TOL) {
            *root = x2;
            return true;
        }

        x0 = x1; f0 = f1;
        x1 = x2; f1 = f(x1);
    }
    return false;
}
