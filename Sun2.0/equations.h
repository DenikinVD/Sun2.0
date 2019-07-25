#ifndef EQUATIONS_H_INCLUDED
#define EQUATIONS_H_INCLUDED

#include "cycle.h"
#include <cmath>

static Pos VDPD(double time, const Pos& pos, double lambda, double mu) {
    Pos res;
    time += 0;
    res.x = pos.y;
    res.y = -mu * (pos.x * pos.x - 1) * pos.y -
            pos.x * (1 + lambda * pos.x * pos.x);
    return res;
}

static Pos SQ_VDPD(double time, const Pos& pos, double lambda, double mu) {
    Pos res;
    time += 0;
    res.x = 1;
    res.y = -mu * (pos.x * pos.x - 1) * sqrt(2 * pos.y) -
            pos.x * (1 + lambda * pos.x * pos.x);
    return res;
}

static Pos mu_derivative(double time, const Pos& pos, double lambda, double mu, const Cycle& sol) {
    Pos sol_pos = sol[time];
    double x = sol_pos.x;
    double y = sol_pos.y;
    Pos res;
    res.x = pos.y;
    res.y = (-2 * mu * x * y - 3 * lambda * x * x - 1) * pos.x
            - mu * (x * x - 1) * pos.y - (x * x - 1) * y;
    return res;
}

static Pos lambda_derivative(double time, const Pos& pos, double lambda, double mu, const Cycle& sol) {
    Pos sol_pos = sol[time];
    double x = sol_pos.x;
    double y = sol_pos.y;
    Pos res;
    res.x = pos.y;
    res.y = (-2 * mu * x * y - 3 * lambda * x * x - 1) * pos.x
            - mu * (x * x - 1) * pos.y - x * x * x;
    return res;
}

static Pos denorm_VDPD(double time, const Pos& pos, double lambda, double mu, double epsilon, double omega) {
    Pos res;
    time += 0;
    res.x = pos.y;
    res.y = - mu * (epsilon * pos.x * pos.x - 1) * pos.y -
            omega * omega * pos.x - lambda * pos.x * pos.x * pos.x;
    return res;
}

#endif // EQUATIONS_H_INCLUDED
