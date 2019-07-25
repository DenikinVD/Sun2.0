#ifndef APPROX_H_INCLUDED
#define APPROX_H_INCLUDED

#include <algorithm>
#include "rk4.h"
#include "equations.h"
#include <set>
#include <tuple>
#include <future>

const double amplitude_approx_step = 1.0 / static_cast<double>(1 << 17);
const double parameter_step = 1.0 / static_cast<double>(1 << 10);
const double approx_time_step = 1.0 / static_cast<double>(1 << 18);
const double regular_time_step = 1.0 / static_cast<double>(1 << 13);
const double close_enough = 1.0 / static_cast<double>(1 << 8);
const double min_amplitude = static_cast<double>(249037) * amplitude_approx_step;
const double max_amplitude = static_cast<double>(275251) * amplitude_approx_step;
// получено из тестирующего запуска с большим шагом

struct Parameters {
    double omega = 1;
    double epsilon = 1;
    double lambda = 0;
    double mu = 0;
    double amplitude = 0;
    double period = 0;
    void normalize() { //нормирует так, что omega и epsilon = 1
        lambda /= (epsilon * omega * omega);
        mu /= omega;
        period *= omega;
        amplitude *= sqrt(epsilon);
        omega = 1;
        epsilon = 1;
    }

    void renormalize(double _omega, double _epsilon) { //меняет нормировку
        this->normalize();
        omega = _omega;
        epsilon = _epsilon;
        lambda *= (epsilon * omega * omega);
        mu *= omega;
        period /= omega;
        amplitude /= sqrt(epsilon);
    }
};

static Pos approx_start_point(double lambda, double mu) {
    Pos left (-max_amplitude, 0);
    Pos right (-min_amplitude, 0);
    Limit positive_velocity {neg_inf, 0, pos_inf, pos_inf, pos_inf};
    while ((right - left).x > amplitude_approx_step) {
        Pos start = (left + right) / 2;
        Cycle cur_res = RK4(amplitude_approx_step, 0, start, SQ_VDPD, positive_velocity, lambda, mu);
        Pos last = cur_res.last_pos();
        if (fabs(last.x) > fabs(start.x)) {
            right = start;
        } else {
            left = start;
        }
    }
    return left;
}

static void calculate_init_positions(std::ofstream& out) {
    Limit positive_velocity {neg_inf, 0, pos_inf, pos_inf, pos_inf};
    for (double mu = 4; mu >= parameter_step; mu -= parameter_step) {
        double lambda = 4 - mu;
        Pos start = approx_start_point(lambda, mu);
        Cycle c = RK4(approx_time_step, 0, start, VDPD, positive_velocity, lambda, mu);
        out << std::fixed << std::setprecision(12) << lambda << " " << mu << " " << start.x << " " << c.period() << "\n";
    }
}

static std::vector<Parameters> extract_init_pos(std::ifstream& in) {
    Parameters p;
    std::vector<Parameters> res;
    while (in >> p.lambda >> p.mu >> p.amplitude >> p.period) {
        p.amplitude = -p.amplitude;
        res.push_back(p);
    }
    return res;
}

//меняет cycle, чтобы был предельный цикл с тем же периодом
static void approx_to_lim_cycle(Cycle& cycle, double& lambda, double& mu) {//предполагается, что cycle не сильно отличается от предельного цикла
    Cycle exept_cycle = cycle;
    double exept_mu = mu;
    double exept_lambda = lambda;
    double lambda_step = 0.5;
    double mu_step = 0.5;

    Pos start(0, 0);
    Limit t_lim {neg_inf, neg_inf, pos_inf, pos_inf, cycle.get_time_finish()};
    Pos approxible_end = -cycle.first_pos();
    size_t max_steps_num = 15;
    size_t counter = 0;
    double dist = distance(cycle.last_pos(), approxible_end);
    while (dist >= close_enough && counter < max_steps_num) {// критерий слабее равенства для экономии времени
        Pos div_mu = RK4(regular_time_step, cycle.get_time_start(), start,
                         mu_derivative, t_lim, lambda, mu, cycle).last_pos();
        Pos div_lambda = RK4(regular_time_step, cycle.get_time_start(), start,
                             lambda_derivative, t_lim, lambda, mu, cycle).last_pos();
        Pos diff = approxible_end - cycle.last_pos();
        Pos koefs = best_kramer(div_mu, div_lambda, diff);

        // защита от резких шагов
        while (lambda < fabs(lambda_step)) {
            lambda_step /= 2;
        }
        while (mu < fabs(mu_step)) {
            mu_step /= 2;
        }
        if (fabs(koefs.x) > mu_step || fabs(koefs.y) > lambda_step) {
            koefs *= std::min(mu_step / fabs(koefs.x), lambda_step / fabs(koefs.y));
        }

        if (lambda < 0.0000001 || mu < 0.0000001) {
            throw "too_small_parameters";
        }

        mu += koefs.x;
        lambda += koefs.y;
        cycle = RK4(regular_time_step, cycle.get_time_start(), cycle.first_pos(), VDPD, t_lim, lambda, mu);
        dist = distance(cycle.last_pos(), approxible_end);
        if (lambda > 5 || mu > 5) {
            lambda = exept_lambda;
            cycle = exept_cycle;
            mu = exept_mu;
            throw "too_big_parameters";
        }
        ++counter;
    }
}

static double integ_difference(const Cycle& a, const Cycle& b) {
    double res = 0;
    for (double time = b.get_time_start(); time < b.get_time_finish(); time += b.get_step()) {
        res += (a[time].x - b[time].x) * (a[time].x - b[time].x) * b.get_step();
    }
    return res;
}

static double dist_difference(const Cycle& a, const Cycle& b) {
    double i = 0;
    double res = 0;
    for (double time = b.get_time_start(); time < b.get_time_finish(); time += b.get_step()) {
        res += distance(a[time], b[time]);
        i += 1;
    }
    return res / i;
}


static bool comparator(const Parameters& a, const Parameters& b) {
    return a.period > b.period;
}

static Cycle find_first_approx(double period, const std::vector<Parameters>& init_par, double& lambda, double& mu) {
    Parameters p;
    p.period = period;
    auto iter = std::upper_bound(init_par.begin(), init_par.end(), p, comparator);
    if (iter == init_par.end()) {
        throw "approximation_error";
    }
    lambda = iter->lambda;
    mu = iter->mu;
    Pos start(-iter->amplitude, 0);
    Limit t_lim{neg_inf, neg_inf, pos_inf, pos_inf, period};
    Cycle first_approx = RK4(regular_time_step, 0, start, VDPD, t_lim, lambda, mu);
    try {
        approx_to_lim_cycle(first_approx, lambda, mu);
    } catch (...) {
        throw "approximation error";
    }
    return first_approx;
}


static Pos find_denorm_start(const Cycle& norm_cycle, const Cycle& approxible_cycle, const Parameters& p) { // находит начало цикла для сравнения
    Pos start = norm_cycle.get_start_with_zero_x(approxible_cycle.half_period_amplitude_sign());
    start.y *= p.omega / sqrt(p.epsilon);
    return start;
}

template <class Function>
static double calculate_mistake(const Pos& denorm_start, const Cycle& approxible_cycle,
                                const Parameters& p, Function difference) {
    Limit t_lim {neg_inf, neg_inf, pos_inf, pos_inf, approxible_cycle.get_time_finish()};
    Cycle denorm_cycle = RK4(approxible_cycle.get_step(), approxible_cycle.get_time_start(), denorm_start, denorm_VDPD,
                           t_lim, p.lambda, p.mu, p.epsilon, p.omega);
    return difference(denorm_cycle, approxible_cycle);
}

template <class Function>
static std::tuple<double, Pos, Parameters> best_in_changing_amplitude(const Cycle& approxible_cycle,
                                                             const Cycle& first_norm_cycle,
                                                             const Parameters& denorm_parameters,
                                                             double step, Function difference) {
    Parameters norm_parameters = denorm_parameters;
    norm_parameters.normalize();
    Parameters best_parameters = denorm_parameters;

    Limit t_lim {neg_inf, neg_inf, pos_inf, pos_inf, first_norm_cycle.period()};
    Cycle norm_cycle = first_norm_cycle;
    Pos best_start = find_denorm_start(norm_cycle, approxible_cycle, denorm_parameters);
    double min_mistake = calculate_mistake(best_start, approxible_cycle, denorm_parameters, difference);
    size_t counter = 0;
    while (norm_parameters.amplitude + step <= max_amplitude &&
           norm_parameters.amplitude + step >= min_amplitude) {
        norm_parameters.amplitude += step;
        Pos norm_start(-norm_parameters.amplitude, 0);
        norm_cycle = RK4(regular_time_step, 0, norm_start, VDPD, t_lim,
                         norm_parameters.lambda, norm_parameters.mu);
        try {
            approx_to_lim_cycle(norm_cycle, norm_parameters.lambda, norm_parameters.mu);
            norm_parameters.renormalize(denorm_parameters.omega, denorm_parameters.epsilon);
            Pos start = find_denorm_start(norm_cycle, approxible_cycle, norm_parameters);
            double mistake = calculate_mistake(start, approxible_cycle, norm_parameters, difference);
            if (mistake < min_mistake) {
                min_mistake = mistake;
                best_start = start;
                best_parameters = norm_parameters;
            }
            norm_parameters.normalize();
        } catch (...) {
            return {min_mistake, best_start, best_parameters};
        }
        ++counter;
    }
    norm_parameters.renormalize(denorm_parameters.omega, denorm_parameters.epsilon);
    return {min_mistake, best_start, best_parameters};
}

template <class Function>
static std::tuple<double, Pos, Parameters> best_with_current_normaliszation(const Cycle& approxible_cycle,
                                                                            const std::vector<Parameters>& init_par,
                                                                            const Pos& omega_epsilon,
                                                                            Function difference) {

    Parameters denorm_parameters;
    Cycle first_approx;
    Pos exept_pos;
    try {
        first_approx = find_first_approx(omega_epsilon.x * approxible_cycle.period(),
                                           init_par, denorm_parameters.lambda, denorm_parameters.mu);
    } catch (...) {
        return {pos_inf, exept_pos, denorm_parameters};
    }
    denorm_parameters.period = first_approx.period();
    denorm_parameters.amplitude = -first_approx.first_pos().x;
    denorm_parameters.renormalize(omega_epsilon.x, omega_epsilon.y);
    auto first_future_res = std::async(best_in_changing_amplitude<Function>, approxible_cycle,
                                       first_approx, denorm_parameters, 0.0005, difference);
    auto second_future_res = std::async(best_in_changing_amplitude<Function>, approxible_cycle,
                                       first_approx, denorm_parameters, -0.0005, difference);

    auto first_res = first_future_res.get();
    auto second_res = second_future_res.get();
    if (std::get<0>(first_res) < std::get<0>(second_res)) {
        return first_res;
    } else {
        return second_res;
    }
}

static std::set<Pos> get_all_moves() {
    std::set<Pos> res;
    for (int i = -1; i != 2; ++i) {
        for (int j = -1; j != 2; ++j) {
            if (i != 0 || j != 0) {
                res.emplace(i, j);
            }
        }
    }
    return res;
}

static std::set<Pos> get_not_visited(const std::set<Pos>& all_moves, const Pos& move) {
    std::set<Pos> res = all_moves;
    for (auto iter = all_moves.begin(); iter != all_moves.end(); ++iter) {
        if (all_moves.find(*iter - move) != all_moves.end()) {
            res.erase(*iter - move);
        }
    }
    return res;
}

template <class Function>
static std::tuple<double, Pos, Parameters> approx_sun_cycle(const Cycle& approxible_cycle,
                                                            const std::vector<Parameters>& init_par,
                                                            const Pos& start_omega_epsilon,
                                                            Function difference) {
    //создание начальной позиции
    Pos omega_epsilon = start_omega_epsilon; // вектор в пространстве omega-epsilon
    auto init_results = best_with_current_normaliszation(approxible_cycle, init_par, omega_epsilon, difference);
    double min_mistake = std::get<0>(init_results);
    Pos best_start = std::get<1>(init_results);
    Parameters best_parameters = std::get<2>(init_results);
    std::set<Pos> all_moves = get_all_moves();
    std::set<Pos> unchecked_moves = all_moves;;
    double step = 1.0 / 4;
    Limit borders{0, 0, 1, 1, pos_inf}; //pos_inf - показывает отсутствие времени
    while (step > 0.00005) {
        Pos best_move(0, 0);
        for (auto iter = unchecked_moves.begin(); iter != unchecked_moves.end(); ++iter) {
            if (borders.fit_lim(omega_epsilon + step * (*iter) , neg_inf)) { //neg_inf - показывает отсутствие времени
                auto res = best_with_current_normaliszation(approxible_cycle, init_par, omega_epsilon + step * (*iter), difference);
                if (std::get<0>(res) < min_mistake) {
                    min_mistake = std::get<0>(res);
                    best_move = step * (*iter);
                    best_start = std::get<1>(res);
                    best_parameters = std::get<2>(res);
                }
            }
        }
        if (best_move.is_not_zero()) {
            omega_epsilon += best_move;
            unchecked_moves = get_not_visited(all_moves, best_move);
        } else {
            unchecked_moves = all_moves;
            step /= 4;
        }
    }
    return {min_mistake, best_start, best_parameters};
}


#endif // APPROX_H_INCLUDED
