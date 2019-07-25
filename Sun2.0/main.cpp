#include <iostream>
#include <string>
#include "equations.h"
#include "rk4.h"
#include "approx.h"

static std::vector<std::vector<Pos>> create_init_normalization() {
    std::vector<std::vector<Pos>> res;
    std::vector<Pos> to_add_1 {Pos(0.244, 0.0431), Pos(0.2755, 0.0319), Pos(0.2872, 0.0355), Pos(0.2761, 0.0237), Pos(0.2885, 0.0206),
                               Pos(0.2771, 0.0164), Pos(0.264, 0.0279), Pos(0.2888, 0.0193), Pos(0.3056, 0.0208), Pos(0.2404, 0.0262)};
    std::vector<Pos> to_add_2 {Pos(0.2458, 0.0436), Pos(0.2655, 0.0316), Pos(0.2866, 0.0349), Pos(0.2771, 0.0235), Pos(0.2893, 0.0204),
                               Pos(0.2911, 0.0175), Pos(0.269, 0.2605), Pos(0.2891, 0.0193), Pos(0.2881, 0.0209), Pos(0.2414, 0.0276)};
    std::vector<Pos> to_add_3 {Pos(0.2463, 0.0426), Pos(0.2038, 0.0309), Pos(0.2876, 0.0352), Pos(0.2735, 0.0232), Pos(0.2864, 0.0205),
                               Pos(0.2765, 0.0169), Pos(0.2448, 0.0262), Pos(0.2212, 0.018), Pos(0.2826, 0.0206), Pos(0.2459, 0.0261)};
    std::vector<Pos> to_add_4 {Pos(0.2452, 0.0428), Pos(0.2027, 0.0317), Pos(0.2876, 0.0349), Pos(0.2725, 0.0236), Pos(0.2876, 0.0207),
                               Pos(0.2755, 0.0175), Pos(0.2477, 0.0263), Pos(0.2925, 0.0192), Pos(0.2826, 0.0212), Pos(0.2449, 0.0274)};
    //получено из предыдущей работы
    res.push_back(to_add_1);
    res.push_back(to_add_2);
    res.push_back(to_add_3);
    res.push_back(to_add_4);

    return res;
}

static std::string approx_dir(size_t type) {
    std::string res;
    if (type == 0) {
        res = "/home/viacheslav/progs/Sun2.0/distance_no_prob_approx";
    }
    if (type == 1) {
        res = "/home/viacheslav/progs/Sun2.0/integral_no_prob_approx";
    }
    if (type == 2) {
        res = "/home/viacheslav/progs/Sun2.0/distance_prob_approx";
    }
    if (type == 3) {
        res = "/home/viacheslav/progs/Sun2.0/integral_prob_approx";
    }
    if (type > 3) throw "wrong type";
    return res;
}

static void execute_approx(size_t type, const std::vector<std::vector<Pos>>& init_normalization,
                    const std::vector<Parameters>& init_parameters) {
    std::string dir = approx_dir(type);
    std::ifstream data;
    data.open(dir + "/approxible_data.txt");
    std::vector<Cycle> portraits;
    portraits.reserve(10);
    for (size_t i = 0; i != 10; ++i) {
        portraits.emplace_back(data);
    }
    data.close();

    std::vector<std::tuple<double, Pos, Parameters>> results(10);
    std::vector<std::future<std::tuple<double, Pos, Parameters>>> future_results(10);
    for (size_t i = 0; i != 10; ++i) {
        if (type % 2) {
            future_results[i] = async(approx_sun_cycle<decltype(integ_difference)>, portraits[i],
                                  init_parameters, init_normalization[type][i], integ_difference);
        } else {
            future_results[i] = async(approx_sun_cycle<decltype(dist_difference)>, portraits[i],
                                  init_parameters, init_normalization[type][i], dist_difference);
        }
    }

    for (size_t i = 0; i != 10; ++i) {
        results[i] = future_results[i].get();
    }

    std::ofstream param_and_mist;
    std::ofstream general_sol;
    param_and_mist.open(dir + "/p&m.txt");
    general_sol.open(dir + "/general_sol.txt");
    for (size_t i = 0; i != 10; ++i) {
        param_and_mist << std::fixed << std::setprecision(12) << std::get<0>(results[i]) << " ";
        Parameters p = std::get<2>(results[i]);
        param_and_mist << p.omega << " " << p.epsilon << " " << p.lambda << " " << p.mu << " "
                        << p.amplitude << " " << p.period << "\n";
        std::ofstream cycle_print;
        cycle_print.open(dir + "/cycle_" + std::to_string(i + 1) + ".txt");
        Limit t_lim{neg_inf, neg_inf, pos_inf, pos_inf, portraits[i].get_time_finish()};
        Cycle c = RK4(portraits[i].get_step(), portraits[i].get_time_start(), std::get<1>(results[i]),
                      denorm_VDPD, t_lim, p.lambda, p.mu, p.epsilon, p.omega);
        for (double time = portraits[i].get_time_start();
             time < portraits[i].get_time_finish() + 0.0000001;
             time += portraits[i].get_step()) {
            general_sol << std::fixed << std::setprecision(12) << time << " " << c[time].x << "\n";
        }
        c.print(cycle_print);
        cycle_print.close();
    }
    general_sol.close();
    param_and_mist.close();
}

int main() {
    // код для генерации вспомогательного файла, работает ~11.5 часов;
    /*
    std::ofstream out;
    out.open("init_positions_2.txt");
    calculate_init_positions(out);
    out.close();
    */

    //приблизительно 4h
    std::ifstream init_pos;
    init_pos.open("init_positions_2.txt");
    std::vector<Parameters> init_par = extract_init_pos(init_pos);
    std::vector<std::vector<Pos>> init_normalization = create_init_normalization();
    for (size_t i = 0; i != 4; ++i) {
        execute_approx(i, init_normalization, init_par);
    }

}
