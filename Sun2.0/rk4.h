#ifndef RK4_H_INCLUDED
#define RK4_H_INCLUDED

#include "cycle.h"

template <class Function, class... Args>
Cycle RK4(const double step, const double time_start, const Pos& start,
           Function&& equation, const Limit& lim, Args&&... args) {
    double time = time_start;
    Pos cur_pos = start;
    Cycle res;
    while (lim.fit_lim(cur_pos, time)) {
        res.add_dot(cur_pos);
        Pos k1 = equation(time, cur_pos, args...);
        Pos inc = step * k1;
        if (lim.fit_lim(cur_pos + (step / 2) * k1, time + step / 2)) {
            Pos k2 = equation(time + step / 2, cur_pos + (step / 2) * k1, args...);
            if (lim.fit_lim(cur_pos + (step / 2) * k2, time + step / 2)) {
                Pos k3 = equation(time + step / 2, cur_pos + (step / 2) * k2, args...);
                if (lim.fit_lim(cur_pos + step * k3, time + step)) {
                    Pos k4 = equation(time + step, cur_pos + step * k3, args...);
                    inc = (step / 6) * (k1 + 2 * k2 + 2 * k3 + k4);
                }
            }
        }
        cur_pos += inc;
        time += step;
    }
    Pos last_pos = create_last_pos(res.last_pos(), cur_pos, time, step, lim);
    res.add_dot(last_pos);
    res.set_time_parameters(time_start, time, step);
    return res;
}

#endif // RK4_H_INCLUDED
