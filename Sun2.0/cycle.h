#ifndef CYCLE_H_INCLUDED
#define CYCLE_H_INCLUDED

#include <vector>
#include <iostream>
#include <fstream>
#include <limits>
#include <initializer_list>
#include <iomanip>
#include <cmath>

const double pos_inf = std::numeric_limits<double>::infinity();
const double neg_inf = -std::numeric_limits<double>::infinity();

class Pos {
private:
    double equality_dist = 0.0000001;

public:
    double x;
    double y;
    Pos() : x(0), y(0) {}
    Pos(double _x, double _y) : x(_x), y(_y) {}
    Pos(const Pos& pos) : x(pos.x), y(pos.y) {}
    friend std::ostream& operator<< (std::ostream &out, const Pos &pos) {
        out << pos.x << " " << pos.y;
        return out;
    }
    Pos& operator= (const Pos& other);
    friend Pos operator+ (const Pos& a, const Pos& b);
    friend Pos operator- (const Pos& a, const Pos& b);
    friend Pos operator* (double s, const Pos& a);
    friend Pos operator/ (const Pos& a, double s);
    Pos operator-();
    Pos& operator/= (const double s);
    Pos& operator+= (const Pos &a);
    Pos& operator-= (const Pos &a);
    Pos& operator*= (const double s);
    friend bool operator< (const Pos& a, const Pos& b);
    friend bool operator== (const Pos& a, const Pos& b);
    friend bool operator!= (const Pos& a, const Pos& b);
    friend Pos best_kramer(Pos a, Pos b, Pos c);
    double length() {
        return sqrt(x * x + y * y);
    }
    bool is_not_zero() {
        return this->length() > equality_dist;
    }
    Pos(std::initializer_list<double> l) : x(*l.begin()), y(*(l.begin() + 1)) {}
    ~Pos() {}
};

Pos operator* (double s, const Pos &a) {
    return {a.x * s, a.y * s};
}

Pos operator/ (const Pos& a, double s) {
    return {a.x / s, a.y / s};
}

Pos operator+ (const Pos& a, const Pos& b) {
    return {a.x + b.x, a.y + b.y};
}

Pos operator- (const Pos& a, const Pos& b) {
    return {a.x - b.x, a.y - b.y};
}

Pos& Pos::operator= (const Pos& other) {
    this->x = other.x;
    this->y = other.y;
    return *this;
}

Pos& Pos::operator+= (const Pos &a) {
    this->x += a.x;
    this->y += a.y;
    return *this;
}

Pos& Pos::operator-= (const Pos &a) {
    this->x -= a.x;
    this->y -= a.y;
    return *this;
}

Pos& Pos::operator*= (const double s) {
    this->x *= s;
    this->y *= s;
    return *this;
}

Pos& Pos::operator/= (const double s) {
    this->x /= s;
    this->y /= s;
    return *this;
}

Pos Pos::operator-() {
    return Pos(-(this->x), -(this->y));
}

static double distance(const Pos& a, const Pos& b) {
    return (a - b).length();
}

bool operator== (const Pos& a, const Pos& b) {
    return distance(a, b) < a.equality_dist;
}

bool operator!= (const Pos& a, const Pos& b) {
    return !(a == b);
}

bool operator< (const Pos& a, const Pos& b) {
    return (a != b) && ((a.x < b.x) || (!(a.x > b.x) && a.y < b.y));
}

static Pos kramer_2d(const Pos& a, const Pos& b, const Pos& c) {
    double det = a.x * b.y - a.y * b.x;
    double det_1 = c.x * b.y - b.x * c.y;
    double det_2 = a.x * c.y - a.y * c.x;
    if (std::isnan(det_1 / det) || std::isnan(det_2 / det)) {
        std::cout << "liniar related\n";
        throw "error";
    }
    return Pos(det_1 / det, det_2 / det);
}

Pos best_kramer(Pos a, Pos b, Pos c) {
    Pos start = kramer_2d(a, b, c);
    double min_len = start.length();
    double step = 2;
    std::vector<Pos> moves;
    for (int i = -1; i != 1; ++i) {
        for (int j = -1; j != 1; ++j) {
            if (i != 0 || j != 0) {
                moves.emplace_back(i, j);
            }
        }
    }
    while (step <= 8) {
        std::vector<Pos> best_moves(3);
        for (Pos move_1 : moves) {
            for (Pos move_2 : moves) {
                for (Pos move_3 : moves) {
                    move_1 *= a.equality_dist / step;
                    move_2 *= a.equality_dist / step;
                    move_3 *= a.equality_dist / step;
                    try {
                        double new_len = kramer_2d(a + move_1,
                                                   b + move_2,
                                                   c + move_3).length();
                        if (new_len < min_len) {
                            min_len = new_len;
                            best_moves[0] = move_1;
                            best_moves[1] = move_2;
                            best_moves[2] = move_3;
                        }
                    } catch(...) {}
                }
            }
        }
        step *= 2;
        a += best_moves[0];
        b += best_moves[1];
        c += best_moves[2];
    }
    return kramer_2d(a, b, c);
}

class Limit {
private:
    Pos dot1;
    Pos dot2;
    double time_limit;
public:
    Limit() : dot1 (neg_inf, neg_inf), dot2(pos_inf, pos_inf), time_limit(pos_inf) {}
    Limit(std::initializer_list<double> l) :
        dot1(*l.begin(), *(l.begin() + 1)),
        dot2(*(l.begin() + 2), *(l.begin() + 3)),
        time_limit(*(l.begin() + 4)) {}
    Limit(const Pos& pos1, const Pos& pos2, double _time_limit) :
    dot1(pos1), dot2(pos2), time_limit(_time_limit) {}
    Limit(const Limit& l) : dot1(l.dot1), dot2(l.dot2), time_limit(l.time_limit) {}
    bool fit_lim(const Pos& pos, double time) const {
        return (pos.x >= dot1.x && pos.x <= dot2.x &&
                pos.y >= dot1.y && pos.y <= dot2.y &&
                time <= time_limit);
    }
    friend Pos create_last_pos(const Pos& fitting_lim, const Pos& not_fitting_lim,
                     double& time, double step, const Limit& lim);
    ~Limit() {}
};

Pos create_last_pos(const Pos& fitting_lim, const Pos& not_fitting_lim,
                     double& time, double step, const Limit& lim) {
    double min_koef = pos_inf;
    time -= step;
    double koef = (lim.dot1.x - fitting_lim.x) / (not_fitting_lim.x - fitting_lim.x);
    if (koef >= 0 && koef < min_koef) min_koef = koef;
    koef = (lim.dot1.y - fitting_lim.y) / (not_fitting_lim.y - fitting_lim.y);
    if (koef >= 0 && koef < min_koef) min_koef = koef;
    koef = (lim.dot2.x - fitting_lim.x) / (not_fitting_lim.x - fitting_lim.x);
    if (koef >= 0 && koef < min_koef) min_koef = koef;
    koef = (lim.dot2.y - fitting_lim.y) / (not_fitting_lim.y - fitting_lim.y);
    if (koef >= 0 && koef < min_koef) min_koef = koef;
    koef = (lim.time_limit - time) / step;
    if (koef >= 0 && koef < min_koef) min_koef = koef;
    time += min_koef * step;
    return fitting_lim + min_koef * (not_fitting_lim - fitting_lim);
}

class Cycle {
private:
    std::vector<Pos> dots;
    double time_start;
    double time_finish;
    double step;
    bool time_paramrters_are_set;

public:

    Cycle() : dots(), time_start(0), time_finish(0),
            step(0), time_paramrters_are_set(false) {}
    Cycle(const std::vector<Pos>& _dots) : dots(_dots), time_start(0), time_finish(0),
                                            step(0), time_paramrters_are_set(false) {}
    Cycle(const Cycle& c) : dots(c.dots), time_start(c.time_start), time_finish(c.time_finish),
                            step(c.step), time_paramrters_are_set(c.time_paramrters_are_set) {}
    Cycle& operator=(const Cycle& c) = default;
    Cycle(std::ifstream &in) : time_paramrters_are_set(true) {
        in >> time_start >> time_finish;
        size_t cycle_size;
        in >> cycle_size;
        dots.reserve(cycle_size);
        for (size_t i = 0; i != cycle_size; ++i) {
            double x, y;
            in >> x >> y;
            dots.emplace_back(x, y);
        }
        step = (time_finish - time_start) / static_cast<double>(cycle_size - 1);
    }

    void print(std::ofstream &out) const {
        for (size_t i = 0; i != dots.size(); ++i) {
            out << std::fixed << std::setprecision(12) << dots[i] << "\n";
        }
    }

    void add_dot(const Pos& pos) {
        dots.push_back(pos);
    }

    Pos first_pos() const {
        return dots[0];
    }

    Pos last_pos() const {
       return dots.back();
    }

    double get_time_start() const {
        return time_start;
    }

    double get_time_finish() const {
        return time_finish;
    }

    double get_step() const {
        return step;
    }

    double period() const {
        return (time_finish - time_start);
    }

    void set_time_parameters(double start, double finish, double _step) {
        if (!time_paramrters_are_set) {
            time_start = start;
            time_finish = finish;
            step = _step;
        } else {
            std::cout << "parameters are already set\n";
            throw "error";
        }
    }

    Pos operator[] (const double time) const {
        if (time > time_finish + 0.0000001) {
            std::cout << "wrong time\n";
            throw "error";
        }

        double unround_index = std::floor((time - time_start) / step);
        double rest = (time - time_start) - unround_index * step;
        size_t index = static_cast<size_t>(unround_index);
        Pos diff = dots[index + 1] - dots[index];
        return dots[index] + rest * diff;
    }

    Pos get_start_with_zero_x(bool with_pos_amplitude) const {
        for (size_t i = 0; i != dots.size() - 1; ++i) {
            if (dots[i].x <= 0 && dots[i + 1].x > 0) {
                Pos dir = dots[i + 1] - dots[i];
                double koef = (0 - dots[i].x) / dir.x;
                Pos res = dots[i] + koef * dir;
                if ((with_pos_amplitude && res.y < 0) ||
                    (!with_pos_amplitude && res.y > 0)) {
                    res = -res;
                }
                return res;
            }
        }
        return Pos(0, 0);
    }

    bool half_period_amplitude_sign() const {
        return (*this)[(time_start + time_finish) / 2].x > 0;
    }

    ~Cycle() {}
};


#endif // CYCLE_H_INCLUDED
