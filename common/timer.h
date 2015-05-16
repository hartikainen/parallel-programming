#ifndef TIMER_H
#define TIMER_H

#include <iomanip>
#include <iostream>
#include <sys/time.h>

class Timer {
public:
    Timer(bool add_tab_=false, std::ostream& f_=std::cout) : start{get_time()}, add_tab{add_tab_}, f{f_}
    {}

    ~Timer() {
        double now = get_time();
        std::ios_base::fmtflags oldf = f.flags(std::ios::right | std::ios::fixed);
        std::streamsize oldp = f.precision(3);
        f << (now - start);
        if (add_tab) {
            f << "\t";
        }
        f << std::flush;
        f.flags(oldf);
        f.precision(oldp);
        f.copyfmt(std::ios(NULL));
    }

private:
    static double get_time() {
        struct timeval tm;
        gettimeofday(&tm, NULL);
        return static_cast<double>(tm.tv_sec) + static_cast<double>(tm.tv_usec) / 1E6;
    }

    double start;
    bool add_tab;
    std::ostream& f;
};

#endif
