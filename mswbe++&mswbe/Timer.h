/*
 * Timer.h
 *
 *  Created on: 5Dec.,2017
 *      Author: Lijun Chang
 *      Email: ljchang@outlook.com
 */

#ifndef TIMER_H_
#define TIMER_H_

#include <cstdlib>
#include <sys/time.h>

class Timer {
public:
    Timer() { m_start = timestamp(); }
    void restart() { m_start = timestamp(); }
    long long elapsed() { return (timestamp() - m_start)/1000000; }
    long long elapsed_in_millisec() { return (timestamp() - m_start)/1000; }
    
private:
    long long m_start;
    
    // Returns a timestamp ('now') in microseconds
    long long timestamp() {
        struct timeval tp;
        gettimeofday(&tp, nullptr);
        return ((long long)(tp.tv_sec))*1000000 + tp.tv_usec;
    }
};

#endif /* TIMER_H_ */

