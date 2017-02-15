#ifndef TIMER_HPP
#define TIMER_HPP

class Timer {
private:
	double prevtime, curtime;

public:
	static unsigned long long rdtsc(void);

	Timer(void);
	double update_time(void);
	double get_time(void);
};

#endif
