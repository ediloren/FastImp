// Template Numerical Toolkit (TNT) for Linear Algebra 
//
// R. Pozo
// Applied and Computational Mathematics Division
// National Institute of Standards and Technology
//
// Modefied by Ben Song and Zhenhai Zhu

#ifndef STPWATCH_H
#define STPWATCH_H

// for clock() and CLOCKS_PER_SEC
#include <unistd.h>
#include <ctime>
#include "sys/times.h"

namespace TNT {

  class stopwatch {

  public:
    // void    start()     : start timing
    // double  stop()      : stop timing
    // void    reset()     : set elapsed time to 0.0
    // double  read()      : read elapsed time (in seconds)

    stopwatch() : running(false), last_utime(0), elapsed(0.0) {
      ticks_per_sec = sysconf(_SC_CLK_TCK);
      usec_per_tick = 1000000 / ticks_per_sec;
      assert(usec_per_tick > 0);
      assert(usec_per_tick < 100000);
    }

    void start() {
      if (! running) { last_utime = second();  running = true;  }
    }

    double stop()  { 
      if (running)  {
	clock_t current_utime = second();
	elapsed += (current_utime - last_utime) / ticks_per_sec
	  +((current_utime - last_utime) % ticks_per_sec) * usec_per_tick*1e-6;
	running = false;
      }
      return elapsed; 
    }

    void reset() {  running = false;  last_utime = 0; elapsed =0.0; }

    double read()   {  
      if (running) {
	clock_t current_utime = second();
	elapsed += (current_utime - last_utime) / ticks_per_sec
	  +((current_utime - last_utime) % ticks_per_sec) * usec_per_tick*1e-6;
	last_utime = second();
      }
      return elapsed;
    }       

  private:
    bool running;
    clock_t last_utime;
    clock_t ticks_per_sec;
    clock_t usec_per_tick;
    double elapsed;

    clock_t second(void) {
      struct tms buf;
      times(&buf);
      return buf.tms_utime;
    }
  };
 
} // namespace TNT


#endif
    

            
