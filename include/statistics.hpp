#ifndef NP_STATISTICS_HPP
#define NP_STATISTICS_HPP
#include <iostream>
#include <ostream>
#include <string>
#include <chrono>

//#define STS(A) A
#define STS(A)

namespace NP {
  struct StatCollect {
    int count[50];
    int total;
    int div;
    int max;
    int rettrue;
    int retfalse;
    std::chrono::time_point<std::chrono::high_resolution_clock> startTime;
    std::chrono::time_point<std::chrono::high_resolution_clock> endTime;
    double sumTime;
    int timeCount;
    int printat;
    bool printit;
    const char *label;

    StatCollect(const char *lab, int printevery=1000)
      : total(0), div(1), max(0),rettrue(0),retfalse(0),printit(false),label(lab),printat(printevery)
    {
      for (int i=0; i<50; i++) count[i]=0;
      startTime = std::chrono::high_resolution_clock::now();
      endTime = startTime;
      sumTime = 0.0;
      timeCount=0;
    }

    void tick(int value)
    {
      if (value>max) {
	printit=true;
	max = value;
      }
      total++;
      if (total%printat==0) printit=true;
      if (value<50) count[value]++;
    }

    void start()
    {
      startTime = std::chrono::high_resolution_clock::now();
    }

    void end()
    {
      endTime = std::chrono::high_resolution_clock::now();
      auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(endTime-startTime);
      sumTime += duration.count();
      timeCount++;
    }

    void tick(bool retvalue)
    {
      if (retvalue)  rettrue++;
      else retfalse++;
    }

    void print()
    {
      if (printit) {
	std::cerr << label << " (" << total << ", " << max << ") : "<< count[0];
	for (int i=1; i<50 && i<=max; i++)
	  std::cerr << ", " << count[i];
	std::cerr << "\n";
	if (rettrue > 0 || retfalse>0)
	  std::cerr << label << " returns true: " << rettrue << ", false: "<< retfalse << "\n";
	if (timeCount > 0)
	  std::cerr << label << " avg. time: " << sumTime/(1.0 * timeCount) << "ns (" << timeCount << " samples)\n"; 
      }
      printit=false;
    }
  };
}
#endif
