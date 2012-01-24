#ifndef PROFILER_HPP_
#define PROFILER_HPP_

#ifdef PROFILER_ON
#include "cycle.h"
#include <iostream>
#include <map>

class ProfileManager {
 public:
  struct ProfilingResult {
    ProfilingResult() : calls(0), time(0) {}
    
    void clear() {
      calls = 0; time = 0;
    }
    size_t calls;
    double time;		
  };

  typedef std::map<const char*, ProfilingResult> Results;

  static ProfileManager* getInstance() {
    static ProfileManager prof = ProfileManager();
    return &prof;
  }

  static void printProfilerData() {
    ProfileManager& prof = *getInstance();
    const Results& results = prof.getResults();
    std::clog << "name,calls,time" << std::endl;
    for(Results::const_iterator it = results.begin(); it != results.end(); ++it) {
      std::clog << it->first << "," << it->second.calls << ","
                << it->second.time << std::endl;
    }
  }
  
 
  void profilingResult(const char* name, double time) {
    ProfilingResult& p = m_results[name];
    p.calls++;
    p.time += time;
  }
  
  void clear() {
    m_results.clear();
  }
  
  const Results& getResults() const {
    return m_results;
  }
  
 private:
  Results m_results;
};

class AutoProfile {
 public:
  AutoProfile(const char* name) : m_name(name), m_start(getticks()) {}

  ~AutoProfile() {
    ticks end = getticks();
    double t = elapsed(end, m_start);
    ProfileManager::getInstance()->profilingResult(m_name, t);
  }

 private:
  const char* m_name;
  ticks m_start;
};


#define PROFILE(name) AutoProfile __scoped_profiler__(name)
#define PRINT_PROFILE_DATA ProfileManager::printProfilerData();
#else
#define PROFILE(name)
#define PRINT_PROFILE_DATA
#endif

#endif
