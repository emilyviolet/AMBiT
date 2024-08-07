#ifndef TIMER_H
#define TIMER_H
#include <array>
#include <chrono>
#include <map>
#include <string>

namespace Ambit 
{
/** Singleton class to keep track of timing data across the calculation */
class Timer
{
public:
    static Timer* Instance();
    /** Mapping from region labels to walltime measurements. Walltimes must be stored as std::chrono
     * durations with "float" representation to make sure we can measure sub-second intervals.*/
    std::map<std::string, std::chrono::duration<float> > walltimes_map;

protected:
    /** Compile-time array of labels for the different regions we're timing */
    std::array<std::string, 8> timer_labels {"HartreeFock", "OneBodyMBPT", "TwoBodyMBPT",
                                             "Brueckner", "SlaterIntegrals", "GenerateMatrix",
                                             "SolveMatrix", "Transitions"};
protected:
    Timer();

public:
    ~Timer();
};
}

#endif
