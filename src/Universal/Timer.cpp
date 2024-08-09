#include "Include.h"
#include "Timer.h"

namespace Ambit
{
    /** Only want one instance of the walltimes per MPI rank, so use a singleton pattern here */
    Timer* Timer::Instance()
    {
        static Timer instance;
        return &instance;
    }

    Timer::Timer()
    {
        // Initialise the mapping between timer labels and walltime measurements
        for(auto label : timer_labels)
        {
            walltimes_map.insert(std::make_pair(label, 0.0f));
        }
    }

    Timer::~Timer()
    {}

    /** Return the current time stamp from chrono's "steady_clock" */
    Timer::timestamp Timer::get_timestamp()
    {
        return(std::chrono::steady_clock::now());
    }

}
