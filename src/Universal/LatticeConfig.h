#ifndef LATTICECONFIG_H
#define LATTICECONFIG_H

#include <variant>

namespace Ambit
{

struct LatticeHybridConfig {
    unsigned num_points = 1000;
    double start_point = 1.0e-6;
    double end_point = 50;
};

struct LatticeExpConfig {
    unsigned num_points = 300;
    double start_point = 1.0e-5;
    double H = 0.05;
};

using LatticeConfig = std::variant<LatticeHybridConfig, LatticeExpConfig>;

} // namespace Ambit


#endif // ndef LATTICECONFIG_H

