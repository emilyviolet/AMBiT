#ifndef SPECIFICATION_CLASS_H
#define SPECIFICATION_CLASS_H

#include <string>
#include <variant>

#include "Universal/LatticeConfig.h"

namespace Ambit
{

struct GlobalSpecification {
    // Lattice parameters.
    //
    // Default values of zero => unset by user configuration and should be replaced by default
    // values that are dependent upon other settings.

    unsigned lattice_num_points = 0;
    double lattice_start_point = 0;
    double lattice_end_point = 0;
    bool lattice_exponential = false;
    double lattice_H = 0.05;

    LatticeConfig getLatticeConfig() const;
};

// On success, return empty string.
// On failure, return (long) error message.

std::string importSpecificationFile(GlobalSpecification&, const std::string& fileName);
std::string importSpecificationKV(GlobalSpecification&, const std::string& assignment);

// Perform global validation of specifications. Return non-empty error message on failure.
std::string validateSpecification(const GlobalSpecification&);

} // namespace Ambit

#endif //ndef SPECIFICATION_CLASS_H
