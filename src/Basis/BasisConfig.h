#ifndef BASIS_CONFIG_H
#define BASIS_CONFIG_H

#include <variant>
#include <optional>
#include <string>
#include <vector>

namespace Ambit
{
// Need to explicitly request the width of this enum to be able to use it in generic classes
enum class SplineType : int;

// Generic parent type for all config types
struct BaseBasisConfig {
    std::string valence_basis;
    std::string frozen_core;
    std::optional<std::string> MBPT_basis; // Technically defined as the MBPT basis, but only gets accessed
                               // when generating the "high" basis orbitals
    std::vector<std::string> include_valence;
    std::vector<std::string> exclude_valence;
    std::string residue;
    std::string inject_orbitals;
    std::string hf_orbitals;
};

// Variant types for each kind of BSpline. We pass this to std::variant to ensure we get exactly
// one version of the BSplines
struct BSplineBasisConfig : BaseBasisConfig {
    double RMax = 50;
    double R0 = 0.0;
    unsigned K = 7;
    unsigned N = 40;
    SplineType spline_type;
};

struct HFBasisConfig : BaseBasisConfig {
};

struct XRBasisConfig : BaseBasisConfig {
    std::string custom_orbitals;
};

using BasisConfig = std::variant<BSplineBasisConfig, HFBasisConfig, XRBasisConfig>;

} // namespace Ambit


#endif // BASISCONFIG_H
