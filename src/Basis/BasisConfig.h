#ifndef BASIS_CONFIG_H
#define BASIS_CONFIG_H

#include <memory>
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
    std::optional<std::string> frozen_core;
    std::optional<std::string> basis_size;
    std::optional<std::string> mbpt_basis; // Technically defined as the MBPT basis, but only gets accessed
                               // when generating the "high" basis orbitals
    std::vector<std::string> include_valence;
    std::vector<std::string> exclude_valence;
    std::optional<std::string> residue;
    std::optional<std::vector<std::string>> inject_orbitals;
    std::optional<std::string> hf_orbitals;
    bool reorthogonalise;
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
    std::optional<std::vector<std::string> > custom_orbitals;
};

typedef std::variant<BSplineBasisConfig, HFBasisConfig, XRBasisConfig> BasisConfig ;
typedef std::unique_ptr<BSplineBasisConfig> pBSplineBasisConfig;
} // namespace Ambit


#endif // BASISCONFIG_H
