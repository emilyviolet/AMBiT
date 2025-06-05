#ifndef HF_CONFIG_H
#define HF_CONFIG_H

#include <string>
#include <optional>

namespace Ambit
{
// Configs for decorators
struct QEDConfig {
    bool uehling;
    bool self_energy;
    bool use_nuclear_density;
    std::optional<double> nuclear_rms_radius;
    bool no_magnetic;
    bool no_electric;
    bool skip_offmass;
    bool use_electron_screening;
};

struct YukawaConfig {
    // Note: we only ever use `mass` when calculating stuff, so we need to eventually convert to
    // hartrees if the user has supplied massEV or rc
    std::optional<double> mass;
    std::optional<double> massEV;
    std::optional<double> rc;
    double scale;
};

struct NuclearPolarisabilityConfig {
   double alphaE;
   double ebarMeV;
};

struct LocalPotentialConfig {
    std::string filename;
    double scale;
};

class HFConfig {
public:
    unsigned Z;
    std::optional<unsigned> N;
    std::optional<int> charge;
    std::string configuration;
    bool breit;
    bool sms;
    bool nms;
    bool only_rel_nms;
    bool nonrel_mass_shift;
    bool include_lower_mass;
    bool local_exchange;
    double xalpha;
    double alpha_squared_variation = 0.0;
    double nuclear_thickness;
    double nuclear_radius;
    double nuclear_inverse_mass;

    std::optional<QEDConfig> qed_config;
    std::optional<YukawaConfig> yukawa_config;
    std::optional<LocalPotentialConfig> local_potential_config;
    std::optional<NuclearPolarisabilityConfig> nuclear_polarisability_config;
};


} // namespace Ambit


#endif // HF_CONFIG_H
