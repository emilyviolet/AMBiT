#ifndef SPECIFICATION_CLASS_H
#define SPECIFICATION_CLASS_H

#include <string>
#include <vector>
#include <optional>
#include <memory>
#include "SpecificationMap.h"

namespace Ambit
{


struct GlobalSpecification {
    friend class parapara::specification<GlobalSpecification>;
    // Lattice parameters.
    //
    // Default values of zero => unset by user configuration and should be replaced by default
    // values that are dependent upon other settings.
    
    // Ungrouped options
    std::string ID;
    unsigned Z = 0;
    std::optional<std::string> level_directory;
    bool no_new_mbpt = false;
    bool check_sizes = false;
    bool clean_run = false;
    bool print_basis = false;
    bool ci_complete = false;
    bool no_ci = false;
    bool configuration_average = false;
    double alpha_squared_variation = 0;
    double nuclear_inverse_mass = 0;
    double nuclear_radius = 0;
    double nuclear_thickness = 0;
    bool m = false;
    // Lattice
    unsigned lattice_num_points = 0;
    double lattice_start_point = 0;
    double lattice_end_point = 0;
    bool lattice_exponential = false;
    double lattice_H = 0.05;
    // HF
    std::optional<unsigned> hf_N;
    std::optional<int> hf_charge;
    std::string hf_configuration;
    bool hf_breit = false;
    bool hf_sms = false;
    bool hf_nms = false;
    bool hf_only_rel_nms = false;
    bool hf_nonrel_mass_shift = false;
    bool hf_include_lower_mass = false;
    bool hf_local_exchange = false;
    double hf_xalpha = 1.0;
    // HF/QED
    bool hf_do_qed = false; // Checks if the [HF/QED] section is defined
    bool hf_qed_uehling = false;
    bool hf_qed_self_energy = false;
    bool hf_qed_use_nuclear_density = false;
    std::optional<double> hf_qed_nuclear_rms_radius;
    bool hf_qed_no_magnetic = false;
    bool hf_qed_no_electric = false;
    bool hf_qed_skip_offmass = false;
    bool hf_qed_use_electron_screening = false;
    // HF/NuclearPolarisability
    bool hf_do_nuclear_polarisability = false;
    double hf_nuclear_polarisability_alpha_e;
    double hf_nuclear_polarisability_ebar_mev;
    // HF/Yukawa
    bool hf_do_yukawa = false;
    std::optional<double> hf_yukawa_mass;
    std::optional<double> hf_yukawa_massev;
    std::optional<double> hf_yukawa_rc;
    double hf_yukawa_scale = 1.0;
    // HF/AddLocalPotential
    bool hf_do_local_potential = false;
    std::string hf_addlocal_filename;
    double hf_addlocal_scale = 1.0;
    // Basis
    std::string basis_valence;
    std::optional<std::string> basis_frozen_core;
    std::optional<std::string> basis_size;
    // TODO EVK: Should this be a std::optional?
    std::vector<std::string> basis_include_valence;
    std::vector<std::string> basis_exclude_valence;
    std::optional<std::string> basis_residue;
    std::optional<std::vector<std::string>> basis_inject_orbitals;
    bool basis_reorthogonalise = false;
    bool basis_hf = false;
    bool basis_bspline = true;
    bool basis_xr = false;
    std::optional<std::string> basis_hf_orbitals;
    std::optional<std::vector<std::string> > basis_xr_custom_orbitals;
    // Basis/BSpline
    std::optional<double> basis_bspline_rmax; 
    double basis_bspline_r0 = 0;
    double basis_bspline_k = 7;
    double basis_bspline_N = 40;
    std::string basis_bspline_splinetype;
    // CI
    std::string ci_leading_configurations;
    std::string ci_leading_rel_configurations;
    std::string ci_extra_configurations;
    std::string ci_extra_rel_configurations;
    // TODO EVK: This is really annoying, because the spec says that it can be *either* an integer
    // (e.g. CI/ElectronExcitations=2) or a string (e.g. ElectronExcitations = '1,5spdf,2,5spd')
    // and this is really annoying to deal with
    unsigned ci_electron_excitations;
    std::optional<std::vector<std::string>> ci_excitation_bounds;
    unsigned ci_hole_excitations = 0;
    std::vector<unsigned> ci_even_parity_twoj;
    std::vector<unsigned> ci_odd_parity_twoj;
    unsigned ci_num_solutions = 6;
    bool ci_all_symmetries = false;
    bool ci_gfactors = false;
    bool ci_no_gfactors = false;
    bool ci_memory_saver = false;
    bool ci_single_configuration_ci = false;
    bool ci_print_configurations = false;
    bool ci_print_rel_configurations = false;
    bool ci_scalapack = false;
    std::optional<double> ci_max_energy;
    std::optional<std::pair<double, double> > ci_configuration_average_energy_range;
    unsigned ci_chunksize = 4;
    bool ci_sort_matrix_by_configuration = false;
    // CI/Output
    bool ci_output_print_hamiltonian = false;
    bool ci_output_write_hamiltonian = false;
    std::optional<double> ci_output_max_displayed_energy;
    std::optional<double> ci_output_min_displayed_percent;
    bool ci_output_print_inline = false;
    std::optional<std::string> ci_output_separator;
    bool ci_output_print_relativistic_configurations = false;
    bool ci_output_no_configs = false;
    // CI/SmallSide
    std::string ci_smallside_leading_configurations;
    unsigned ci_smallside_electron_excitations;
    std::optional<std::vector<std::string>> ci_smallside_excitation_bounds;
    unsigned ci_smallside_hole_excitations = 0;
    bool ci_smallside_print_configurations = false;
    bool ci_smallside_print_rel_configurations = false;
    std::optional<std::pair<double, double> > ci_smallside_configuration_average_energy_range;
    // MBPT
    std::optional<std::string> mbpt_basis;
    std::optional<std::string> mbpt_energy_denom_orbitals;
    bool mbpt_use_valence;
    bool mbpt_no_core;
    bool mbpt_use_subtraction;
    bool mbpt_no_subtraction;
    bool mbpt_no_extra_box;
    double mbpt_energy_denom_floor = 0.01;
    double mbpt_delta = 0.0;
    std::optional<std::vector<unsigned> > mbpt_twobody_storage_limits;
    std::optional<std::vector<double> > mbpt_onebody_scaling;
    bool mbpt_brueckner = false;
    // NOTE: Always use these when checking which MBPT diagrams to include
    bool mbpt_one_body = false;
    bool mbpt_two_body = false;
    bool mbpt_three_body = false;
    // MBPT/Brueckner
    double mbpt_brueckner_startpoint = 4.35e-5;
    double mbpt_brueckner_endpoint = 8.0;
    unsigned mbpt_brueckner_stride = 4;
    std::vector<double> mbpt_brueckner_scaling;
    std::vector<double> mbpt_brueckner_energy_scaling;
    bool mbpt_brueckner_use_lower = false;
    bool mbpt_brueckner_use_lower_lower = false;
    bool mbpt_brueckner_excited = false;
    // We need to explicitly enumerate the options for MBPT diagrams: s1, s2, s3, s12, s123, s13
    // and s23. It's kind of clunky to use these directly in the main body of the code, so we
    // obfuscate them here and then, during the spec validation, set the more general attributes:
    // mbpt_one_body, mbpt_two_body, mbpt_three_body as appropriate
    // NOTE: Don't use these directly in AMBiT, use GlobalSpecification::mbpt_one_body and friends
    // TODO EVK: Really want to be able to make these private, but GlobalSpecification is declared
    // globally, so it doesn't work. Probably need to move the parsing stuff into a separate class
    // so I can make it a friend
    bool _s1 = false;
    bool _s2 = false;
    bool _s3 = false;
    bool _s123 = false;
    bool _s12 = false;
    bool _s13 = false;
    bool _s23 = false;
};

// On success, return empty string.
// On failure, return (long) error message.

std::string importSpecificationFile(GlobalSpecification&, const std::string& fileName);
std::string importSpecificationKV(GlobalSpecification&, const std::string& assignment);

// Perform global validation of specifications. Return non-empty error message on failure.
std::string validateAndNormaliseSpecification(GlobalSpecification&);

// Return a map-like interface to the configuration, which lets us do nice things like 
// config_map["key"] = "value"
typedef keyed_record_view<GlobalSpecification> SpecificationMap;
typedef std::shared_ptr<SpecificationMap> pSpecificationMap;
SpecificationMap get_config_map_view(GlobalSpecification& gs);

} // namespace Ambit

#endif //ndef SPECIFICATION_CLASS_H
