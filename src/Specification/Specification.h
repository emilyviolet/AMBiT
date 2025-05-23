#ifndef SPECIFICATION_CLASS_H
#define SPECIFICATION_CLASS_H

#include <string>
#include <vector>
#include <variant>
#include <optional>

#include "Basis/BasisConfig.h"
#include "HartreeFock/HFConfig.h"
#include "Universal/LatticeConfig.h"

namespace Ambit
{

struct GlobalSpecification {
    // Lattice parameters.
    //
    // Default values of zero => unset by user configuration and should be replaced by default
    // values that are dependent upon other settings.
    
    // Ungrouped options
    std::string ID;
    unsigned Z = 0;
    std::string level_directory;
    bool s1 = false;
    bool s2 = false;
    bool s3 = false;
    bool no_new_mbpt = false;
    bool check_sizes = false;
    bool clean_run = false;
    bool print_basis = false;
    bool ci_complete = false;
    bool no_ci = false;
    bool configuration_average = false;
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
    // NOTE: Moving AlphaSquared into HF specification since that's where it's used
    double hf_alpha_squared_variation = 0;
    // Also moving nuclear parameters
    double hf_nuclear_inverse_mass = 0;
    double hf_nuclear_radius = 0;
    double hf_nuclear_thickness = 0;
    // HF/QED
    bool hf_qed_uehling = false;
    bool hf_qed_self_energy = false;
    bool hf_qed_use_nuclear_density = false;
    double hf_qed_nuclear_rms_radius = 0;
    bool hf_qed_no_magnetic = false;
    bool hf_qed_no_electric = false;
    bool hf_qed_skip_offmass = false;
    bool hf_qed_use_electron_screening = false;
    // HF/NuclearPolarisability
    bool hf_nuclear_polarisability_alpha_e = false;
    bool hf_nuclear_polarisability_ebar_mev = false;
    // HF/Yukawa
    double hf_yukawa_mass = 1.0;
    double hf_yukawa_massev = 1.0;
    double hf_yukawa_rc = 1.0;
    double hf_yukawa_scale = 1.0;
    // HF/AddLocalPotential
    std::string hf_addlocal_filename;
    double hf_addlocal_scale = 1.0;
    // Basis
    std::string basis_valence;
    std::string basis_frozen_core;
    std::string basis_include_valence;
    std::string basis_exclude_valence;
    std::string basis_residue;
    std::string basis_inject_orbitals;
    bool basis_reorthogonalise;
    bool basis_hf;
    bool basis_bspline;
    bool basis_xr;
    std::string basis_hf_orbitals;
    std::string basis_custom_orbitals;
    // Basis/BSpline
    double basis_bspline_rmax = 0; 
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
    double ci_max_energy = 0.0;
    std::pair<double, double> ci_configuration_average_energy_range;
    unsigned ci_chunksize = 4;
    bool ci_sort_matrix_by_configuration = false;
    // CI/Output
    bool ci_output_print_hamiltonian = false;
    bool ci_output_write_hamiltonian = false;
    double ci_output_max_displayed_energy = 0.0;
    double ci_min_displayed_percent = 0.0;
    bool ci_output_display_print_inline = false;
    std::string ci_output_string_separator = " ";
    bool ci_output_print_relativistic_configurations = false;
    // CI/SmallSide
    std::string ci_smallside_leading_configurations;
    unsigned ci_smallside_electron_excitations;
    std::optional<std::vector<std::string>> ci_smallside_excitation_bounds;
    unsigned ci_smallside_hole_excitations = 0;
    bool ci_smallside_print_configurations = false;
    bool ci_smallside_print_rel_configurations = false;
    std::pair<double, double> ci_smallside_configuration_average_energy_range;
    // MBPT
    std::string mbpt_basis;
    std::string mbpt_energy_denom_orbitals;
    bool mbpt_use_valence;
    bool mbpt_no_core;
    bool mbpt_use_subtraction;
    bool mbpt_no_subtraction;
    bool mbpt_no_extra_box;
    double mbpt_energy_denom_floor = 0.01;
    double mbpt_delta = 0.0;
    std::vector<unsigned> mbpt_twobody_storage_limits;
    std::vector<double> mbpt_onebody_scaling;
    bool mbpt_brueckner = false;
    // MBPT/Brueckner
    double mbpt_brueckner_startpoint = 4.35e-5;
    double mbpt_brueckner_endpoint = 8.0;
    unsigned mbpt_brueckner_stride = 4;
    std::vector<double> mbpt_brueckner_scaling;
    std::vector<double> mbpt_brueckner_energy_scaling;
    bool mbpt_brueckner_use_lower = false;
    bool mbpt_brueckner_use_lower_lower = false;


    LatticeConfig getLatticeConfig() const;
    BasisConfig getBasisConfig() const;
    HFConfig getHFConfig() const;
};

// On success, return empty string.
// On failure, return (long) error message.

std::string importSpecificationFile(GlobalSpecification&, const std::string& fileName);
std::string importSpecificationKV(GlobalSpecification&, const std::string& assignment);

// Perform global validation of specifications. Return non-empty error message on failure.
std::string validateSpecification(const GlobalSpecification&);

} // namespace Ambit

#endif //ndef SPECIFICATION_CLASS_H
