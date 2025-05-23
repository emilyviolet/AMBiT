#include <algorithm>
#include <fstream>
#include <istream>
#include <string>
#include <string_view>

#include "Specification.h"
#include "parapara/parapara.h"
#include "Include.h"

namespace P = parapara;

namespace Ambit
{

template <typename Record>
P::hopefully<void> custom_import_ini(Record& rec, const P::specification_set<Record>& specs, std::istream& in);

inline auto positive = P::greater_than(0, "must be positive");

P::specification<GlobalSpecification> global_specifications[] = {
    // Ungrouped options
    {"ID",                      &GlobalSpecification::ID, P::nonempty()},
    {"Z",                       &GlobalSpecification::Z, positive},
    {"LevelDirectory",          &GlobalSpecification::level_directory},
    {"-s1",                     &GlobalSpecification::s1},
    {"-s2",                     &GlobalSpecification::s2},
    {"-s3",                     &GlobalSpecification::s3},
    {"--no-new-mbpt",           &GlobalSpecification::no_new_mbpt},
    {"--check-sizes",           &GlobalSpecification::check_sizes},
    {"-c",                      &GlobalSpecification::clean_run},
    {"-p",                      &GlobalSpecification::print_basis},
    {"--ci-complete",           &GlobalSpecification::ci_complete},
    {"--no-ci",                 &GlobalSpecification::no_ci},
    {"--configuration-average", &GlobalSpecification::configuration_average},
    // Lattice
    {"Lattice/NumPoints",           &GlobalSpecification::lattice_num_points,    P::nonzero()},
    {"Lattice/StartPoint",          &GlobalSpecification::lattice_start_point,   positive},
    {"Lattice/EndPoint",            &GlobalSpecification::lattice_end_point,     positive},
    {"Lattice/--exp-lattice",       &GlobalSpecification::lattice_exponential},
    {"Lattice/H",                   &GlobalSpecification::lattice_H,             positive},
    // HF
    {"HF/N",                            &GlobalSpecification::hf_N, positive},
    {"HF/Charge",                       &GlobalSpecification::hf_charge},
    {"HF/Configuration",                &GlobalSpecification::hf_configuration},
    {"HF/--breit",                      &GlobalSpecification::hf_breit},
    {"HF/--sms",                        &GlobalSpecification::hf_sms},
    {"HF/--nms",                        &GlobalSpecification::hf_nms},
    {"HF/--only-relativistic-nms",      &GlobalSpecification::hf_only_rel_nms},
    {"HF/--nonrelativistic-mass-shift", &GlobalSpecification::hf_nonrel_mass_shift},
    {"HF/--include-lower-mass",         &GlobalSpecification::hf_include_lower_mass},
    {"HF/--local-exchange",             &GlobalSpecification::hf_local_exchange},
    {"HF/Xalpha",                       &GlobalSpecification::hf_xalpha},
    {"HF/AlphaSquaredVariation",        &GlobalSpecification::hf_alpha_squared_variation},
    {"HF/NuclearInverseMass",           &GlobalSpecification::hf_nuclear_inverse_mass},
    {"HF/NuclearRadius",                &GlobalSpecification::hf_nuclear_radius},
    {"HF/NuclearThickness",             &GlobalSpecification::hf_nuclear_thickness},
    // HF/QED
    {"HF/QED/--uehling",                &GlobalSpecification::hf_qed_uehling},
    {"HF/QED/--self-energy",            &GlobalSpecification::hf_qed_self_energy},
    {"HF/QED/--use-nuclear-density",    &GlobalSpecification::hf_qed_use_nuclear_density},
    {"HF/QED/NuclearRMSRadius",         &GlobalSpecification::hf_qed_nuclear_rms_radius, P::nonzero()},
    {"HF/QED/--no-magnetic",            &GlobalSpecification::hf_qed_no_magnetic},
    {"HF/QED/--no-electric",            &GlobalSpecification::hf_qed_no_electric},
    {"HF/QED/--skip-offmass",           &GlobalSpecification::hf_qed_skip_offmass},
    {"HF/QED/--use-electron-screening", &GlobalSpecification::hf_qed_use_electron_screening},
    // HF/NuclearPolarisability
    {"HF/NuclearPolarisability/AlphaE",     &GlobalSpecification::hf_nuclear_polarisability_alpha_e},
    {"HF/NuclearPolarisability/EbarMeV",    &GlobalSpecification::hf_nuclear_polarisability_ebar_mev},
    // HF/Yukawa
    {"HF/Yukawa/Mass",      &GlobalSpecification::hf_yukawa_mass},
    {"HF/Yukawa/MassEV",    &GlobalSpecification::hf_yukawa_massev},
    {"HF/Yukawa/Rc",        &GlobalSpecification::hf_yukawa_rc},
    {"HF/Yukawa/Scale",     &GlobalSpecification::hf_yukawa_scale},
    // HF/AddLocalPotential
    {"HF/AddLocalPotential/Filename",   &GlobalSpecification::hf_addlocal_filename},
    {"HF/AddLocalPotential/Scale",      &GlobalSpecification::hf_addlocal_scale},
    // Basis
    {"Basis/Valence",           &GlobalSpecification::basis_valence},
    {"Basis/FrozenCore",        &GlobalSpecification::basis_frozen_core},
    {"Basis/IncludeValence",    &GlobalSpecification::basis_include_valence},
    {"Basis/ExcludeValence",    &GlobalSpecification::basis_exclude_valence},
    {"Basis/Residue",           &GlobalSpecification::basis_residue},
    {"Basis/InjectOrbitals",    &GlobalSpecification::basis_inject_orbitals},
    {"Basis/--reorthogonalise", &GlobalSpecification::basis_reorthogonalise},
    {"Basis/--hf-basis",        &GlobalSpecification::basis_hf},
    {"Basis/--bspline-basis",   &GlobalSpecification::basis_bspline},
    {"Basis/--xr-basis",        &GlobalSpecification::basis_xr},
    {"Basis/HFOrbitals",        &GlobalSpecification::basis_hf_orbitals},
    {"Basis/CustomOrbitals",    &GlobalSpecification::basis_custom_orbitals},
    // Basis/BSpline
    {"Basis/BSpline/Rmax",          &GlobalSpecification::basis_bspline_rmax},
    {"Basis/BSpline/R0",            &GlobalSpecification::basis_bspline_r0},
    {"Basis/BSpline/K",             &GlobalSpecification::basis_bspline_k},
    {"Basis/BSpline/N",             &GlobalSpecification::basis_bspline_N},
    {"Basis/BSpline/SplineType",    &GlobalSpecification::basis_bspline_splinetype},
    // CI
    {"CI/LeadingConfigurations",                &GlobalSpecification::ci_leading_configurations},
    {"CI/LeadingRelativisticConfigurations",    &GlobalSpecification:: ci_leading_rel_configurations},
    {"CI/ExtraConfigurations",                  &GlobalSpecification::ci_extra_configurations},
    {"CI/ExtraRelativisticConfigurations",      &GlobalSpecification::ci_extra_rel_configurations},
    {"CI/ElectronExcitations",                  &GlobalSpecification::ci_electron_excitations},
    {"CI/ExcitationBounds",                     &GlobalSpecification::ci_excitation_bounds},
    {"CI/HoleExcitations",                      &GlobalSpecification::ci_hole_excitations},
    {"CI/EvenParityTwoJ",                       &GlobalSpecification::ci_even_parity_twoj},
    {"CI/OddParityTwoJ",                        &GlobalSpecification::ci_odd_parity_twoj},
    {"CI/NumSolutions",                         &GlobalSpecification::ci_num_solutions},
    {"CI/--all-symmetries",                     &GlobalSpecification::ci_all_symmetries},
    {"CI/--gfactors",                           &GlobalSpecification::ci_gfactors},
    {"CI/--no-gfactors",                        &GlobalSpecification::ci_no_gfactors},
    {"CI/--memory-saver",                       &GlobalSpecification::ci_memory_saver},
    {"CI/--single-configuration-CI",            &GlobalSpecification::ci_single_configuration_ci},
    {"CI/--print-configurations",               &GlobalSpecification::ci_print_configurations},
    {"CI/--print-relativistic-configurations",  &GlobalSpecification:: ci_print_rel_configurations},
    {"CI/--scalapack",                          &GlobalSpecification::ci_scalapack},
    {"CI/MaxEnergy",                            &GlobalSpecification::ci_max_energy},
    {"CI/ConfigurationAverageEnergyRange",      &GlobalSpecification::ci_configuration_average_energy_range},
    {"CI/ChunkSize",                            &GlobalSpecification::ci_chunksize},
    {"CI/--sort-matrix-by-configuration",       &GlobalSpecification::ci_sort_matrix_by_configuration}

};

P::specification_set<GlobalSpecification> global_specifications_dict(global_specifications, P::keys_lc_nows);

LatticeConfig GlobalSpecification::getLatticeConfig() const {
    if (lattice_exponential) {
        LatticeExpConfig config;

        config.H = lattice_H;
        if (lattice_num_points>0) config.num_points = lattice_num_points;
        if (lattice_start_point>0) config.start_point = lattice_start_point;
        return config;
    }
    else {
        LatticeHybridConfig config;

        if (lattice_num_points>0) config.num_points = lattice_num_points;
        if (lattice_start_point>0) config.start_point = lattice_start_point;
        if (lattice_end_point>0) config.end_point = lattice_end_point;
        return config;
    }
}

BasisConfig GlobalSpecification::getBasisConfig() const {
    if(basis_xr) {
        XRBasisConfig config;
        config.frozen_core = basis_frozen_core;
        config.valence_basis = basis_valence;
        config.include_valence = basis_include_valence;
        config.exclude_valence = basis_exclude_valence;
        config.residue = basis_residue;
        config.inject_orbitals = basis_inject_orbitals;
        config.hf_orbitals = basis_hf_orbitals;
        config.custom_orbitals = basis_custom_orbitals;

        return(config);
    } else if (basis_hf) {
        HFBasisConfig config;
        config.frozen_core = basis_frozen_core;
        config.valence_basis = basis_valence;
        config.include_valence = basis_include_valence;
        config.exclude_valence = basis_exclude_valence;
        config.residue = basis_residue;
        config.inject_orbitals = basis_inject_orbitals;
        config.hf_orbitals = basis_hf_orbitals;

        return(config);
    } else {
        BSplineBasisConfig config;
        if(basis_bspline_rmax > 0) {
            config.RMax = basis_bspline_rmax;
        // If Basis/BSpline/RMax is unset, then it defaults to the value of Lattice/EndPoint,
        // otherwise it defaults to 40.0;
        } else if(lattice_end_point > 0) {
            config.RMax = lattice_end_point;
        } else {
            config.RMax = 40.0;
        }
        config.R0 = basis_bspline_r0;
        config.K = basis_bspline_k;
        
        // Finally, get the kind of spline to use. Convert eveything to lower case, since we don't
        // really care about capitalisation here. Default to Reno type
        std::string spline_type = basis_bspline_splinetype;
        std::transform(spline_type.begin(), spline_type.end(), spline_type.begin(), ::toupper);
        if(spline_type == "vanderbilt") {
            config.spline_type = SplineType::Vanderbilt;
        } else if (spline_type == "notredame") {
            config.spline_type = SplineType::NotreDame;
        } else {
            config.spline_type = SplineType::Reno;
        }
        config.frozen_core = basis_frozen_core;
        config.valence_basis = basis_valence;
        config.include_valence = basis_include_valence;
        config.exclude_valence = basis_exclude_valence;
        config.residue = basis_residue;
        config.inject_orbitals = basis_inject_orbitals;
        config.hf_orbitals = basis_hf_orbitals;

        return(config);
    }
}

HFConfig GlobalSpecification::getHFConfig() const {
    HFConfig config;

    // TODO: Need some options to go through all the different possible decorators here
    config.Z = Z;
    config.charge = hf_charge;
    config.N = hf_N;
    config.configuration = hf_configuration;
    config.breit = hf_breit;
    config.sms = hf_sms;
    config.nms = hf_nms;
    config.only_rel_nms = hf_only_rel_nms;
    config.nonrel_mass_shift = hf_nonrel_mass_shift;
    config.include_lower_mass = hf_include_lower_mass;
    config.local_exchange = hf_local_exchange;
    config.xalpha = hf_xalpha;
    config.alpha_squared_variation = hf_alpha_squared_variation;
    config.nuclear_thickness = hf_nuclear_thickness;
    config.nuclear_radius = hf_nuclear_radius;
    config.nuclear_inverse_mass = hf_nuclear_inverse_mass;

    return(config);
}

std::string importSpecificationFile(GlobalSpecification& gs, const std::string& fileName) {
    std::ifstream in(fileName);
    if (!in) return "unable to open input file '"+fileName+"'";

    auto h = custom_import_ini(gs, global_specifications_dict, in);
    if (!h) return P::explain(h.error(), true);

    return "";
}

std::string importSpecificationKV(GlobalSpecification& gs, const std::string& assignment) {
    auto h = P::import_k_eq_v(gs, global_specifications_dict, assignment);
    if (!h) return P::explain(h.error(), true);
    else return "";
}

// Perform global validation of specifications. Return non-empty error message on failure.
std::string validateSpecification(const GlobalSpecification& gs) {
    // Check for consistent Lattice settings:

    if (gs.lattice_exponential && gs.lattice_end_point>0)
        return "Lattice/EndPoint is ignored if Lattice/--exp-lattice is set";

    // Check for consistency in HF decorators (e.g. cannot request both nonrelativistic and
    // relativistic-only constraints)
    if (gs.hf_only_rel_nms && gs.hf_nonrel_mass_shift)
        return "HF/--only-relativistic-nms and HF/--nonrelativistic-mass-shift cannot both be set at the same time";
    // Must have at least one of HF/N and HF/Charge, so error out if these are missing
    if((!gs.hf_charge) && (!gs.hf_N))
        return "Must specify at least one of HF/N or HF/Charge";

    return "";
}

// Custom INI import implementation:

using P::ini_record;
using P::ini_record_kind;

ini_record custom_ini_parser(std::string_view v) {
    using token = ini_record::token;
    using size_type = std::string_view::size_type;

    constexpr size_type npos = std::string_view::npos;
    constexpr std::string_view ws{" \t\f\v\r\n"};

    auto comment_at = [](std::string_view v, size_type p) { return v.substr(p, 2)=="//"; };
    size_type b = v.find_first_not_of(ws);

    // empty or comment?
    if (b==npos || comment_at(v, b)) return ini_record{ini_record_kind::empty};

    // section heading?
    if (v[b]=='[') {
        size_type e = v.find(']');

        // check for malformed heading
        if (e==npos) return {ini_record_kind::syntax_error, token("", b+1)};

        if (e+1<v.length()) {
            auto epilogue = v.find_first_not_of(ws, e+1);
            if (epilogue!=npos && !comment_at(v, epilogue)) {
                return {ini_record_kind::syntax_error, token("", epilogue)};
            }
        }

        b = v.find_first_not_of(ws, b+1);
        e = v.find_last_not_of(ws, e-1);

        return {ini_record_kind::section, token(v.substr(b, e>=b? e+1-b: 0), b+1)};
    }

    // expect key first, followed by ws and eol, =, or //.
    size_type j = std::min(v.find('=', b), v.find("//", b));
    token key_token{v.substr(b, j==b? 0: v.find_last_not_of(ws, j-1)+1-b), b+1};

    // key without value?
    if (j==npos || v[j]!='=') {
        return {ini_record_kind::key, key_token};
    }

    // skip to text after =, look for value
    size_type eq = j;
    size_type value_cindex = eq;

    if (j<v.length()) {
        j = v.find_first_not_of(ws, j+1);
        if (j!=npos && !comment_at(v, j)) {
            value_cindex = j+1;

            // if value is not quoted, take text up to eol or first eol, discarding trailing ws
            if (v[j]!='\'') {
                size_type end = v.find("//", j);
                if (end!=npos) --end;

                return {ini_record_kind::key_value, key_token,
                        token{v.substr(j, v.find_last_not_of(ws, end)-j+1), value_cindex}};
            }
            else {
                // quoted value; take until next unescaped '
                std::string value;
                size_type epilogue = npos;
                bool esc = false;

                for (size_type i = j+1; i<v.length(); ++i) {
                    if (esc) {
                        value += v[i];
                        esc = false;
                    }
                    else if (v[i]=='\'') {
                        epilogue = i+1;
                        break;
                    }
                    else if (v[i]=='\\') {
                        esc = true;
                    }
                    else {
                        value += v[i];
                    }
                }

                // unterminated quoted value or escaped eol?
                if (epilogue==npos || esc) {
                    return {ini_record_kind::syntax_error, token("", j)};
                }

                // extra stuff following value that is not a comment?
                if (epilogue<v.length()) {
                    epilogue = v.find_first_not_of(ws, epilogue);
                    if (epilogue!=npos && !comment_at(v, epilogue)) {
                        return {ini_record_kind::syntax_error, token("", epilogue)};
                    }
                }

                return {ini_record_kind::key_value, key_token, token{value, value_cindex}};
            }
        }
    }
    // key with empty value
    return {ini_record_kind::key_value, key_token, token{"", eq}};
}

// Use a custom line-by-line ini importer to handle relative section headings
template <typename Record>
P::hopefully<void> custom_import_ini(Record& rec, const P::specification_set<Record>& specs, std::istream& in)
{
    constexpr auto npos = std::string_view::npos;

    P::ini_style_importer importer(custom_ini_parser, in);
    while (importer) {
        std::string prev_sec{importer.section()};
        auto h = importer.run_one(rec, specs, P::default_reader(), "/");

        // Ignore unrecognized keys until specification coverage is complete
        if (!h) {
             if (h.error().error!=P::failure::unrecognized_key) return P::unexpected(h.error());
             else continue; 
        }

        if (h.value() != P::ini_record_kind::section) continue;

        std::string_view new_sec = importer.section();
        if (new_sec.substr(0, 2)=="./") {
            std::string section{prev_sec};
            section += new_sec.substr(1);
            importer.section(section);
        }
        else if (new_sec.substr(0, 3)=="../") {
            std::string_view prefix(prev_sec);
            do {
                new_sec.remove_prefix(3);
                if (auto tail = prefix.rfind('/'); tail != npos) {
                    prefix = prefix.substr(0, tail);
                }
            } while (new_sec.substr(0, 3)=="../");

            std::string section{prefix};
            section += "/";
            section += new_sec;
            importer.section(section);
        }
    }
    return {};
}

} // namespace Ambit
