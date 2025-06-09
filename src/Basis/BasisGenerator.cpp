#include "BasisGenerator.h"
#include "Include.h"
#include "HartreeFock/ConfigurationParser.h"
#include "HartreeFock/Integrator.h"
#include "HartreeFock/HartreeFocker.h"
#include "HartreeFock/NucleusDecorator.h"
#include "ExternalField/NormalMassShiftDecorator.h"
#include "ExternalField/SpecificMassShiftDecorator.h"
#include "ExternalField/TwoBodySMSOperator.h"
#include "ExternalField/BreitHFDecorator.h"
#include "ExternalField/RadiativePotential.h"
#include "ExternalField/NuclearPolarisability.h"
#include "ExternalField/YukawaPotential.h"
#include "Universal/Lattice.h"
#include "Universal/LatticeConfig.h"
#include "Specification/Specification.h"

namespace Ambit
{
BasisGenerator::BasisGenerator(pLattice lat, MultirunOptions& userInput, GlobalSpecification specification, pPhysicalConstant physical_constant):
    lattice(lat), user_input(userInput), physical_constant(physical_constant), 
    specification(std::move(specification)),
    open_core(nullptr)
{
    orbitals = pOrbitalManager(new OrbitalManager(lattice));

    // User input configuration objects
    basis_config = specification.getBasisConfig();
    hf_config = specification.getHFConfig();
    lattice_config = specification.getLatticeConfig();
}

BasisGenerator::~BasisGenerator()
{}

void BasisGenerator::InitialiseHF(pHFOperator& undressed_hf)
{
    unsigned int Z = hf_config.Z;

    // HF/Charge and HF/N may or may not be present in the input file, so either grab the value if
    // it exists, or calculate it based on the electronic parameters
    int Charge;
    if(hf_config.charge)
    {
        Charge = hf_config.charge.value();
    } else 
    {   if(hf_config.N && Z >= hf_config.N.value())
        {
            Charge = Z - hf_config.N.value();
        } else {
            Charge = 0;
        }
    }

    //TODO: Error message if Charge or N is missing or incorrect.
    std::string config = hf_config.configuration;

    // Get orbitals and occupancies
    std::string open_shell_string;
    std::string closed_shell_string;
    size_t colon_pos = config.find(':');
    if(colon_pos == std::string::npos)
    {   open_shell_string = config;
        closed_shell_string = config;
    }
    else
    {   open_shell_string = config;
        open_shell_string.erase(colon_pos, 1);
        closed_shell_string = config.substr(0,colon_pos);
    }

    OccupationMap open_shell_occupations = ConfigurationParser::ParseFractionalConfiguration(open_shell_string);

    // Set open_core occupancies
    open_core->SetOccupancies(open_shell_occupations);
    if(open_core->NumElectrons() != Z - Charge)
    {   *errstream << "Core::BuildFirstApproximation: Incorrect electron count in configuration." << std::endl;
        exit(1);
    }

    pIntegrator integrator(new SimpsonsIntegrator(lattice));
    pODESolver ode_solver(new AdamsSolver(integrator));
    pCoulombOperator coulomb(new CoulombOperator(lattice, ode_solver));

    if(physical_constant == nullptr)
    {
        physical_constant = pPhysicalConstant(new PhysicalConstant());
        double alpha_variation = hf_config.alpha_squared_variation;
        if(alpha_variation)
            physical_constant->SetAlphaSquaredIncreaseRatio(alpha_variation);
    }

    undressed_hf = pHFOperator(new HFOperator(Z, open_core, physical_constant, integrator, coulomb));
    hf = undressed_hf;

    // Add nuclear potential
    double nuclear_radius = hf_config.nuclear_radius;
    if(nuclear_radius)
    {
        nucleus = std::make_shared<NucleusDecorator>(hf, coulomb, integrator);
        double nuclear_thickness = hf_config.nuclear_thickness;
        nucleus->SetFermiParameters(nuclear_radius, nuclear_thickness);
        nucleus->SetCore(open_core);
        *outstream << "Nuclear RMS radius = " << nucleus->CalculateNuclearRMSRadius() << std::endl;
        hf = nucleus;
        undressed_hf = hf;
    }

    // Hartree operator
    hartreeY = pHartreeY(new HartreeY(integrator, coulomb));

    // Add additional operators
    double NuclearInverseMass = hf_config.nuclear_inverse_mass;
    if(NuclearInverseMass)
    {
        bool do_nms = hf_config.nms;
        bool do_sms = hf_config.sms;
        bool nonrel_ms = hf_config.nonrel_mass_shift;
        bool relativistic_nms = hf_config.only_rel_nms;
        bool lower_sms = hf_config.include_lower_mass;

        // Default: do specific mass shift
        if(!do_nms && !do_sms && !relativistic_nms)
            do_sms = true;

        if(do_nms)
        {
            pNormalMassShiftDecorator nms_op = std::make_shared<NormalMassShiftDecorator>(hf, relativistic_nms, nonrel_ms);
            nms_op->SetInverseMass(NuclearInverseMass);
            nms_op->SetCore(open_core);
            hf = nms_op;
        }

        if(do_sms)
        {
            // HF decorator
            pSpecificMassShiftDecorator sms_op = std::make_shared<SpecificMassShiftDecorator>(hf, nonrel_ms, lower_sms);
            sms_op->SetInverseMass(NuclearInverseMass);
            sms_op->SetCore(open_core);
            hf = sms_op;

            // HartreeY decorator
            pSMSOperator Ysms;
            if(nonrel_ms)
            {
                Ysms = std::make_shared<TwoBodySMSOperator>(hartreeY, lower_sms);
            }
            else
            {
                double Zalpha = Z * physical_constant->GetAlpha();
                Ysms = std::make_shared<TwoBodySMSOperator>(hartreeY, Zalpha);
            }

            Ysms->SetInverseMass(NuclearInverseMass);
            hartreeY = Ysms;
        }
    }

    if(hf_config.breit)
    {
        pHartreeY breit = std::make_shared<BreitZero>(std::make_shared<HartreeYBase>(), integrator, coulomb);
        pHFOperator breit_hf = std::make_shared<BreitHFDecorator>(hf, breit);
        hf = breit_hf;

        // Decorate HartreeY function
        hartreeY = std::make_shared<BreitZero>(hartreeY, integrator, coulomb);
    }

    // QED options
    if(hf_config.qed_config)
    {
        // Get the nuclear RMS radius to use for QED calculations. This will be either a user
        // specified value (QED/NuclearRMSRadius) or the value used in the rest of the calculation 
        double nuclear_rms_radius;
        if(hf_config.qed_config->nuclear_rms_radius)
            nuclear_rms_radius = hf_config.qed_config->nuclear_rms_radius.value();
        else
            nuclear_rms_radius = GetNuclearRMSRadius();

        // Uehling options
        if(hf_config.qed_config->uehling)
        {
            pUehlingDecorator uehling;

            if(hf_config.qed_config->use_nuclear_density)
            {   
                uehling.reset(new UehlingDecorator(hf, nucleus->GetNuclearDensity()));
            }
            else
            {   
                uehling.reset(new UehlingDecorator(hf, nuclear_rms_radius));
            }

            hf = uehling;
        }

        // Self-energy options
        if(hf_config.qed_config->self_energy)
        {
            pElectricSelfEnergyDecorator electricQED;
            pMagneticSelfEnergyDecorator magneticQED;

            // Use the nuclear density from the rest of the calculation
            if(hf_config.qed_config->use_nuclear_density)
            {
                if(!hf_config.qed_config->no_magnetic)
                {   
                    magneticQED.reset(new MagneticSelfEnergyDecorator(hf, nucleus->GetNuclearDensity()));
                    hf = magneticQED;
                }
                if(!hf_config.qed_config->no_electric)
                {   
                    electricQED.reset(new ElectricSelfEnergyDecorator(hf, nucleus->GetNuclearDensity()));
                    hf = electricQED;
                }
            }
            // Use a user-specified value for nuclear density/RMS
            else
            {                   
                if(!hf_config.qed_config->no_magnetic)
                {   
                    magneticQED.reset(new MagneticSelfEnergyDecorator(hf, nuclear_rms_radius));
                    hf = magneticQED;
                }
                if(!hf_config.qed_config->no_electric)
                {   
                    bool skip_offmass = hf_config.qed_config->skip_offmass; 
                    electricQED.reset(new ElectricSelfEnergyDecorator(hf, nuclear_rms_radius, !skip_offmass));
                    hf = electricQED;
                }
            }
        }
    }

    // Nuclear polarisability options
    if(hf_config.nuclear_polarisability_config)
    {
        double alphaE = hf_config.nuclear_polarisability_config->alphaE; 
        double Ebar = hf_config.nuclear_polarisability_config->ebarMeV; 

        hf = std::make_shared<NuclearPolarisability>(hf, alphaE, Ebar);
    }

    // Yukawa options
    if(hf_config.yukawa_config)
    {
        double mass = 1.0;
        // Note that there are multiple different, equivalent ways of specifying the mass. The
        // specification guarantees that exactly one of these is set, to avoid conflicting values
        if(hf_config.yukawa_config->mass)
            mass = hf_config.yukawa_config->mass.value();
        else if (hf_config.yukawa_config->massEV)
            mass = hf_config.yukawa_config->massEV.value()/MathConstant::Instance()->ElectronMassInEV;
        else if(hf_config.yukawa_config->rc)
            mass = 1./(physical_constant->GetAlpha() * hf_config.yukawa_config->rc.value());

        double scale = hf_config.yukawa_config->scale;
        hf = std::make_shared<YukawaDecorator>(hf, mass, scale);
    }

    if(hf_config.local_exchange)
    {
        double xalpha = hf_config.xalpha;
        pHFOperator localexch = std::make_shared<LocalExchangeApproximation>(hf, coulomb, xalpha);
        localexch->SetCore(open_core);
        hf = localexch;
        hf->IncludeExchange(false);
        undressed_hf = hf;
    }

    // Local potential decorator options
    if(hf_config.local_potential_config)
    {
        std::string filename = hf_config.local_potential_config->filename;
        double scale = hf_config.local_potential_config->scale;
        pImportedPotentialDecorator loc(new ImportedPotentialDecorator(hf, filename));
        loc->SetScale(scale);
        hf = loc;
    }

    // Set closed core occupancies
    OccupationMap closed_shell_occupations = ConfigurationParser::ParseFractionalConfiguration(closed_shell_string);

    // Make closed shell core. Ensure that all shells are completely filled.
    for(OccupationMap::iterator it = closed_shell_occupations.begin(); it != closed_shell_occupations.end(); it++)
        it->second = 2. * abs(it->first.Kappa());

    // Create closed core with empty pointers for all occupied orbitals
    closed_core = pCore(new Core(lattice));
    closed_core->SetOccupancies(closed_shell_occupations);
}

void BasisGenerator::SetOrbitalMaps()
{
    // Transfer from all to closed core
    OrbitalMap& all = *orbitals->all;
    for(auto core_occupation: closed_core->GetOccupancies())
    {
        closed_core->AddState(all.GetState(core_occupation.first));
    }
    orbitals->core = closed_core;

    // Hole and deep states.
    // Easiest to start with all core states in deep and modify from there.
    orbitals->deep = std::make_shared<OrbitalMap>(lattice);
    orbitals->hole = std::make_shared<OrbitalMap>(lattice);
    *orbitals->deep = *orbitals->core;

    OrbitalMap& deep = *orbitals->deep;
    OrbitalMap& hole = *orbitals->hole;

    // std::visitor to deal with variant types of the basis config (e.g. different kinds of basis
    // functions)
    // TODO EVK: Not sure if I like having to do this every time I access the base class...
    std::string deep_states = std::visit([](auto &&var) {return var.frozen_core;},basis_config);
    if(deep_states.length())
    {
        std::vector<int> max_deep_pqns = ConfigurationParser::ParseBasisSize(deep_states);
        auto it = deep.begin();
        while(it != deep.end())
        {
            // Not deep
            if(it->first.L() >= max_deep_pqns.size() ||
               it->first.PQN() > max_deep_pqns[it->first.L()])
            {
                hole.AddState(it->second);
                it = deep.erase(it);
            }
            else
                it++;
        }
    }

    // IncludeValence moves deep orbitals into valence holes
    int num_unfrozen = std::visit([](auto&& var) -> int {return(var.include_valence.size());}, basis_config);
    for(int i = 0; i < num_unfrozen; i++)
    {
        NonRelInfo nrorb = std::visit([&](auto&& var) -> NonRelInfo {
                return(ConfigurationParser::ParseOrbital(var.include_valence[i]));
                }, basis_config);
        for(auto& orbinfo: nrorb.GetRelativisticInfos())
        {
            auto it = deep.find(orbinfo);
            if(it != deep.end())
            {
                hole.AddState(it->second);
                deep.erase(it);
            }
            else
            {   *errstream << "BasisGenerator::SetOrbitalMaps: Basis/IncludeValence "
                           << orbinfo.Name() << " not found in frozen core." << std::endl;
            }
        }
    }

    // Transfer from all to excited states
    std::string valence_states = std::visit([](auto&& var) -> std::string {return var.valence_basis;},basis_config);
    std::vector<int> max_pqn_per_l = ConfigurationParser::ParseBasisSize(valence_states);

    orbitals->particle = std::make_shared<OrbitalMap>(lattice);
    OrbitalMap& particle = *orbitals->particle;

    for(auto& orbital: all)
    {
        if(orbital.first.L() < max_pqn_per_l.size()
           && orbital.first.PQN() <= max_pqn_per_l[orbital.first.L()]
           && closed_core->GetState(orbital.first) == nullptr)
        {
            particle.AddState(orbital.second);
        }
    }

    // high (virtual) states.
    // Two type magic things happening here: std::visit to concretize the basis variant type, then
    // a value_or since the MBOPT basis might not exist
    std::string virtual_states = std::visit([](auto&& var) -> std::string {
            return(var.MBPT_basis.value_or(""));
            }, basis_config);
    orbitals->excited = std::make_shared<OrbitalMap>(lattice);
    orbitals->high = std::make_shared<OrbitalMap>(lattice);

    OrbitalMap& excited = *orbitals->excited;
    OrbitalMap& high = *orbitals->high;

    if(virtual_states.size())
    {
        max_pqn_per_l = ConfigurationParser::ParseBasisSize(virtual_states);

        for(auto& orbital: all)
        {
            if(orbital.first.L() < max_pqn_per_l.size()
               && orbital.first.PQN() <= max_pqn_per_l[orbital.first.L()]
               && closed_core->GetState(orbital.first) == nullptr)
            {
                excited.AddState(orbital.second);
                if(particle.GetState(orbital.first) == nullptr)
                    high.AddState(orbital.second);
            }
        }
    }
    else
    {   // high is empty, excited is just particles
        *orbitals->excited = *orbitals->particle;
    }

    // ExcludeValence moves particle orbitals into high states
    int num_excluded = std::visit([](auto&& var) -> int {return(var.exclude_valence.size());}, basis_config);
    for(int i = 0; i < num_excluded; i++)
    {
        NonRelInfo nrorb = std::visit([&](auto&& var) -> NonRelInfo {
                return(ConfigurationParser::ParseOrbital(var.exclude_valence[i]));
                }, basis_config);
        for(auto& orbinfo: nrorb.GetRelativisticInfos())
        {
            auto it = particle.find(orbinfo);
            if(it != particle.end())
            {
                high.AddState(it->second);
                particle.erase(it);
            }
            else
            {   *errstream << "BasisGenerator::SetOrbitalMaps: Basis/ExcludeValence "
                           << orbinfo.Name() << " not found in valence particle set." << std::endl;
            }
        }
    }

    // Make valence orbitals
    orbitals->valence = std::make_shared<OrbitalMap>(lattice);
    orbitals->valence->AddStates(*orbitals->particle);
    orbitals->valence->AddStates(*orbitals->hole);
}

void BasisGenerator::UpdateNonSelfConsistentOperators()
{
    if(hf_config.qed_config->use_electron_screening)
    {
        if(nucleus == nullptr || !hf_config.qed_config->use_nuclear_density)
        {
            *logstream << "Cannot have screened Uehling without finite sized nucleus." << std::endl;
            return;
        }

        RadialFunction density(nucleus->GetNuclearDensity());
        for(const auto& orb: *open_core)
        {
            density -= orb.second->GetDensity() * open_core->GetOccupancy(orb.first);
        }

        pUehlingDecorator uehling;
        pMagneticSelfEnergyDecorator magneticQED;
        pElectricSelfEnergyDecorator electricQED;

        // Traverse HFOperatorDecorator stack in hf to find QED decorators.
        std::shared_ptr<HFBasicDecorator> hfdecorator = std::dynamic_pointer_cast<HFBasicDecorator>(hf);
        while(hfdecorator)
        {
            uehling = std::dynamic_pointer_cast<UehlingDecorator>(hfdecorator);
            magneticQED = std::dynamic_pointer_cast<MagneticSelfEnergyDecorator>(hfdecorator);
            electricQED = std::dynamic_pointer_cast<ElectricSelfEnergyDecorator>(hfdecorator);

            hfdecorator = std::dynamic_pointer_cast<HFBasicDecorator>(hfdecorator->GetWrapped());
        }

        if(uehling)
            uehling->GenerateUehling(density);
        if(magneticQED)
            magneticQED->GenerateMagnetic(density);
        if(electricQED)
        {   electricQED->GenerateEhigh(density);
            electricQED->GenerateElow(density);
        }
    }
}

pCore BasisGenerator::GenerateHFCore(pCoreConst open_shell_core)
{
    open_core = pCore(new Core(lattice));
    hf = nullptr;
    hartreeY = nullptr;

    if(open_shell_core)
    {   // Copy, use same lattice
        open_core.reset(open_shell_core->Clone());
        lattice = open_core->GetLattice();
        orbitals = pOrbitalManager(new OrbitalManager(lattice));
    }

    pHFOperator undressed_hf;
    InitialiseHF(undressed_hf);

    // Create Hartree-Fock solver; define integrators.
    pIntegrator integrator(new SimpsonsIntegrator(lattice));
    pODESolver ode_solver(new AdamsSolver(integrator));
    HartreeFocker HF_Solver(ode_solver);

    // TODO: Check occupancies match
    if(!open_shell_core)
    {   HF_Solver.StartCore(open_core, undressed_hf);
        HF_Solver.SolveCore(open_core, undressed_hf);
    }

    // Update any non-self-consistent screening operators (e.g. radiative potentials)
    UpdateNonSelfConsistentOperators();
    HF_Solver.SolveCore(open_core, hf);

    // Resize lattice according to larger of core or user input.
    unsigned int core_size = open_core->LargestOrbitalSize();
    unsigned int original_lattice_size = std::visit([](auto&& var) -> unsigned int {
            return(var.num_points);
            }, lattice_config);

    lattice->resize(mmax(core_size, original_lattice_size));

    return open_core;
}

pHFOperator BasisGenerator::RecreateBasis(pOrbitalManager orbital_manager)
{
    lattice = orbital_manager->GetLattice();
    open_core = pCore(new Core(lattice));

    pHFOperator undressed_hf;
    InitialiseHF(undressed_hf);

    // Copy orbitals from orbital_manager to open_core
    for(auto pair: *open_core)
    {
        pOrbital state = orbital_manager->all->GetState(pair.first);
        if(state == nullptr)
        {   *errstream << "BasisGenerator::CreateHFOperator(): orbital " << pair.first.Name() << " not found." << std::endl;
            exit(1);
        }
        *pair.second = *state;
    }

    hf->SetCore(open_core);
    UpdateNonSelfConsistentOperators();

    // Modify orbital manager maps according to input file
    orbitals = orbital_manager;
    SetOrbitalMaps();

    return hf;
}

pOrbitalManagerConst BasisGenerator::GenerateBasis()
{
    // Make sure hf is correct
    std::string residue ;//= user_input("Basis/Residue", "");
    std::visit([](auto &&var) -> std::string {
            return(var.residue);
            }, basis_config);

    if(residue.empty())
    {
        hf->SetCore(open_core);
    }
    else
    {   // Strip any errant colon
        size_t colon_pos = residue.find(':');
        if(colon_pos != std::string::npos)
            residue.erase(colon_pos, 1);

        // No need to clone, since we are not changing the core orbitals
        pCore residual_core = std::make_shared<Core>(*open_core);

        OccupationMap residual_occupations = ConfigurationParser::ParseFractionalConfiguration(residue);
        residual_core->SetOccupancies(residual_occupations);

        hf->SetCore(residual_core);
    }

    // Generate excited states
    std::string all_states = user_input("Basis/BasisSize", "");
    if(all_states.empty())
        all_states = user_input("MBPT/Basis", "");
    if(all_states.empty())
        all_states = user_input("Basis/ValenceBasis", "");

    bool reorth = user_input.search("Basis/--reorthogonalise"); // Perform extra orthogonalisation

    std::vector<int> max_pqn_per_l = ConfigurationParser::ParseBasisSize(all_states);
    pOrbitalMap excited;

    if(user_input.search("Basis/--hf-basis"))
    {   excited = GenerateHFExcited(max_pqn_per_l);
    }
    else if(user_input.search("Basis/--xr-basis"))
    {   excited = GenerateXRExcited(max_pqn_per_l);
    }
    else // default "Basis/--bspline-basis"
    {   user_input.search("Basis/--bspline-basis"); // Just to clear UFO from user_input.
        excited = GenerateBSplines(max_pqn_per_l);

        // Replace requested valence orbitals with HF orbitals
        std::string hf_valence_states = user_input("Basis/HFOrbitals", "");
        if(!hf_valence_states.empty())
        {
            UpdateHFOrbitals(ConfigurationParser::ParseBasisSize(hf_valence_states), excited);
        }
    }

    // Inject any special orbitals from another basis, and push the old ones to higher pqn
    int number_injected = user_input.vector_variable_size("Basis/InjectOrbitals");
    for(int i = 0; i < number_injected; i++)
    {
        std::string inject_string = user_input("Basis/InjectOrbitals", "", i);
        InjectOrbitals(inject_string, excited);
        reorth = true;
    }

    // Place all orbitals in orbitals->all.
    // Finally create orbitals->all and the state index
    orbitals->all = pOrbitalMap(new OrbitalMap(lattice));
    orbitals->all->AddStates(*open_core);
    orbitals->all->AddStates(*excited);

    orbitals->MakeStateIndexes();

    // Organise orbitals
    SetOrbitalMaps();

    if(reorth)
    {
        for(auto excited_orbital_pair: *orbitals->excited)
            Orthogonalise(excited_orbital_pair.second);
    }

    if(DebugOptions.OutputHFExcited())
    {   OrbitalInfo max_i(-1, 1), max_j(-1, 1);
        double orth = TestOrthogonality(max_i, max_j);
        *logstream << "<" << max_i.Name() << " | " << max_j.Name() << "> = " << orth << std::endl;
    }

    return orbitals;
}

void BasisGenerator::InjectOrbitals(const std::string& input, pOrbitalMap excited) const
{
    std::string inputfile;
    std::string origin;
    std::string target;
    size_t colon_pos = input.find(':');
    size_t arrow_pos = input.find("->");

    if(colon_pos == std::string::npos || arrow_pos == std::string::npos)
    {
        *errstream << "Basis/InjectOrbitals incorrectly specified." << std::endl;
        exit(1);
    }

    inputfile = input.substr(0, colon_pos);
    origin = input.substr(colon_pos+1, arrow_pos);
    target = input.substr(arrow_pos+2, input.size());

    // Remove whitespace from inputfile; import orbitals
    inputfile.erase(std::remove_if(inputfile.begin(), inputfile.end(), isspace), inputfile.end());
    pOrbitalManager imported_orbitals = std::make_shared<OrbitalManager>(inputfile);

    NonRelInfo nonrelorigin = ConfigurationParser::ParseOrbital(origin);
    NonRelInfo nonreltarget = ConfigurationParser::ParseOrbital(target);

    for(int step = 1; step <= (nonreltarget.L()? 2: 1); step++)
    {
        // Get injected orbital
        pOrbital inject;
        if(step == 1)
        {   inject = imported_orbitals->all->GetState(nonrelorigin.GetFirstRelativisticInfo());
            inject->SetKappa(nonreltarget.GetFirstRelativisticInfo().Kappa());
        }
        else
        {   inject = imported_orbitals->all->GetState(nonrelorigin.GetSecondRelativisticInfo());
            inject->SetKappa(nonreltarget.GetSecondRelativisticInfo().Kappa());
        }
        inject->SetPQN(nonreltarget.PQN());

        // Move all orbitals with higher PQN out
        std::vector<pOrbital> moved_orbitals;
        int max_moved_pqn = 0;
        auto it = excited->begin();
        while(it != excited->end())
        {
            if((it->first.Kappa() == inject->Kappa()) &&
               (it->first.PQN() >= inject->PQN()))
            {
                pOrbital moved = it->second;

                // Increment the PQN
                moved->SetPQN(it->first.PQN() + 1);
                max_moved_pqn = mmax(max_moved_pqn, moved->PQN());
                moved_orbitals.push_back(moved);

                it = excited->erase(it);
            }
            else
                ++it;
        }

        // Inject new orbital
        excited->AddState(inject);

        // Move other orbitals back
        for(auto orbital: moved_orbitals)
        {
            if(orbital->PQN() != max_moved_pqn)
                excited->AddState(orbital);
        }
    }
}

void BasisGenerator::CreateBruecknerOrbitals(pBruecknerDecorator brueckner)
{
    // Set hf operator to brueckner for the rest of the calculation
    hf = brueckner;

    pOrbitalMap orbitals_to_update = orbitals->valence;
    if(user_input.search("MBPT/Brueckner/--excited"))
        orbitals_to_update = orbitals->excited;

    // Get max PQN for l
    std::vector<int> max_pqn;
    for(auto& pair: *orbitals_to_update)
    {
        int l = pair.first.L();
        if(l+1 > max_pqn.size())
            max_pqn.resize(l+1);

        max_pqn[l] = mmax(max_pqn[l], pair.first.PQN());
    }

    pOrbitalMap brueckner_orbitals = GenerateBSplines(max_pqn);

    for (auto &pair: *orbitals_to_update)
    {
        pOrbital brueckner_orbital = brueckner_orbitals->GetState(pair.first);

        // Copy back to orbital manager
        if (brueckner_orbital)
            *pair.second = *brueckner_orbital;
    }

    // Update HF orbitals
    std::string hf_valence_states = user_input("Basis/HFOrbitals", "");
    if(!hf_valence_states.empty())
    {
        UpdateHFOrbitals(ConfigurationParser::ParseBasisSize(hf_valence_states), orbitals_to_update);
    }
}

void BasisGenerator::Orthogonalise(pOrbital current) const
{
    pIntegrator integrator(hf->GetIntegrator());
    current->ReNormalise(integrator);

    // Orthogonalise to core
    if(orbitals->core)
    {
        auto it = orbitals->core->begin();
        while(it != orbitals->core->end())
        {
            pOrbitalConst other = it->second;
            if((other->Kappa() == current->Kappa()) && (other->PQN() < current->PQN()))
            {
                double S = integrator->GetInnerProduct(*other, *current);
                (*current) -= (*other) * S;

                current->ReNormalise(integrator);
            }
            it++;
        }
    }

    // Orthogonalise to other excited states.
    if(orbitals->excited)
    {
        auto it = orbitals->excited->begin();
        while(it != orbitals->excited->end())
        {
            pOrbitalConst other = it->second;
            if((other->Kappa() == current->Kappa()) && (other->PQN() < current->PQN()))
            {
                double S = integrator->GetInnerProduct(*other, *current);
                (*current) -= (*other) * S;

                current->ReNormalise(integrator);
            }
            it++;
        }
    }

    current->SetEnergy(hf->GetMatrixElement(*current, *current));
}

void BasisGenerator::Orthogonalise(pOrbital current, pOrbitalMapConst orbitals) const
{
    pIntegrator integrator(hf->GetIntegrator());
    current->ReNormalise(integrator);

    // Orthogonalise to core
    for(auto it: *orbitals)
    {
        pOrbitalConst other = it.second;
        if((other->Kappa() == current->Kappa()) && (other->PQN() < current->PQN()))
        {
            double S = integrator->GetInnerProduct(*other, *current);
            (*current) -= (*other) * S;

            current->ReNormalise(integrator);
        }
    }

    current->SetEnergy(hf->GetMatrixElement(*current, *current));
}

double BasisGenerator::TestOrthogonality(OrbitalInfo& max_i, OrbitalInfo& max_j) const
{
    double max_orth = 0.;
    pIntegrator integrator = hf->GetIntegrator();

    pOrbitalMap all_states;
    if(orbitals->all)
        all_states = orbitals->all;
    else if(open_core)
        all_states = open_core;
    else
        return max_orth;

    auto it = all_states->begin();
    while(it != all_states->end())
    {
        auto jt = all_states->begin();
        while(jt != all_states->end() && (it->first != jt->first))
        {
            if(it->first.Kappa() == jt->first.Kappa())
            {
                double orth = fabs(integrator->GetInnerProduct(*it->second, *jt->second));
                if(orth > max_orth)
                {   max_orth = orth;
                    max_i = it->first;
                    max_j = jt->first;
                }
            }
            jt++;
        }

        it++;
    }

    return max_orth;
}

double BasisGenerator::GetNuclearRMSRadius() const
{
    if(nucleus)
        return nucleus->CalculateNuclearRMSRadius();
    else
        return 0.0;
}
}
