#include "Include.h"

#ifdef AMBIT_USE_OPENMP
    #include <omp.h>
#endif

// Below purposely not included: this file is for a template class and should be included in the header.
// #include "SlaterIntegrals.h"

namespace Ambit
{
template <class MapType>
SlaterIntegrals<MapType>::SlaterIntegrals(pOrbitalManagerConst orbitals, bool two_body_reverse_symmetry_exists):
    SlaterIntegralsInterface(orbitals, two_body_reverse_symmetry_exists)
{
    NumStates = orbitals->size();
    SetUpMap();
}

template <class MapType>
SlaterIntegrals<MapType>::SlaterIntegrals(pOrbitalManagerConst orbitals, pHartreeY hartreeY_op, bool two_body_reverse_symmetry_exists):
    SlaterIntegralsInterface(orbitals, two_body_reverse_symmetry_exists), hartreeY_operator(hartreeY_op)
{
    NumStates = orbitals->size();
    SetUpMap();
}

template <class MapType>
SlaterIntegrals<MapType>::SlaterIntegrals(pOrbitalManagerConst orbitals, pHartreeY hartreeY_op):
    SlaterIntegralsInterface(orbitals, hartreeY_op->ReverseSymmetryExists()), hartreeY_operator(hartreeY_op)
{
    NumStates = orbitals->size();
    SetUpMap();
}

template <class MapType>
unsigned int SlaterIntegrals<MapType>::CalculateTwoElectronIntegrals(pOrbitalMapConst orbital_map_1, pOrbitalMapConst orbital_map_2, pOrbitalMapConst orbital_map_3, pOrbitalMapConst orbital_map_4, bool check_size_only)
{
    // NOTE: For each set of orbitals, we actually calculate
    //       R^k(12,34) = < 4 | Y^k_{31} | 2 >
    // on the assumption that i1 and i2 are smaller.

    unsigned int i1, i2, i3, i4;
    int k;
    pOrbitalConst s1, s2, s3, s4;

    // Get Y^k_{31}
    auto it_1 = orbital_map_1->begin();

    // First, get the number of integrals we need to store. This can be done in parallel, although 
    // there will be some overhead, since we need to keep found_keys sync'd between threads.

#ifdef AMBIT_USE_OPENMP
    // The HartreeY operator is not thread-safe, so make a separate clone for each thread
    std::vector<pHartreeY> hartreeY_operators;
    for(int ii = 0; ii < omp_get_max_threads(); ++ii){
        hartreeY_operators.emplace_back(hartreeY_operator->Clone());
    }
#endif

    std::set<KeyType> found_keys;   // Keep track of valid keys
    hartreeY_operator->SetLightWeightMode(true); // Makes check_sizes faster

    while(it_1 != orbital_map_1->end())
    {
        i1 = orbitals->state_index.at(it_1->first);
        s1 = it_1->second;

        auto it_3 = orbital_map_3->begin();
        if(two_body_reverse_symmetry && orbital_map_1 == orbital_map_3)
        {   it_3 = it_1;
            i3 = i1;
        }

        while(it_3 != orbital_map_3->end())
        {
            i3 = orbitals->state_index.at(it_3->first);
            s3 = it_3->second;

            // Limits on k. This is the expensive part to calculate
#ifdef AMBIT_USE_OPENMP
            k = hartreeY_operators[omp_get_thread_num()]->SetOrbitals(s3, s1);
#else
            k = hartreeY_operator->SetOrbitals(s3, s1);
#endif
            while(k != -1)
            {
                auto it_2 = orbital_map_2->begin();
                while(it_2 != orbital_map_2->end())
                {
                    i2 = orbitals->state_index.at(it_2->first);
                    s2 = it_2->second;

                    auto it_4 = orbital_map_4->begin();
                    while(it_4 != orbital_map_4->end())
                    {
                        i4 = orbitals->state_index.at(it_4->first);
                        s4 = it_4->second;

                        // Check max_pqn conditions and k conditions
                        if(((s1->L() + s2->L() + s3->L() + s4->L())%2 == 0) &&
                           (2 * k >= abs(s2->TwoJ() - s4->TwoJ())) &&
                           (2 * k <= s2->TwoJ() + s4->TwoJ()))
                        {
                            KeyType key = GetKey(k, i1, i2, i3, i4);
                            found_keys.insert(key);
                        }
                        it_4++;
                    }
                    it_2++;
                }
#ifdef AMBIT_USE_OPENMP
                k = hartreeY_operators[omp_get_thread_num()]->NextK();
#else
                k = hartreeY_operator->NextK();
#endif
            } // K loop
            it_3++;
        }
        it_1++;
    }
    hartreeY_operator->SetLightWeightMode(false);

    // If we're running with check-sizes then we're finished, so return the number of keys
    if(check_size_only)
    {   
        return found_keys.size();
    }
    
    // Now create a vector containing all the valid keys, and another to hold the corresponding integrals
    std::vector<KeyType> keys;

    for(auto key : found_keys)
    {
        keys.push_back(key);
    }
    std::vector<double> values(keys.size());

    // TODO: Loop over keys goes here!
#ifdef AMBIT_USE_OPENMP
    #pragma omp parallel for private(k, i1, i2, i3, i4) shared(keys, values, hartreeY_operators, orbitals)
#endif
    for(int ii = 0; ii < keys.size(); ii++)
    {
        auto key = keys[ii];

        // Expand the current key into a set of k and orbital indices i1 to i4
        auto expanded_key = ReverseKey(orbitals->size(), key);
        k = std::get<0>(expanded_key);
        i1 = std::get<1>(expanded_key);
        i2 = std::get<2>(expanded_key);
        i3 = std::get<3>(expanded_key);
        i4 = std::get<4>(expanded_key);

        // Now get pointers to the orbitals corresponding to our indices
        // This is slightly convoluted, but the logic is this: 
        // The reverse_state_index gives us an OrbitalInfo, which we pass to the various orbital_maps
        // to get a key-value pair (since the underlying storage is a std::map), the second element
        // of which is a pointer to an orbital (as required).
        auto s1 = orbital_map_1->find(orbitals->reverse_state_index.at(i1))->second;
        auto s2 = orbital_map_2->find(orbitals->reverse_state_index.at(i2))->second;
        auto s3 = orbital_map_3->find(orbitals->reverse_state_index.at(i3))->second;
        auto s4 = orbital_map_4->find(orbitals->reverse_state_index.at(i4))->second;

        // Now actually calculate the matrix elements
#ifdef AMBIT_USE_OPENMP
        hartreeY_operators[omp_get_thread_num()]->SetOrbitals(s3, s1);
        hartreeY_operators[omp_get_thread_num()]->SetK(k);

        double radial = hartreeY_operators[omp_get_thread_num()]->GetMatrixElement(*s4, *s2);
        values[ii] = radial;
#else
        hartreeY_operator->SetOrbitals(s3, s1);
        hartreeY_operator->SetK(k);
        double radial = hartreeY_operator->GetMatrixElement(s4, s2);
        values[ii] = radial;
#endif

    }

    // Gather all the threads' local copies of the integrals into the global TwoElectronIntegrals
    for(int ii = 0; ii < keys.size(); ii++)
    {
        TwoElectronIntegrals.insert(std::make_pair(keys[ii], values[ii]));
    }
    return TwoElectronIntegrals.size();
}

template <class MapType>
auto SlaterIntegrals<MapType>::GetKey(unsigned int k, unsigned int i1, unsigned int i2, unsigned int i3, unsigned int i4) const -> KeyType
{
    if(two_body_reverse_symmetry)
    {   // Ordering of indices:
        // (i1 <= i3) && (i2 <= i4) && (i1 <= i2) && (if i1 == i2, then (i3 <= i4))
        // therefore (i1 <= i2 <= i4) and (i1 <= i3)
        if(i3 < i1)
            std::swap(i3, i1);
        if(i4 < i2)
            std::swap(i4, i2);
        if(i2 < i1)
        {   std::swap(i2, i1);
            std::swap(i3, i4);
        }
        if((i1 == i2) && (i4 < i3))
            std::swap(i3, i4);
    }
    else
    {   // Ordering of indices:
        // i1 is smallest && (if i1 == i2, then (i3 <= i4))
        //                && (if i1 == i3, then (i2 <= i4))
        //                && (if i1 == i4, then (i2 <= i3))

        // Assert one of i1, i3 is smallest
        if(mmin(i1, i3) > mmin(i2, i4))
        {   std::swap(i1, i2);
            std::swap(i3, i4);
        }
        // Assert i1 <= i3
        if(i1 > i3)
        {   std::swap(i1, i3);
            std::swap(i2, i4);
        }

        if((i1 == i2) && (i4 < i3))
            std::swap(i3, i4);
        if((i1 == i3) && (i4 < i2))
            std::swap(i2, i4);
        if((i1 == i4) && (i3 < i2))
            std::swap(i2, i3);
    }

    KeyType key = k  * NumStates*NumStates*NumStates*NumStates +
                  i1 * NumStates*NumStates*NumStates +
                  i2 * NumStates*NumStates +
                  i3 * NumStates +
                  i4;
    return key;
}

template <class MapType>
auto SlaterIntegrals<MapType>::GetKey(ExpandedKeyType expanded_key) const -> KeyType
{
    return GetKey(std::get<0>(expanded_key), std::get<1>(expanded_key), std::get<2>(expanded_key), std::get<3>(expanded_key), std::get<4>(expanded_key));
}

template <class MapType>
double SlaterIntegrals<MapType>::GetTwoElectronIntegral(unsigned int k, const OrbitalInfo& s1, const OrbitalInfo& s2, const OrbitalInfo& s3, const OrbitalInfo& s4) const
{
    unsigned int i1 = orbitals->state_index.at(s1);
    unsigned int i2 = orbitals->state_index.at(s2);
    unsigned int i3 = orbitals->state_index.at(s3);
    unsigned int i4 = orbitals->state_index.at(s4);

    KeyType key = GetKey(k, i1, i2, i3, i4);
    double radial = 0.;

    auto it = TwoElectronIntegrals.find(key);
    if(it != TwoElectronIntegrals.end())
    {
        radial = it->second;
    }
    else if((s1.L() + s3.L() + k)%2 == 0 && (s2.L() + s4.L() + k)%2 == 0)
    {   // Only print error if requested integral has correct parity rules
#ifdef AMBIT_USE_OPENMP
        #pragma omp critical(ERRSTREAM)
#endif
        *errstream << "SlaterIntegrals::GetTwoElectronIntegral() failed to find integral."
                   << "\n  R^" << k << " ( " << s1.Name() << " " << s2.Name()
                   << ", " << s3.Name() << " " << s4.Name() << "):  key = "
                   << key << "  num_states = " << NumStates << "\n";
    }

    return radial;
}

template <class MapType>
void SlaterIntegrals<MapType>::Read(const std::string& filename)
{
    FILE* fp = file_err_handler->fopen(filename.c_str(), "rb");
    if(!fp)
        return;

    OrbitalIndex old_state_index;
    ReadOrbitalIndexes(old_state_index, fp);

    unsigned int old_key_size;
    file_err_handler->fread(&old_key_size, sizeof(unsigned int), 1, fp);

    switch(old_key_size)
    {
        case sizeof(unsigned long long int):
        {
            unsigned int num_integrals;
            unsigned long long int old_key;
            double value;

            file_err_handler->fread(&num_integrals, sizeof(unsigned int), 1, fp);

            unsigned long long int old_num_states = old_state_index.size();

            for(unsigned int i = 0; i < num_integrals; i++)
            {
                file_err_handler->fread(&old_key, sizeof(unsigned long long int), 1, fp);
                file_err_handler->fread(&value, sizeof(double), 1, fp);

                ExpandedKeyType temp_expanded = ReverseKey(old_num_states, old_key);
                KeyType new_key = GetKey(temp_expanded);

                auto it = TwoElectronIntegrals.find(new_key);
                if(it == TwoElectronIntegrals.end())
                    TwoElectronIntegrals[new_key] = value;
                else
                    it->second += value;
            }
            break;
        }

        case sizeof(unsigned int):
        {
            unsigned int num_integrals;
            unsigned int old_key;
            double value;

            file_err_handler->fread(&num_integrals, sizeof(unsigned int), 1, fp);

            unsigned long long int old_num_states = old_state_index.size();

            for(unsigned int i = 0; i < num_integrals; i++)
            {
                file_err_handler->fread(&old_key, sizeof(unsigned int), 1, fp);
                file_err_handler->fread(&value, sizeof(double), 1, fp);

                ExpandedKeyType temp_expanded = ReverseKey(old_num_states, old_key);
                KeyType new_key = GetKey(temp_expanded);

                auto it = TwoElectronIntegrals.find(new_key);
                if(it == TwoElectronIntegrals.end())
                    TwoElectronIntegrals[new_key] = value;
                else
                    it->second += value;
            }
            break;
        }
    }

    file_err_handler->fclose(fp);
}

template <class MapType>
auto SlaterIntegrals<MapType>::ReverseKey(KeyType num_states, KeyType key) -> ExpandedKeyType
{
    KeyType running_power = num_states * num_states * num_states * num_states;
    KeyType remainder = key;
    ExpandedKeyType expanded_key;

    std::get<0>(expanded_key) = remainder/running_power;
    remainder -= std::get<0>(expanded_key) * running_power;
    running_power = running_power/num_states;

    std::get<1>(expanded_key) = remainder/running_power;
    remainder -= std::get<1>(expanded_key) * running_power;
    running_power = running_power/num_states;

    std::get<2>(expanded_key) = remainder/running_power;
    remainder -= std::get<2>(expanded_key) * running_power;
    running_power = running_power/num_states;

    std::get<3>(expanded_key) = remainder/running_power;
    remainder -= std::get<3>(expanded_key) * running_power;
    running_power = running_power/num_states;

    std::get<4>(expanded_key) = remainder/running_power;
    remainder -= std::get<4>(expanded_key) * running_power;

    return expanded_key;
}

template <class MapType>
void SlaterIntegrals<MapType>::Write(const std::string& filename) const
{
    if(ProcessorRank == 0)
    {
        FILE* fp = file_err_handler->fopen(filename.c_str(), "wb");

        // Write state index
        WriteOrbitalIndexes(orbitals->state_index, fp);

        unsigned int KeyType_size = sizeof(KeyType);
        file_err_handler->fwrite(&KeyType_size, sizeof(unsigned int), 1, fp);

        unsigned int num_integrals = size();
        file_err_handler->fwrite(&num_integrals, sizeof(unsigned int), 1, fp);

        for(auto& pair: TwoElectronIntegrals)
        {
            const double value = pair.second;   // Convert to double
            file_err_handler->fwrite(&pair.first, sizeof(KeyType), 1, fp);
            file_err_handler->fwrite(&value, sizeof(double), 1, fp);
        }

        file_err_handler->fclose(fp);
    }
}
}
