#ifndef EXCITED_STATES_H
#define EXCITED_STATES_H

#include "HartreeFock/StateManager.h"
#include "HartreeFock/StateIterator.h"
#include "HartreeFock/Core.h"
#include "HartreeFock/SigmaPotential.h"
#include <set>
#include <map>
#include <vector>

typedef std::map<int, SigmaPotential*> SigmaMap;
typedef std::map<StateInfo, double> SigmaAmount;

class ExcitedStates : public StateManager
{
public:     // Methods for StateManager role
    ExcitedStates(Lattice* lattice, const Core* atom_core);
    virtual ~ExcitedStates();

    virtual void AddState(DiscreteState* s);

    /** Get a copy of the state, including iterating sigma (if available). */
    virtual DiscreteState GetStateWithSigma(const StateInfo& info) const;

public:     // Methods for managing excited single particle basis states
    /** Create virtual states above the core.
        num_states_per_l is simply a vector to say how many states
        should be included above the core for each angular momentum, L.
        eg: {3, 2, 1} would make 3 s-wave, 2 p-wave, 1 d-wave, etc.
        If a state already exists, then this routine just skips it.
     */
    virtual void CreateExcitedStates(const std::vector<unsigned int>& num_states_per_l) = 0;

    /** Update all of the excited states because the core has changed. */
    virtual void Update() = 0;

public:     // Methods for managing sigma
    void SetIdentifier(const std::string* s) { identifier = s; }

    /** Create a sigma operator for the given state to second order.
        This operator will henceforth be included in the exchange potential of
        the state (to turn this off set SigmaAmount to zero).
        Return value is the energy of the state including second order matrix elements.
     */
    double GetSecondOrderSigma(const StateInfo& info);

    /** Retreive a sigma operator that already exists. Return true if successful. */
    bool RetrieveSecondOrderSigma(const StateInfo& info);

    /** Create sigma operator if one doesn't already exist.
    	Set SigmaAmount to give correct energy (in Hartree units).
     */
    void SetEnergyViaSigma(const StateInfo& info, double energy);

    /** Get/Set the amount of sigma to be used when calculating state for which sigma exists. */
    void SetSigmaAmount(const StateInfo& info, double amount);
    double GetSigmaAmount(const StateInfo& info) const;

    const Core* GetCore() const { return core; }

protected:
    /** current.f = r * previous.f
        PRE: previous.kappa = current.kappa
     */
    void MultiplyByR(const DiscreteState* previous, DiscreteState* current) const;

    /** current.f = sin(kr) * previous.f, where k = Pi/R_max
        PRE: previous.kappa = current.kappa
     */
    void MultiplyBySinR(const DiscreteState* previous, DiscreteState* current) const;

    /** current.f = r * sin(kr) * previous.f
        PRE: current.l = previous.l + 1
     */
    void MultiplyByRSinR(const DiscreteState* previous, DiscreteState* current) const;

    /** Orthogonalise to all states that have the same kappa and principal quantum number
        less than current (both in core and in excited states).
        PRE: all states with smaller pqn must already be orthogonal
     */
    void Orthogonalise(DiscreteState* current) const;

protected:
    const Core* core;

    const std::string* identifier;  // needed to store sigmas
    SigmaMap SecondOrderSigma;
    SigmaAmount SecondOrderAmount;
};

#endif
