#ifndef LOCAL_POTENTIAL_DECORATOR
#define LOCAL_POTENTIAL_DECORATOR

#include "HFOperator.h"

/** Generic (abstract-ish) class to add an extra local potential to a HF operator.
    Local member extraLocalPotential must be set by subclasses via SetCore() and/or SetODEParameters().
 */
class LocalPotentialDecorator : public HFOperatorDecorator
{
public:
    LocalPotentialDecorator(pHFOperator wrapped_hf, pOPIntegrator integration_strategy = pOPIntegrator());

public:
    virtual void GetODEFunction(unsigned int latticepoint, const SpinorFunction& fg, double* w) const;
    virtual void GetODECoefficients(unsigned int latticepoint, const SpinorFunction& fg, double* w_f, double* w_g, double* w_const) const;
    virtual void GetODEJacobian(unsigned int latticepoint, const SpinorFunction& fg, double** jacobian, double* dwdr) const;

public:
    virtual SpinorFunction ApplyTo(const SpinorFunction& a) const;

protected:
    RadialFunction extraLocalPotential;
};

/** Decorate HF operator with a local approximation to the exchange potential. Good for first approximations. */
class LocalExchangeApproximation : public LocalPotentialDecorator
{
public:
    LocalExchangeApproximation(pHFOperator wrapped_hf, pOPIntegrator integration_strategy = pOPIntegrator());

    virtual void SetCore(const Core* hf_core);
};

#endif
