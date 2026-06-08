#ifndef _FLUX_FUNCTIONS
#define _FLUX_FUNCTIONS

#include "DataStructs.h"

// base class for flux functions
template<class T>
class FluxFunction
{
  public:
    FluxFunction();

    // in U; out F
    virtual void computeFlux(DataStruct<T> &U, DataStruct<T> &F) = 0;

    // same as above but at the node level
    virtual T computeFlux(const T &Ui) = 0;
};

template<class T>
class LinearFlux : public FluxFunction<T>
{
  private:
    T c;

  public:
    LinearFlux();

    // in U; out F
    virtual void computeFlux(DataStruct<T> &U, DataStruct<T> &F);

    // same as above but at the node level
    virtual T computeFlux(const T &Ui);
};

// ============================================================
// EULER 1D FLUX FUNCTIONS
// ============================================================

// Structure for Euler state vector: (rho, rho*u, rho*E)
template<class T>
struct EulerState
{
  T rho;      // density
  T rho_u;    // momentum (rho * u)
  T rho_E;    // total energy (rho * E)
  
  // Equation of state: p = (gamma - 1) * rho * (E - u^2/2)
  // where gamma = 1.4 for air
  T computePressure(T gamma = 1.4) const
  {
    T u = rho_u / rho;
    T E = rho_E / rho;
    return (gamma - 1.0) * rho * (E - 0.5 * u * u);
  }
};

// Euler flux function for 1D
template<class T>
class EulerFlux : public FluxFunction<T>
{
  private:
    T gamma;  // heat capacity ratio (default 1.4 for air)

  public:
    EulerFlux();
    
    // These methods are not used for Euler (uses the vector versions below)
    virtual void computeFlux(DataStruct<T> &U, DataStruct<T> &F) {}
    virtual T computeFlux(const T &Ui) { return 0; }

    // Compute flux for all points
    // Input: rho, rho_u, rho_E (density, momentum, total energy)
    // Output: f_rho, f_rho_u, f_rho_E
    void computeFlux(DataStruct<T> &rho, DataStruct<T> &rho_u, DataStruct<T> &rho_E,
                     DataStruct<T> &f_rho, DataStruct<T> &f_rho_u, DataStruct<T> &f_rho_E);
    
    // Compute flux at a single point
    void computeFluxAtPoint(const EulerState<T> &state, EulerState<T> &flux);
};

#endif // _FLUX_FUNCTIONS