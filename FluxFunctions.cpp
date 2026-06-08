#include "FluxFunctions.h"

template<class T>
FluxFunction<T>::FluxFunction()
{
  
};

template<class T>
LinearFlux<T>::LinearFlux()
{
  c = 1.;
};

template<class T>
void LinearFlux<T>::computeFlux(DataStruct<T> &U, DataStruct<T> &F)
{
  T *dataU = U.getData();
  T *dataF = F.getData();

  for(int n = 0; n < U.getSize(); n++)
  {
    dataF[n] = c * dataU[n];
  };
};

template<class T>
T LinearFlux<T>::computeFlux(const T &Ui)
{
  return c * Ui;
};

// ============================================================
// EULER FLUX IMPLEMENTATION
// ============================================================

template<class T>
EulerFlux<T>::EulerFlux()
{
  gamma = 1.4;  // air
};

template<class T>
void EulerFlux<T>::computeFluxAtPoint(const EulerState<T> &state, EulerState<T> &flux)
{
  // Compute velocity and pressure
  T u = state.rho_u / state.rho;
  T p = state.computePressure(gamma);
  T E = state.rho_E / state.rho;
  
  // Flux: f = [rho*u, rho*u^2 + p, u*(rho*E + p)]
  flux.rho = state.rho_u;
  flux.rho_u = state.rho_u * u + p;
  flux.rho_E = u * (state.rho_E + p);
}

template<class T>
void EulerFlux<T>::computeFlux(DataStruct<T> &rho, DataStruct<T> &rho_u, DataStruct<T> &rho_E,
                                DataStruct<T> &f_rho, DataStruct<T> &f_rho_u, DataStruct<T> &f_rho_E)
{
  T *data_rho = rho.getData();
  T *data_rho_u = rho_u.getData();
  T *data_rho_E = rho_E.getData();
  
  T *data_f_rho = f_rho.getData();
  T *data_f_rho_u = f_rho_u.getData();
  T *data_f_rho_E = f_rho_E.getData();
  
  for(int n = 0; n < rho.getSize(); n++)
  {
    EulerState<T> state;
    state.rho = data_rho[n];
    state.rho_u = data_rho_u[n];
    state.rho_E = data_rho_E[n];
    
    EulerState<T> flux;
    computeFluxAtPoint(state, flux);
    
    data_f_rho[n] = flux.rho;
    data_f_rho_u[n] = flux.rho_u;
    data_f_rho_E[n] = flux.rho_E;
  }
}

template class FluxFunction<float>;
template class FluxFunction<double>;

template class LinearFlux<float>;
template class LinearFlux<double>;

template class EulerFlux<float>;
template class EulerFlux<double>;