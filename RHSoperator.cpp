#include "RHSoperator.h"

template<class T>
RHSOperator<T>::RHSOperator()
{

}


template<class T>
RHSOperator<T>::~RHSOperator()
{

}

template<class T>
Central1D<T>::Central1D(DataStruct<T> &_U, 
                     DataStruct<T> &_mesh, 
                     FluxFunction<T> &_F):
U(_U), mesh(_mesh), F(_F)
{
  RHS.setSize(_U.getSize());
}

template<class T>
Central1D<T>::~Central1D()
{

}

template<class T>
void Central1D<T>::evalRHS(DataStruct<T> &Uin)
{
  // the BC should be included in the mesh
  // momentarily done here by hand
  T *dataRHS = RHS.getData();
  const T *dataU = Uin.getData();
  const T *dataMesh = mesh.getData();
  const int len = U.getSize();
  
  for(int j = 0; j < len; j++)
  {
    T dx;
    if(j == 0)
    {
      dx = dataMesh[len-1] - dataMesh[len-2];
      dx += dataMesh[1] - dataMesh[0];
      dataRHS[0] = -(F.computeFlux(dataU[1]) - F.computeFlux(dataU[len-2]))/dx;
    }
    else
    {
      dx = dataMesh[j+1] - dataMesh[j-1];
      dataRHS[j] = -(F.computeFlux(dataU[j+1]) - F.computeFlux(dataU[j-1]))/dx;
    }
  }

  dataRHS[len-1] = dataRHS[0];
}

template<class T>
void Central1D<T>::eval()
{
  evalRHS(U);
}

template<class T>
void Central1D<T>::eval(DataStruct<T> &Uin)
{
  evalRHS(Uin);
}

template<class T>
DataStruct<T>& Central1D<T>::ref2RHS()
{
  return RHS;
}


template class RHSOperator<float>;
template class RHSOperator<double>;

template class Central1D<float>;
template class Central1D<double>;

// ============================================================
// EULER 1D RHS OPERATOR IMPLEMENTATION
// ============================================================

template<class T>
Central1DEuler<T>::Central1DEuler(DataStruct<T> &_rho, DataStruct<T> &_rho_u, DataStruct<T> &_rho_E,
                                   DataStruct<T> &_mesh, EulerFlux<T> &_F):
rho(_rho), rho_u(_rho_u), rho_E(_rho_E), mesh(_mesh), F(_F)
{
  int n = _rho.getSize();
  RHS_rho.setSize(n);
  RHS_rho_u.setSize(n);
  RHS_rho_E.setSize(n);
  
  f_rho.setSize(n);
  f_rho_u.setSize(n);
  f_rho_E.setSize(n);
}

template<class T>
Central1DEuler<T>::~Central1DEuler()
{
}

template<class T>
void Central1DEuler<T>::evalRHS(DataStruct<T> &rho_in, DataStruct<T> &rho_u_in, DataStruct<T> &rho_E_in)
{
  // Compute fluxes
  F.computeFlux(rho_in, rho_u_in, rho_E_in, f_rho, f_rho_u, f_rho_E);
  
  T *dataRHS_rho = RHS_rho.getData();
  T *dataRHS_rho_u = RHS_rho_u.getData();
  T *dataRHS_rho_E = RHS_rho_E.getData();
  
  const T *dataF_rho = f_rho.getData();
  const T *dataF_rho_u = f_rho_u.getData();
  const T *dataF_rho_E = f_rho_E.getData();
  
  const T *dataMesh = mesh.getData();
  const int len = rho_in.getSize();
  
  // Apply central differences with periodic boundary conditions
  for(int j = 0; j < len; j++)
  {
    T dx;
    int j_plus = (j + 1) % len;
    int j_minus = (j - 1 + len) % len;
    
    dx = dataMesh[j_plus] - dataMesh[j_minus];
    
    dataRHS_rho[j] = -(dataF_rho[j_plus] - dataF_rho[j_minus]) / dx;
    dataRHS_rho_u[j] = -(dataF_rho_u[j_plus] - dataF_rho_u[j_minus]) / dx;
    dataRHS_rho_E[j] = -(dataF_rho_E[j_plus] - dataF_rho_E[j_minus]) / dx;
  }
}

template<class T>
void Central1DEuler<T>::eval()
{
  evalRHS(rho, rho_u, rho_E);
}

template<class T>
void Central1DEuler<T>::eval(DataStruct<T> &Uin)
{
  // dummy implementation
}

template<class T>
void Central1DEuler<T>::eval(DataStruct<T> &rho_in, DataStruct<T> &rho_u_in, DataStruct<T> &rho_E_in)
{
  evalRHS(rho_in, rho_u_in, rho_E_in);
}

template class Central1DEuler<float>;
template class Central1DEuler<double>;

template class Central1D<float>;
template class Central1D<double>;