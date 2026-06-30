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
  ghostLeft = T(0.);
  ghostRight = T(0.);
  dx = mesh.getSize() > 1 ? mesh.getData()[1] - mesh.getData()[0] : T(1.);
}

template<class T>
Central1D<T>::Central1D(DataStruct<T> &_U, 
                     DataStruct<T> &_mesh, 
                     FluxFunction<T> &_F,
                     T _dx):
U(_U), mesh(_mesh), F(_F)
{
  RHS.setSize(_U.getSize());
  ghostLeft = T(0.);
  ghostRight = T(0.);
  dx = _dx;
}

template<class T>
Central1D<T>::~Central1D()
{

}

template<class T>
void Central1D<T>::setGhostValues(T left, T right)
{
  ghostLeft = left;
  ghostRight = right;
}

template<class T>
void Central1D<T>::evalRHS(DataStruct<T> &Uin)
{
  T *dataRHS = RHS.getData();
  const T *dataU = Uin.getData();
  const int len = U.getSize();
  const T oneDivTwoDx = T(1.) / (T(2.) * dx);
  
  for(int j = 0; j < len; j++)
  {
    const T uLeft = j == 0 ? ghostLeft : dataU[j-1];
    const T uRight = j == len-1 ? ghostRight : dataU[j+1];
    dataRHS[j] = -(F.computeFlux(uRight) - F.computeFlux(uLeft))*oneDivTwoDx;
  }
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
