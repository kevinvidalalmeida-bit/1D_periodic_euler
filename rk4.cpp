#include "rk4.h"
#include <cassert>

template<class T>
RungeKutta4<T>::RungeKutta4(DataStruct<T> &_Un):
Un(_Un)
{
  nSteps = 4;
  currentStep = 0;

  coeffsA = new T[4];
  coeffsB = new T[4];
  coeffsA[0] = 0.;
  coeffsA[1] = 0.5;
  coeffsA[2] = 0.5;
  coeffsA[3] = 1.;
  coeffsB[0] = 1.;
  coeffsB[1] = 2.;
  coeffsB[2] = 2.;
  coeffsB[3] = 1.;

  fi = new DataStruct<T>[nSteps];

  fi[0].setSize(Un.getSize());
  fi[1].setSize(Un.getSize());
  fi[2].setSize(Un.getSize());
  fi[3].setSize(Un.getSize());

  Ui.setSize(Un.getSize());
};

template<class T>
RungeKutta4<T>::~RungeKutta4()
{
  delete[] fi;
  delete[] coeffsA;
  delete[] coeffsB;
};

template<class T>
int RungeKutta4<T>::getNumSteps()
{
  return nSteps;
};

template<class T>
void RungeKutta4<T>::initRK()
{
  currentStep = 0;
};

template<class T>
void RungeKutta4<T>::stepUi(T dt)
{
  assert(currentStep < nSteps);

  if(currentStep == 0)
  {
    T *dataUi = Ui.getData();
    const T *dataU  = Un.getData();

    for(int n = 0; n < Ui.getSize(); n++)
    {
      dataUi[n] = dataU[n];
    }
  }
  else
  {
    T *datafi = fi[currentStep-1].getData();
    T *dataUi = Ui.getData();
    const T *dataU  = Un.getData();

    for(int n = 0; n < Ui.getSize(); n++)
    {
      dataUi[n] = dataU[n] + coeffsA[currentStep]*dt* datafi[n];
    }
  }
};

template<class T>
void RungeKutta4<T>::finalizeRK(const T dt)
{
  T *dataUn = Un.getData();
  T *dataUi = Ui.getData();

  // set Ui to 0
  for(int n = 0; n < Ui.getSize(); n++)
  {
    dataUi[n] = 0.;
  }
  
  for(int s = 0; s < nSteps; s++)
  {
    const T *dataFi = fi[s].getData();
    const T b = coeffsB[s];

    for(int n = 0; n < Ui.getSize(); n++)
    {
      dataUi[n] += b * dataFi[n];
    }
  }

  const T oneDiv6 = 1. / 6.;
  for(int n = 0; n < Ui.getSize(); n++)
  {
    dataUn[n] += dt * oneDiv6 * dataUi[n];
  }
};

template<class T>
void RungeKutta4<T>::setFi(DataStruct<T> &_F)
{
  T *dataFi = fi[currentStep].getData();
  const T *dataF  = _F.getData();

  for(int n = 0; n < Ui.getSize(); n++)
  {
    dataFi[n] = dataF[n];
  }

  currentStep++;
};

template<class T>
DataStruct<T> * RungeKutta4<T>::currentU()
{
  return &Ui;
};


template class RungeKutta4<float>;
template class RungeKutta4<double>;

// ============================================================
// RUNGE KUTTA 4 FOR EULER 1D IMPLEMENTATION
// ============================================================

template<class T>
RungeKutta4Euler<T>::RungeKutta4Euler(DataStruct<T> &_Un_rho, DataStruct<T> &_Un_rho_u, DataStruct<T> &_Un_rho_E):
Un_rho(_Un_rho), Un_rho_u(_Un_rho_u), Un_rho_E(_Un_rho_E)
{
  nSteps = 4;
  currentStep = 0;

  coeffsA = new T[4];
  coeffsB = new T[4];
  coeffsA[0] = 0.;
  coeffsA[1] = 0.5;
  coeffsA[2] = 0.5;
  coeffsA[3] = 1.;
  coeffsB[0] = 1.;
  coeffsB[1] = 2.;
  coeffsB[2] = 2.;
  coeffsB[3] = 1.;

  int size = _Un_rho.getSize();
  
  fi_rho = new DataStruct<T>[nSteps];
  fi_rho_u = new DataStruct<T>[nSteps];
  fi_rho_E = new DataStruct<T>[nSteps];

  for(int i = 0; i < nSteps; i++)
  {
    fi_rho[i].setSize(size);
    fi_rho_u[i].setSize(size);
    fi_rho_E[i].setSize(size);
  }

  Ui_rho.setSize(size);
  Ui_rho_u.setSize(size);
  Ui_rho_E.setSize(size);
  
  RHS_rho.setSize(size);
  RHS_rho_u.setSize(size);
  RHS_rho_E.setSize(size);
};

template<class T>
RungeKutta4Euler<T>::~RungeKutta4Euler()
{
  delete[] fi_rho;
  delete[] fi_rho_u;
  delete[] fi_rho_E;
  delete[] coeffsA;
  delete[] coeffsB;
};

template<class T>
int RungeKutta4Euler<T>::getNumSteps()
{
  return nSteps;
};

template<class T>
void RungeKutta4Euler<T>::initRK()
{
  currentStep = 0;
};

template<class T>
void RungeKutta4Euler<T>::stepUi(T dt)
{
  int size = Un_rho.getSize();
  
  if(currentStep == 0)
  {
    T *dataUi_rho = Ui_rho.getData();
    T *dataUi_rho_u = Ui_rho_u.getData();
    T *dataUi_rho_E = Ui_rho_E.getData();
    
    const T *dataU_rho = Un_rho.getData();
    const T *dataU_rho_u = Un_rho_u.getData();
    const T *dataU_rho_E = Un_rho_E.getData();

    for(int n = 0; n < size; n++)
    {
      dataUi_rho[n] = dataU_rho[n];
      dataUi_rho_u[n] = dataU_rho_u[n];
      dataUi_rho_E[n] = dataU_rho_E[n];
    }
  }
  else
  {
    T *dataUi_rho = Ui_rho.getData();
    T *dataUi_rho_u = Ui_rho_u.getData();
    T *dataUi_rho_E = Ui_rho_E.getData();
    
    const T *dataU_rho = Un_rho.getData();
    const T *dataU_rho_u = Un_rho_u.getData();
    const T *dataU_rho_E = Un_rho_E.getData();
    
    const T *datafi_rho = fi_rho[currentStep-1].getData();
    const T *datafi_rho_u = fi_rho_u[currentStep-1].getData();
    const T *datafi_rho_E = fi_rho_E[currentStep-1].getData();

    for(int n = 0; n < size; n++)
    {
      dataUi_rho[n] = dataU_rho[n] + coeffsA[currentStep]*dt*datafi_rho[n];
      dataUi_rho_u[n] = dataU_rho_u[n] + coeffsA[currentStep]*dt*datafi_rho_u[n];
      dataUi_rho_E[n] = dataU_rho_E[n] + coeffsA[currentStep]*dt*datafi_rho_E[n];
    }
  }
};

template<class T>
void RungeKutta4Euler<T>::finalizeRK(const T dt)
{
  int size = Un_rho.getSize();
  T *dataUn_rho = Un_rho.getData();
  T *dataUn_rho_u = Un_rho_u.getData();
  T *dataUn_rho_E = Un_rho_E.getData();
  
  T *dataRHS_rho = RHS_rho.getData();
  T *dataRHS_rho_u = RHS_rho_u.getData();
  T *dataRHS_rho_E = RHS_rho_E.getData();

  // Initialize RHS to 0
  for(int n = 0; n < size; n++)
  {
    dataRHS_rho[n] = 0.;
    dataRHS_rho_u[n] = 0.;
    dataRHS_rho_E[n] = 0.;
  }
  
  // Accumulate RK stages
  for(int s = 0; s < nSteps; s++)
  {
    const T *dataFi_rho = fi_rho[s].getData();
    const T *dataFi_rho_u = fi_rho_u[s].getData();
    const T *dataFi_rho_E = fi_rho_E[s].getData();
    const T b = coeffsB[s];

    for(int n = 0; n < size; n++)
    {
      dataRHS_rho[n] += b * dataFi_rho[n];
      dataRHS_rho_u[n] += b * dataFi_rho_u[n];
      dataRHS_rho_E[n] += b * dataFi_rho_E[n];
    }
  }

  // Update solution
  const T oneDiv6 = 1. / 6.;
  for(int n = 0; n < size; n++)
  {
    dataUn_rho[n] += dt * oneDiv6 * dataRHS_rho[n];
    dataUn_rho_u[n] += dt * oneDiv6 * dataRHS_rho_u[n];
    dataUn_rho_E[n] += dt * oneDiv6 * dataRHS_rho_E[n];
  }
};

template<class T>
void RungeKutta4Euler<T>::setFi(DataStruct<T> &_F_rho, DataStruct<T> &_F_rho_u, DataStruct<T> &_F_rho_E)
{
  int size = Ui_rho.getSize();
  T *dataFi_rho = fi_rho[currentStep].getData();
  T *dataFi_rho_u = fi_rho_u[currentStep].getData();
  T *dataFi_rho_E = fi_rho_E[currentStep].getData();
  
  const T *dataF_rho = _F_rho.getData();
  const T *dataF_rho_u = _F_rho_u.getData();
  const T *dataF_rho_E = _F_rho_E.getData();

  for(int n = 0; n < size; n++)
  {
    dataFi_rho[n] = dataF_rho[n];
    dataFi_rho_u[n] = dataF_rho_u[n];
    dataFi_rho_E[n] = dataF_rho_E[n];
  }

  currentStep++;
};

template<class T>
void RungeKutta4Euler<T>::currentU(DataStruct<T> &out_rho, DataStruct<T> &out_rho_u, DataStruct<T> &out_rho_E)
{
  out_rho = Ui_rho;
  out_rho_u = Ui_rho_u;
  out_rho_E = Ui_rho_E;
};

template class RungeKutta4Euler<float>;
template class RungeKutta4Euler<double>;