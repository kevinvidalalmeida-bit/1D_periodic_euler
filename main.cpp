#include <iostream> // std::cout
#include <fstream> // ofile
#include <string>
#include <math.h>
#include <iomanip> // set precision
#include <mpi.h>

#include "DataStructs.h"
#include "rk4.h"
#include "FluxFunctions.h"
#include "RHSoperator.h"

#ifdef _DOUBLE_
#define FLOATTYPE double
#else
#define FLOATTYPE float
#endif

// declare supporting functions
void write2FileEuler(DataStruct<FLOATTYPE> &X, DataStruct<FLOATTYPE> &rho, 
                     DataStruct<FLOATTYPE> &rho_u, DataStruct<FLOATTYPE> &rho_E, std::string name);
FLOATTYPE calcL2normEuler(DataStruct<FLOATTYPE> &rho, DataStruct<FLOATTYPE> &rho_init);


int main(int narg, char **argv)
{
  int numPoints =  80;
  FLOATTYPE dummy = 0.; // not used for Euler

  if(narg != 2)
  {
    std::cout<< "Wrong number of arguments. You should include:" << std::endl;
    std::cout<< "    Num points" << std::endl;
    return 1;
  }else
  {
    numPoints = std::stoi(argv[1]);
  }

  // Solution data: rho, rho*u, rho*E
  DataStruct<FLOATTYPE> rho(numPoints), rho_u(numPoints), rho_E(numPoints), xj(numPoints);

  // Flux function for Euler 1D
  EulerFlux<FLOATTYPE> euler_flux;

  // Time solver for Euler
  RungeKutta4Euler<FLOATTYPE> rk(rho, rho_u, rho_E);

  // Initial Condition: density profile with small perturbation
  FLOATTYPE *datax = xj.getData();
  FLOATTYPE *data_rho = rho.getData();
  FLOATTYPE *data_rho_u = rho_u.getData();
  FLOATTYPE *data_rho_E = rho_E.getData();
  
  FLOATTYPE gamma = 1.4;
  FLOATTYPE rho0 = 1.0;
  FLOATTYPE u0 = 0.0;
  FLOATTYPE p0 = 1.0;
  
  for(int j = 0; j < numPoints; j++)
  {
    // Position
    datax[j] = FLOATTYPE(j)/FLOATTYPE(numPoints-1);
    
    // Add small density perturbation
    FLOATTYPE pert = 0.1 * sin(2.0 * M_PI * datax[j]);
    data_rho[j] = rho0 + pert;
    
    // Momentum
    data_rho_u[j] = data_rho[j] * u0;
    
    // Total energy: E = p / ((gamma-1)*rho) + u^2/2
    FLOATTYPE E = p0 / ((gamma - 1.0) * data_rho[j]) + 0.5 * u0 * u0;
    data_rho_E[j] = data_rho[j] * E;
  }

  // Store initial condition
  DataStruct<FLOATTYPE> rho_init;
  rho_init = rho;

  // RHS Operator
  Central1DEuler<FLOATTYPE> rhs(rho, rho_u, rho_E, xj, euler_flux);

  // CFL condition (need to be conservative for nonlinear systems)
  FLOATTYPE CFL = 0.5;
  FLOATTYPE dx = datax[1] - datax[0];
  FLOATTYPE dt = CFL * dx;

  // Output Initial Condition
  write2FileEuler(xj, rho, rho_u, rho_E, "initialCondition.csv");

  FLOATTYPE t_final = 0.1;
  FLOATTYPE time = 0.;
  
  // Temporary storage for intermediate Ui
  DataStruct<FLOATTYPE> Ui_rho(numPoints), Ui_rho_u(numPoints), Ui_rho_E(numPoints);

  // init timer
  double compTime = MPI_Wtime();

  // main loop
  while(time < t_final)
  {
    if(time + dt >= t_final) dt = t_final - time;

    // take RK step
    rk.initRK();
    for(int s = 0; s < rk.getNumSteps(); s++)
    {
      rk.stepUi(dt);
      rk.currentU(Ui_rho, Ui_rho_u, Ui_rho_E);
      
      // Evaluate RHS at intermediate values
      rhs.eval(Ui_rho, Ui_rho_u, Ui_rho_E);
      
      // Set Fi
      rk.setFi(rhs.ref2RHS_rho(), rhs.ref2RHS_rho_u(), rhs.ref2RHS_rho_E());
    }
    rk.finalizeRK(dt);
    time += dt;
  }

  // finish timer
  compTime = MPI_Wtime() - compTime;

  write2FileEuler(xj, rho, rho_u, rho_E, "final.csv");

  // L2 norm error (on density)
  FLOATTYPE err = calcL2normEuler(rho, rho_init);
  std::cout << std::setprecision(4) << "Comp. time: " << compTime;
  std::cout << " sec. Error: " << err;
  std::cout << std::endl;

  return 0;
}


// ==================================================================
// AUXILIARY FUNCTIONS
// ==================================================================
void write2FileEuler(DataStruct<FLOATTYPE> &X, DataStruct<FLOATTYPE> &rho, 
                     DataStruct<FLOATTYPE> &rho_u, DataStruct<FLOATTYPE> &rho_E, std::string name)
{
  std::ofstream file;
  file.open(name, std::ios_base::trunc);
  if(!file.is_open()) 
  {
    std::cout << "Couldn't open file: " << name << std::endl;
    exit(1);
  }
  
  const FLOATTYPE *datax = X.getData();
  const FLOATTYPE *data_rho = rho.getData();
  const FLOATTYPE *data_rho_u = rho_u.getData();
  const FLOATTYPE *data_rho_E = rho_E.getData();
  
  FLOATTYPE gamma = 1.4;
  
  file << "x,rho,u,p,E" << std::endl;
  
  for(int j = 0; j < rho.getSize(); j++)
  {
    FLOATTYPE u = data_rho_u[j] / data_rho[j];
    FLOATTYPE E = data_rho_E[j] / data_rho[j];
    FLOATTYPE p = (gamma - 1.0) * data_rho[j] * (E - 0.5 * u * u);
    
    file << datax[j] << "," << data_rho[j] << "," << u << "," << p << "," << E << std::endl;
  }

  file.close();
}

FLOATTYPE calcL2normEuler(DataStruct<FLOATTYPE> &rho, DataStruct<FLOATTYPE> &rho_init)
{
  FLOATTYPE err = 0.;
  const FLOATTYPE *data_rho = rho.getData();
  const FLOATTYPE *data_init = rho_init.getData();

  for(int n = 0; n < rho.getSize(); n++)
  {
    err += (data_rho[n] - data_init[n])*(data_rho[n] - data_init[n]);
  }

  return sqrt( err );
}