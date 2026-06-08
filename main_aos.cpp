#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include <iomanip>
#include <chrono>
#include <cstring>

#ifdef _DOUBLE_
#define FLOATTYPE double
#else
#define FLOATTYPE float
#endif

struct ConservativeVars {
  FLOATTYPE rho, rho_u, rho_E;
};

int main(int narg, char **argv) {
  int numPoints = 80;
  if(narg != 2) {
    std::cout << "Usage: " << argv[0] << " <numPoints>" << std::endl;
    return 1;
  }
  numPoints = std::stoi(argv[1]);

  ConservativeVars *state = new ConservativeVars[numPoints];
  ConservativeVars *state_init = new ConservativeVars[numPoints];
  ConservativeVars *flux = new ConservativeVars[numPoints];
  ConservativeVars *RHS = new ConservativeVars[numPoints];
  ConservativeVars *Ui = new ConservativeVars[numPoints];
  FLOATTYPE *xj = new FLOATTYPE[numPoints];

  FLOATTYPE gamma = 1.4, rho0 = 1.0, u0 = 0.0, p0 = 1.0;
  
  for(int j = 0; j < numPoints; j++) {
    xj[j] = FLOATTYPE(j) / FLOATTYPE(numPoints - 1);
    FLOATTYPE pert = 0.1 * sin(2.0 * M_PI * xj[j]);
    state[j].rho = rho0 + pert;
    state[j].rho_u = state[j].rho * u0;
    FLOATTYPE E = p0 / ((gamma - 1.0) * state[j].rho) + 0.5 * u0 * u0;
    state[j].rho_E = state[j].rho * E;
    state_init[j] = state[j];
  }

  FLOATTYPE CFL = 0.5;
  FLOATTYPE dx = xj[1] - xj[0];
  FLOATTYPE dt = CFL * dx;
  FLOATTYPE t_final = 0.1;
  FLOATTYPE time = 0.;

  auto startTime = std::chrono::high_resolution_clock::now();

  while(time < t_final) {
    if(time + dt >= t_final) dt = t_final - time;
    
    for(int j = 0; j < numPoints; j++) {
      FLOATTYPE u = state[j].rho_u / state[j].rho;
      FLOATTYPE E = state[j].rho_E / state[j].rho;
      FLOATTYPE p = (gamma - 1.0) * state[j].rho * (E - 0.5 * u * u);
      flux[j].rho = state[j].rho_u;
      flux[j].rho_u = state[j].rho_u * u + p;
      flux[j].rho_E = u * (state[j].rho_E + p);
    }
    
    for(int j = 0; j < numPoints; j++) {
      int jp = (j + 1) % numPoints;
      int jm = (j - 1 + numPoints) % numPoints;
      FLOATTYPE dxinv = 1.0 / (xj[jp] - xj[jm]);
      RHS[j].rho = -(flux[jp].rho - flux[jm].rho) * dxinv;
      RHS[j].rho_u = -(flux[jp].rho_u - flux[jm].rho_u) * dxinv;
      RHS[j].rho_E = -(flux[jp].rho_E - flux[jm].rho_E) * dxinv;
    }
    
    for(int j = 0; j < numPoints; j++) {
      state[j].rho += dt / 6.0 * RHS[j].rho;
      state[j].rho_u += dt / 6.0 * RHS[j].rho_u;
      state[j].rho_E += dt / 6.0 * RHS[j].rho_E;
    }
    
    time += dt;
  }

  auto endTime = std::chrono::high_resolution_clock::now();
  double compTime = std::chrono::duration<double>(endTime - startTime).count();

  FLOATTYPE err = 0.;
  for(int n = 0; n < numPoints; n++) {
    FLOATTYPE diff = state[n].rho - state_init[n].rho;
    err += diff * diff;
  }
  err = sqrt(err);

  std::cout << std::setprecision(4) << "Comp. time: " << compTime;
  std::cout << " sec. Error: " << err << std::endl;

  delete[] state;
  delete[] state_init;
  delete[] flux;
  delete[] RHS;
  delete[] Ui;
  delete[] xj;

  return 0;
}
