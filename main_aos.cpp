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

void computeRHS(const ConservativeVars *state, ConservativeVars *rhs,
                ConservativeVars *flux, const FLOATTYPE *xj,
                int numPoints, FLOATTYPE gamma);

int main(int narg, char **argv) {
  int numPoints = 80;
  if(narg != 2) {
    std::cout << "Usage: " << argv[0] << " <numPoints>" << std::endl;
    return 1;
  }
  numPoints = std::stoi(argv[1]);
  if(numPoints < 3) {
    std::cout << "numPoints must be at least 3" << std::endl;
    return 1;
  }

  ConservativeVars *state = new ConservativeVars[numPoints];
  ConservativeVars *state_init = new ConservativeVars[numPoints];
  ConservativeVars *flux = new ConservativeVars[numPoints];
  ConservativeVars *prev_rhs = new ConservativeVars[numPoints];
  ConservativeVars *accum_rhs = new ConservativeVars[numPoints];
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
  const FLOATTYPE coeffsA[4] = {0.0, 0.5, 0.5, 1.0};
  const FLOATTYPE coeffsB[4] = {1.0, 2.0, 2.0, 1.0};

  auto startTime = std::chrono::high_resolution_clock::now();

  while(time < t_final) {
    if(time + dt >= t_final) dt = t_final - time;

    for(int j = 0; j < numPoints; j++) {
      accum_rhs[j].rho = 0.0;
      accum_rhs[j].rho_u = 0.0;
      accum_rhs[j].rho_E = 0.0;
    }

    for(int s = 0; s < 4; s++) {
      const ConservativeVars *stage_state = (s == 0) ? state : Ui;
      computeRHS(stage_state, prev_rhs, flux, xj, numPoints, gamma);

      const FLOATTYPE b = coeffsB[s];
      if(s < 3) {
        const FLOATTYPE a = coeffsA[s + 1] * dt;
        for(int j = 0; j < numPoints; j++) {
          accum_rhs[j].rho += b * prev_rhs[j].rho;
          accum_rhs[j].rho_u += b * prev_rhs[j].rho_u;
          accum_rhs[j].rho_E += b * prev_rhs[j].rho_E;

          Ui[j].rho = state[j].rho + a * prev_rhs[j].rho;
          Ui[j].rho_u = state[j].rho_u + a * prev_rhs[j].rho_u;
          Ui[j].rho_E = state[j].rho_E + a * prev_rhs[j].rho_E;
        }
      } else {
        for(int j = 0; j < numPoints; j++) {
          accum_rhs[j].rho += b * prev_rhs[j].rho;
          accum_rhs[j].rho_u += b * prev_rhs[j].rho_u;
          accum_rhs[j].rho_E += b * prev_rhs[j].rho_E;
        }
      }
    }

    for(int j = 0; j < numPoints; j++) {
      state[j].rho += dt / 6.0 * accum_rhs[j].rho;
      state[j].rho_u += dt / 6.0 * accum_rhs[j].rho_u;
      state[j].rho_E += dt / 6.0 * accum_rhs[j].rho_E;
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
  delete[] prev_rhs;
  delete[] accum_rhs;
  delete[] Ui;
  delete[] xj;

  return 0;
}

void computeRHS(const ConservativeVars *state, ConservativeVars *rhs,
                ConservativeVars *flux, const FLOATTYPE *xj,
                int numPoints, FLOATTYPE gamma) {
  for(int j = 0; j < numPoints; j++) {
    FLOATTYPE u = state[j].rho_u / state[j].rho;
    FLOATTYPE E = state[j].rho_E / state[j].rho;
    FLOATTYPE p = (gamma - 1.0) * state[j].rho * (E - 0.5 * u * u);
    flux[j].rho = state[j].rho_u;
    flux[j].rho_u = state[j].rho_u * u + p;
    flux[j].rho_E = u * (state[j].rho_E + p);
  }

  FLOATTYPE dx0 = (xj[numPoints - 1] - xj[numPoints - 2]) + (xj[1] - xj[0]);
  rhs[0].rho = -(flux[1].rho - flux[numPoints - 2].rho) / dx0;
  rhs[0].rho_u = -(flux[1].rho_u - flux[numPoints - 2].rho_u) / dx0;
  rhs[0].rho_E = -(flux[1].rho_E - flux[numPoints - 2].rho_E) / dx0;

  for(int j = 1; j < numPoints - 1; j++) {
    FLOATTYPE dxinv = 1.0 / (xj[j + 1] - xj[j - 1]);
    rhs[j].rho = -(flux[j + 1].rho - flux[j - 1].rho) * dxinv;
    rhs[j].rho_u = -(flux[j + 1].rho_u - flux[j - 1].rho_u) * dxinv;
    rhs[j].rho_E = -(flux[j + 1].rho_E - flux[j - 1].rho_E) * dxinv;
  }

  rhs[numPoints - 1] = rhs[0];
}
