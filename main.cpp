#include <iostream> // std::cout
#include <fstream> // ofile
#include <string>
#include <math.h>
#include <iomanip> // set precision
#include <vector>
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

#ifdef _DOUBLE_
#define MPI_FLOATTYPE MPI_DOUBLE
#else
#define MPI_FLOATTYPE MPI_FLOAT
#endif

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// declare supporting functions
void write2File(DataStruct<FLOATTYPE> &X, DataStruct<FLOATTYPE> &U, std::string name);
FLOATTYPE calcL2normSquaredExact(DataStruct<FLOATTYPE> &x, DataStruct<FLOATTYPE> &u, FLOATTYPE k, FLOATTYPE time);
void buildPartition(int numPoints, int worldSize, int rank, int &localNumPoints, int &offset);
void buildCountsDisplacements(int numPoints, int worldSize, std::vector<int> &counts, std::vector<int> &displacements);
void exchangeGhosts(DataStruct<FLOATTYPE> &U, FLOATTYPE ghosts[2], int rank, int worldSize);
void gatherToRoot(DataStruct<FLOATTYPE> &local, DataStruct<FLOATTYPE> &global, int numPoints, int rank, int worldSize);


int main(int narg, char **argv)
{
  MPI_Init(&narg, &argv);

  int rank, worldSize;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &worldSize);

  int numPoints =  80;
  FLOATTYPE k = 2.; // wave number
  FLOATTYPE t_final = 1.;

  if(narg != 3 && narg != 4)
  {
    if(rank == 0)
    {
      std::cout<< "Wrong number of arguments. You should include:" << std::endl;
      std::cout<< "    Num points" << std::endl;
      std::cout<< "    Wave number" << std::endl;
      std::cout<< "    Final time (optional)" << std::endl;
    }
    MPI_Finalize();
    return 1;
  }else
  {
    numPoints = std::stoi(argv[1]);
    k         = std::stod(argv[2]);
    if(narg == 4) t_final = std::stod(argv[3]);
  }

  if(numPoints < worldSize)
  {
    if(rank == 0)
    {
      std::cout << "The number of points must be at least the number of MPI processes." << std::endl;
    }
    MPI_Abort(MPI_COMM_WORLD, 2);
  }

  int localNumPoints, offset;
  buildPartition(numPoints, worldSize, rank, localNumPoints, offset);

  // solution data
  DataStruct<FLOATTYPE> u(localNumPoints), xj(localNumPoints);

  // flux function
  LinearFlux<FLOATTYPE> lf;

  // time solver
  RungeKutta4<FLOATTYPE> rk(u);

  // Initial Condition
  FLOATTYPE *datax = xj.getData();
  FLOATTYPE *dataU = u.getData();
  const FLOATTYPE dx = FLOATTYPE(1.) / FLOATTYPE(numPoints);
  for(int j = 0; j < localNumPoints; j++)
  {
    const int globalIndex = offset + j;

    // xj
    datax[j] = FLOATTYPE(globalIndex) * dx;

    // init Uj
    dataU[j] = sin(k*2. * M_PI * datax[j]);
  }

  // Operator
  Central1D<FLOATTYPE> rhs(u,xj,lf,dx);

  FLOATTYPE CFL = 2.4;
  FLOATTYPE dt = CFL*dx;

  // Output Initial Condition
  DataStruct<FLOATTYPE> globalX, globalU;
  gatherToRoot(xj, globalX, numPoints, rank, worldSize);
  gatherToRoot(u, globalU, numPoints, rank, worldSize);
  if(rank == 0) write2File(globalX, globalU, "initialCondition.csv");

  FLOATTYPE time = 0.;

  // init timer
  MPI_Barrier(MPI_COMM_WORLD);
  double compTime = MPI_Wtime();

  // main loop
  while(time < t_final)
  {
    if(time+dt >= t_final) dt = t_final - time;

    // take RK step
    rk.initRK();
    for(int s = 0; s < rk.getNumSteps(); s++)
    {
      rk.stepUi(dt);

      FLOATTYPE ghosts[2];
      DataStruct<FLOATTYPE> *Ui = rk.currentU();
      exchangeGhosts(*Ui, ghosts, rank, worldSize);
      rhs.setGhostValues(ghosts[0], ghosts[1]);
      rhs.eval(*Ui);

      rk.setFi(rhs.ref2RHS());
    }
    rk.finalizeRK(dt);
    time += dt;
  }

  // finishe timer
  compTime = MPI_Wtime() - compTime;
  double globalCompTime;
  MPI_Reduce(&compTime, &globalCompTime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

  gatherToRoot(u, globalU, numPoints, rank, worldSize);
  if(rank == 0) write2File(globalX, globalU, "final.csv");

  // L2 norm
  FLOATTYPE localErr2 = calcL2normSquaredExact(xj, u, k, t_final);
  FLOATTYPE globalErr2;
  MPI_Reduce(&localErr2, &globalErr2, 1, MPI_FLOATTYPE, MPI_SUM, 0, MPI_COMM_WORLD);

  if(rank == 0)
  {
    const FLOATTYPE err = sqrt(globalErr2);
    std::cout << std::setprecision(4) << "MPI processes: " << worldSize;
    std::cout << ". Comp. time: " << globalCompTime;
    std::cout << " sec. Error: " << err/k;
    std::cout << " kdx: " << k*dx*2.*M_PI;
    std::cout << std::endl;
  }

  MPI_Finalize();
  return 0;
}


// ==================================================================
// AUXILIARY FUNCTIONS
// ==================================================================
void write2File(DataStruct<FLOATTYPE> &X, DataStruct<FLOATTYPE> &U, std::string name)
{
  std::ofstream file;
  file.open(name,std::ios_base::trunc);
  if(!file.is_open()) 
  {
    std::cout << "Couldn't open file for Initial Condition" << std::endl;
    exit(1);
  }
  
  for(int j = 0; j < U.getSize(); j++)
  {
    file << X.getData()[j] << " ," << U.getData()[j] << std::endl;
  }

  file.close();
}

FLOATTYPE calcL2normSquaredExact(DataStruct<FLOATTYPE> &x, DataStruct<FLOATTYPE> &u, FLOATTYPE k, FLOATTYPE time)
{
  FLOATTYPE err = 0.;
  const FLOATTYPE *dataX = x.getData();
  const FLOATTYPE *dataU = u.getData();

  for(int n = 0; n < u.getSize(); n++)
  {
    const FLOATTYPE exact = sin(k*2. * M_PI * (dataX[n] - time));
    err += (dataU[n] - exact)*(dataU[n] - exact);
  }

  return err;
}

void buildPartition(int numPoints, int worldSize, int rank, int &localNumPoints, int &offset)
{
  const int base = numPoints / worldSize;
  const int remainder = numPoints % worldSize;

  localNumPoints = base;
  if(rank == worldSize-1) localNumPoints += remainder;

  offset = rank * base;
}

void buildCountsDisplacements(int numPoints, int worldSize, std::vector<int> &counts, std::vector<int> &displacements)
{
  counts.resize(worldSize);
  displacements.resize(worldSize);

  const int base = numPoints / worldSize;
  const int remainder = numPoints % worldSize;

  for(int r = 0; r < worldSize; r++)
  {
    counts[r] = base;
    if(r == worldSize-1) counts[r] += remainder;
    displacements[r] = r * base;
  }
}

void exchangeGhosts(DataStruct<FLOATTYPE> &U, FLOATTYPE ghosts[2], int rank, int worldSize)
{
  FLOATTYPE *dataU = U.getData();
  const int last = U.getSize()-1;
  const int leftRank = (rank - 1 + worldSize) % worldSize;
  const int rightRank = (rank + 1) % worldSize;

  MPI_Sendrecv(&dataU[0], 1, MPI_FLOATTYPE, leftRank, 10,
               &ghosts[1], 1, MPI_FLOATTYPE, rightRank, 10,
               MPI_COMM_WORLD, MPI_STATUS_IGNORE);

  MPI_Sendrecv(&dataU[last], 1, MPI_FLOATTYPE, rightRank, 20,
               &ghosts[0], 1, MPI_FLOATTYPE, leftRank, 20,
               MPI_COMM_WORLD, MPI_STATUS_IGNORE);
}

void gatherToRoot(DataStruct<FLOATTYPE> &local, DataStruct<FLOATTYPE> &global, int numPoints, int rank, int worldSize)
{
  std::vector<int> counts, displacements;
  buildCountsDisplacements(numPoints, worldSize, counts, displacements);

  if(rank == 0 && global.getSize() == 0) global.setSize(numPoints);

  MPI_Gatherv(local.getData(), local.getSize(), MPI_FLOATTYPE,
              rank == 0 ? global.getData() : NULL, counts.data(), displacements.data(),
              MPI_FLOATTYPE, 0, MPI_COMM_WORLD);
}
