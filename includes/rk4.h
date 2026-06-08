#ifndef _RUNGE_KUTTA
#define _RUNGE_KUTTA

/*
RUNGE KUTTA 4 CLASS

This class allows time stepping. For flexibility each RK step is taken as follows:
  initRK();
  for(s = 0; s<rk.getNumSteps(); s++)
    1. stepUi(dt)
    2. impose BCs over the current Ui which can be obtained through currentU()
    3. compute current F (i.e., Fi = F(Ui) )
    5. setFi(Fi)
  finalizeRK(dt);

*/

#include "DataStructs.h"
#include "FluxFunctions.h"

template<class T>
class RungeKutta4 
{
  private:
    int nSteps;
    int currentStep;

    T *coeffsA, *coeffsB;

    // reference to solution (Un)
    DataStruct<T> &Un;

    // intermediate solution
    DataStruct<T> Ui;

    // RHS 
    DataStruct<T> *fi;

  public:

    // default constructor
    RungeKutta4(DataStruct<T> &_Un);

    // default destructor
    ~RungeKutta4();

    int getNumSteps();

    // initialize the RK
    void initRK();

    // finalizes the RK (updates Un)
    void finalizeRK(const T dt);

    /*
    For the step to work properly, the user must provide the appropriate F for the current Ui.
    This is done this way becase the user might want to modify the Ui or Fi so that 
    Boundary conditions can be imposed
    */
   void stepUi(T dt);
   void setFi(DataStruct<T> &_F);

   // current Ui
   DataStruct<T>* currentU();
};

// ============================================================
// RUNGE KUTTA 4 FOR EULER 1D (SYSTEMS)
// ============================================================

template<class T>
class RungeKutta4Euler
{
  private:
    int nSteps;
    int currentStep;

    T *coeffsA, *coeffsB;

    // References to solutions (Un for each variable)
    DataStruct<T> &Un_rho;
    DataStruct<T> &Un_rho_u;
    DataStruct<T> &Un_rho_E;

    // Intermediate solutions
    DataStruct<T> Ui_rho;
    DataStruct<T> Ui_rho_u;
    DataStruct<T> Ui_rho_E;

    // RHS arrays (4 for RK4)
    DataStruct<T> *fi_rho;
    DataStruct<T> *fi_rho_u;
    DataStruct<T> *fi_rho_E;
    
    // Accumulator for final RHS
    DataStruct<T> RHS_rho;
    DataStruct<T> RHS_rho_u;
    DataStruct<T> RHS_rho_E;

  public:
    // Constructor
    RungeKutta4Euler(DataStruct<T> &_Un_rho, DataStruct<T> &_Un_rho_u, DataStruct<T> &_Un_rho_E);

    // Destructor
    ~RungeKutta4Euler();

    int getNumSteps();

    // Initialize the RK
    void initRK();

    // Finalize the RK (updates Un)
    void finalizeRK(const T dt);

    // Step Ui
    void stepUi(T dt);
    
    // Set Fi for current step
    void setFi(DataStruct<T> &_F_rho, DataStruct<T> &_F_rho_u, DataStruct<T> &_F_rho_E);

    // Get current Ui
    void currentU(DataStruct<T> &out_rho, DataStruct<T> &out_rho_u, DataStruct<T> &out_rho_E);
};

#endif // _RUNGE_KUTTA