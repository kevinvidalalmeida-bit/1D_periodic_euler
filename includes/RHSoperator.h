#ifndef _RHS_OPERATOR
#define _RHS_OPERATOR

#include "DataStructs.h"
#include "FluxFunctions.h"

template<class T>
class RHSOperator
{
public:
  RHSOperator();
  ~RHSOperator();

  virtual void eval() = 0;
  virtual void eval(DataStruct<T> &Uin) = 0;

  virtual DataStruct<T>& ref2RHS() = 0;
};

template<class T>
class Central1D : public RHSOperator<T>
{
private:

  // structure containing the RHS values
  DataStruct<T> RHS;
  
  // reference to current solution
  DataStruct<T> &U;

  // reference to mesh 
  // TODO: change to a mesh structure
  DataStruct<T> &mesh;

  // reference to flux function
  FluxFunction<T> &F;

  void evalRHS(DataStruct<T> &Uin);

public:
  Central1D(DataStruct<T> &_U, DataStruct<T> &_mesh, FluxFunction<T> &_F);
  ~Central1D();

  virtual void eval();
  virtual void eval(DataStruct<T> &Uin);

  virtual DataStruct<T>& ref2RHS();

};

// ============================================================
// EULER 1D RHS OPERATOR
// ============================================================

template<class T>
class Central1DEuler : public RHSOperator<T>
{
private:
  // RHS for each variable (rho, rho_u, rho_E)
  DataStruct<T> RHS_rho;
  DataStruct<T> RHS_rho_u;
  DataStruct<T> RHS_rho_E;
  
  // References to solution variables
  DataStruct<T> &rho;
  DataStruct<T> &rho_u;
  DataStruct<T> &rho_E;
  
  // Reference to mesh
  DataStruct<T> &mesh;
  
  // Reference to flux function
  EulerFlux<T> &F;
  
  // Temporary flux storage
  DataStruct<T> f_rho;
  DataStruct<T> f_rho_u;
  DataStruct<T> f_rho_E;
  
  void evalRHS(DataStruct<T> &rho_in, DataStruct<T> &rho_u_in, DataStruct<T> &rho_E_in);

public:
  Central1DEuler(DataStruct<T> &_rho, DataStruct<T> &_rho_u, DataStruct<T> &_rho_E,
                 DataStruct<T> &_mesh, EulerFlux<T> &_F);
  ~Central1DEuler();

  virtual void eval();
  virtual void eval(DataStruct<T> &Uin);  // dummy implementation
  
  // Evaluate RHS at intermediate values (for RK stages)
  void eval(DataStruct<T> &rho_in, DataStruct<T> &rho_u_in, DataStruct<T> &rho_E_in);
  
  virtual DataStruct<T>& ref2RHS() { return RHS_rho; }  // dummy
  
  // Get references to RHS for each variable
  DataStruct<T>& ref2RHS_rho() { return RHS_rho; }
  DataStruct<T>& ref2RHS_rho_u() { return RHS_rho_u; }
  DataStruct<T>& ref2RHS_rho_E() { return RHS_rho_E; }

};

#endif // _RHS_OPERATOR