/*
*******************************************************
******* Author : Chloe Audebert
******* Last mod : 17/10/2018
*******************************************************
******* Description ***********
ICG pharmaco model : all body circulation taken into account 
Liver with 3 compartments : Sinusoid, hepatocytes and bile canaliculi
Initial value in the blood is assumed non-zero (just after injection)
Model is considered to fit El-Desoky et al. 1999  data of ICG fluo in rabbits, 
concentrations are in arbitrary units
*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include <ida/ida.h>
#include <ida/ida_dense.h>
#include <nvector/nvector_serial.h>
#include <sundials/sundials_math.h>
#include <sundials/sundials_types.h>

#define NEQ 4

#define delta_t RCONST(0.1)
#define T0 RCONST(0.0)
#define Tf RCONST(1500.0) //time in second
//solution file name
static char WriteIN[] = "./Res/data.txt";

#define Ksh RCONST(20.3)/*sinusoides <--> hepatocytes (volume/time)*/
#define Qhb RCONST(0.46) /*bile excression volume/time*/
#define S RCONST(0.0006) /*saturation parameter 1/concentration  */
#define Fbc RCONST(0.16) /*mL/s bile canaliculi duct flow*/

//Ha Pv flow from El-Desoky measurements
#define Fha RCONST(0.25) /*ml/s HA flow*/
#define Fpv RCONST(1.75) /*ml/s PV flow*/
//Control HA+PV : 1.75
//10% occl HA+PH: 1.58
//25% occl HA+PH: 1.31
//40% occl HA+PH:1.05
//70% occ HA+PH: 0.53
//90% occl HA+PH:0.175
//intial concentration in the blood compartment
#define Y0 RCONST(83000.0)

/*Volume in mL*/
#define Vblood RCONST(160.0) // "total" blood volume 56ml/kg and rabbit weight = 2.9kg
// Liver volume assumed 80mL
#define Vs RCONST(12.3) //sinosoids volume, about 15.3% of the estimated  Liver volume 
#define Vbc RCONST(2.7) //Bile canaliculi volume in mL 3.4% of the estimated liver volume
#define Vh RCONST(65.0) //hepatocytes volume, the rest: 81.3% of the estimated liver volume


/* Macro to define dense matrix elements, indexed from 1. */
#define IJth(A,i,j) DENSE_ELEM(A,i-1,j-1)

/* Prototypes of functions called by IDA */
int resrob(realtype tres, N_Vector yy, N_Vector yp,
   N_Vector resval, void *user_data);

int jacrob(long int Neq, realtype tt, realtype cj,
   N_Vector yy, N_Vector yp, N_Vector resvec,
   DlsMat JJ, void *user_data,
   N_Vector tempv1, N_Vector tempv2, N_Vector tempv3);

/* Prototypes of private functions */
static void PrintHeader(realtype rtol, N_Vector avtol, N_Vector y);

static void PrintOutput(void *mem, realtype t, N_Vector y, FILE* dataOutput);

static void PrintRootInfo(int root_f1, int root_f2);

static void PrintFinalStats(void *mem);

static int check_flag(void *flagvalue, char *funcname, int opt);

/* Global variables */
int i;

/*
*--------------------------------------------------------------------
* Main Program
*--------------------------------------------------------------------
*/
int main(void)
{  
 /*def outputfile and ida variables*/
  FILE* dataOutput = NULL;
  dataOutput = fopen(WriteIN, "w+");
  
  void *mem;
  N_Vector yy, yp, avtol;
  realtype rtol, *yval, *ypval, *atval;
  realtype t0, tfinal,tout, tout1, tret, time;
  realtype abs_tolerance;
  int retval, retvalr;

  mem = NULL;
  yy = yp = avtol = NULL;
  yval = ypval = atval = NULL;

  /* Allocate N-vectors. */
  yy = N_VNew_Serial(NEQ);
  if(check_flag((void *)yy, "N_VNew_Serial", 0)) return(1);

  yp = N_VNew_Serial(NEQ);
  if(check_flag((void *)yp, "N_VNew_Serial", 0)) return(1);

  avtol = N_VNew_Serial(NEQ);
  if(check_flag((void *)avtol, "N_VNew_Serial", 0)) return(1);

  /* Create and initialize y, y'=yp, and absolute tolerance vectors. */
  time = T0;
  yval = NV_DATA_S(yy);
  ypval = NV_DATA_S(yp);
  
  /*initial condition --> supose no ICG is present in the body*/
  for(i=0;i<NEQ;i++)
    yval[i] = 0.0;
  //Except in the blood circulation, t=0 is just after total bolus injection
  yval[0] = Y0;

  //initial y and yp have to verify the equations:
  /*initial condition yp : compute from equations*/
  ypval[0] = (Fha+Fpv)*yval[1]/Vs - (Fha+Fpv)*yval[0]/Vblood;

  ypval[1] = (Fha+Fpv)*yval[0]/Vblood - (Fha+Fpv)*yval[1]/Vs - Ksh*(yval[1]/Vs - yval[2]/Vh);
  
  ypval[2] = Ksh*(yval[1]/Vs - yval[2]/Vh) - ( Qhb/(1.0 + (S*yval[2]/Vh) ) )*yval[2]/Vh;
  
  ypval[3] = ( Qhb/(1.0 + (S*yval[2]/Vh) ) )*yval[2]/Vh - (Fbc/Vbc)*yval[3];

 
  /*define tolerance*/
  rtol = RCONST(1.0e-06);
  abs_tolerance=1.0e-6;
  atval = NV_DATA_S(avtol);
  for(i=0;i<NEQ;i++)
    atval[i] = RCONST(abs_tolerance);

  PrintHeader(rtol, avtol, yy);

  /* Call IDACreate and IDAInit to initialize IDA memory */
  mem = IDACreate();
  if(check_flag((void *)mem, "IDACreate", 0)) return(1);

  retval = IDAInit(mem, resrob, T0, yy, yp);
  if(check_flag(&retval, "IDAInit", 1)) return(1);

  /* Call IDASVtolerances to set tolerances */
  retval = IDASVtolerances(mem, rtol, avtol);
  if(check_flag(&retval, "IDASVtolerances", 1)) return(1);
  /* Free avtol */
  N_VDestroy_Serial(avtol);

  /* Call IDADense and set up the linear solver. */
  
  retval = IDADense(mem, NEQ);
  if(check_flag(&retval, "IDADense", 1)) return(1);
 
  retval = IDADlsSetDenseJacFn(mem, jacrob);
  if(check_flag(&retval, "IDADlsSetDenseJacFn", 1)) return(1);
  
  tout = T0 + delta_t;
  while(1) {

    retval = IDASolve(mem, tout, &tret, yy, yp, IDA_NORMAL);
    PrintOutput(mem, tret, yy, dataOutput);
  
    if(check_flag(&retval, "IDASolve", 1)) return(1);

    if (retval == IDA_SUCCESS) {
      tout += delta_t;
    }
    if (tout>Tf) break;
  }
  PrintFinalStats(mem);

  /*Close outputData file*/
  fclose(dataOutput);
  /* Free memory */
  IDAFree(&mem);
  N_VDestroy_Serial(yy);
  N_VDestroy_Serial(yp);
  
  return(0);

}
/*
*--------------------------------------------------------------------
* Functions called by IDA
*--------------------------------------------------------------------
*/

/*
* Define the system residual function.
*/
int resrob(realtype time, N_Vector yy, N_Vector yp, N_Vector rr, void *user_data)
{
  realtype *yval, *ypval, *rval;
  int i;
 
  yval = NV_DATA_S(yy);
  ypval = NV_DATA_S(yp);
  rval = NV_DATA_S(rr);

  rval[0] = ypval[0] - (Fha+Fpv)*yval[1]/Vs + (Fha+Fpv)*yval[0]/Vblood;

  rval[1] = ypval[1] - (Fha+Fpv)*yval[0]/Vblood + (Fha+Fpv)*yval[1]/Vs + Ksh*(yval[1]/Vs - yval[2]/Vh);

  rval[2] = ypval[2] - Ksh*(yval[1]/Vs - yval[2]/Vh) + ( Qhb/(1.0 + (S*yval[2]/Vh) ) )*yval[2]/Vh;

  rval[3] = ypval[3] - ( Qhb/(1.0 + (S*yval[2]/Vh) ) )*yval[2]/Vh + (Fbc/Vbc)*yval[3];
   
  return(0);
}

/*

* Define the Jacobian function.

*/
int jacrob(long int Neq, realtype tt, realtype cj,
   N_Vector yy, N_Vector yp, N_Vector resvec,
   DlsMat JJ, void *user_data,
   N_Vector tempv1, N_Vector tempv2, N_Vector tempv3)

{
  int i,j;
  realtype *yval;
  yval = NV_DATA_S(yy);
  for(i=1;i<NEQ+1;i++){
    for(j=1; j<NEQ+1;j++){
      IJth(JJ,i,j) = 0.0;
    }
  }
  //Change non-zero values : 
  IJth(JJ,1,1) = cj + (Fha+Fpv)/Vblood;
  IJth(JJ,2,1) = -(Fha+Fpv)/Vblood;

  IJth(JJ,1,2) = -(Fha+Fpv)/Vs;
  IJth(JJ,2,2) = cj + (Ksh/Vs) + ((Fha+Fpv)/Vs);
  IJth(JJ,3,2) = -Ksh/Vs;
  
  IJth(JJ,2,3) = -Ksh/Vh;
  IJth(JJ,3,3) = cj + (Ksh/Vh) + ( (Qhb/Vh) / ( (1+(S*yval[2]/Vh))*(1+(S*yval[2]/Vh)) ) );
  IJth(JJ,4,3) = - ( (Qhb/Vh) / ( (1+(S*yval[2]/Vh))*(1+(S*yval[2]/Vh)) ) );

  IJth(JJ,4,4) = cj + Fbc/Vbc;

  return(0);
}
/*
*--------------------------------------------------------------------
* Private functions
*--------------------------------------------------------------------
*/

 

/*
* Print first lines of output (problem description)
*/
static void PrintHeader(realtype rtol, N_Vector avtol, N_Vector y)

{
  realtype *atval, *yval;

  atval = NV_DATA_S(avtol);
  yval = NV_DATA_S(y);

  printf("ICG model solution with blood circulation, Sinusoids, hepatocytes and bile canaliculi compartments \n");
  printf("Tolerance parameters: rtol = %1.4e atol = %1.4e\n", rtol, atval[0]);
  printf("Initial conditions y0 = (");
  for(i=0;i<NEQ;i++){
    printf("%1.4e ", yval[i]); 
  }
  printf(")\n");
 
}
/*
* Print Output
*/
static void PrintOutput(void *mem, realtype t, N_Vector y, FILE* dataOutput)
{
  realtype *yval;
  int retval, kused;
  long int nst;
  realtype hused;

  yval = NV_DATA_S(y);

  retval = IDAGetLastOrder(mem, &kused);
  check_flag(&retval, "IDAGetLastOrder", 1);

  retval = IDAGetNumSteps(mem, &nst);
  check_flag(&retval, "IDAGetNumSteps", 1);

  retval = IDAGetLastStep(mem, &hused);
  check_flag(&retval, "IDAGetLastStep", 1);

  fprintf(dataOutput,"%10.10e ", t);
  for(i=0;i<NEQ;i++)
    fprintf(dataOutput, "%12.10e ", yval[i]); 

  fprintf(dataOutput, " %3ld %1d %12.4e\n", nst, kused, hused );

}
/*
* Print final integrator statistics
*/
static void PrintFinalStats(void *mem)
{

 int retval;

 long int nst, nni, nje, nre, nreLS, netf, ncfn, nge;

 retval = IDAGetNumSteps(mem, &nst);
 check_flag(&retval, "IDAGetNumSteps", 1);

 retval = IDAGetNumResEvals(mem, &nre);
 check_flag(&retval, "IDAGetNumResEvals", 1);

 retval = IDADlsGetNumJacEvals(mem, &nje);
 check_flag(&retval, "IDADlsGetNumJacEvals", 1);

 retval = IDAGetNumNonlinSolvIters(mem, &nni);
 check_flag(&retval, "IDAGetNumNonlinSolvIters", 1);

 retval = IDAGetNumErrTestFails(mem, &netf);
 check_flag(&retval, "IDAGetNumErrTestFails", 1);

 retval = IDAGetNumNonlinSolvConvFails(mem, &ncfn);
 check_flag(&retval, "IDAGetNumNonlinSolvConvFails", 1);

 retval = IDADlsGetNumResEvals(mem, &nreLS);
 check_flag(&retval, "IDADlsGetNumResEvals", 1);

 retval = IDAGetNumGEvals(mem, &nge);
 check_flag(&retval, "IDAGetNumGEvals", 1);

 printf("\nFinal Run Statistics: \n\n");

 printf("Number of steps = %ld\n", nst);

 printf("Number of residual evaluations = %ld\n", nre+nreLS);

 printf("Number of Jacobian evaluations = %ld\n", nje);

 printf("Number of nonlinear iterations = %ld\n", nni);

 printf("Number of error test failures = %ld\n", netf);

 printf("Number of nonlinear conv. failures = %ld\n", ncfn);

 printf("Number of root fn. evaluations = %ld\n", nge);

}

/*

* Check function return value...

* opt == 0 means SUNDIALS function allocates memory so check if

* returned NULL pointer

* opt == 1 means SUNDIALS function returns a flag so check if

* flag >= 0

* opt == 2 means function allocates memory so check if returned

* NULL pointer

*/

static int check_flag(void *flagvalue, char *funcname, int opt)
{

 int *errflag;
 /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
 if (opt == 0 && flagvalue == NULL) {
   fprintf(stderr,

       "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",

       funcname);

   return(1);

 } else if (opt == 1) {
   /* Check if flag < 0 */
   errflag = (int *) flagvalue;
   if (*errflag < 0) {

     fprintf(stderr,

         "\nSUNDIALS_ERROR: %s() failed with flag = %d\n\n",

         funcname, *errflag);

     return(1);

   }

 } else if (opt == 2 && flagvalue == NULL) {
   /* Check if function returned NULL pointer - no memory allocated */
   fprintf(stderr,

       "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",

       funcname);

   return(1);

 }
 return(0);
 
}
