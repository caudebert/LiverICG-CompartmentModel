/*
Author : Chloe Audebert
******* Description ***********
ICG pharmaco model : all body circulation taken into account 
Liver with 3 compartments : Sinusoid, hepatocytes and bile canaliculi
Initial value in the blood is assumed non-zero, t0 is juste after injection
Model considere to fit ElDesoky et al. paper on ICG in rabbits, concentration are in arbitrary units
Sensitivity analysis : 
Observation : Liver amount
Parameter of interest : Ksh Qhb S and Fb

*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include <idas/idas.h>
#include <idas/idas_dense.h>
#include <nvector/nvector_serial.h>
#include <sundials/sundials_math.h>
#include <sundials/sundials_types.h>

/*Variables*/
#define NEQ 5
#define NS 4
#define NP 4

#define delta_t RCONST(0.5)
#define T0 RCONST(0.0)
#define Tf RCONST(1500.0)

#define Ksh RCONST(1.0) /*sinusoides <--> hepatocytes*/
#define Qhb RCONST(6.50) /*bile excression rate in 1/s*/
#define S RCONST(0.001)

/*Volume in mL*/
#define Vblood RCONST(160.0) // "total" blood volume
// Liver volume assumed 78.3mL~80mL 2.7% body weight (from litterature)
#define Vs RCONST(12.3) //sinosoids volume, about 15.3% of the estimated  Liver volume 
#define Vbc RCONST(2.7) //Bile canaliculi volume in mL 3.4% of the estimated liver volume
#define Vh RCONST(65.0) //hepatocytes volume, the rest: 81.3% of the estimated liver volume

//Ha Pv flow measured
#define Fha RCONST(0.25) /*mL/s HA flow*/
#define Fpv RCONST(1.5) /*mL/s PV flow*/

#define Fbc RCONST(0.05) /*mL/s bile flow*/


#define Y0 RCONST(83000.0) //init qty in Blood

/*the variable from which the sensitivity is computed*/
/*Liver: 5*/
#define varSens 5 
//solution file name
static char WriteIN[] = "data.txt";
static char SensFile[] = "Sens_L.txt";
static char NormFile[] = "D_L.txt";

/* Type : UserData */
typedef struct {
  realtype p[NS];           /* problem parameters */
  realtype coef;
} *UserData;

/* Macro to define dense matrix elements, indexed from 1. */
#define IJth(A,i,j) DENSE_ELEM(A,i-1,j-1)
#define Ith(v,i)    NV_Ith_S(v,i-1)       /* i-th vector component i=1..NEQ */

/* Prototypes of functions called by IDA */
int resrob(realtype tres, N_Vector yy, N_Vector yp,
   N_Vector resval, void *user_data);

int jacrob(long int Neq, realtype tt, realtype cj,
           N_Vector yy, N_Vector yp, N_Vector resvec,
           DlsMat JJ, void *user_data,
           N_Vector tempv1, N_Vector tempv2, N_Vector tempv3);

int resS(int Ns, realtype time, 
	 N_Vector yy, N_Vector yp, N_Vector resval,
	 N_Vector *yyS, N_Vector *ypS, N_Vector *resvalS,
	 void *user_data, 
	 N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

/* Prototypes of private functions */
static void PrintHeader(realtype rtol, N_Vector avtol, N_Vector y);

static void PrintOutput(void *mem, realtype t, N_Vector y, FILE* dataOutput, N_Vector *yS, FILE *file_s);

static void PrintRootInfo(int root_f1, int root_f2);

static void PrintFinalStats(void *mem);

static int check_flag(void *flagvalue, char *funcname, int opt);

int Update_D(N_Vector D, N_Vector *yS);

/*Global variable*/
int i;

/*
*--------------------------------------------------------------------
* Main Program
*--------------------------------------------------------------------
*/
int main(void)
{  
  /*sensi files*/
  FILE* file_S = NULL; /*sensitivity file*/
  FILE* file_D = NULL; /*D file scaling values*/
  file_S = fopen(SensFile, "w+");
  file_D = fopen(NormFile, "w+");

 /*def outputfile and ida variables*/
  FILE* dataOutput = NULL;
  dataOutput = fopen(WriteIN, "w+");

  void *mem;
  N_Vector yy, yp, avtol;
  realtype rtol, *yval, *ypval, *atval;
  realtype t0, tfinal,tout, tout1, tret, time;
  realtype abs_tolerance;
  int iout, retval, retvalr;

  realtype pbar[NS];
  int is, sensi_meth;
  booleantype err_con;  // Flag for error control
  UserData data; /*parameter for sensitivity*/
  N_Vector *yS, *ypS;
  N_Vector id, D; // The diagonal vector of the scaling matrix: scaled_sens = diag(D) * sens
  realtype *dval; 
  int retvals;

  mem = NULL;
  yy = yp = avtol = NULL;
  yval = ypval = atval = NULL;
  yS = ypS = NULL;
  D = id = NULL;
  dval = NULL;

  /* User data structure */
  data = (UserData) malloc(sizeof *data);
  if (check_flag((void *)data, "malloc", 2)) return(1);
  data->p[0] = Ksh;
  data->p[1] = Qhb;
  data->p[2] = Fbc;
  data->p[3] = S;
  data->coef = 0.5;
  /* Allocate N-vectors. */
  yy = N_VNew_Serial(NEQ);
  if(check_flag((void *)yy, "N_VNew_Serial", 0)) return(1);

  yp = N_VNew_Serial(NEQ);
  if(check_flag((void *)yp, "N_VNew_Serial", 0)) return(1);

  avtol = N_VNew_Serial(NEQ);
  if(check_flag((void *)avtol, "N_VNew_Serial", 0)) return(1);

  D = N_VNew_Serial(NS);
  if(check_flag((void *)D, "N_VNew_Serial", 0)) return(1);
  /* Create and initialize y, y', and absolute tolerance vectors. */
  time = T0;
  yval = NV_DATA_S(yy);
  ypval = NV_DATA_S(yp);
  
  /*initial condition --> supose no ICG is present in the body*/
  for(i=0;i<NEQ;i++)
    yval[i] = 0.0;
  
  yval[0] = Y0;
  yval[4] = yval[1] + yval[2] + yval[3]; 

  /*initial condition yp : compute from equations*/
  ypval[0] = (Fha+Fpv)*yval[1]/Vs - (Fha+Fpv)*yval[0]/Vblood;
  
  ypval[1] = (Fha+Fpv)*yval[0]/Vblood - (Fha+Fpv)*yval[1]/Vs - (Ksh/Vs)*yval[1] + (Ksh/Vh)*yval[2];
  
  ypval[2] = (Ksh/Vs)*yval[1] - (Ksh/Vh)*yval[2] - ( Qhb/(1.0 + (S*yval[2]/Vh) ) )*(yval[2]/Vh);
  
  ypval[3] = ( Qhb/(1.0 + (S*yval[2]/Vh) ) )*(yval[2]/Vh) - (Fbc/Vbc)*yval[3];
  
  ypval[4] = ypval[3] + ypval[1] + ypval[2];
  
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

  /* Attach user data */
  retval = IDASetUserData(mem, data);
  if (check_flag(&retval, "IDASetUserData", 1)) return(1);
  
  /* Call IDADense and set up the linear solver. */
  retval = IDADense(mem, NEQ);
  if(check_flag(&retval, "IDADense", 1)) return(1);

  retval = IDADlsSetDenseJacFn(mem, jacrob);
  if(check_flag(&retval, "IDADlsSetDenseJacFn", 1)) return(1);
  /* Sensitivity-related settings */
  sensi_meth = IDA_STAGGERED;
  err_con = FALSE;

  pbar[0] = data->p[0];
  pbar[1] = data->p[1];
  pbar[2] = data->p[2];
  pbar[3] = data->p[3];
  
  /* Create and initialize to zero the sensitivity variables and their derivatives */
  yS = N_VCloneVectorArray_Serial(NS, yy);
  if (check_flag((void *)yS, "N_VCloneVectorArray_Serial", 0)) return(1);
  for (is=0;is<NS;is++) N_VConst(0.0, yS[is]);
  
  ypS = N_VCloneVectorArray_Serial(NS, yp);
  if (check_flag((void *)ypS, "N_VCloneVectorArray_Serial", 0)) return(1);
  for (is=0;is<NS;is++) N_VConst(0.0, ypS[is]);

  /*Ksh*/
  Ith(ypS[0], 2) = (-yval[1]/Vs) + (yval[2]/Vh);
  Ith(ypS[0], 3) = (yval[1]/Vs) - (yval[2]/Vh);
  /*Qhb*/
  Ith(ypS[1], 3) = - (1.0/(1+ (S*yval[2]/Vh) ) )*yval[2]/Vh;
  Ith(ypS[1], 4) = (1.0/(1+ (S*yval[2]/Vh)) )*yval[2]/Vh;
  /*Fbc*/
  Ith(ypS[2], 4) = -yval[3]/Vbc;
  /*S */
  Ith(ypS[3], 3) = (1.0/( (1+(S*yval[2]/Vh))*(1+(S*yval[2]/Vh)) ))*(Qhb/Vh)*yval[2]*yval[2]/Vh;
  Ith(ypS[3], 4) = -(1.0/( (1+(S*yval[2]/Vh))*(1+(S*yval[2]/Vh)) ))*(Qhb/Vh)*yval[2]*yval[2]/Vh;
  
  retval = IDASensInit(mem, NS, sensi_meth, resS, yS, ypS);
  if(check_flag(&retval, "IDASensInit", 1)) return(1);
  
  retval = IDASensEEtolerances(mem);
  if(check_flag(&retval, "IDASensEEtolerances", 1)) return(1);
  
  retval = IDASetSensErrCon(mem, err_con);
  if (check_flag(&retval, "IDASetSensErrCon", 1)) return(1);

  retval = IDASetSensParams(mem, data->p, pbar, NULL);
  if (check_flag(&retval, "IDASetSensParams", 1)) return(1);
 
  /*loop to resolve*/
  tout = T0 + delta_t;
  while(1) {
    /*solve eq*/
    retval = IDASolve(mem, tout, &tret, yy, yp, IDA_NORMAL);
    if(check_flag(&retval, "IDASolve", 1)) return(1);
    /*solve sensi residual eq*/
    retvals = IDAGetSens(mem, &tret, yS);
    if(check_flag(&retvals, "IDAGetSens", 1)) return(1);
   
    PrintOutput(mem, tret, yy, dataOutput, yS, file_S);
    /*update D*/
    Update_D(D,yS);
 
    if (retval == IDA_SUCCESS) {
      tout += delta_t;
    }
    if (tout>Tf) break;
  }
  PrintFinalStats(mem);

  /* Write the updated D vector*/
  dval = NV_DATA_S(D);
  for(i=0;i<NP;i++)
    fprintf(file_D, "%12.12e ", 1/sqrt(dval[i]));

  fprintf(file_D, "\n ");
  /*Close outputData file*/
  fclose(dataOutput);
  fclose(file_S);
  fclose(file_D);

  /* Free memory */
  IDAFree(&mem);
  N_VDestroy_Serial(yy);
  N_VDestroy_Serial(yp);
  N_VDestroyVectorArray_Serial(yS, NS);
  N_VDestroyVectorArray_Serial(ypS, NS);
  free(data);
 
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
  /*parameters*/
  UserData data;
  data = (UserData) user_data;
  double Ksh1 = data->p[0];
  double Qhb1 = data->p[1];
  double Fbc1 = data->p[2];
  double S1 = data->p[3];
  
  yval = NV_DATA_S(yy);
  ypval = NV_DATA_S(yp);
  rval = NV_DATA_S(rr);

  rval[0] = ypval[0] - (Fha+Fpv)*yval[1]/Vs + (Fha+Fpv)*yval[0]/Vblood;

  rval[1] = ypval[1] - (Fha+Fpv)*yval[0]/Vblood + (Fha+Fpv)*yval[1]/Vs + Ksh1*(yval[1]/Vs - yval[2]/Vh);

  rval[2] = ypval[2] - Ksh1*(yval[1]/Vs - yval[2]/Vh) + ( Qhb1/(1.0 + (S1*yval[2]/Vh) ) )*yval[2]/Vh;

  rval[3] = ypval[3] - ( Qhb1/(1.0 + (S1*yval[2]/Vh) ) )*yval[2]/Vh + (Fbc1/Vbc)*yval[3];
  
  rval[4] = yval[4] - yval[1] - yval[2] - yval[3];

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
  /*parameters*/
  UserData data;
  data = (UserData) user_data;
  double Ksh1 = data->p[0];
  double Qhb1 = data->p[1];
  double Fbc1 = data->p[2];
  double S1 = data->p[3];

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
  IJth(JJ,2,2) = cj + (Ksh1/Vs) + ((Fha+Fpv)/Vs);
  IJth(JJ,3,2) = -Ksh1/Vs;
  IJth(JJ,5,2) = -1.0;
  
  IJth(JJ,2,3) = -Ksh1/Vh;
  IJth(JJ,3,3) = cj + (Ksh1/Vh) + ( (Qhb1/Vh) / ( (1+(S1*yval[2]/Vh))*(1+(S1*yval[2]/Vh)) ) );
  IJth(JJ,4,3) = - ( (Qhb1/Vh) / ( (1+(S1*yval[2]/Vh))*(1+(S1*yval[2]/Vh)) ) );
  IJth(JJ,5,3) = -1.0;
 
  IJth(JJ,4,4) = cj + Fbc1/Vbc;
  IJth(JJ,5,4) = -1.0;
 
  IJth(JJ,5,5) = 1.0;

  return(0);

}

/*

* Define the system residual function for sensitivity equations.

*/
int resS(int Ns, realtype time, 
	 N_Vector yy, N_Vector yp, N_Vector resval,
	 N_Vector *yyS, N_Vector *ypS, N_Vector *resvalS,
	 void *user_data, 
	 N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  UserData data;
  realtype s1, s2, s3, s4, s5;
  realtype sp1, sp2, sp3, sp4, sp5;
  realtype *yval;  

  yval = NV_DATA_S(yy); 
  int is;

  /*parameters*/
  data = (UserData) user_data;
  double Ksh1 = data->p[0];
  double Qhb1 = data->p[1];
  double Fbc1 = data->p[2];
  double S1 = data->p[3];
  
  for (is=0; is<NS; is++) {
    s1 = Ith(yyS[is],1);
    s2 = Ith(yyS[is],2);
    s3 = Ith(yyS[is],3);
    s4 = Ith(yyS[is],4);
    s5 = Ith(yyS[is],5);
   
    sp1 = Ith(ypS[is],1);
    sp2 = Ith(ypS[is],2);
    sp3 = Ith(ypS[is],3);
    sp4 = Ith(ypS[is],4);
    sp5 = Ith(ypS[is],5);
    //senstivity eq : dyF si + dypF sip + dpF
    //commmon equations = dyF si + dypF sip, dypF = id(5), dyF = jac   
    Ith(resvalS[is],1) = sp1 + ((Fha+Fpv)/Vblood)*s1 - ((Fha+Fpv)/Vs)*s2;

    Ith(resvalS[is],2) = sp2 - ((Fha+Fpv)/Vblood)*s1 + (Ksh1/Vs)*s2 + ((Fha+Fpv)/Vs)*s2 - (Ksh1/Vh)*s3;

    Ith(resvalS[is],3) = sp3 - (Ksh1/Vs)*s2 + (Ksh1/Vh)*s3 +  ((Qhb1/Vh) /( (1+(S1*yval[2]/Vh))*(1+(S1*yval[2]/Vh)) ) )*s3;

    Ith(resvalS[is],4) = sp4 - ((Qhb1/Vh)/( (1+(S1*yval[2]/Vh))*(1+(S1*yval[2]/Vh)) ) ) * s3 + (Fbc1/Vbc)*s4;

    Ith(resvalS[is],5) = -s2 - s3 - s4 + s5;
    //specific correction dpF
    switch (is) {
    case 0:/*Ksh*/
      Ith(resvalS[is],2) = Ith(resvalS[is],2) + yval[1]/Vs - yval[2]/Vh;
      Ith(resvalS[is],3) = Ith(resvalS[is],3) - yval[1]/Vs + yval[2]/Vh;
    case 1:/*Qhb*/
      Ith(resvalS[is],3) = Ith(resvalS[is],3) + (1.0/Vh)*yval[2]/( 1+ (S1*yval[2]/Vh));
      Ith(resvalS[is],4) = Ith(resvalS[is],4) - (1.0/Vh)*yval[2]/( 1+ (S1*yval[2]/Vh));
      break;
    case 2 : /*Fbc*/ 
      Ith(resvalS[is],4) = Ith(resvalS[is],4) + yval[3]/Vbc;
      break;
    case 3 : /*S*/
      Ith(resvalS[is],3) = Ith(resvalS[is],3) - (Qhb1/Vh)*yval[2]*(yval[2]/Vh)/( (1+(S1*yval[2]/Vh))*(1+(S1*yval[2]/Vh)) );
      Ith(resvalS[is],4) = Ith(resvalS[is],4) + (Qhb1/Vh)*yval[2]*(yval[2]/Vh)/( (1+(S1*yval[2]/Vh))*(1+(S1*yval[2]/Vh)) );
      break;
    
    }
  }
  
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

  printf("ICG model solution with blood circulation, Liver(3 compart) compartments \n");

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
static void PrintOutput(void *mem, realtype t, N_Vector y, FILE* dataOutput, N_Vector *yS, FILE *file_s)

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

  /*the sensitivity of y[varSens]*/
  double dys1 = Ith(yS[0],varSens);
  double dys2 = Ith(yS[1],varSens);
  double dys3 = Ith(yS[2],varSens); 
  double dys4 = Ith(yS[3],varSens); 

  fprintf(file_s, "%12.12e %12.12e %12.12e %12.12e \n", dys1, dys2, dys3, dys4);
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

int Update_D(N_Vector D, N_Vector *yS)
{
  realtype dys1, dys2, dys3, dys4;
  realtype *dval;
  /*the sensitivity of y VarSens*/
  dys1 = Ith(yS[0],varSens);
  dys2 = Ith(yS[1],varSens);
  dys3 = Ith(yS[2],varSens);
  dys4 = Ith(yS[3],varSens);

  dval = NV_DATA_S(D);
  
  dval[0] = dval[0] + dys1*dys1;
  dval[1] = dval[1] + dys2*dys2;
  dval[2] = dval[2] + dys3*dys3;
  dval[3] = dval[3] + dys4*dys4;
  
  return(0);
}

