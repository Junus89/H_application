#ifndef FDSolver_H_
#include<complex.h> 
#include<math.h>
#define FDSolver_H_
//#define PI atan(1,1);
#define PI	3.14159265358979323846


//double Tint,TR,fR,OmegaR;
double DT;
int DSNum, BNum, FNum, TNum;
double R,Rou_0,C_0,RMaTip;
double VectorX,VectorY,VectorZ,VuN,VvN,VUN;
double OX,OY,OZ, MaX, MaY, MaZ;
double Theta=0;
double DX,DY,DZ,DR,DVX,DVY,DVZ,Dp,Drou,Du,Dv,Dw,Gu,Gv,Gw;
double Time[1440];





double A =0;
double *OTime;
double OmegaM, OmegaR, Gamma,  DX, DY, DZ,  DR, OX, OY,\
	 OZ, Theta, *DataXR, *DataYR, *DataZR, *DOrX, *DOrY,  *DOrZ, *DOr, *DORStar, *DOR,\
		 *DORStarX, *DORStarY, *DORStarZ, *DORX, *DORY,*DORZ, *RGamma, *Vx, *Vy, *Vz, *Q,\
			 *Lx, *Ly, *Lz, *FxM, *FyM, *FzM, *FxP, *FyP, *FzP;
double *FRStarM, *FRStarP, *FRM, *FRP;
double complex z1 = 0.0 + 1.0*I;
double complex Inicomp = 0.0 + 0.0*I;

double complex SP1;
double complex SP2;
double complex SP3;




double **make_dmatrix(size_t m, size_t n);// construct 2D array type of double
void free_dmatrix(double **a, size_t n);// freeing it

double ***make_3Ddmatrix(size_t p, size_t q, size_t r);// construct 3D array type of double
void free_3Ddmatrix(double ***a, size_t q, size_t r); // freeing it

double ****make_4Ddmatrix(size_t o, size_t p, size_t q, size_t r); // construct 4D array type of double
void free_4Ddmatrix(double ****a, size_t p, size_t q, size_t r);


// ----fifthloop prepreation-----



//void fifthLoop();
//void fifthLoop(double OmegaR, double Omega, double MaX, double MaY, double MaZ, double ka, double **pF,double complex SP1,double complex SP2,double complex SP3);
void fifthLoop(double OmegaR, double Omega, double MaX, double MaY, double MaZ, double ka);


double *Ones(int n);// Ones function


#endif


