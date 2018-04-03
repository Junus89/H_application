#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<complex.h> /* Standard Library of Complex Numbers */

#include"FDSolver.h"

/* ---------- function definitions --------- */

/* Memory allocation function for different dimension arrays */
//***************************1D*******************************************/
float *vector(long nl,long nh)
/*	
 */
{
    float *outputAr; 

    outputAr = (float *)malloc((nh+1)*sizeof(float));
    if (outputAr == NULL)
    {
        printf("No enough memory space.");
        exit(1);
    }
    return outputAr;
}

void free_vector(float *inputAr,long nl,long nh)
{
    free(inputAr);
    return;
}

/*----------------------------------------------------------- */
double *dvector(long nl,long nh)
/*	
 */
{
    double *outputAr; 

    outputAr = (double *)malloc((nh+1)*sizeof(double));
    if (outputAr == NULL)
    {
        printf("No enough memory space.");
        exit(1);
    }
    return outputAr;
}


void free_dvector(double *inputAr,long nl,long nh)
{
    free(inputAr);
    return;
}

//***************************2D*******************************************/
float **matrix(long nrl,long nrh, long ncl,long nch)
/*	
 */
{
    int i;
    float **outputAr;
    outputAr = (float **)malloc((nrh+1)*sizeof(float *));
    if (outputAr == NULL)
    {
        printf("No enough memory space.");
        exit(1);
    }
	
    for (i = 0; i < nrh+1; i++)
    {
        outputAr[i] = (float *)malloc((nch+1)*sizeof(float));
        if (outputAr[i] == NULL)
		{
	       printf("No enough memory space.");
	       exit(1);
		}
    }		
    return outputAr;
}

void free_matrix(float **inputAr,long nrl,long nrh,long ncl,long nch)
/*	
 */
{
    int i;
    for (i = 0; i < nrh+1; i++)	free(inputAr[i]);
    free(inputAr);		
    return;
}


/*-----------------------------------------------------------------------*/

double **dmatrix(long nrl,long nrh, long ncl,long nch)
/*	
 */
{
    int i;
    double **outputAr;
    outputAr = (double **)malloc((nrh+1)*sizeof(double *));
    if (outputAr == NULL)
    {
        printf("No enough memory space.");
        exit(1);
    }
	
    for (i = 0; i < nrh+1; i++)
    {
        outputAr[i] = (double *)malloc((nch+1)*sizeof(double));
        if (outputAr[i] == NULL)
		{
	        printf("No enough memory space.");
	        exit(1);
		}
	}		
    return outputAr;
}

void free_dmatrix(double **inputAr,long nrl,long nrh,long ncl,long nch)
/*	
 */
{
    int i;
    for (i = 0; i <= nrh; i++)  free(inputAr[i]);
    free(inputAr);		
    return;
}

//***************************3D*******************************************/

float ***f3dmatrix(long nrl,long nrh, long ncl,long nch,long ndl,long ndh)
/*
 */
{
    int i,j;
    float ***outputAr;
    outputAr = (float ***)malloc((nrh+1)*sizeof(float**));
    if (outputAr == NULL)
    {
        printf("No enough memory space.");
        exit(1);
    }
	
    for (i = 0; i < nrh+1; i++)
    {
        outputAr[i] = (float **)malloc((nch+1)*sizeof(float*));
        if (outputAr[i] == NULL)
		{
	        printf("No enough memory space.");
	        exit(1);
		}
	}

    for (i = 0; i < nrh+1; i++)
    {
        for (j = 0; j < nch+1; j++)
		{

	        outputAr[i][j] = (float*)malloc((ndh+1)*sizeof(float));
	        if (outputAr[i][j] == NULL)
			{
	            printf("No enough memory space.");
	            exit(1);
			}
		}
    }	
    return outputAr;
}



void free_f3dmatrix(float ***inputAr,long nrl,long nrh,long ncl,long nch,long ndl,long ndh)
/*	
 */
{
    int i,j;
	
    for (i = 0; i < nrh+1; i++)
    {
        for (j = 0; j < nch+1; j++)
		{
	        free(inputAr[i][j]);
		}
	}
    for (i = 0; i < nrh+1; i++)	free(inputAr[i]);
    free(inputAr);		
    return;
}


/*-----------------------------------------------------------------------*/

double ***d3dmatrix(long nrl,long nrh, long ncl,long nch,long ndl,long ndh)
/*
 */
{
    int i,j;
    double ***outputAr;
    outputAr = (double ***)malloc((nrh+1)*sizeof(double **));
    if (outputAr == NULL)
    {
        printf("No enough memory space.");
        exit(1);
    }
	
    for (i = 0; i < nrh+1; i++)
    {
        outputAr[i] = (double **)malloc((nch+1)*sizeof(double*));
        if (outputAr[i] == NULL)
		{
	        printf("No enough memory space.");
	        exit(1);
		}
	}

    for (i = 0; i < nrh+1; i++)
    {
        for (j = 0; j < nch+1; j++)
		{

	        outputAr[i][j] = (double *)malloc((ndh+1)*sizeof(double));
	        if (outputAr[i][j] == NULL)
			{
	            printf("No enough memory space.");
	            exit(1);
			}
		}
	}	
    return outputAr;
}

void free_d3dmatrix(double ***inputAr,long nrl,long nrh,long ncl,long nch,long ndl,long ndh)
/*	
 */
{
    int i,j;
	
    for (i = 0; i < nrh+1; i++)	
    {
        for (j = 0; j < nch+1; j++)
		{
	        free(inputAr[i][j]);
		}
	}
    for (i = 0; i < nrh+1; i++)	free(inputAr[i]);
    free(inputAr);		
    return;
}


//***************************4D*******************************************/


float ****f4dmatrix(long nrl, long nrh, long ncl, long nch, long ndl, long ndh, long nfl, long nfh)
{
	int i,j, k;
	float ****outputAr;
	outputAr = (float ****)malloc((nrh+1)*sizeof(float ***));
	if(outputAr == NULL)
	{
		printf("No enough memory space !");
		exit(1);
	}
	for(i=0;i< nrh+1;i++)
	{
		outputAr[i] = (float ***)malloc((nch+1)*sizeof(float **));
		if(outputAr == NULL)
		{
			printf("No enough memory space!");
			exit(1);
		}
	}
	for (i=0;i<nrh+1;i++)
	{
		for(j=0;j<nch+1;j++)
		{
			outputAr[i][j] = (float **)malloc((ndh+1)*sizeof(float*));
			if(outputAr[i][j] == NULL)
			{
				printf("No enough memory space!");
				exit(1);
			}
		}
	}
	
	for(i = 0;i<nrh+1;i++)
	{
		for(j=0;j<nch+1;j++)
		{
			for(k=0;k<ndh+1;k++)
			{
				outputAr[i][j][k] = (float *)malloc((nfh+1)*sizeof(float));
				if(outputAr[i][j][k]==NULL)
				{
					printf("No enough memory space!");
					exit(1);
				}
			}
		}	
	}
	
	return outputAr;
}

void free_f4dmatrix(float ****inputAr, long nrl,long nrh, long ncl, long nch, long ndl, long ndh,long nfl, long nfh)
{
	int i,j,k;

	for(i=0;i<nrh+1;i++)
		{
		for(j=0;j<nch+1;j++)
			{
				for(k=0;k<ndh+1;k++)
					{
						free(inputAr[i][j][k]);
					}
					free(inputAr[i][j]);
			}
			free(inputAr[i]);
		}
/*
	for(i=0;i<nrh+1;i++)
		{
			for(j=0;j<nch+1;j++)
				{
					free(inputAr[i][j]);
				}
		}

	for(i=0;i<nrh+1;i++) free(inputAr[i]);*/
	free(inputAr);

}

/*-----------------------------------------------------------------------*/

double ****d4dmatrix(long nrl, long nrh, long ncl, long nch, long ndl, long ndh, long nfl, long nfh)
{
	int i,j, k;
	double ****outputAr;
	outputAr = (double ****)malloc((nrh+1)*sizeof(double ***));
	if(outputAr == NULL)
	{
		printf("No enough memory space !");
		exit(1);
	}
	for(i=0;i< nrh+1;i++)
	{
		outputAr[i] = (double ***)malloc((nch+1)*sizeof(double **));
		if(outputAr == NULL)
		{
			printf("No enough memory space!");
			exit(1);
		}
	}
	for (i=0;i<nrh+1;i++)
	{
		for(j=0;j<nch+1;j++)
		{
			outputAr[i][j] = (double **)malloc((ndh+1)*sizeof(double*));
			if(outputAr[i][j] == NULL)
			{
				printf("No enough memory space!");
				exit(1);
			}
		}
	}
	
	for(i = 0;i<nrh+1;i++)
	{
		for(j=0;j<nch+1;j++)
		{
			for(k=0;k<ndh+1;k++)
			{
				outputAr[i][j][k] = (double *)malloc((nfh+1)*sizeof(double));
				if(outputAr[i][j][k]==NULL)
				{
					printf("No enough memory space!");
					exit(1);
				}
			}
		}	
	}
	
	return outputAr;
}

void free_d4dmatrix(double ****inputAr, long nrl,long nrh, long ncl, long nch, long ndl, long ndh,long nfl, long nfh)
{
	int i,j,k;

	for(i=0;i<nrh+1;i++)
		{
		for(j=0;j<nch+1;j++)
			{
				for(k=0;k<ndh+1;k++)
					{
						free(inputAr[i][j][k]);
					}
					free(inputAr[i][j]);
			}
			free(inputAr[i]);
		}
/*
	for(i=0;i<nrh+1;i++)
		{
			for(j=0;j<nch+1;j++)
				{
					free(inputAr[i][j]);
				}
		}

	for(i=0;i<nrh+1;i++) free(inputAr[i]);*/
	free(inputAr);

}



double *Ones(int Tnum)
{
  double *One;
  One = dvector(1,Tnum);
  for(int i=0;i<Tnum;i++)
  {
	  One[i] = 1.0;
  }
  return One;
}


//void fifthLoop(double OmegaR,double Omega, double MaX, double MaY, double MaZ, double ka, double **pF,double complex SP1,double complex SP2,double complex SP3)
void fifthLoop(double OmegaR,double Omega, double MaX, double MaY, double MaZ, double ka)
{
  double *DataXR, *DataYR, *DataZR, *DOrX, *DOrY,  *DOrZ, *DOr, *DORStar, *DOR,\
 		 *DORStarX, *DORStarY, *DORStarZ, *DORX, *DORY,*DORZ, *RGamma, *Vx, *Vy, *Vz, *Q,\
 			 *Lx, *Ly, *Lz, *FxM, *FyM, *FzM, *FxP, *FyP, *FzP, *FRStarM, *FRStarP, *FRM, *FRP;
  DataXR = dvector(1,TNum);
	DataYR = dvector(1,TNum);
	DataZR = dvector(1,TNum);

  
  DOrX=dvector(1,TNum);
  DOrY=dvector(1,TNum);
  DOrZ=dvector(1,TNum);
  DOr=dvector(1,TNum);
	

  
  DORStar=dvector(1,TNum);
  DOR=dvector(1,TNum);
  DORStarX=dvector(1,TNum);
  DORStarY=dvector(1,TNum);
  DORStarZ=dvector(1,TNum);
	
  
  DORX=dvector(1,TNum);
  DORY=dvector(1,TNum);
  DORZ=dvector(1,TNum);
  RGamma=dvector(1,TNum);
  
  Vx=dvector(1,TNum);
  Vy=dvector(1,TNum);
  Vz=dvector(1,TNum);
  Q =dvector(1,TNum);
  
  Lx=dvector(1,TNum);
  Ly=dvector(1,TNum);
  Lz=dvector(1,TNum);
	
  FxM=dvector(1,TNum);
  FyM=dvector(1,TNum);
  FzM=dvector(1,TNum);
	
  FxP=dvector(1,TNum); 
  FyP=dvector(1,TNum);
  FzP=dvector(1,TNum);
	
  FRStarM=dvector(1,TNum);
  FRStarP=dvector(1,TNum);
	
  FRM=dvector(1,TNum);
  FRP=dvector(1,TNum);
	
	double complex SP1T; double complex SP1B=0.0+0.0*I; double complex SP1E=0.0+0.0*I;
	double complex SP2T=0.0+0.0*I; double complex SP2B=0.0+0.0*I; double complex SP2E=0.0+0.0*I;
	double complex SP3T=0.0+0.0*I; double complex SP3B=0.0+0.0*I; double complex SP3E=0.0+0.0*I;
	SP1T=0.0+0.0*I;
	SP1=0.0 + 0.0*I;
	SP2=0.0 + 0.0*I;
	SP3=0.0 + 0.0*I;

	//double *One,VectorX,VectorY,VectorZ,VuN,VvN,VUN;
	double *One,VectorX,VectorY,VectorZ,VuN,VvN,VUN;
  One = Ones(TNum);
	
	
  for(int i=0;i<TNum;i++)
	{
		
		Dp = Dp*One[i];
		Drou = Drou*One[i];
		
    Gu = Du*cos(Time[i]*OmegaR+Theta)-Dv*sin(Time[i]*OmegaR+Theta);
   	Gv = Du*sin(Time[i]*OmegaR+Theta)+Dv*cos(Time[i]*OmegaR+Theta);
    Gw = Dw*One[i];
		
		
		DataXR[i] = DR*cos(OmegaR*Time[i]+atan2(DY,DX)+Theta);
    DataYR[i] = DR*sin(OmegaR*Time[i]+atan2(DY,DX)+Theta);
   	DataZR[i] = DZ*One[i];
		
		DOrX[i] = OX-DataXR[i];
		DOrY[i] = OY-DataYR[i];
		DOrZ[i] = OZ-DataZR[i];
		DOr[i] = sqrt(pow(DOrX[i],2)+pow(DOrY[i],2)+pow(DOrZ[i],2));

		//printf("DOR[%i] = %4.9f\n",i,DOR[i]);
		DORStar[i] = sqrt(pow(DOr[i],2)+pow(Gamma,2)*pow((MaX*DOrX[i]+MaY*DOrY[i]+MaZ*DOrZ[i]),2))/Gamma;
		DOR[i] = pow(Gamma,2)*(DORStar[i]-(MaX*DOrX[i]+MaY*DOrY[i]+MaZ*DOrZ[i]));
		DORStarX[i] = (DOrX[i]+pow(Gamma,2)*(MaX*DOrX[i]+MaY*DOrY[i]+MaZ*DOrZ[i])*MaX)/(pow(Gamma,2)*DORStar[i]);
		DORStarY[i] = (DOrY[i]+pow(Gamma,2)*(MaX*DOrX[i]+MaY*DOrY[i]+MaZ*DOrZ[i])*MaY)/(pow(Gamma,2)*DORStar[i]);
		DORStarZ[i] = (DOrZ[i]+pow(Gamma,2)*(MaX*DOrX[i]+MaY*DOrY[i]+MaZ*DOrZ[i])*MaZ)/(pow(Gamma,2)*DORStar[i]);
		
		
		DORX[i] = pow(Gamma,2)*(DORStarX[i]-MaX);
		DORY[i] = pow(Gamma,2)*(DORStarY[i]-MaY);
		DORZ[i] = pow(Gamma,2)*(DORStarZ[i]-MaZ);
		
		
		RGamma[i] = Time[i]+DOR[i]/C_0;
		//printf("RGamma[%i] = %4.9f\n",i,RGamma[i]);
		
		Vx[i] = -DR*OmegaR*sin(OmegaR*Time[i]+atan2(DY,DX)+Theta);
		Vy[i] = DR*OmegaR*cos(OmegaR*Time[i]+atan2(DY,DX)+Theta);
    Vz[i] = 0.0*One[i];
		
		VectorX = DVX*cos(OmegaR*Time[i]+Theta)-DVY*sin(OmegaR*Time[i]+Theta);
		VectorY = DVX*sin(OmegaR*Time[i]+Theta)+DVY*cos(OmegaR*Time[i]+Theta);
		VectorZ = DVZ*One[i];
		
		VuN = Gu*VectorX+Gv*VectorY+Gw*VectorZ;
   	VvN = Vx[i]*VectorX+Vy[i]*VectorY+Vz[i]*VectorZ;
    VUN = C_0*(MaX*VectorX+MaY*VectorY+MaZ*VectorZ);
		
		
    Q[i] = A*Drou*(VuN-(VvN-VUN))+Rou_0*A*(VvN-VUN);
    
   	Lx[i] = A*Drou*Gu*(VuN-(VvN-VUN))+A*Dp*VectorX;
   	Ly[i] = A*Drou*Gv*(VuN-(VvN-VUN))+A*Dp*VectorY;
   	Lz[i] = A*Drou*Gw*(VuN-(VvN-VUN))+A*Dp*VectorZ;
		
		/*
		Q[i] = cos(OmegaM*Time[i])*A;
		
		Lx[i] = cos(OmegaM*Time[i])*cos(OmegaR*Time[i]+atan2(DY,DX)+Theta)*A;
		Ly[i] = cos(OmegaM*Time[i])*sin(OmegaR*Time[i]+atan2(DY,DX)+Theta)*A;
		Lz[i] = A*One[i]; */

			
		FxM[i] = Lx[i];
		FyM[i] = Ly[i];
		FzM[i] = Lz[i];
		
		FxP[i] = -C_0*MaX*Q[i]; // here MaX, MaY, and MaZ are zero, is given from the file so the remaining result is zero
		FyP[i] = -C_0*MaY*Q[i];
		FzP[i] = -C_0*MaZ*Q[i];
		
		FRStarM[i] = FxM[i]*DORStarX[i]+FyM[i]*DORStarY[i]+FzM[i]*DORStarZ[i];
		FRStarP[i] = FxP[i]*DORStarX[i]+FyP[i]*DORStarY[i]+FzP[i]*DORStarZ[i];
		
		FRM[i] = FxM[i]*DORX[i]+FyM[i]*DORY[i]+FzM[i]*DORZ[i];
		FRP[i] = FxP[i]*DORX[i]+FyP[i]*DORY[i]+FzP[i]*DORZ[i];
		
		//--------------following the time integration --------------//
		

		
		SP1T += DT*(z1*Omega*Q[i]/DORStar[i])*cexp(-1*z1*Omega*RGamma[i]);
		//SP1T += cexp(-1*z1*RGamma[i]);
		SP1B = (DT/2)*(z1*Omega*Q[0]*cexp(-1*z1*Omega*RGamma[0])/DORStar[0]);
		SP1E = (DT/2)*(z1*Omega*Q[TNum-1]*cexp(-1*z1*Omega*RGamma[TNum-1])/DORStar[TNum-1]);
		SP1 = SP1T-SP1B-SP1E;
		
		SP2T += DT*((z1*ka*FRM[i]/DORStar[i]+FRStarM[i]/pow(DORStar[i],2))*cexp(-1*z1*Omega*RGamma[i]));
		SP2B = (DT/2)*(z1*ka*FRM[0]/DORStar[0]+FRStarM[0]/pow(DORStar[0],2))*cexp(-1*z1*Omega*RGamma[0]);
		SP2E = (DT/2)*(z1*ka*FRM[TNum-1]/DORStar[TNum-1]+FRStarM[TNum-1]/pow(DORStar[TNum-1],2))*cexp(-1*z1*Omega*RGamma[TNum-1]);
		SP2 = SP2T-SP2B-SP2E;
		
		SP3T += DT*(z1*ka*FRP[i]/DORStar[i]+FRStarP[i]/pow(DORStar[i],2))*cexp(-1*z1*Omega*RGamma[i]);
		SP3B = (DT/2)*(z1*ka*FRP[0]/DORStar[0]+FRStarP[0]/pow(DORStar[0],2))*cexp(-1*z1*Omega*RGamma[0]);
		SP3E = (DT/2)*(z1*ka*FRP[TNum-1]/DORStar[TNum-1]+FRStarP[TNum-1]/pow(DORStar[TNum-1],2))*cexp(-1*z1*Omega*RGamma[TNum-1]);
		SP3 = SP3T-SP3B-SP3E;
		
	  
    }
    free_dvector(One,1,TNum);
		
		free_dvector(DataXR,1,TNum);
		free_dvector(DataYR,1,TNum);
		free_dvector(DataZR,1,TNum);
  
		free_dvector(DOrX,1,TNum);
		free_dvector(DOrY,1,TNum);
		free_dvector(DOrZ,1,TNum);
		free_dvector(DOr,1,TNum);

		free_dvector(DORStar,1,TNum);
		free_dvector(DOR,1,TNum);
		free_dvector(DORStarX,1,TNum);
		free_dvector(DORStarY,1,TNum);
		free_dvector(DORStarZ,1,TNum);
	

	  free_dvector(DORX,1,TNum);
	  free_dvector(DORY,1,TNum);
	  free_dvector(DORZ,1,TNum);
	  free_dvector(RGamma,1,TNum);
	  free_dvector(Vx,1,TNum); 
	  free_dvector(Vy,1,TNum); 
	  free_dvector(Vz,1,TNum); 
	  free_dvector(Q,1,TNum); 
	  free_dvector(Lx,1,TNum);
	  free_dvector(Ly,1,TNum);
	  free_dvector(Lz,1,TNum);
	  free_dvector(FxM,1,TNum);
	  free_dvector(FyM,1,TNum); 
	  free_dvector(FzM,1,TNum); 
	  free_dvector(FxP,1,TNum); 
	  free_dvector(FyP,1,TNum); 
	  free_dvector(FzP,1,TNum);
	  free_dvector(FRStarM,1,TNum); 
	  free_dvector(FRStarP,1,TNum);
	  free_dvector(FRM,1,TNum); 
	  free_dvector(FRP,1,TNum);
}


