#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "printMat.h"
#include "findCovar.h"
double eVals[3];


void cubic(double* A, double* B, double* C, double* D)
{
	double* f = NULL;
	double* g = NULL;
	double* h = NULL;
	double* x1 = NULL;
	double* x2 = NULL;
	double* x3 = NULL;
	double* i = NULL;
	double* j = NULL;
	double* k = NULL;
	double* L = NULL;
	double* m = NULL;
	double* n = NULL;
	double* p = NULL;
	double* r = NULL;
	double* s = NULL;
	double* t = NULL;
	double* x2_i = NULL;
	double* x3_i = NULL;
	double* u = NULL;
	
	f = (double*) malloc(sizeof(double));
	g = (double*) malloc(sizeof(double));
	h = (double*) malloc(sizeof(double));
	
	if( f == NULL || g == NULL || h == NULL )
	{
		printf("Problem with memory allocation!");
	}
	
	*f = ( ( ( 3.0 * (*C / *A) ) ) - ( ( *B * (*B ) ) / ( *A * ( *A ) ) ) ) / 3.0;
	*g = ( ( ( 2.0 * *B * *B * *B) / ( *A * *A * *A ) ) - ( ( 9.0 * *B * *C) / ( *A * *A ) ) + ( ( 27.0 * *D) / *A ) ) / 27.0;
	*h = ( ( *g * *g ) / 4.0 ) + ( ( *f * *f * *f ) / 27.0 );
	
	if( *h == 0 && *f == 0 && *g == 0 )
	{
		x1 = (double*) malloc(sizeof(double));
		
		if ( x1 == NULL )
		{
			printf("Problem with memory allocation");
		}
		
		*x1 = pow( ( *D / *A ), 1.0 / 3.0) * ( -1.0 );
    eVals[0]=*x1;
    eVals[1]=*x1;
    eVals[2]=*x1;
		
		printf("All 3 Roots are: %.4lf\n", *x1);
		
		free(f);
		free(g);
		free(h);
		free(x1);
		
	}
	
	else if( *h <= 0 )
	{
		
		i = (double*) malloc(sizeof(double));
		j = (double*) malloc(sizeof(double));
		k = (double*) malloc(sizeof(double));
		L = (double*) malloc(sizeof(double));
		m = (double*) malloc(sizeof(double));
		n = (double*) malloc(sizeof(double));
		p = (double*) malloc(sizeof(double));
		x1 = (double*) malloc(sizeof(double));
		x2 = (double*) malloc(sizeof(double));
		x3 = (double*) malloc(sizeof(double));
		
		if ( i == NULL || j == NULL || k == NULL || L == NULL || m == NULL || n == NULL || p == NULL || x1 == NULL || x2 == NULL || x3 == NULL )
		{
			printf("Problem with memory allocation");
		}
		
		*i = pow( ( ( ( *g * *g ) / 4.0 ) - *h ), 0.5 );
		*j = pow( *i, 1.0 / 3.0 );
		*k = acos( - ( *g / ( 2.0 * *i ) ) );
		*L = *j * ( -1.0 );
		*m = cos( *k / 3.0 );
		*n = sqrt( 3 ) * sin( *k / 3.0 );
		*p = ( *B / ( 3 * *A ) ) * ( -1.0 );
		*x1 = 2.0 * *j * cos( *k / 3.0 ) - ( *B / ( 3.0 * *A ) );
		*x2 = *L * ( *m + *n ) + *p;
		*x3 = *L * ( *m - *n ) + *p;
		printf("Roots are: x1 = %.4lf, x2 = %.4lf, x3 = %.4lf\n", *x1, *x2, *x3);

    eVals[0]=*x1;
    eVals[1]=*x2;
    eVals[2]=*x3;
		//*x2 = *x1;
		free(f);
		free(g);
		free(h);
		free(i);
		free(j);
		free(k);
		free(L);
		free(m);
		free(n);
		free(p);
		free(x1);
		free(x2);
		free(x3);	
	}
	
}	
int main(void)
{

  int i,j,k;
  double arr1[10];
  double arr2[10];
  double arr3[10];
  double covar_mat[3][3]={0.0};
  double temp[4]={0};
  double determinant=0;
  double PCA[10][2]={0.0};
  double a,b,c;


  double data_set[10][3]={ {7,4,6},
                          {4,1,8},
                          {6,3,5},
                          {8,6,1},
                          {8,5,7},
                          {7,2,9},
                          {5,3,3},
                          {9,5,8},
                          {7,4,5},
                          {8,2,2} 
                          };

//Spliting the Data Set into 3 matrix column wise
 for ( i=0;i<10;i++){
      arr1[i]=data_set[i][0];
 }
 
 for ( i=0;i<10;i++){
      arr2[i]=data_set[i][1];
 }
 for ( i=0;i<10;i++){
      arr3[i]=data_set[i][2];
 }
//find size of each column matrix
int m = sizeof(arr1) / sizeof(arr1[0]); 
int n = sizeof(arr2) / sizeof(arr2[0]); 
int o = sizeof(arr3) / sizeof(arr3[0]);
  
//Populating Covariance matrix
   covar_mat[0][0]=covariance(arr1,arr1,m);
   covar_mat[1][1]=covariance(arr2,arr2,m);
   covar_mat[2][2]=covariance(arr3,arr3,m);

   a=covariance(arr2,arr1,m);
   covar_mat[1][0]=a;
   covar_mat[0][1]=a;

   b=covariance(arr3,arr1,m);
   covar_mat[2][0]=b;
   covar_mat[0][2]=b;
   c=covariance(arr3,arr2,m);
   covar_mat[2][1]=c;
   covar_mat[1][2]=c;
  
  printf("\nCorrelation covar_mat\n");

  print3x3(covar_mat);


  
  float d,g,e;
  double md1[2][2]={0};
  double md2[2][2]={0};
  double md3[2][2]={0};
  double eVector[3][2];
  double y=0;
  double x=0;
  double dm1=0,dm2=0,dm3=0;

              
  md1[0][0]=covar_mat[1][1];
  md1[0][1]=covar_mat[1][2];
  md1[1][0]=covar_mat[2][1];
  md1[1][1]=covar_mat[2][2];

  md2[0][0]=covar_mat[0][0];
  md2[0][1]=covar_mat[0][2];
  md2[1][0]=covar_mat[2][0];
  md2[1][1]=covar_mat[2][2];

  md3[0][0]=covar_mat[0][0];
  md3[0][1]=covar_mat[0][1];
  md3[1][0]=covar_mat[1][0];
  md3[1][1]=covar_mat[1][1];
                
  printf("\n\nCovarriance Matrix is:\n");
 for(i=0;i<3;i++){
   for(j=0;j<3;j++){
     printf("%.2f ",covar_mat[i][j]);
   }
   printf("\n");
 }

 determinant = covar_mat[0][0] * ((covar_mat[1][1]*covar_mat[2][2]) - (covar_mat[2][1]*covar_mat[1][2])) -covar_mat[0][1] * (covar_mat[1][0]
   * covar_mat[2][2] - covar_mat[2][0] * covar_mat[1][2]) + covar_mat[0][2] * (covar_mat[1][0] * covar_mat[2][1] - covar_mat[2][0] * covar_mat[1][1]);

   dm1= (md1[0][0] * md1[1][1]) - (md1[1][0] * md1[0][1]);
   //printf("\n%d",dm1);
   dm2= (md2[0][0] * md2[1][1]) - (md2[0][1] * md2[1][0]);
    //printf("\n%d",dm2);
   dm3= (md3[0][0] * md3[1][1]) - (md3[0][1] * md3[1][0]);
    //printf("\n%d",dm3);
   y=dm1+dm2+dm3;
    for(i=0;i<3;i++){
      for(j=0;j<3;j++){
        if(i==j)
            {
                x=x+covar_mat[i][j];
            }
      }
    }
  printf("\nDeterminant of 3X3 covar_mat: %.4f", determinant);

  printf("\nThe Characteristics equation is: lamda^3 - %.4f lamda^2 + %.4f lamda - %.4f = 0 \n",x,y,determinant);
  double A=1;
  double B =x*(-1);
  double C=y;
  double D =determinant*(-1);
  //finding roots of the Characteristics equation
  cubic(&A,&B,&C,&D);

 double xyz;
 for (i = 0; i < 3; ++i) 
        {
            for (j = i + 1; j < 3; ++j)
            {
                if (eVals[i] < eVals[j]) 
                {
                    xyz =  eVals[i];
                    eVals[i] = eVals[j];
                    eVals[j] = xyz;
                }
            }
    }
 printf("Smallest eigen Value will be discarted: ");
  for(i=0;i<3;i++){
    printf("%.4f ",eVals[i]);
  }

   for(int i=0;i<2;i++){

            printf("\nThe eigenvector corresponding to eigenvalue %f is:\n\n",eVals[i]);
            d=(covar_mat[0][1]*covar_mat[1][2])-((covar_mat[1][1]-eVals[i])*covar_mat[0][2]);
            e=(covar_mat[0][2]*covar_mat[1][0])-((covar_mat[0][0]-eVals[i])*covar_mat[1][2]);
            g=((covar_mat[0][0]-eVals[i])*(covar_mat[1][1]-eVals[i]))-(covar_mat[1][0]*covar_mat[0][1]);
            printf("\t|\t%f\t|\n\t|\t%f\t|\n\t|\t%f\t|\n",d,e,g);

            eVector[0][i]=d;
            eVector[1][i]=e;
            eVector[2][i]=g;
   }
   printf("\n");
   printf("\nEigen Vectors\n\n");
   print3x2(eVector);
//Computing PC
for(i=0; i<10; i++){
        for(j=0; j<2; j++){
            for(k=0; k<3; k++)
            {
                PCA[i][j]+=data_set[i][k]*eVector[k][j];
            }
        }
}
//printing PC
printf("\n\nPCA\n");
    for(i=0;i<10;i++){
      for(j=0;j<2;j++){
          printf("%.3f ",PCA[i][j]);
      }
      printf("\n");
    }
}
