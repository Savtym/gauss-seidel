/*
 * www.twitter.com/torayeff
 * 
 * Gauss-Seidel method:
 * http://math.semestr.ru/optim/methodzeidel.php
 * http://en.wikipedia.org/wiki/Gauss–Seidel_method
 * http://ru.wikipedia.org/wiki/Метод_Гаусса_—_Зейделя
 * 
 * To solve A*x = b transform to 
 * A(t)*A*x = A(t)*b <---> C*x = d 
 * (A(t) transpose of matrix A)
 * (in this state, Gauss-Seidel method always congerges)
*/
#include <iostream>
#include <cstdio>
#include <cmath>
using namespace std;

#define MAXN 100

int N; //matrix dimension
double eps = 0.001;
double err = eps;
double A[MAXN][MAXN];
double b[MAXN][MAXN];
double At[MAXN][MAXN];
double C[MAXN][MAXN];
double d[MAXN][MAXN];

//convergence check
bool converge(double *xk, double *xkp)
{
    double norm = 0;
    for(int i=0; i<N; i++) 
    {
      norm += (xk[i] - xkp[i])*(xk[i] - xkp[i]);
    }
    if(sqrt(norm)>=eps)
      return false;
    err = eps;
    return true;
}
 

int main()
{
	printf("Input dimension N: ");
	scanf("%d", &N);
	
	for(int i=0; i<N; i++)
	{
		for(int j=0; j<N; j++)
		{
			double tmp;
			printf("A[%d][%d] = ", i+1, j+1);
			scanf("%lf", &tmp);
			A[i][j] = tmp;
		}
	}
	
	for(int i=0; i<N; i++)
	{
		double tmp;
		printf("b[%d] = ", i+1);
		scanf("%lf", &tmp);
		b[i][0] = tmp;
	}
	

	//transpose A
	for(int i=0; i<N; i++)
	{
		for(int j=0; j<N; j++)
		{
			At[i][j] = A[j][i];
		}
	}
	
	//multiply At*A = C
	for(int i=0; i<N; i++)
	{
		for(int j=0; j<N; j++)
		{
			C[i][j] = 0;
			for(int k=0; k<N; k++)
			{
				C[i][j] = C[i][j] + At[i][k]*A[k][j];
			}
		}
	}
	
	//multiply At*b = d
	for(int i=0; i<N; i++)
	{
		for(int j=0; j<1; j++)
		{
			d[i][j] = 0;
			for(int k=0; k<N; k++)
			{
				d[i][j] = d[i][j] + At[i][k]*b[k][j];
			}
		}
	}
	
	//Solve Gauss-Seidel method for C*x = d
	double x[N]; //current values
	double p[N]; //previous values
	
	for(int i=0; i<N; i++)
	{
		x[i] = 0;
		p[i] = 0;
	}
	
	do
	{
		for(int i = 0; i<N; i++)
			p[i] = x[i];
			
		for(int i=0; i<N; i++)
		{
			double v = 0.0;
			for(int j=0; j<i; j++)
			{
				double cij;
				if(i==j) cij = 0;
				else cij = -1.0*C[i][j]/C[i][i];
				v += cij*x[j];
			}
			
			for(int j=i; j<N; j++)
			{
				double cij;
				if(i==j) cij = 0;
				else cij = -1.0*C[i][j]/C[i][i];
				v += cij*p[j];
			}
			v +=  1.0*d[i][0]/C[i][i];
			x[i] = v;
		}
		
	}while (!converge(x, p));
	
	//print solution
	printf("Err. val.: %lf\n", err);
	for(int i=0; i<N; i++)
	{
		printf("%lf ", x[i]);
	}
	
	return 0;
}
