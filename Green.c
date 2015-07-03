/*
 ============================================================================
 Name        : Green.c
 Author      : Ivan Syzonenko
 Version     : 1.0
 Copyright   : Free to copy for everyone
 Description : Hello World in C, Ansi-style
 ============================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

// #define _X_YES_

void allocate(	/*IN  =>*/const double a, const double b, const double *h,
				/*OUT =>*/double **x, double **u, double **u_prime, unsigned int *totiter);

void rungekutta(/*IN  =>*/const double a, const double b, const double h,
				/*OUT =>*/double *x, double *u, double *u_prime, const unsigned int totiter);

// #ifdef _X_YES_
	double p2 (const double x);
	double p1 (const double x);
	double p0 (const double x);
// #else
// 	double p2 ();
// 	double p1 ();
// 	double p0 ();
// #endif
	double RHS_fx_func(const double x);
	double f1(const double x, const double u, const double v);
	double f2(const double x, const double u, const double v);

void matrINV(double  in_matr[][2], double  out_matr[][2], const char debug);

int main(void)
{
	clock_t begin, end;
	double 	a = 0, // lower integration limit
			b = 1, // upper integration limit
			h = 0.0000001, // runge-kutta step
	// below are our boundary coefficients
			alpha0 = 1, // for lower limit near first derivative
			alpha1 = -1, // for lower limit near zero derivative
			beta0 = 0, // for upper limit near first derivative
			beta1 = 1; // for upper limit near zero derivative

	unsigned int totiter = 0;
	double *x, *u, *u_prime;
	allocate(a, b, &h, &x, &u, &u_prime, &totiter);

	double *C1 = (double *)malloc(totiter * sizeof(double));
	double *C2 = (double *)malloc(totiter * sizeof(double));
	double *P  = (double *)malloc(totiter * sizeof(double));
	double *W  = (double *)malloc(totiter * sizeof(double));

	double Am[2][2], AmINV[2][2], Rm[2], H[2];

	double partA;
	double partB;

	double *solution = (double *)malloc(totiter * sizeof(double));

	begin = clock();
	rungekutta(a, b, h, x, u, u_prime, totiter);

	// printf("u[0] = %f\n", u[0]);
	// printf("u[totiter-1] = %f\n", u[totiter-1]);
	// printf("u[totiter] = %f\n", u[totiter]);
	// printf("u[totiter + totiter - 1] = %f\n\n", u[totiter + totiter-1]);

	// printf("u_prime[0] = %f\n", u_prime[0]);
	// printf("u_prime[totiter-1] = %f\n", u_prime[totiter-1]);
	// printf("u_prime[totiter] = %f\n", u_prime[totiter]);
	// printf("u_prime[totiter + totiter - 1] = %f\n\n", u_prime[totiter + totiter-1]);


	for(int i=0; i < totiter; ++i)
		W[i] = u[i] * u_prime[totiter + i] - u[totiter + i] * u_prime[i];

	partA = beta0 * u_prime[totiter - 1] + beta1 * u[totiter - 1] ;
	partB = beta0 * u_prime[totiter + totiter - 1] + beta1 * u[totiter + totiter - 1];

	C1[0] = 0;
	C2[0] = 0;
	P[0] = 0;

	for(int i=0; i < totiter-1; ++i)
	{
		C1[i+1] = C1[i] - 0.5 * h * (   u[totiter+i]   * RHS_fx_func(x[i] )   / ( p0( x[i] )   * W[i] )
									+   u[totiter+i+1] * RHS_fx_func(x[i+1] ) / ( p0( x[i+1] ) * W[i+1] ) );

		C2[i+1] = C2[i] + 0.5 * h * (   u[i]           * RHS_fx_func(x[i] )   / ( p0( x[i] )   * W[i] )
									+   u[i+1]         * RHS_fx_func(x[i+1] ) / ( p0( x[i+1] ) * W[i+1] ) );

		P[i+1]  = P[i]  + 0.5 * h * ( ( u[totiter+i]     * partA - u[i] * partB) *  RHS_fx_func(x[i]  ) / ( p0( x[i]   ) * W[i]   )
									+ ( u[totiter + i+1] * partA - u[i+1] * partB) *  RHS_fx_func(x[i+1]) / ( p0( x[i+1] ) * W[i+1] ) );
	}

//	printf("C1[0] = %5.10f\n", C1[0]);
//	printf("C1[totiter-1] = %5.10f\n", C1[totiter-1]);
//	printf("C2[0] = %5.10f\n", C2[0]);
//	printf("C2[totiter-1] = %5.10f\n", C2[totiter-1]);
//	printf("P[0] = %5.10f\n", P[0]);
//	printf("P[totiter-1] = %5.10f\n", P[totiter-1]);




	Am[0][0] = alpha0 * u_prime[0]                   + alpha1 * u[0];
	Am[0][1] = alpha0 * u_prime[totiter]             + alpha1 * u[totiter];
	Am[1][0] = beta0  * u_prime[totiter-1]           + beta1  * u[totiter-1];
	Am[1][1] = beta0  * u_prime[totiter + totiter-1] + beta1  * u[totiter + totiter-1];

	Rm[0] = 0;
	Rm[1] = P[totiter-1];

	matrINV(Am, AmINV, 0);

//	printf("%f\t%f\n", AmINV[0][0], AmINV[0][1]);
//	printf("%f\t%f\n", AmINV[1][0], AmINV[1][1]);

	H[0] = AmINV[0][0] * Rm[0] + AmINV[0][1] * Rm[1];
	H[1] = AmINV[1][0] * Rm[0] + AmINV[1][1] * Rm[1];


	for(int i=0; i < totiter; ++i)
	     solution[i] = ( H[0] + C1[i] ) * u[i] + ( H[1] + C2[i] ) * u[totiter+i];
	end = clock();
	printf("Total time %f\n", (double)(end - begin) / CLOCKS_PER_SEC);
	// DONE WITH COMPUTATION. GOOD TIME TO CHECK THE RESULT

	printf("Now checking result against function: %s\n\
			recommended input is :\n\
			p2 = %s\n\
			p1 = %s\n\
			p0 = %s\n\
			RHS_fx_func = %s\n",
			"0.25 * (2 * x[i] * x[i] - x[i] - 1)",
			"1", "4", "8", "6.75 + (3.75 + 0.5 * x) * x");
	double *check = (double*)malloc(totiter * sizeof(double));
	for (int i = 0; i < totiter; ++i)
	{
		check[i] = 0.25 * (2 * x[i] * x[i] - x[i] - 1);
	}

	double total_error = 0;

	for (int i = 0; i < totiter; ++i)
	{
		total_error += fabs(check[i] - solution[i]);
		//printf("%5.10f - %5.10f = %5.10f\n", check[i], solution[i], check[i] - solution[i]);
	}

	printf("Total error is : %5.5f\n", total_error);

	printf("Done with all calculations. Exiting...\n");

	return EXIT_SUCCESS;
}


void rungekutta(const double a, const double b, const double h,
			/*OUT ->*/	 double *x, double *u,  double *u_prime, const unsigned int totiter)
{

	double *v = u_prime;
	x[0] = a;
	u[0] = 0;
	v[0] = 1;
	u[totiter] = 1;
	v[totiter] = 0;

	double 	k[4], l[4];

	for(int j = 0; j < 2; ++j)
	{
		int cur_idx_1 = j*totiter; //will be optimized out
		for(int i=0; i < totiter - 1; ++i)
		{
			int cur_idx_2 = cur_idx_1 + i;//will be optimized out
			k[0] = h * f1(x[i], u[cur_idx_2] , v[cur_idx_2]);
			l[0] = h * f2(x[i], u[cur_idx_2] , v[cur_idx_2]);

			k[1] = h * f1(x[i] + h/2.0, u[cur_idx_2] + k[0]/2.0 , v[cur_idx_2] + l[0]/2.0);
			l[1] = h * f2(x[i] + h/2.0, u[cur_idx_2] + k[0]/2.0 , v[cur_idx_2] + l[0]/2.0);

			k[2] = h * f1(x[i] + h/2.0, u[cur_idx_2] + k[1]/2.0 , v[cur_idx_2] + l[1]/2.0);
			l[2] = h * f2(x[i] + h/2.0, u[cur_idx_2] + k[1]/2.0 , v[cur_idx_2] + l[1]/2.0);

			k[3] = h * f1(x[i] + h, u[cur_idx_2] + k[2], v[cur_idx_2] + l[2]);
			l[3] = h * f2(x[i] + h, u[cur_idx_2] + k[2], v[cur_idx_2] + l[2]);

			u[cur_idx_2+1] = u[cur_idx_2] + (k[0] + 2*(k[1] + k[2]) + k[3])/6.0;
			v[cur_idx_2+1] = v[cur_idx_2] + (l[0] + 2*(l[1] + l[2]) + l[3])/6.0;
		}
	}
}


void allocate(const double a, const double b, const double *h, double **x, double **u, double **u_prime, unsigned int *totiter)
{
	if (*h != 0)
	{
		*totiter = (unsigned int)round(abs(b - a)/ *h);
	}
	else
	{
		printf("Attempted to divide by zero\nExiting...\n\n");
		exit(EXIT_FAILURE);
	}

	*x		 = (double*)malloc( *totiter * sizeof(double));
	*u		 = (double*)malloc( *totiter * 2 * sizeof(double));
	*u_prime = (double*)malloc( *totiter * 2 * sizeof(double));

	for (int i = 0; i < *totiter-1; ++i)
		(*x)[i+1] = (*x)[i] + *h;
}

double f1(const double x, const double u, const double v)
{
	return v;
}


double f2(const double x, const double u, const double v)
{
    return -( p1(x) * v + p2(x) * u )/p0(x);
}

// #ifdef _X_YES_
double p2 (const double x)
// #else
// double p2 ()
// #endif
{
	return 1;
}

// #ifdef _X_YES_
double p1 (const double x)
// #else
// double p1 ()
// #endif
{
	return 4;
}

// #ifdef _X_YES_
double p0 (const double x)
// #else
// double p0 ()
// #endif
{
	return 8;
}

// #ifdef _X_YES_
double RHS_fx_func(const double x)
// #else
// double RHS_fx_func()
// #endif
{
	return 6.75 + (3.75 + 0.5 * x) * x;
}

void matrINV(double  in_matr[][2], double  out_matr[][2], const char debug)
// void matrINV(const double* const * in_matr, double** out_matr, const char debug)
{
	
	double D = in_matr[0][0]*in_matr[1][1] - in_matr[0][1]*in_matr[1][0];

	if ( D == 0)
	{
		printf("Division by zero during inverse matrix calculations...\n");
		exit(EXIT_FAILURE);
	}
	
	out_matr[0][0] = in_matr[1][1] / D;
	out_matr[0][1] = -in_matr[0][1] / D;
	out_matr[1][0] = -in_matr[1][0] / D;
	out_matr[1][1] = in_matr[0][0] / D;


	if (debug == 1)
	{
		printf("DEBUG print of matrINV:\n");

		for (int i = 0; i < 2; i++)
			for (int j = 0; j < 2; j++)
			{
				if (j == 0)
					printf("\n");
				printf("%f  \t", out_matr[i][j]);
			}

		printf("\nThis is the last line(matrINV).\n");
	}
}
