#include "math.h"

#include <gsl/gsl_math.h>

#include <gsl/gsl_multimin.h>

#include <gsl/gsl_math.h>

#include <stdlib.h>

#include <stdio.h>


#define M 300



#define pop_dim  100

//int num_times = pop_dim;

#define num_bins 100000


#define my_inf 1e20



static double time[pop_dim];
static long double q[M][M][pop_dim];
static long double r[M][pop_dim];
static long double product[M][M][pop_dim];
static long double product2[M][M][pop_dim];
static long double product3[M][pop_dim];
static long double product4[M][M][pop_dim];

static long double C1[M][M][pop_dim];
static long double Zfirst1[M][M][pop_dim];
static long double Zsecond1[M][M][pop_dim];

static long double C2;
static long double Zfirst2[pop_dim];
static long double Zsecond2[pop_dim];

static long double H1[M][pop_dim];
static long double Yfirst1[M][pop_dim];
static long double Ysecond1[M][pop_dim];

static long double H2;
static long double Yfirst2[pop_dim];
static long double Ysecond2[pop_dim];

static long double I1[M][pop_dim];
static long double Wfirst1[M][pop_dim];
static long double Wsecond1[M][pop_dim];

static long double I2;
static long double Wfirst2[pop_dim];
static long double Wsecond2[pop_dim];

static long double G_exp[M][pop_dim];
static long double G_exp_2[M][pop_dim];

static long double B1[M][pop_dim];
static long double B2[M][pop_dim];
static long double Z1[M][pop_dim];
static long double Z2[M][pop_dim];

static long double exp_part_L[M];
static long double exp_part_other[M][M][pop_dim];



/* Getting longest match vector from a file */




// mu is the per site per generation mutation rate

//   static double mu = .000000025;
static long double mu = .000000025;

// large stiffness biases toward more constant population size
// see my_f and my_df

static double stiffness = 0.001;

// complexity is the number of terms k>j for which we compute 
// Prob( <L at kth coal given = L at jth coal) accurately

static int complexity = 1;

static long double w[M][M];
static long double new_w[M][M][M];
static long double A[M];
//static int num_haps = M;
//static int pairs = M*(M-1) / 2;


int num_loci;
static int first[10000]; // ={ [ 0 ... 9999 ] = 1 };
static long double rho[10000];









//static int num_bins = 1000;
static double bin_bound[num_bins]; // = { [0 ... 999] = 0.0};
static double bin_center[num_bins]; // = { [0 ... 999] = 0.0};
static int bin_counter[num_bins]; // = { [0 ... 999] = 0};

static double quant_list[1000];
static double quant_values[1000];
static int num_quantiles;







void make_times(int num_haps, int num_times)
{


double pairs = num_haps * (num_haps -1) / 2.0;




    int j = 0;

    time[0] = 0.0;

//printf("time(%d) = %f\n", j, time[0]);

    for (j=1; j<num_times; j++)

    {
time[j] = 1.0 *  j / num_times * 2.0 *  (10000 / pairs); 
//time[j] = (pow(j, 2) / pow(num_times, 2)) * 2.0 * (20000 / pairs); 
//    time[j] =  (1.0 / pairs) * 5000 * pow(1.5, j);

//time[1] = 300.0;

 	printf("time(%d) = %f \n", j, time[j]);

    }

if(num_times == 70)
{
time[1] = 10;
time[2] = 15;
time[3] = 20;
time[4] = 25;
time[5] = 30;
time[6] = 40;
}

if(num_times == 120)
{
time[1] = 8;
time[2] = 12;
time[3] = 16;
time[4] = 20;
time[5] = 24;
time[6] = 28;
time[7] = 32;
time[8] = 36;
time[9] = 40;
time[10] = 46;
time[11] = 52;
}



//time[1] = 2000.0;

}


void make_bins(int num_haps)
{
int k =0;
double pairs = num_haps * (num_haps -1) / 2.0;


for(k=0; k < 10; k++)
{

bin_bound[k] =  pairs * 1.0 * (k+1) ;

bin_center[k] = bin_bound[k];
//printf("bin_center[%d] = %f \n", k, bin_center[k]);
}



for(k=10; k < num_bins; k++)

	{

		bin_bound[k] = bin_bound[k-1] * 1.12 ;//1.048808848;

			if( bin_bound[k] > pow(10,10) ){return;}
		bin_center[k] = bin_bound[k];
		//bin_center[k] = pairs * bin_bound[k-1] * 1.01 ; //* 1.048808848;
		//if (k < 500){printf("bin_center[%d] = %f \n", k, bin_center[k]);}
	}

}



void make_custom_bins(int num_haps)
{
int k =0;
double pairs = num_haps * (num_haps -1) / 2.0;

//printf("first(%d) = %d", k, first[k]);
for(k=0; k < num_loci; k++)
{



bin_center[k] = 1.0 * first[k];
if(bin_center[k] > 0) {bin_counter[k] = 1;}
else{bin_counter[k] = 0;}
//printf("bin_center[%d] = %f \n", k, bin_center[k]);
}

}

void precompute(int num_haps, int num_times)
{

int i;
int j; 
int k;
int l;
//int m;

for(i=1; i<num_haps; i++)
{
	A[i] = (num_haps - i + 1) * (num_haps - i) / 2;

}

for(j=1; j<num_haps; j++)
{
    for(i=1; i<=j; i++)
    {
        
        long double temp =1;
        for (k =1; k <=j; k++)
        {
            if (k != i) {temp = temp * A[k] / (A[k] - A[i]) ; }
        }
       w[j][i] = temp;
	
    }

}


for(j=1; j<num_haps; j++)
{
long double test = 0;
	 for(i=1; i<=j; i++)
	{
	test = test + w[j][i];
	}

}





int t;
for(j=1; j<num_haps; j++)
{
	for (k =j+1; k < num_haps; k++)
	{
		for(l = j+1; l <=k; l++)
		{	
			long double new_temp = 1;
			for(t = j; t<= k; t++)
			{
				if(t != i) {new_temp = new_temp * A[t] / (A[t] - A[i]) ;}
			}

	
			new_w[j][k][l] = new_temp;
		}
	}

}



}


void precompute_L_ind_quantities
(int num_haps, int num_times, double T[pop_dim],  double pop_size[num_times] )
{
int m;
double lambda[num_times];
for (m=0; m< num_times; m++)  {lambda[m] = (1 / pop_size[m]);}

//  Should move q and r outside of PHI because they do not depend on L
// This could be expensive to recompute for every L
int i;
int l;
int j;

for(i = 1; i < num_haps; i++)
{
	for(l = 1; l < num_haps; l++)
	{
		double temp_q = 1;
		q[i][l][0] = 1;
		for(j = 1; j < num_times; j++)
		{
			temp_q = temp_q * exp(-(A[i]-A[l]) * lambda[j-1] * ( T[j] - T[j-1] ) );
			q[i][l][j] = temp_q * exp( (A[i]-A[l]) * lambda[j] * T[j]);
		}
		
		double temp_r = 1;
		r[l][0] = 1;
		for(j = 1; j < num_times; j++)
		{
			temp_r = temp_r * exp(-A[l] * lambda[j-1] * ( T[j] - T[j-1] ) );
			r[l][j] = temp_r * exp( A[l] * lambda[j] * T[j]);
		}
	}
}

for(i = 1; i < num_haps; i++)
{
	for(l = 1; l < num_haps; l++)
	{
		for (m=0; m< num_times; m++)
		{

//  exp_part_other and product are just here to remove some baggage from inside PHI
// exp_part_other removes unneeded exponentiation, product is less significant

		exp_part_other[i][l][m] = exp(-T[m] * lambda[m] * (A[i] - A[l]) );

		product[i][l][m] = lambda[m] * q[i][l][m];
		product2[i][l][m] =  product[i][l][m] * exp_part_other[i][l][m];
		product3[l][m] = lambda[m] * r[l][m];
	//	product4[i][l][m] = product[i][l][m] * Z1[i][m];


// C1 is the L-independent part of C in PHI
// we precompute to remove all exponentiation from PHI
// the use of C1 inside an exponential is hidden in Z which is in PHI


		C1[i][l][m] = lambda[m] * (A[i] - A[l]);
		Zfirst1[i][l][m] = exp( -T[m] * C1[i][l][m] );
		if(m == num_times -1){ Zsecond1[i][l][m] =0; }
		else{  Zsecond1[i][l][m] = exp( -T[m+1] * C1[i][l][m] );  }
		}
	}
}

for(i = 1; i < num_haps; i++)
{
	for (m=0; m< num_times; m++)
		{

		H1[i][m] = lambda[m] * A[i] ;
		Yfirst1[i][m] = exp( -T[m] * H1[i][m] );
		if(m == num_times -1){ Ysecond1[i][m] =0; }
		else{  Ysecond1[i][m] = exp( -T[m+1] * H1[i][m] );  }

		I1[i][m] = lambda[m] * A[i] ;
		Wfirst1[i][m] = exp( -T[m] * I1[i][m] );
		if(m == num_times -1){ Wsecond1[i][m] =0; }
		else{  Wsecond1[i][m] = exp( -T[m+1] * I1[i][m] );  }

		if(m>0){G_exp[i][m] = exp( - A[i] * lambda[m-1] * ( T[m]- T[m-1] )   );}
		G_exp_2[i][m] = exp ( A[i] * T[m] * lambda[m] ) ;
	

		}
}


}



long double G(int num_haps, int num_times,int i, double L, double T[pop_dim], double lambda[num_times], long double var1 )
{
//printf("var1 = %.10Lf\n", var1);
// G is like F except it is for prob( < L given T ) whereas F is for prob( = L given T)

double sum =0;
double prefix[num_times];
prefix[0] = 1;
//printf("pop 0 = %f\n", lambda[0]);
double C[num_times];
double new_mess[num_times];
double term[num_times];
// The main computation
int j=0;
for (j=0; j< num_times; j++)
{

double temp_j;
double temp_j_plus_1;

    C[j] = 2*var1*L + A[i] * lambda[j];

 //if(j==3){printf("\t\t\tC(%d) = %f\t%f\t\t%f\n", j, C[j], (1 / pop_size[j]), pop_size[j] );}
	if(j==0)
	{
	//printf("TESTING %f\n", 1.9);
	prefix[j] = 1 ;
	}

    	else
	{
	//printf("PREFIX(%d) =  %f\n",j,prefix[j-1]);
	// prefix[j] = prefix[j-1] * exp( - A[i] * lambda[j-1] * ( T[j]- T[j-1] )   );
	prefix[j] = prefix[j-1] *G_exp[i][j];
	}

    if(j == num_times - 1)
    {// printf("Hello!!! %f\n", 2.8);
        new_mess[j] = exp(-C[j]*T[j]) * (1 / C[j] );
	//new_mess[j] = Yfirst1[i][j] * Yfirst2[j]* (1 / C[j] );
    }
    else
    {//printf("temp = %f\n", (-C[j]*T[j]));
        new_mess[j] = exp(-C[j]*T[j]) * (1 / C[j] )   - exp(-C[j]*T[j+1]) * (1 / C[j] ) ;

	  //new_mess[j] = (Yfirst1[i][j] * Yfirst2[j] - Ysecond1[i][j] * Ysecond2[j] ) * (1 / C[j] ) ;

    }

// the only difference between G and F is the lack of a 2*mu*T, which shows up in the next line 
// with the 2*mu removed, and the new_mess instead of mess

    //	term[j] =    lambda[j] * exp ( A[i] * T[j] * lambda[j] ) * prefix[j] * new_mess[j];
	term[j] =    lambda[j] * G_exp_2[i][j] * prefix[j] * new_mess[j];
    sum = sum + term[j];
//printf("sum_G = %f\n", sum);
}  // end loop in j 
//printf("temp = %f\t%f\t%f\n",  C[4], 2*mu*L,  A * (1.0 / pop_size[4]) );
//printf("mess(%d)=%f\n",j, new_mess[j]);


double new_sum = A[i]*sum;

//printf(" new_sum = %f\n", new_sum);

return new_sum;
}














long double F(int num_haps, int num_times, int i, double L, double T[pop_dim], double lambda[num_times], long double var1, long double var2)
{

double sum =0;
double prefix[num_times];
prefix[0] = 1;

double C[num_times];
double mess[num_times];
double term[num_times];


// The main computation
int j=0;
for (j=0; j< num_times; j++)
{

double temp_j;
double temp_j_plus_1;



	//double lambda = (1 / pop_size[j]);
    C[j] = 2 * var1 * L + A[i] *  lambda[j];
//printf("A = %Lf\n", A[i]);
	//double C_squared_inverse = (1 /  C[j]*C[j] );
	//double ct = C[j]*T[j];
//printf("L = %f\tC(%d) = %f\t\t pop_size = %f\n", L, j, C[j], pop_size[j]);
	if(j==0)
	{
	//printf("TESTING %f\n", 1.9);
	prefix[j] = 1 ;
	}

    	else
	{
	//printf("PREFIX(%d) =  %f\n",j,prefix[j-1]);
	 //prefix[j] = prefix[j-1] * exp( - A * lambda[j-1] * ( T[j]- T[j-1] )   );
	prefix[j] = prefix[j-1] * G_exp[i][j];
	}

    if(j == num_times - 1)
    {// printf("Hello!!! %f\n", 2.8);
        mess[j] = exp(-C[j]*T[j]) * ( C[j]*T[j] + 1 ) * (1 /  (C[j]*C[j]) )  ;
 //    mess[j] = Yfirst1[i][j] * Yfirst2[j] * ( C[j]*T[j] + 1 ) * (1 /  (C[j]*C[j]) )  ;
    }
    else
    {//printf("temp = %f\n", (-C[j]*T[j]));
        mess[j] = (exp(-C[j]*T[j]) * ( C[j]*T[j] + 1 ) - exp(-C[j]*T[j+1]) * ( C[j]*T[j+1] + 1 ) )
		* (1 /  (C[j]*C[j]) ) ;

//  mess[j] = ( Yfirst1[i][j] * Yfirst2[j] * ( C[j]*T[j] + 1 ) - Ysecond1[i][j] * Ysecond2[j] * ( C[j]*T[j+1] + 1 ) )	* (1 /  (C[j]*C[j]) ) ;
    }

 	//term[j] =   2*mu * A[i] *   lambda[j] * exp ( A * T[j] * lambda[j] ) * prefix[j] * mess[j];
 	term[j] =   2*var2 * A[i] *   lambda[j] * G_exp_2[i][j] * prefix[j] * mess[j];
	//printf("mess(%d)=%f\n",j, mess[j]);
    sum = sum + term[j];
	

}  // end loop in j 

//double new_sum = sum *  2*mu * A;

//printf("new_sum = %f\n", sum);
//printf("old_sum = %f\n", compute_likelihood( L, T, pop_size));

return sum;
}







double compute_prob(int num_haps, int num_times, double L, double T[pop_dim], double pop_size[num_times], double lambda[num_times])
{

long double RHO = .000000025;

int marker = 0;
int i=0;
// int k = 0;
// int j = 0;
// int l = 0;
for(i=0; i<num_times; i++)
{
//printf("pop_size[%d] = %f\n",i, pop_size[i]);
if (pop_size[i] < 100)
{
 marker = 1;
//	printf("MARKER = %f\n", 1.0);
return(10000000000000);
}
}

//printf("lambda = %f\n", lambda[0]);

//printf("mu = %.10Lf\tRHO = %.10Lf\n", mu, RHO);

long double f = F(num_haps, num_times,  1, L, T, lambda, mu + RHO, mu  );
long double g = G(num_haps, num_times,  1, L, T, lambda, mu + RHO);
//num_steps

//long double varphi[


int num_steps = 10;
double delta = L / num_steps;
long double new_total = 0;
long double new_new_total = 0;
double x = 0;
int j =0;
int m = 0;
long double F_vector[num_steps];
long double G_vector[num_steps];
for(j =1; j< num_steps; j++)
{

F_vector[j] = F(2, num_times, 1, j * delta + 1,  T, lambda, mu+RHO, mu);
G_vector[j] = G(2, num_times, 1, j * delta + 1,  T, lambda, mu+RHO);

}
//long double total = 0;
long double total_F_F = 0;
long double total_F_G= 0;
long double R_1_total = 0;
long double S_1_total = 0;
long double F_conv_F_vector[num_steps];
long double F_conv_G_vector[num_steps];

for(j =1; j< num_steps; j++)
        {

	// one copy of F changed to G to compute P(>L) rather than P(=L)

        //R_1_total = R_1_total + delta * F_vector[j] * F_vector[num_steps - j];
 	R_1_total = R_1_total + delta * F_vector[j] * G_vector[num_steps - j];
        }
//total = 0 ;


for(j =1; j< num_steps; j++)
{ 
	total_F_F = 0 ;
	total_F_G = 0 ;
	for(m =1; m< j; m++)
	{
		total_F_F = total_F_F + delta * F_vector[m] * F_vector[j  - m];
		total_F_G = total_F_G + delta * F_vector[m] * G_vector[j  - m];
		//printf("FconvF (%d) = %f\n", 
	}
	F_conv_F_vector[j] = total_F_F;
	F_conv_G_vector[j] = total_F_G;
}

long double R_2_total = 0;
for(j =1; j< num_steps; j++)
        {
	//printf("now  = %Lf\n", R_2_total);
        //R_2_total = R_2_total + delta* F_conv_F_vector[j] * F_vector[num_steps  - j];
	R_2_total = R_2_total + delta* F_conv_F_vector[j] * G_vector[num_steps  - j];
        }

long double R_3_total = 0;
for(j =1; j< num_steps; j++)
        {
	//R_2_total = R_2_total;
        //R_3_total = R_3_total + delta* F_conv_F_vector[j] * F_conv_F_vector[num_steps  - j];
	R_3_total = R_3_total + delta* F_conv_F_vector[j] * F_conv_G_vector[num_steps  - j];
        }



//printf("g = %Lf\n", g);
//return f  + R_1_total + R_2_total + R_3_total;
//return g;
return g  + R_1_total + R_2_total + R_3_total;
 
//return prob_L_max;




}




double Joint_like(int num_times, double T[pop_dim], double pop_size[num_times], double p[2])
{
int num_haps = p[0];
//int num_times = p[1];
double lambda[num_times];
int marker = 0;
int i=0;


for(i=0; i<num_times; i++)
{

//printf("pop_size[%d] = %f\n",i, pop_size[i]);
 if (pop_size[i] < 100) //&& pop_size[i] < 100000000){}
// else{ }//break;   }
{
 marker = 1;
//	printf("MARKER = %f\n", 1.0);
return(10000000000000);
}

lambda[i] = (1 / pop_size[i]);
}


precompute_L_ind_quantities( num_haps, num_times, T, pop_size );



double joint = 0;
// int l=0;
//     for(l=0; l<num_bins; l++)
//     {
// 	int m;
// 	for(m = 0; m < num_times; m++)
// 	{	
// 			exp_part_L[m] = exp( T[m] * mu * bin_center[l] );
// 	}

int l=0;
    for(l=0; l< num_loci; l++)
    {
	if(bin_center[l] > 1)
	{

	//	printf("asdasdfbin_center = %f", bin_center[l]);
		int m;
		for(m = 0; m < num_times; m++)
		{	
				exp_part_L[m] = exp( T[m] * mu * bin_center[l] );
		}
	}
    }

double temp = 0;
double tile = 0;
double total = 0;
for(i=0; i< num_quantiles; i++)
{
	temp =compute_prob(num_haps, num_times, quant_values[i], T, pop_size, lambda);
	tile = quant_list[i];
	total = total +  pow( (1 - temp) - tile, 2);
	//if (i % 100 == 0) { printf("computed_quant_values(%d) = %f \t tile = %f\n", i, 1- temp, tile); }
}
total = sqrt(total / num_quantiles);


 //double temp1 = compute_prob(num_haps, num_times, 236.0, T, pop_size, lambda);
//double temp2 = compute_prob(num_haps, num_times, 658.0, T, pop_size, lambda);
// double temp3 = compute_prob(num_haps, num_times, 1640.0, T, pop_size, lambda);
// //printf("temp2 = %f\n", temp2);
// 
//total = pow(fabs(temp3- .25),2);

//double temp4 = compute_prob(num_haps, num_times, bin_center[l], T, pop_size, lambda);
printf("total = %f\n", total);
if (total >-100 && total < 100){}
else{abort();}
return ( total );
}



double first_deriv(int num_times, double x[num_times])
{

int j=0;
double ans = 0;
for (j=1; j< num_times; j++)
{
if (num_times !=1) {ans = ans + (1.0 / (num_times -1) ) * fabs(x[j] - x[j-1]);}
else{ans = 0;}
}

return ans;
}





double my_f (const gsl_vector *v, void *params)
{
double *p = (double *)params;  // only here to make params used, for debugging/optimizing purposes
//p[0] = 2.0;

int num_times = p[1];
// double *dp = (double *)params;  //commented b/c unused in formulaAug18.c
double x[num_times];
double variation = 0;
double second_deriv = 0;
double third_deriv = 0;

// double negative_penalty = 0;//commented b/c unused in formulaAug18.c

int j;
for (j=0; j<num_times; j++)
{
x[j] = gsl_vector_get(v,j);

if ( x[j] < 100 )
{
return 100000000000;
}
}
//printf("pop 0 = %f\n", x[0]);
double temp1 = Joint_like( num_times,time, x, p);// + stiffness * first_deriv(num_times,x);
//printf("join_like = %f\n", temp1);
return( temp1 );  //
//return(Joint_like( num_times,time, x, p) + stiffness * first_deriv(num_times,x) );  //
}









void my_df (const gsl_vector *v, void *params,
       gsl_vector *df)
{
double *p = (double *)params;  // only here to make params used, for debugging/optimizing purposes
//p[0] = 2.0;

int num_times = p[1];

  double x[num_times], y[num_times];
//  double *dp = (double *)params;  //commented b/c unused in formulaAug18.c

int j;

// IF delta is NOT 1 then need to reinstate inv_delta and include it below in gsl_vector_set...
double delta = 1;
// double inv_delta = 1 / delta ;  //commented b/c unused in formulaAug18.c
for (j=0; j<num_times; j++)
{
	x[j] = gsl_vector_get(v, j);
	y[j] = gsl_vector_get(v, j);
}

for (j=0; j<num_times; j++)
{
	y[j] = y[j]+delta;
/*
gsl_vector_set(df, j, (1.0 / delta) * ( Joint_like(num_times, time, y, p) + stiffness * first_deriv(num_times, y) - Joint_like(num_times, time, x, p) - stiffness * first_deriv(num_times, x) )  );*/
gsl_vector_set(df, j, (1.0 / delta) * ( Joint_like(num_times, time, y, p)  - Joint_like(num_times, time, x, p) ));
	y[j] = x[j];
}


}



/* Compute both f and df together. */


void my_fdf (const gsl_vector *x, void *params,

        double *f, gsl_vector *df)

{

  *f = my_f(x, params);

  my_df(x, params, df);

}




int read_binned_file(char filename[])
{

char file_name[25];
FILE *fp;

fp = fopen(filename,"r"); // read mode
if( fp == NULL )
{
      perror("Error while opening the data file.\n");
      exit(EXIT_FAILURE);
}

int temp;
double temp_rho;
int i=0;

while(fscanf(fp, "%d\t%lf\n", &temp, &temp_rho) != EOF && i < 5005)
{ 
	first[i]=temp;
	bin_center[i] = temp;
	bin_counter[i] = 1;//temp_rho;
	rho[i] = 0;
	num_loci = i;
	i++;
}
printf("NUMLOCI = %d\n", num_loci);
fclose(fp);
return 0;
}



int read_quant_file(char filename[])
{
printf("Entering read_quant_file... %f\n", 1.0);
char file_name[25];
FILE *fp;
fp = fopen(filename,"r"); // read mode

if( fp == NULL )
{
	perror("Error while opening the data file.\n");
	exit(EXIT_FAILURE);
}

double temp_A;
double temp_B;
int i=0;

while(fscanf(fp, "%lf\t%lf\n", &temp_A, &temp_B) != EOF && i < 5005)
{ 
	quant_list[i] = .01 * temp_A;
	quant_values[i] = temp_B;
	printf("%f -th percentile = %f\n", temp_A, quant_values[i]);
	i++;
}
num_quantiles = i;
fclose(fp);
return 0;
}



int main (int argc, char *argv[])
{



printf("Number of haplotypes: %d\n", atoi(argv[1]) );
int temp =0;
temp = atoi(argv[1]);
printf("Number of time intervals:  %d\n", atoi(argv[2]) );

printf("File name:  %s\n", argv[3] );

 int num_haps = atoi(argv[1]);
 int num_times = atoi(argv[2]);


//make_bins(num_haps);
//readfile(argv[3]);
read_binned_file(argv[3]);
//make_custom_bins(num_haps); // one bin per match length
// printf("num_loci = %d\n", num_loci);
// printf("marker %f\n", 1.7);
// printf("a;lskdjf %d\n", 3);


make_times(num_haps, num_times);
read_quant_file(argv[4]);
precompute(num_haps, num_times);




/*
 struct scenario the_scenario;
the_scenario.num_haps = argv[1];
the_scenario.num_times = argv[2];*/



  size_t iter = 0;
//  size_t i; // commented out to fix the problem labelled XYRQ below
		// there was a comparison between signed/unsiged int
  
int status;

  const gsl_multimin_fdfminimizer_type *T;
  gsl_multimin_fdfminimizer *s;

  /* Position of the minimum (1,2). */

int j=0;

double par[2];

par[0] = num_haps;
par[1] = num_times;

  gsl_vector *x;
  gsl_multimin_function_fdf my_func;

  my_func.f = &my_f;
  my_func.df = &my_df;
  my_func.fdf = &my_fdf;
  my_func.n = num_times;
  my_func.params = &par;


  /* Starting point, x = (5,7) */

  x = gsl_vector_alloc (num_times);

for (j=0; j<num_times; j++)
{
if(j<=50){gsl_vector_set (x, j, 30000  );}
else {gsl_vector_set (x, j, 20000   );}
}

 T = gsl_multimin_fdfminimizer_conjugate_fr;
//   T = gsl_multimin_fdfminimizer_vector_bfgs2;

  s = gsl_multimin_fdfminimizer_alloc (T, num_times);

// likelihood setting
  gsl_multimin_fdfminimizer_set (s, &my_func, x, 10000, .000001);

//  quantile setting
  //gsl_multimin_fdfminimizer_set (s, &my_func, x, 3000.0, .100);


  do
    {

      iter++;
      status = gsl_multimin_fdfminimizer_iterate (s);

	//printf("help %d\n", 3);

      if (status){
printf ("Iteration terminated at:\nSTART HERE\n");

	//printf ("%5zu ", iter);  // using %zu instead of %d b/c iter is type size_t not type int
			int k=0;
			for (k = 0; k < num_times; k++)  
			{
			printf ("%f\t%10.3f\n", time[k], gsl_vector_get (s->x, k));  // problem XYRQ
						// using k here renders size_t i above unused
			}
	printf ("f() = %7.3f \n", s->f);
//  q
        break;}

//printf("help %d\n", 4);


 //.likelihood setting
	status = gsl_multimin_test_gradient (s->gradient, 1e-10 );

// quantile setting
//	status = gsl_multimin_test_gradient (s->gradient, .010 );


      if (status == GSL_SUCCESS)
 { printf ("Minimum found at:\nSTART HERE\n");}
			int k=0;
 			for (k = 0; k < num_times; k++)  
 			{
 			printf ("%f\t%10.3f\n", time[k], gsl_vector_get (s->x, k)); // problem XYRQ
 						// using k here renders size_t i above unused
 			}
	//printf ("f() = %7.3f \n", s->f);
      
// 	  printf ("Minimum found at:\n");
// 
// 	printf ("%5zu ", iter);  // using %zu instead of %d b/c iter is type size_t not type int
// 			int k=0;
// 			for (k = 0; k < num_times; k++)  
// 			{
// 			printf ("%10.3e ", gsl_vector_get (s->x, k)); // problem XYRQ
// 						using k here renders size_t i above unused
// 			}
 	

	
//       printf ("%5d %.5f %.5f %10.5f\n", iter,
//               gsl_vector_get (s->x, 0),
//               gsl_vector_get (s->x, 1),
//               s->f);

    }

  while (status == GSL_CONTINUE && iter < 1000);

  gsl_multimin_fdfminimizer_free (s);
  gsl_vector_free (x);


  return 0;

}













