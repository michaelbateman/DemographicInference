#include "math.h"

#include <gsl/gsl_math.h>

#include <gsl/gsl_multimin.h>

#include <gsl/gsl_math.h>

#include <stdlib.h>

#include <stdio.h>

#include "expfirenov28.h"

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
//static long double mu = .000000025;
static long double mu = 2 * .000000025;

double scale_factor = 1e-15;

int num_pop_sizes = 2;
int num_growth_rates = 2;



// large stiffness biases toward more constant population size
// see my_f and my_df

static double stiffness = 0.01;

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









void make_times(int num_haps, int num_times)
{


double pairs = num_haps * (num_haps -1) / 2.0;




    int j = 0;

    time[0] = 0.0;

//printf("time(%d) = %f\n", j, time[0]);

double gap = (1.0  / num_times) * 2.0 *  (40000 / pairs); 

time[1] = 3.0 * gap;
    for (j=2; j<num_times - 2; j++)

    {


time[j] = time[j-1] + gap;
//time[j] = (pow(j, 2) / pow(num_times, 2)) * 2.0 * (20000 / pairs); 
//    time[j] =  (1.0 / pairs) * 5000 * pow(1.5, j);

//time[1] = 300.0;

 	printf("time(%d) = %f \n", j, time[j]);

    }
time[num_times - 2] = time[num_times - 3] + 3.0 * gap;
time[num_times - 1] = time[num_times - 2] + 3.0 * gap;



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
//printf("A(%d) = %f\n", i, A[i]);
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
	//if(j< 40 && j > 20) {printf("w(%d, %d) = %Lf\n", j, i, w[j][i]);}
    }

}


for(j=1; j<num_haps; j++)
{
long double test = 0;
	 for(i=1; i<=j; i++)
	{
	test = test + w[j][i];
	}
//printf("sum w[%d][.] = %.20Lf\n", j, test);
}



/*
for(i=1; i<num_haps; i++)
{
	for(l=1; l<num_haps; l++)
	{
		for(m=1; m<num_haps; m++)
		{
			pre[i][l][m] = 0;
		}
	}
}

*/

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

void precompute_other_quantities( int num_haps, int num_times, double L, double T[pop_dim], double pop_size[num_times])
{

// C2 is the L-dependent but (i,l,m)-independent part of C in PHI
// we precompute to remove all exponentiation from PHI
// the use of C2 inside an exponential is hidden in Z which is in PHI 

C2 = - mu * L;  // negative sign

H2 = 2 * mu * L;  // positive sign

I2 =  mu * L;  // positive sign

int m;
double lambda[num_times];
for (m=0; m< num_times; m++)  
{
	lambda[m] = (1 / pop_size[m]);
// 	C2 = - mu * L;
Zfirst2[m] = exp( -T[m] * C2);
if(m == num_times -1){ Zsecond2[m] = 0; }
else{Zsecond2[m] = exp( -T[m+1] * H2);}


Yfirst2[m] = exp( -T[m] * H2);
if(m == num_times -1){ Ysecond2[m] = 0; }
else{Ysecond2[m] = exp( -T[m+1] * H2);}


Wfirst2[m] = exp( -T[m] * I2);
if(m == num_times -1){ Wsecond2[m] = 0; }
else{Wsecond2[m] = exp( -T[m+1] * I2);}


}




int i;

for(i = 1; i < num_haps; i++)
{
	for (m=0; m< num_times; m++) 
	{
	//	B1[i][m] = mu * L + A[i] * lambda[m];
	//	B2[i][m] = 2 * mu * L + A[i] * lambda[m];

	//	Z1[i][m] = Z(m, B1[i][m]);
	//	Z2[i][m] = Z(m, B2[i][m]);
	}
}


}


// There is some redundancy in computing F and G separately, since they 
// rely on many of the same sub-calculations
// a more elegant solution would be to have one function H that returns 
// both F,G in a vector but this seems a bit complicated for now
long double F(int num_haps, int num_times, int i, double L, double T[pop_dim], double lambda[num_times] )
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
    C[j] = 2*mu*L + A[i] *  lambda[j];
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
 	term[j] =   2*mu * A[i] *   lambda[j] * G_exp_2[i][j] * prefix[j] * mess[j];
	//printf("mess(%d)=%f\n",j, mess[j]);
    sum = sum + term[j];
	

}  // end loop in j 

//double new_sum = sum *  2*mu * A;

//printf("new_sum = %f\n", sum);
//printf("old_sum = %f\n", compute_likelihood( L, T, pop_size));

return sum;
}

long double G(int num_haps, int num_times,int i, double L, double T[pop_dim], double lambda[num_times] )
{
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

    C[j] = 2*mu*L + A[i] * lambda[j];
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
      //  new_mess[j] = exp(-C[j]*T[j]) * (1 / C[j] );
	new_mess[j] = Yfirst1[i][j] * Yfirst2[j]* (1 / C[j] );
    }
    else
    {//printf("temp = %f\n", (-C[j]*T[j]));
   //     new_mess[j] = exp(-C[j]*T[j]) * (1 / C[j] )             - exp(-C[j]*T[j+1]) * (1 / C[j] ) ;

	  new_mess[j] = (Yfirst1[i][j] * Yfirst2[j] - Ysecond1[i][j] * Ysecond2[j] ) * (1 / C[j] ) ;

    }

// the only difference between G and F is the lack of a 2*mu*T, which shows up in the next line 
// with the 2*mu removed, and the new_mess instead of mess

    //	term[j] =    lambda[j] * exp ( A[i] * T[j] * lambda[j] ) * prefix[j] * new_mess[j];
	term[j] =    lambda[j] * G_exp_2[i][j] * prefix[j] * new_mess[j];
    sum = sum + term[j];

}  // end loop in j 
//printf("temp = %f\t%f\t%f\n",  C[4], 2*mu*L,  A * (1.0 / pop_size[4]) );
//printf("mess(%d)=%f\n",j, new_mess[j]);
//printf("sum_G = %f\n", sum);

double new_sum = A[i]*sum;

return new_sum;
}


long double PHI(int num_haps, int num_times, int i, int l, double L, double T[pop_dim], double lambda[num_times])
{

long double C;

long double sum = 0;
long double X[num_times];
long double interior_X[num_times];
long double boundary_X[num_times];
int m=0;
int j =0;

long double new_Z;
long double new_Y;
long double new_W;
for (m=0; m< num_times; m++) 
{
C = lambda[m] *( A[i] - A[l] ) - mu * L;
new_Z = (1 / C) * ( Zfirst1[i][l][m] * Zfirst2[m] - Zsecond1[i][l][m] * Zsecond2[m] );

//printf("C = %f\t\t m = %d\n", C, m);
//if(m != num_times -1){interior_X[m] = lambda[m] * q[i][l][m] *Z(m, C) ;}

//if(m != num_times -1){interior_X[m] = product[i][l][m] *Z(m, C) ;}
if(m != num_times -1){interior_X[m] = product[i][l][m] * new_Z ;}

//boundary_X[m] = lambda[m] * q[i][l][m] * exp(-T[m] * C) / C ;
//boundary_X[m] = lambda[m] * q[i][l][m] * exp_part_other[i][l][m] * exp_part_L[m] / C ;
boundary_X[m] = product2[i][l][m] * exp_part_L[m] / C ;

/*
if( exp(-T[m+1] * C) < 1e20 && m < num_times - 1 )  
{
if(m==0){X[m] = boundary_X[m];}
else{X[m] = X[m-1] + interior_X[m-1] - boundary_X[m-1] + boundary_X[m];}
//printf("int_X = %f\n", interior_X[m]);
}

else if( exp(-T[m] * C) < 1e20 && m ==  - 1) 
{
if(m==0){X[m] = boundary_X[m];}
else{X[m] = X[m-1] + interior_X[m-1] - boundary_X[m-1] + boundary_X[m];}
}
else
{
X[m] = 1e20;
}



The reason for these cases is that Z(m,C) and exp(-T[m] * C) and 
consequently X[m] can be quite large
however if they are quite large then B2 below is even larger and
consequently term1 below is very very small despite X[m] being large

*/

if(m==0){X[m] = boundary_X[m];}
else{X[m] = X[m-1] + interior_X[m-1] - boundary_X[m-1] + boundary_X[m];}

// if(X[m] > - 1e20 && X[m] < 1e20){}
// else if (X[m] <= - 1e20){ X[m] = - 1e20;}
// else { X[m] =  1e20;}

if(X[m] > - my_inf && X[m] < my_inf){}
else if (X[m] <= - my_inf)
{ 
X[m] = - my_inf; //printf("L = %f\n",L);
}
else 
{ 
X[m] =  my_inf;//printf("L = %f\n",L);
}


//printf("bdry_X = %f\t\t%f\t\tL = %f\n", boundary_X[m], exp(-T[m] * C), L);
//printf("X = %f\n", X[m]);


//}




// for (m=0; m< num_times; m++) 
// {


//C = lambda[m] *( A[i] - A[l] ) - mu * L;
//double B1 = mu * L + A[i] * lambda[m];
//double B2 = 2 * mu * L + A[l] * lambda[m];
//printf("hellow\n");



//term2 = lambda[m] * q[i][l][m] * Z1[i][m] / C;

//double temp_term2 = product4[i][l][m]  / C;
//product4[i][l][m] = product[i][l][m] * Z1[i][m];

//term1 = X[m] * Z( m, B2 );

long double H = H1[l][m] + H2;
new_Y = (1 / H) * ( Yfirst1[l][m] * Yfirst2[m] - Ysecond1[l][m] * Ysecond2[m] );

//double term1 = X[m] * Z2[l][m];
long double term1 = X[m] * new_Y;


long double I = I1[i][m] + I2;
new_W = (1 / I) * ( Wfirst1[i][m] * Wfirst2[m] - Wsecond1[i][m] * Wsecond2[m] );

//double term2 = product[i][l][m] * Z1[i][m] / C;
long double term2 = product[i][l][m] * new_W / C;



//term2 = lambda[m] * q[i][l][m] * Z(m, B1 ) / C;

//psi = lambda[m] * r[l][m] * ( term1 - term2);
long double psi = product3[l][m] * ( term1 - term2);

//if(term1 > -1e-20 && term1 < 1e-20 ){ abort();}
//printf("term1 = %.20f\n", term1);
//printf("term2= %.20f\n", term2);
//printf("temp_term2= %.20f\n", temp_term2);
sum = sum + psi;
//if (psi <=2  && psi> -5){}//printf("\t\t\t\t\t\tFOUND IT %f\n", 0.0); break;}
//else{printf("psi = %f\n", Z(1,C)); printf("C = %f\t\t L = %f\n", X[1], L);abort();}
//printf("sum = %.20f\n",term1);
}
//if (L < 3000){printf("L = %f\t\t sum = %.20Lf\n", L, sum);}

return sum; 
}

double compute_prob(int num_haps, int num_times, double L, double T[pop_dim], double pop_size[num_times], double lambda[num_times])
{




int marker = 0;
int i=0;
int k = 0;
int j = 0;
int l = 0;
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
 

long double  prob_L_from_coalescent_time[num_haps];
long double  prob_less_than_L_from_coalescent_time[num_haps];
// double total_D = 1;  //commented b/c unused in formulaAug18.c

prob_L_from_coalescent_time[0] = 0;
prob_less_than_L_from_coalescent_time[0] = 0;

// first we precompute G( A[i], L, T, pop_size) and F( A[i], L, T, pop_size) 
// for i = 1, ..., num_haps - 1 that way we don't compute them in the full
// double loop over i and j below


// Code for computing F and G simultaneously exists in formulaJuly29.c


long double f[num_haps];
long double g[num_haps];
f[0] = 0;
g[0] = 1;

for(i = 1; i < num_haps; i++)
{
f[i] = F(num_haps, num_times,  i, L, T, lambda );
g[i] =  1.0- G(num_haps, num_times,  i, L, T, lambda );
//printf("f(%d) = %Lf\t\t g = %Lf\t\tL = %f\n", i, f[i], g[i], L);
}


long double phi[num_haps][num_haps];
phi[0][0] = 1.0;
for(i = 1; i < num_haps; i++)
{
	for(l = i+ 1; l < num_haps; l++)
	{
		 phi[i][l]= PHI(num_haps, num_times, i,l, L, T, lambda);

		//printf("{phi[%d][%d] = %.20f\t\t%f\n", i,l,phi[i][l], L);

		if (phi[i][l] <=2  && phi[i][l] > -5){}//printf("\t\t\t\t\t\tFOUND IT %f\n", 0.0); break;}
		//else{printf("phi = %f\n", phi[i][l]); }
	}
}








/*
This double loop makes prob_L_from_coalescent_time[j] and 
prob_less_than_L_from_coalescent_time[j]
This is the weakest point of the algorithm, as it takes time  
(1/2) * num_haps * (num_haps - 1) 

the essential reason for this is that the matrix w is square of width (num_haps - 1), 
and we need half of its entries to compute all of the numbers prob_L_from_coalescent_time[j]
*/

for(j=1; j < num_haps; j++)
{
long double temp_C = 0;
long double temp_D = 0;

	for(i=1; i <=j; i++)
		{
		temp_C = temp_C + w[j][i] * f[i];
		temp_D = temp_D + w[j][i] * g[i];
	//	printf("i = %d\tw = %Lf\n", i, w[j][i]);
		}
	prob_L_from_coalescent_time[j] = temp_C;
	prob_less_than_L_from_coalescent_time[j] = temp_D;

if (temp_D <=50&& temp_D > -50){}//printf("\t\t\t\t\t\tFOUND IT %f\n", 0.0); break;}
else{
//temp_D = 1.0;
printf("temp_D = %Lf\n", temp_D); 
abort();
}
}

long double new_D[num_haps][num_haps];
long double new_d[num_haps][num_haps];
long double temp_new_d = 0.0;
for(j=1; j < num_haps; j++)
{
	for(k=j+1; k <= j + complexity; k++)
	{
		temp_new_d = 0.0;

		for( i = 1; i <=j; i++)
		{
			for( l = j+1; l <=k; l++)
			{
			temp_new_d = temp_new_d + w[j][i] * new_w[j+1][k][l] * A[i] * A[l] * phi[i][l];
			}
		}
	

		new_d[j][k] = temp_new_d;
		new_D[j][k] = 1 - new_d[j][k];
	}
//printf("new_D(%d,%d) = %.20f\n", j,j+1, new_D[j][j+1]);

}

long double E[num_haps];


long double E_total=1;
for(j=1; j < num_haps; j++)
{
	E_total = E_total * prob_less_than_L_from_coalescent_time[j];
}
//printf("E_total = %f\n", E_total);


/*
for(j=1; j < num_haps; j++)
{
double temp_E = 1;
temp_E = E_total  / prob_less_than_L_from_coalescent_time[j] ;
 
	for(k=j+1; k <= j+complexity; k++)
	{
		temp_E = temp_E * new_D[j][k] / prob_less_than_L_from_coalescent_time[k] ;
	}


if (temp_E <=2  && temp_E > -5){}//printf("\t\t\t\t\t\tFOUND IT %f\n", 0.0); break;}
else
{
 temp_E = 0.0;
//printf("temp_E = %f\n", temp_E); 
//abort();
}
E[j] = temp_E;
}

*/

// old version of computing E[j]
// new version has only one loop from j = 1 ... num_haps, old version has two loops
// This made no huge difference in running time *facepalm* when implemented b/c PHI is bottleneck
// and PHI is already computed by this point
// Could still be important though if num_haps gets very large???


for(j=1; j < num_haps; j++)
{
long double temp_E = 1;

	for(k=1; k < num_haps; k++)
	{
		if(k > j && k <= j + complexity )
		{

//new_D[j][k] = 1 - new_d[j][k];
			temp_E = temp_E * new_D[j][k];	
		//	printf("first_temp_E = %f\n", temp_E);	
		}
		else if(k != j)
		{	
			     temp_E = temp_E * prob_less_than_L_from_coalescent_time[k];
		//	printf("second_temp_E = %f\t\t%f\n", temp_E, L);	
		}
	}
	
E[j] = temp_E;

if (temp_E <=2  && temp_E > -5){}//printf("\t\t\t\t\t\tFOUND IT %f\n", 0.0); break;}
else
{
 temp_E = 1.0;
//printf("temp_E = %f\n", temp_E); 

//abort();
}

}






long double prob_L_max = 0;
//double prob_L_max_2 = 0;  //commented b/c unused in formulaAug18.c
for(j=1; j<num_haps; j++)
{
//if (X[j] == 1.0) { printf("YES = %f\t\t\t %d\n", 1.0, j);}
 	prob_L_max = prob_L_max + prob_L_from_coalescent_time[j] * E[j];
//printf("E = %Lf\n", E[j]);
	//prob_L_max = prob_L_max + prob_L_from_coalescent_time[j] * X[j];
/*printf(" \t\t\tNEW = %f\n", prob_L_max);
printf(" exactly L = %f\n", prob_L_from_coalescent_time[j]);
printf(" help = %f\n", X[j]);
printf(" test = %f\n", prob_L_from_coalescent_time[j]* X[j] );*/
}
//printf("pop_size = %f\n", pop_size[0] ) ;
//printf("prob L max = %.20Lf \n", prob_L_max);
return prob_L_max;


}

int diff(double *x, int k, double delta)
{

x[k] = x[k] + delta;
return;
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

int problem_site_counter = 0;

for(l=0; l< num_loci; l++)
	{
	if(bin_center[l] >0)
		{
		//mu = .000000025;// + rho[l];

//printf("help %f\n", 1.0);
		precompute_other_quantities( num_haps, num_times, bin_center[l], T, pop_size);
		double temp4 = compute_prob(num_haps, num_times, bin_center[l], T, pop_size, lambda);
		//if(l == 348){printf("%.20f\t\t\n", temp4);}
		if(temp4 <=0)
		{
			//problem_site_counter = problem_site_counter + 1;
			printf("%.20f\t L(%d) = %f \t mu = %.14Lf\n", temp4, l, bin_center[l], mu);
		
		//mu = .0000000125;
		double temp5 = compute_prob(num_haps, num_times, bin_center[l], T, pop_size, lambda);
		printf("%.20f\t REPRINT \t\t mu = %.14Lf\n", temp5, l, bin_center[l], mu);

		//mu = .0000000125 + 2 * rho[l];
		double temp6 = compute_prob(num_haps, num_times, bin_center[l], T, pop_size, lambda);
		printf("%.20f\t REPRINT \t\t mu = %.14Lf\n", temp6, l, bin_center[l], mu);

		}
			else{
				//bin_counter[l] = bin_counter[l] * bin_center[l];
				joint = joint + log( temp4 ) * bin_counter[l] ;//* bin_center[l]/ 1000.0  ; // extra factor of bin_center[l] in place to account for the collection of all matches rather than only s
				//printf("joint = %f\n", joint);

				//printf("L = %f\t, counter = %d\n", bin_center[l], bin_counter[l]);
			} 
	//	printf("temp4 = %.20f\n", temp4);
	    	//joint = joint + log( compute_prob(num_haps, num_times, bin_center[l], T, pop_size, lambda) ) * bin_counter[l];
	//	printf("temp_JOINT = %f\t\t L = %f\n", joint, bin_center[l]);
		if (joint < 1e20 && joint > -1e20){}
		else{ abort() ;}
		}


	}

	
//printf("hellow!!! %f", 1.0);    
//joint =  log( compute_prob( .00003, T, pop_size) );
//printf("JOINT = %f", joint); 	
return ( - joint );
}


/*
double first_deriv(int num_times, double x[num_times])
{

int j=0;
double ans = 0;
for (j=1; j< num_times; j++)
{
if (num_times !=1) {ans = ans + (1.0 / (num_times -1) ) * fabs(x[j] - x[j-1]);}
else{ans = 0;}
}
printf("ans = %f\n", ans);
return ans;
}*/


double my_f (const gsl_vector *v, void *params)
{
double *p = (double *)params;  

int num_times = p[1];

double x[num_times];
double variation = 0;
double second_deriv = 0;
double third_deriv = 0;



int j;
int breakpoint = num_times / 2 - 1;
int special[4];
special[0] = 0;
special[1] = breakpoint;
special[2] = breakpoint + 1;
special[3] = num_times - 1;
for (j=0; j<4; j++)
{
//printf("special %d = %d\n",j, special[j]);
}

x[0] = gsl_vector_get(v,0);
x[breakpoint] = gsl_vector_get(v,1);
x[breakpoint + 1] = gsl_vector_get(v,2);
x[num_times - 1] = gsl_vector_get(v,3);
double alpha =  (1 / (time[breakpoint] - time[1]) ) * log(x[0] / x[breakpoint]);
double alpha2 = (1 / ( time[num_times - 1] - time[breakpoint + 2 ] ) ) * log(x[breakpoint + 1] / x[num_times - 1]);


for (j=1; j< breakpoint; j++)
{
        x[j] = x[j-1] * exp(-alpha * ( time[j+1] - time[j] ));
}

for (j=breakpoint +  2; j < num_times - 1; j++)
{
        x[j] = x[j-1] * exp(-alpha2 * ( time[j+1] - time[j] ));
}
for (j=0; j< num_times; j++)
{
//printf("my x(%d) = %f, %.10f, %.10f\n", j, x[j], alpha, alpha2);
}



double temp1 = Joint_like( num_times,time, x, p) + stiffness * first_deriv(num_times,x);
return( temp1 );  

}



void my_df (const gsl_vector *v, void *params,
       gsl_vector *df)
{
double *p = (double *)params;  
int num_times = p[1];
double x[num_times], y[num_times];
int j;
double delta = 1;
double pop_delta = 1000;
int breakpoint = num_times / 2 - 1;
int special[4];
special[0] = 0;
special[1] = breakpoint;
special[2] = breakpoint + 1;
special[3] = num_times - 1;

for (j=0; j<4; j++)
{
	//printf("special %d = %d\n",j, special[j]);
}

x[0] = gsl_vector_get(v,0);
x[breakpoint] = gsl_vector_get(v,1);
x[breakpoint + 1] = gsl_vector_get(v,2);
x[num_times - 1] = gsl_vector_get(v,3);
double alpha =  (1 / (time[breakpoint] - time[1]) ) * log(x[0] / x[breakpoint]);
double alpha2 = (1 / ( time[num_times - 1] - time[breakpoint + 2 ] ) ) * log(x[breakpoint + 1] / x[num_times - 1]);

for (j=1; j< breakpoint; j++)
{
	x[j] = x[j-1] * exp(-alpha * ( time[j+1] - time[j] ));
}

for (j=breakpoint +  2; j < num_times - 1; j++)
{
	x[j] = x[j-1] * exp(-alpha2 * ( time[j+1] - time[j] ));
}

for (j = 0; j < num_times; j++)
{//	printf("alpha 2 = %f\n", alpha2);
//	printf("%d   check %f\n", j, x[j]);
	y[j] = x[j];
}

int k ;
for ( k = 0; k < 4; k++)
{
//printf("checka;sdlfj%f\n", 2.0);
//printf("%d  %d\n",k,  special[k]);
y[special[k]] = x[special[k]] + pop_delta;

gsl_vector_set(df, k, (1.0 / pop_delta) * ( Joint_like(num_times, time, y, p) - Joint_like(num_times, time, x, p) ) );
y[special[k]] = x[special[k]];

}





}



/* Compute both f and df together. */


void my_fdf (const gsl_vector *x, void *params, double *f, gsl_vector *df)
{
	*f = my_f(x, params);
	my_df(x, params, df);
}



int readfile(char filename[])
{
   //char ch;  //commented b/c unused in formulaAug18.c
   char file_name[25];
   FILE *fp;
   //printf("Enter the name of file you wish to see\n");
   //gets(file_name);




//  fp = fopen("recover10.max","r"); // read mode

//fp = fopen("double10.max","r"); // read mode
  fp = fopen(filename,"r"); // read mode
//  fp = fopen("test30.checksep.max","r"); // read mode
//  fp = fopen("outputlongunsorted.txt","r"); // read mode
// fp = fopen("constant30.txt","r"); // read mode
   if( fp == NULL )
   {
      perror("Error while opening the data file.\n");
      exit(EXIT_FAILURE);
   }


/*  Printing to screen */

 //  printf("The contents of your file are :%s\n", file_name);
//    while( ( ch = fgetc(fp) ) != EOF )
//       printf("%c",ch);

int temp;
double temp_rho;
int i=0;

//printf("test1");
	while(fscanf(fp, "%d\t%lf\n", &temp, &temp_rho) != EOF && i < 50005)
	{ //if(i> 0){ printf("i = %d \t  hello:  %d\t rho = %f\n",i, temp, temp_rho);}
		first[i]=temp;
		//rho[i] = temp_rho * 1e-8;
		rho[i] = 0;
/*
	int l = 0;
	for(l=0; l< num_bins; l++)
		{
			if( first[i] >= bin_bound[l] && first[i] < bin_bound[l+1] )
				{
				//printf(" success %d", l );
		bin_counter[l] = bin_counter[l] + 1;}
		//else{printf("fail %f", 1.0);}
		}
*/

num_loci = i;

	i++;
	}
printf("NUMLOCI = %d\n", num_loci);
/*
int check =0;
int l=0;
for(l=0; l< num_bins; l++)
 {
//if(l < 5000){printf("counter(%d) = %d\t\tbin_bound = %f\n", l, bin_counter[l], bin_bound[l]);}
check = check + bin_counter[l];
}
//printf("CHECK = %d", check);
*/



   fclose(fp);
   return 0;

}



int read_binned_file(char filename[])
{
   //char ch;  //commented b/c unused in formulaAug18.c
   char file_name[25];
   FILE *fp;
   //printf("Enter the name of file you wish to see\n");
   //gets(file_name);




//  fp = fopen("recover10.max","r"); // read mode

//fp = fopen("double10.max","r"); // read mode
  fp = fopen(filename,"r"); // read mode
//  fp = fopen("test30.checksep.max","r"); // read mode
//  fp = fopen("outputlongunsorted.txt","r"); // read mode
// fp = fopen("constant30.txt","r"); // read mode
   if( fp == NULL )
   {
      perror("Error while opening the data file.\n");
      exit(EXIT_FAILURE);
   }


/*  Printing to screen */

 //  printf("The contents of your file are :%s\n", file_name);
//    while( ( ch = fgetc(fp) ) != EOF )
//       printf("%c",ch);

int temp;
double temp_rho;


int i=0;

//printf("test1");
	while(fscanf(fp, "%d\n", &temp) != EOF && i < 5005)
	{ //if(i> 0){ printf("i = %d \t  hello:  %d\t rho = %f\n",i, temp, temp_rho);}
		first[i]=temp;

		bin_center[i] = temp;
		bin_counter[i] = 1;//temp_rho;
		//printf(" L = %d \t counter = %f\n", bin_center[i], bin_counter[i]);
		//rho[i] = temp_rho * 1e-8;
		rho[i] = 0;
/*
	int l = 0;
	for(l=0; l< num_bins; l++)
		{
			if( first[i] >= bin_bound[l] && first[i] < bin_bound[l+1] )
				{
				//printf(" success %d", l );
		bin_counter[l] = bin_counter[l] + 1;}
		//else{printf("fail %f", 1.0);}
		}
*/

num_loci = i;

	i++;
	}
printf("NUMLOCI = %d\n", num_loci);
/*
int check =0;
int l=0;
for(l=0; l< num_bins; l++)
 {
//if(l < 5000){printf("counter(%d) = %d\t\tbin_bound = %f\n", l, bin_counter[l], bin_bound[l]);}
check = check + bin_counter[l];
}
//printf("CHECK = %d", check);
*/



   fclose(fp);
   return 0;

}


void read_parameters(char filename[], double *par)
{
char file_name[25];
FILE *fp;
fp = fopen(filename,"r"); // read mode

if( fp == NULL )
{
	perror("Error while opening the data file.\n");
	exit(EXIT_FAILURE);
}
double temp;

int i=0;

while (fscanf(fp,"%lf",&temp)!=EOF)
{
	
	par[i] = temp;
	printf("%.10f\n",par[i]);
	i++;
}
fclose(fp);

// par[0] = model.num_haps;
// par[1] = model.num_times;
// par[2] = num_eras;
// par[3] = num_growth_rates;
// par[4] = mu;
// par[5] = RHO;
// 
// par[6] = complexity;
// par[7] = stiffness;
// par[8] = scale_factor;

//return 0;
}


int main (int argc, char *argv[])
{



printf("Number of haplotypes: %d\n", atoi(argv[1]) );
int temp =0;
temp = atoi(argv[1]);
printf("Number of time intervals:  %d\n", atoi(argv[2]) );

printf("File name:  %s\n", argv[4] );

// int num_haps = atoi(argv[1]);
 int num_times = atoi(argv[2]);
int breakpoint = atoi(argv[3]);

// struct Model model;
// model.num_haps = atoi(argv[1]);
// model.num_times = atoi(argv[2]);

read_binned_file(argv[4]);
//make_custom_bins(num_haps); // one bin per match length
// printf("num_loci = %d\n", num_loci);
// printf("marker %f\n", 1.7);
// printf("a;lskdjf %d\n", 3);
double par[9];
read_parameters("temp.params", par);
int i;
for(i =0; i < 9; i++)
{
printf("par[%d] = %.10f\n", i,  par[i]);
}
// par[0] = model.num_haps;
// par[1] = model.num_times;
// par[2] = num_size_choices;
// par[3] = num_growth_rates;
// par[4] = mu;
// par[5] = RHO;
// 
// par[6] = complexity;
// par[7] = stiffness;
// par[8] = scale_factor;
make_times(par[0], par[1]);
precompute(par[0], par[1]);


  size_t iter = 0;
//  size_t i; // commented out to fix the problem labelled XYRQ below
		// there was a comparison between signed/unsiged int
  
int status;

  const gsl_multimin_fdfminimizer_type *T;
  gsl_multimin_fdfminimizer *s;

  /* Position of the minimum (1,2). */

int j=0;


/*
par[0] = model.num_haps;
par[1] = model.num_times;
par[2] = num_size_choices;
par[3] = num_growth_rates;*/

//int num_size_choices = 2;
//int num_growth_rates = 2;
long double mu;
long double RHO;
int complexity;
double stiffness;
double scale_factor;
// 	par[0] = model.num_haps;
// 	par[1] = model.num_times;
// 	par[2] = num_eras;
// 	par[3] = num_growth_rates; # parameter obsolete
// 	par[4] = mu;
// 	par[5] = RHO;
// 	
// 	par[6] = complexity;
// 	par[7] = stiffness;
// 	par[8] = scale_factor;


int num_vars =  2 * par[2];
// num_vars = num_vars;


printf("hellow1  %f\n", 2.0);
  gsl_vector *x;
  gsl_multimin_function_fdf my_func;

  my_func.f = &my_f;
  my_func.df = &my_df;
  my_func.fdf = &my_fdf;
  my_func.n = num_vars;
  my_func.params = &par;


  /* Starting point, x = (5,7) */

  x = gsl_vector_alloc (num_vars);
/*
gsl_vector_set (x, 0, 24000  );
gsl_vector_set (x, 1, 24000);
gsl_vector_set (x, 2, 24000  );
gsl_vector_set (x, 3, 24000);*/

for (j=0; j<num_vars; j++)
{printf("j = %d\n", j);
if(j<=10){gsl_vector_set (x, j, 30000  );}
else {gsl_vector_set (x, j, 30000   );}
}
printf("waiting %f\n", 1.0);

 T = gsl_multimin_fdfminimizer_conjugate_pr;
   //T = gsl_multimin_fdfminimizer_vector_bfgs2;

  s = gsl_multimin_fdfminimizer_alloc (T, num_vars);

// likelihood setting
  gsl_multimin_fdfminimizer_set (s, &my_func, x, 20000, .01);

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
				for (k = 0; k < num_vars; k++)  
 			{
 			printf ("%10.15f\n", gsl_vector_get (s->x, k)); // problem XYRQ
 						// using k here renders size_t i above unused
 			}
			
	printf ("f() = %7.3f \n", s->f);
//  q
        break;}

//printf("help %d\n", 4);


 //.likelihood setting
	status = gsl_multimin_test_gradient (s->gradient, .00001 );

// quantile setting
//	status = gsl_multimin_test_gradient (s->gradient, .010 );


      if (status == GSL_SUCCESS)
 { printf ("Minimum found at:\nSTART HERE\n");}
			int k=0;
 			for (k = 0; k < num_vars; k++)  
 			{
 			printf ("%10.15f\n", gsl_vector_get (s->x, k)); // problem XYRQ
 						// using k here renders size_t i above unused
 			}
			
	printf ("f() = %7.3f \n", s->f);
      
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













