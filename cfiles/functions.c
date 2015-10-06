#include "math.h"

#include <gsl/gsl_math.h>

#include <gsl/gsl_multimin.h>

#include <gsl/gsl_math.h>

#include <stdlib.h>

#include <stdio.h>

#include "expfirenov28.h"

double first_deriv(int num_times, double x[num_times])
{

int j=0;
double ans = 0;
for (j=1; j< num_times; j++)
{
if (num_times !=1) {ans = ans + (1.0 / (num_times -1) ) * fabs(x[j] - x[j-1]);}
else{ans = 0;}
}
//printf("ans = %f\n", ans);
return ans;
}
