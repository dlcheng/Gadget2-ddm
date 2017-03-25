#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>  /*random number generator */

#include "allvars.h"
#include "proto.h"
#include "define.h"

#ifdef HALO_PART_RANDOM_KICK            /* -dlcheng */

double diff_equation(double *param, double t)
{
  return 4.0 * param[0] * pow(t,3) + 3.0 * param[1] * pow(t,2) + 2.0 * param[2] * t + param[3];
}

double motion_equation(double *param, double t)
{
  return param[0] * pow(t,4) + param[1] * pow(t, 3) + param[2] * pow(t, 2) +  param[3] * t + param[4];
}

double find_timestep_random_kick(FLOAT *vel, FLOAT *accel, double l_kick)
{
  double x_1,x_2;
  double param[5];

  param[0] = (pow(accel[0], 2) + pow(accel[1], 2) + pow(accel[2], 2)) / 4.0; /* index of the forth power */
  param[1] = accel[0] * vel[0] + accel[1] * vel[1] + accel[2] * vel[2];
  param[2] = pow(vel[0], 2) + pow(vel[1], 2) + pow(vel[2], 2);
  param[3] = 0.0;
  param[4] = -1.0 * l_kick * l_kick;

  x_1 = 0.0;    /* starts from 0  */
  x_2 = x_1 - motion_equation(param,x_1) / diff_equation(param,x_1);
  
  while(fabs(motion_equation(param,x_2)) > 1e-4)                 /* demand f(x_2) < 1e-4  */
   {
     x_1 = x_2;
     x_2 = x_1 - motion_equation(param, x_1) / diff_equation(param, x_1);
   }

  return x_2;

}
/* calculate the l_eff increase for this step */
double update_l_eff(FLOAT *pos, FLOAT *vel, double t)
{
  double w_kick;
  double v = 0.0;
  int i;

  for(i=0; i<3; i++)
    v += vel[i] * vel[i];  
  v = sqrt(v);

  w_kick = kick_weight(pos);

  return w_kick * v * t; 
}


double kick_weight(FLOAT *pos)
{
  double r = 0.0;
  int i;
  double w_kick;

  for(i=0; i<3; i++)
    r += pos[i] * pos[i];
  r = sqrt(r);

  w_kick = REF_WEIGHT_KICK * pow((*extra_potential_rho)(r) / (*extra_potential_rho)(1.0), SLOPE_KICK);

  return w_kick;
}

#endif

void vel_kick_copy(FLOAT *a, FLOAT *b)
{
  int i;

  for(i=0; i<3; i++)
    b[i] = a[i];

}                                           /*! copy b from a */

void accel_kick_add(FLOAT *a,FLOAT *b, FLOAT *c)
{
  int i;

  for(i=0; i<3; i++)
    c[i] = a[i] + b[i];
}                                            /*! c = a+b */


void do_half_kick(FLOAT *vel, FLOAT *accel, double t)
{
  int i;

  for(i=0; i<3; i++)
    vel[i] += accel[i] * t;

}                                          /* do the half kick to the particle velocity  */

void do_random_kick(FLOAT *vel, double V_kick)
{ 
  random_kick_tab_update();                /* refresh the random number table firstly */

  double a1 = random_kick_tab[0];         /* [0.1) */
  double a2 = random_kick_tab[1];
  
  a1 = 2.0 * a1 - 1;
  a1 = acos(a1);                           /* this is the theta */

  a2 = 2.0 * a2 * PI; 

  vel[0] += V_kick * sin(a1) * cos(a2); 
  vel[1] += V_kick * sin(a1) * sin(a2);
  vel[2] += V_kick * cos(a1);
}                                         /* ! apply V_kick to vel */

void random_kick_tab_update()
{
  int i;

  for(i = 0; i < 2; i++)
    random_kick_tab[i] = gsl_rng_uniform(random_kick_g);
}                                         /*! end random_kick_num_update */
