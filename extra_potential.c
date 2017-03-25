#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>

#include "allvars.h"
#include "proto.h"
#include "define.h"

/*! 
 *! This routine calculate the a 3-pieces power law bulge density profile potential, 
 *! and the acceleration due to this potential to 
 *! the dark matter particles. 
 *! -dlcheng
 */
#ifdef EXTRA_POTENTIAL           /* -dlcheng */

void init_extra_potential()
{
#if (EXTRA_POTENTIAL == 0)
  extra_potential_accel = mw_bulge_potential_accel;
  extra_potential_accel_ad = mw_bulge_potential_accel_ad;
  extra_potential_rho = mw_rho_bulge;
#endif  

#if (EXTRA_POTENTIAL == 1)
  extra_potential_accel = hernquist_potential_accel;
  extra_potential_accel_ad = hernquist_potential_accel_ad;
  extra_potential_rho = hernquist_rho;
#endif
}

#if  (EXTRA_POTENTIAL == 0)                   /* Milky way bulge potential */
static double center[3] = {0.0,0.0,0.0};   /*! x,y,z of the center in unit of kpc/h , usually set h=1 */
static double mbh = 3e-4;                  /*! black hole mass in internal unit 1e10 solar mass /h    */          
static double r_0 = 0.1;                   /*! radius of the flat core                                */
static double r_1 = 1.992;                 /*! the turning radius of the density slope                */
static double r_s1 = 1.0;                  /*! the scale length of the first density power law        */
static double r_s2 = 4.0;                  /*! the scale length of the second density power law       */
static double rho_s1 = 4e-2;               /*! the scale density of the first density power law, 
                                                 in unit of 10^10 solar mass/h                          */
static double rho_s2 = 3e-4;               /*! the scale density of the second density power law, 
                                                 in unit of 10^10 solar mass/h                          */
static double slope_1 = -2.8 ;             /*!  the slope of the first density power law              */
static double slope_2 = -4.25;             /*! the slope of the second density power law              */
static double G_fact = 43007.1;            /*! gravitational constant in internal unit                */


void mw_bulge_potential_accel(FLOAT *pos, FLOAT *bulge_accel)
{
  double r = 0.0;
  int i;
  double mass,drec_r[3];

  for(i=0; i<3; i++)
    {
      r += (pos[i] - center[i]) * (pos[i] - center[i]);
      drec_r[i] = center[i] - pos[i];               /* point to the center */
    }

  r = sqrt(r);
  mass = mw_bulge_mass(r);
  
  for(i=0; i<3; i++)
     bulge_accel[i] = G_fact * mass * drec_r[i] / r / r / r;


}                                                  /* end mw_bulge_potential_accel */

void mw_bulge_potential_accel_ad(FLOAT *pos, FLOAT *bulge_accel, double ad_time)
{
  double r = 0.0;
  int i;
  double mass,drec_r[3];
  double red;

#ifdef AD_EXTRA_POTENTIAL_GROWTH
  red = All.Ti_Current * All.Timebase_interval / ad_time;
#else
#ifdef AD_EXTRA_POTENTIAL_RELEASE
  red = 1.0 - All.Ti_Current * All.Timebase_interval / ad_time;
#endif
#endif

  for(i=0; i<3; i++)
    {
      r += (pos[i] - center[i]) * (pos[i] - center[i]);
      drec_r[i] = center[i] - pos[i];             /* point to the center */
    }

  r = sqrt(r);
  mass = mw_bulge_mass(r);
  
  for(i=0; i<3; i++)
     bulge_accel[i] = G_fact * mass * drec_r[i] * red / r / r / r ;


}                                                /* end mw_bulge_potential_accel_ad */

/* the bulge mass inside radius r*/
double mw_bulge_mass(double r)  
{               
  double rho_core =  rho_s1 * pow(r_0/r_s1,slope_1);

  if(r < r_0)
    return mbh + 4.0 * PI / 3 * rho_core * pow(r,3.0);
  else
    if(r < r_1)
      return mbh + 4.0 * PI / 3 * rho_core * pow(r_0,3.0) + 4.0 * PI * rho_s1 * pow(r_s1, -1.0 * slope_1) / (3.0 + slope_1) * (pow(r,3+ slope_1) - pow(r_0,3 + slope_1)) ;
    else
      return mbh + 4.0 * PI / 3 * rho_core * pow(r_0,3.0) + 4.0 * PI * rho_s1 * pow(r_s1, -1.0 * slope_1) / (3.0 + slope_1) * (pow(r_1,3+ slope_1) - pow(r_0,3 + slope_1)) +  4.0 * PI * rho_s2 * pow(r_s2, -1.0 * slope_2) / (3.0 + slope_2) * (pow(r,3+ slope_2) - pow(r_1,3 + slope_2));

}                                             /* end mw_bulge_mass*/
/* the density at the radius */
double mw_rho_bulge(double r)
{
  
  if(r < r_0)
    return rho_s1 * pow(r_0/r_s1,slope_1);
  else
    if(r < r_1)
      return rho_s1 * pow(r/r_s1,slope_1);
    else
      return rho_s2 * pow(r/r_s2,slope_2);
  
}                                           /* end mw_rho_bulge */

#endif

#if (EXTRA_POTENTIAL == 1)                 /* Hernquist profile for high redshift dwarf galaxy */

static double center[3] = {0.0,0.0,0.0}; /* center in unit of kpc/h , usually set h=1 */
static double Mbt = 1.7e-2;              /* total baryon matter mass in unit of 10^-10 solar mass /h */
static double Rb = 0.01;                 /* scale length of Hernquist model of baryon in unit of kpc */
static double G_fact = 43007.1;          /* gravitational constant in internal unit */

void hernquist_potential_accel(FLOAT *pos,FLOAT *bulge_accel)
{
  double r = 0.0;
  int i;
  double drec_r[3];

  for(i=0; i<3; i++)
    {
      r += (pos[i] - center[i]) * (pos[i] - center[i]);
      drec_r[i] = center[i] - pos[i];   
    }

  r = sqrt(r);
  
  for(i=0; i<3; i++)
     bulge_accel[i] = G_fact * Mbt / r / pow((r+Rb),2.0) * drec_r[i];


}

void hernquist_potential_accel_ad(FLOAT *pos, FLOAT *bulge_accel, double ad_time)
{
  double r = 0.0;
  int i;
  double drec_r[3];
  double red;
#ifdef AD_EXTRA_POTENTIAL_GROWTH
  red = All.Ti_Current * All.Timebase_interval / ad_time;
#else
#ifdef AD_EXTRA_POTENTIAL_RELEASE
  red = 1.0 - All.Ti_Current * All.Timebase_interval / ad_time;
#endif
#endif
  for(i=0; i<3; i++)
    {
      r += (pos[i] - center[i]) * (pos[i] - center[i]);
      drec_r[i] = center[i] - pos[i]; 
    }

  r = sqrt(r);
 
  for(i=0; i<3; i++)
     bulge_accel[i] = G_fact * Mbt / r / pow((r+Rb),2.0) * drec_r[i] * red;
}

double hernquist_rho(double r)
{
  return 1.0 / 2.0 / PI * Mbt * Rb / r / pow((r+Rb),3.0);  
} 
#endif

#endif
