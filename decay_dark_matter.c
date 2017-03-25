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
  * this file contains funtions to handle the creation of daughters and modification of the mother particle, 
  * such as its mass,
  * softening length. The velocity of the mother particle is half-kicked as the reference velocity.
  */
#ifdef DECAY_DARK_MATTER  /* -dlcheng */

void create_daught_and_kick(int i, int local_decay_num, int ti_step)
{
 int daught_i = NumPart + (local_decay_num - 1) * DECAY_DAUGHT_PART_NUM;
 int j, tstart, tend;
 FLOAT half_kick_v[3];
 double dt_gravkick1, dt_gravkick2;
 double mass_ratio = (1.0 - exp(-1.0 * P[i].decay_factor)) * (1.0 - P[i].initial_daugther_mass / P[i].Mass);      
                                                                /*! the decayed fraction w.r.p to the total mass */
 double mass_remain_ratio = 1.0 - mass_ratio;

 if(All.ComovingIntegrationOn)
   {
     tstart = (P[i].Ti_begstep + P[i].Ti_endstep) / 2;	        /*! midpoint of old step */
     tend = P[i].Ti_endstep + ti_step / 2;	                    /*! midpoint of new step */
     dt_gravkick1 = get_gravkick_factor(tstart, P[i].Ti_endstep);
     dt_gravkick2 = get_gravkick_factor(P[i].Ti_endstep, tend);
   }
 else
   {
     tstart = (P[i].Ti_begstep + P[i].Ti_endstep) / 2;	        /*! midpoint of old step */
     tend = P[i].Ti_endstep + ti_step / 2;	                    /*! midpoint of new step */
     dt_gravkick1 = (P[i].Ti_endstep - tstart) * All.Timebase_interval;
     dt_gravkick2 = (tend - P[i].Ti_endstep) * All.Timebase_interval;
   }


 for(j = 0; j < 3; j++)
   half_kick_v[j] = P[i].Vel[j];
 
 do_half_kick(half_kick_v, P[i].GravAccel, dt_gravkick1);
 
 for(j = 0; j < DECAY_DAUGHT_PART_NUM; j++)
   {
     set_part_info_for_daught(i, daught_i, half_kick_v, ti_step, mass_ratio);
     do_half_kick(P[daught_i].Vel, P[daught_i].GravAccel, dt_gravkick2);
#ifdef EXTRA_POTENTIAL
     do_half_kick(P[daught_i].Vel, ex_pot_accel, dt_gravkick2);
#endif
     daught_i++;
   } 

/*! modify parent particle */
  P[i].Mass *= mass_remain_ratio;          
  P[i].decay_factor = 0.0;          
  P[i].Potential *= mass_remain_ratio;
}                                       /*! end create_deacy_daught_part_and_kick */

void set_part_info_for_daught(int i, int daught_i, FLOAT *half_v, int ti_step, double mass_ratio)
{ 
  int j;
  double atime;

  if(All.ComovingIntegrationOn)
    atime = All.Time;
  else
    atime = 1;
  
  for(j = 0; j < 3; j++)
     {
       P[daught_i].Pos[j] = P[i].Pos[j];
       P[daught_i].GravAccel[j] = P[i].GravAccel[j];
#ifdef PMGRID
       P[daught_i].GravPM[j] = P[i].GravPM[j];
#endif
#ifdef FORCETEST
       P[daught_i].GravAccelDirect[j] = P[i].GravAccelDirect[j];
#endif
       P[daught_i].Vel[j] = half_v[j];
     }	

  P[daught_i].Potential = P[i].Potential * mass_ratio / (double) DECAY_DAUGHT_PART_NUM;		
  P[daught_i].OldAcc = P[i].OldAcc;			
  P[daught_i].Type = P[i].Type;	
  	                  
#ifdef FLEXSTEPS
  P[daught_i].FlexStepGrp = P[i].FlexStepGrp;		
#endif
  P[daught_i].GravCost = P[i].GravCost;		
#ifdef PSEUDOSYMMETRIC
  P[daught_i].AphysOld = P[i].AphysOld;  
#endif   

  P[daught_i].Ti_endstep = P[i].Ti_endstep + ti_step;                
  P[daught_i].Ti_begstep = P[i].Ti_endstep;     

  do_random_kick(P[daught_i].Pos, 0); /* particle born at the center */

  P[daught_i].Mass = P[i].Mass * mass_ratio / (double) DECAY_DAUGHT_PART_NUM;

  do_random_kick(P[daught_i].Vel, atime * DECAY_KICK_VEL);
					
        
#ifdef HALO_PART_RANDOM_KICK
  P[daught_i].l_eff = 0.0;                    
#endif
  P[daught_i].decay_flag = 0;           /*! no more decay for the daught particle */                
  P[daught_i].decay_factor = 0.0;   
  
  if(All.ComovingIntegrationOn)        
    P[daught_i].a_born = All.Time;
  else 
    P[daught_i].a_born = 1;   
}                                         /*! end copy_info_for_daught */


double age_of_the_system()
{
  gsl_function F;
  gsl_integration_workspace *workspace;
  double result0,result1,abserr;
  double a0 = 1e-8;
  double b0 = All.TimeBegin;
  double b1 = All.TimeMax;

  if(All.ComovingIntegrationOn)
  {
    workspace = gsl_integration_workspace_alloc(100000);
    F.function = &interg_a;
  
    gsl_integration_qag(&F,a0, b0, 0, 1.0e-8, 100000, GSL_INTEG_GAUSS41, workspace, &result0, &abserr);
    gsl_integration_qag(&F,a0, b1, 0, 1.0e-8, 100000, GSL_INTEG_GAUSS41, workspace, &result1, &abserr);
  
    gsl_integration_workspace_free(workspace);

    return (result1 - result0) / All.Hubble;
  }
  else
    return (b1 - b0); 
}

double interg_a(double a, void *param)
{
  return 1.0 / sqrt(All.Omega0 / a + (1 - All.Omega0) * a * a);
}

void update_decay_softening()
{
  int i;  
  double time_factor = -1 * 0.69314718056 * system_decay_time / DECAY_HALF_LIFE;
  double m_1 = undecay_part_mass * (1.0 - INITIAL_DM_FRAC) * exp(time_factor) + undecay_part_mass * INITIAL_DM_FRAC;
	                                                               /* the mean mass of the big simulation particles */
  double decay_gs_factor;
  
  if(All.ComovingIntegrationOn)
   decay_gs_factor = 2 * All.G * m_1 * All.Time / DECAY_KICK_VEL / DECAY_KICK_VEL;
  else
   decay_gs_factor = 2 * All.G * m_1 / DECAY_KICK_VEL / DECAY_KICK_VEL;
   
  for(i=0; i<NumPart; i++)
    {
	if(P[i].decay_flag == 0)	
	  {
	  P[i].soften_length = decay_gs_factor / P[i].a_born / P[i].a_born;	 
	  
	  if(P[i].soften_length < All.SofteningTable[1])    /* the table is refreshed in gravity_tree(void) before tree reconstructure. -dlcheng */
	    P[i].soften_length = All.SofteningTable[1];
	    		  
	  P[i].soften_length = P[i].soften_length * 2.8;		/*! the soften_length is defined now as the force softening -dlcheng */		    	  	 
      }
    else
      P[i].soften_length = All.ForceSoftening[1];		    
    }  
}	
#endif
