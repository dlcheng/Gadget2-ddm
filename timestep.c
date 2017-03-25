#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include "allvars.h"
#include "proto.h"

/*! \file timestep.c 
 *  \brief routines for 'kicking' particles in momentum space and assigning new timesteps
 */

static double fac1, fac2, fac3, hubble_a, atime, a3inv;
static double dt_displacement = 0;

#ifdef DECAY_DARK_MATTER                /* -dlcheng */
static double ln_2 = 0.69314718056;
#endif

#ifdef RANDOM_KICK_TIMESTEP             /* -dlcheng */
static  FLOAT  vel_t[3];
static  FLOAT  accel_t[3];
#endif

/*! This function advances the system in momentum space, i.e. it does apply
 *  the 'kick' operation after the forces have been computed. Additionally, it
 *  assigns new timesteps to particles. At start-up, a half-timestep is
 *  carried out, as well as at the end of the simulation. In between, the
 *  half-step kick that ends the previous timestep and the half-step kick for
 *  the new timestep are combined into one operation.
 */
void advance_and_find_timesteps(void)
{
  int i, j, no, ti_step, ti_min, tend, tstart;
  double dt_entr, dt_entr2, dt_gravkick, dt_hydrokick, dt_gravkick2, dt_hydrokick2, t0, t1;
  double minentropy, aphys;
  FLOAT dv[3];

#ifdef DECAY_DARK_MATTER               /* -dlcheng */
  int local_decay_num = 0;             /*!  number of particles decayed in this run -dlcheng          */
  unsigned int new_part_num, id_start;
  int temp_flag = 0;
  
  new_part_num = 0;
  id_start = 0;

  for(i = 0; i < NTask; i++)
    {
      decay_tab_in[i] = 0;   /*! initialize the table, used to set up the ID for New particles   -dlcheng */
      decay_tab_out[i] = 0;
    }
#endif

#ifdef PETER_HALO_DECAY
  double decay_p;
  unsigned int peter_local_decay_num = 0;
  unsigned int peter_decay_num;
#endif

#ifdef HALO_PART_RANDOM_KICK  /* -dlcheng */
  int kick_num = 0;
#ifdef RANDOM_KICK_TIMESTEP   /* -dlcheng */
  double half_t;
#endif
#endif

#ifdef FLEXSTEPS
  int ti_grp;
#endif
#if defined(PSEUDOSYMMETRIC) && !defined(FLEXSTEPS)
  double apred, prob;
  int ti_step2;
#endif
#ifdef PMGRID
  double dt_gravkickA, dt_gravkickB;
#endif
#ifdef MAKEGLASS
  double disp, dispmax, globmax, dmean, fac, disp2sum, globdisp2sum;
#endif

  t0 = second();

  if(All.ComovingIntegrationOn)
    {
      fac1 = 1 / (All.Time * All.Time);
      fac2 = 1 / pow(All.Time, 3 * GAMMA - 2);
      fac3 = pow(All.Time, 3 * (1 - GAMMA) / 2.0);
      hubble_a = All.Omega0 / (All.Time * All.Time * All.Time)
	+ (1 - All.Omega0 - All.OmegaLambda) / (All.Time * All.Time) + All.OmegaLambda;

      hubble_a = All.Hubble * sqrt(hubble_a);
      a3inv = 1 / (All.Time * All.Time * All.Time);
      atime = All.Time;
    }
  else
    fac1 = fac2 = fac3 = hubble_a = a3inv = atime = 1;

#ifdef DECAY_DARK_MATTER
   if(All.ComovingIntegrationOn)
    {
    system_decay_time += All.TimeStep / All.Time / hubble_a;
    }
  else
    system_decay_time += All.TimeStep;
#endif

#ifdef NOPMSTEPADJUSTMENT
  dt_displacement = All.MaxSizeTimestep;
#else
  if(Flag_FullStep || dt_displacement == 0)
    find_dt_displacement_constraint(hubble_a * atime * atime);
#endif

#ifdef PMGRID
  if(All.ComovingIntegrationOn)
    dt_gravkickB = get_gravkick_factor(All.PM_Ti_begstep, All.Ti_Current) -
      get_gravkick_factor(All.PM_Ti_begstep, (All.PM_Ti_begstep + All.PM_Ti_endstep) / 2);
  else
    dt_gravkickB = (All.Ti_Current - (All.PM_Ti_begstep + All.PM_Ti_endstep) / 2) * All.Timebase_interval;

  if(All.PM_Ti_endstep == All.Ti_Current)	/* need to do long-range kick */
    {
      /* make sure that we reconstruct the domain/tree next time because we don't kick the tree nodes in this case */
      All.NumForcesSinceLastDomainDecomp = 1 + All.TotNumPart * All.TreeDomainUpdateFrequency;
    }
#endif


#ifdef MAKEGLASS
  for(i = 0, dispmax = 0, disp2sum = 0; i < NumPart; i++)
    {
      for(j = 0; j < 3; j++)
	{
	  P[i].GravPM[j] *= -1;
	  P[i].GravAccel[j] *= -1;
	  P[i].GravAccel[j] += P[i].GravPM[j];
	  P[i].GravPM[j] = 0;
	}

      disp = sqrt(P[i].GravAccel[0] * P[i].GravAccel[0] +
		  P[i].GravAccel[1] * P[i].GravAccel[1] + P[i].GravAccel[2] * P[i].GravAccel[2]);

      disp *= 2.0 / (3 * All.Hubble * All.Hubble);

      disp2sum += disp * disp;

      if(disp > dispmax)
	dispmax = disp;
    }

  MPI_Allreduce(&dispmax, &globmax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(&disp2sum, &globdisp2sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  dmean = pow(P[0].Mass / (All.Omega0 * 3 * All.Hubble * All.Hubble / (8 * M_PI * All.G)), 1.0 / 3);

  if(globmax > dmean)
    fac = dmean / globmax;
  else
    fac = 1.0;

  if(ThisTask == 0)
    {
      printf("\nglass-making:  dmean= %g  global disp-maximum= %g  rms= %g\n\n",
	     dmean, globmax, sqrt(globdisp2sum / All.TotNumPart));
      fflush(stdout);
    }

  for(i = 0, dispmax = 0; i < NumPart; i++)
    {
      for(j = 0; j < 3; j++)
	{
	  P[i].Vel[j] = 0;
	  P[i].Pos[j] += fac * P[i].GravAccel[j] * 2.0 / (3 * All.Hubble * All.Hubble);
	  P[i].GravAccel[j] = 0;
	}
    }
#endif




  /* Now assign new timesteps and kick */

#ifdef FLEXSTEPS
  if((All.Ti_Current % (4 * All.PresentMinStep)) == 0)
    if(All.PresentMinStep < TIMEBASE)
      All.PresentMinStep *= 2;

  for(i = 0; i < NumPart; i++)
    {
      if(P[i].Ti_endstep == All.Ti_Current)
	{
	  ti_step = get_timestep(i, &aphys, 0);

	  /* make it a power 2 subdivision */
	  ti_min = TIMEBASE;
	  while(ti_min > ti_step)
	    ti_min >>= 1;
	  ti_step = ti_min;

	  if(ti_step < All.PresentMinStep)
	    All.PresentMinStep = ti_step;
	}
    }

  ti_step = All.PresentMinStep;
  MPI_Allreduce(&ti_step, &All.PresentMinStep, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);

  if(dt_displacement < All.MaxSizeTimestep)
    ti_step = (int) (dt_displacement / All.Timebase_interval);
  else
    ti_step = (int) (All.MaxSizeTimestep / All.Timebase_interval);

  /* make it a power 2 subdivision */
  ti_min = TIMEBASE;
  while(ti_min > ti_step)
    ti_min >>= 1;
  All.PresentMaxStep = ti_min;


  if(ThisTask == 0)
    printf("Syn Range = %g  PresentMinStep = %d  PresentMaxStep = %d \n",
	   (double) All.PresentMaxStep / All.PresentMinStep, All.PresentMinStep, All.PresentMaxStep);

#endif


  for(i = 0; i < NumPart; i++)
    {
      if(P[i].Ti_endstep == All.Ti_Current)
	{
#ifdef EXTRA_POTENTIAL    /* initialize the potential acceleration  -dlcheng */
#ifdef AD_EXTRA_POTENTIAL_GROWTH
          if(All.Ti_Current * All.Timebase_interval > AD_TIME)
             (*extra_potential_accel)(P[i].Pos, ex_pot_accel);
          else
             (*extra_potential_accel_ad)(P[i].Pos, ex_pot_accel, AD_TIME);
#else
#ifdef AD_EXTRA_POTENTIAL_RELEASE
          if(All.Ti_Current * All.Timebase_interval > AD_TIME)
             ex_pot_accel[0] = ex_pot_accel[1] = ex_pot_accel[2] = 0;
          else
             (*extra_potential_accel_ad)(P[i].Pos, ex_pot_accel, AD_TIME);
#else
          (*extra_potential_accel)(P[i].Pos, ex_pot_accel);
#endif
#endif       
#endif   


#ifdef HALO_PART_RANDOM_KICK   /* -dlcheng */
          if(P[i].l_eff > L_KICK) 
            { 
              do_random_kick(P[i].Vel, VMAX_KICK); /* apply the random kick  */
              P[i].l_eff = 0.0;                    /* set l_eff to zero      */
              kick_num++;         
            } 
#ifdef RANDOM_KICK_TIMESTEP
            vel_kick_copy(P[i].Vel, vel_t);                                             /* copy vel_t from Vel */
#ifdef EXTRA_POTENTIAL
            accel_kick_add(P[i].GravAccel, ex_pot_accel, accel_t);                      /* accel_t = GravAccel + ex_pot_accel */
#else
            accel_t[0] = P[i].GravAccel[0];                                             /* else accel_t = P[i].GravAccel */
            accel_t[1] = P[i].GravAccel[1];
            accel_t[2] = P[i].GravAccel[2];
#endif
            half_t = (P[i].Ti_endstep - P[i].Ti_begstep) * All.Timebase_interval / 2.0; /* do half kick         */
            do_half_kick(vel_t, accel_t, half_t);                                       /* change vel_t                   */
#endif
#endif
	  ti_step = get_timestep(i, &aphys, 0);

	  /* make it a power 2 subdivision */
	  ti_min = TIMEBASE;
	  while(ti_min > ti_step)
	    ti_min >>= 1;
	  ti_step = ti_min;

#ifdef FLEXSTEPS
	  ti_grp = P[i].FlexStepGrp % All.PresentMaxStep;
	  ti_grp = (ti_grp / All.PresentMinStep) * All.PresentMinStep;
	  ti_step = ((P[i].Ti_endstep + ti_grp + ti_step) / ti_step) * ti_step - (P[i].Ti_endstep + ti_grp);
#else

#ifdef PSEUDOSYMMETRIC
	  if(P[i].Type != 0)
	    {
	      if(P[i].Ti_endstep > P[i].Ti_begstep)
		{
		  apred = aphys + ((aphys - P[i].AphysOld) / (P[i].Ti_endstep - P[i].Ti_begstep)) * ti_step;
		  if(fabs(apred - aphys) < 0.5 * aphys)
		    {
		      ti_step2 = get_timestep(i, &apred, -1);
		      ti_min = TIMEBASE;
		      while(ti_min > ti_step2)
			ti_min >>= 1;
		      ti_step2 = ti_min;

		      if(ti_step2 < ti_step)
			{
			  get_timestep(i, &apred, ti_step);
			  prob =
			    ((apred - aphys) / (aphys - P[i].AphysOld) * (P[i].Ti_endstep -
									  P[i].Ti_begstep)) / ti_step;
			  if(prob < get_random_number(P[i].ID))
			    ti_step /= 2;
			}
		      else if(ti_step2 > ti_step)
			{
			  get_timestep(i, &apred, 2 * ti_step);
			  prob =
			    ((apred - aphys) / (aphys - P[i].AphysOld) * (P[i].Ti_endstep -
									  P[i].Ti_begstep)) / ti_step;
			  if(prob < get_random_number(P[i].ID + 1))
			    ti_step *= 2;
			}
		    }
		}
	      P[i].AphysOld = aphys;
	    }
#endif

#ifdef SYNCHRONIZATION
	  if(ti_step > (P[i].Ti_endstep - P[i].Ti_begstep))	/* timestep wants to increase */
	    {
	      if(((TIMEBASE - P[i].Ti_endstep) % ti_step) > 0)
		ti_step = P[i].Ti_endstep - P[i].Ti_begstep;	/* leave at old step */
	    }
#endif
#endif /* end of FLEXSTEPS */

	  if(All.Ti_Current == TIMEBASE)	/* we here finish the last timestep. */
	    ti_step = 0;

	  if((TIMEBASE - All.Ti_Current) < ti_step)	/* check that we don't run beyond the end */
	    ti_step = TIMEBASE - All.Ti_Current;

	  tstart = (P[i].Ti_begstep + P[i].Ti_endstep) / 2;	/* midpoint of old step */
	  tend = P[i].Ti_endstep + ti_step / 2;	/* midpoint of new step */

	  if(All.ComovingIntegrationOn)
	    {
	      dt_entr = (tend - tstart) * All.Timebase_interval;
	      dt_entr2 = (tend - P[i].Ti_endstep) * All.Timebase_interval;
	      dt_gravkick = get_gravkick_factor(tstart, tend);
	      dt_hydrokick = get_hydrokick_factor(tstart, tend);
	      dt_gravkick2 = get_gravkick_factor(P[i].Ti_endstep, tend);
	      dt_hydrokick2 = get_hydrokick_factor(P[i].Ti_endstep, tend);
	    }
	  else
	    {
	      dt_entr = dt_gravkick = dt_hydrokick = (tend - tstart) * All.Timebase_interval;
	      dt_gravkick2 = dt_hydrokick2 = dt_entr2 = (tend - P[i].Ti_endstep) * All.Timebase_interval;
	    }

#ifdef DECAY_DARK_MATTER /* -dlcheng */
      if(P[i].Type == 1 && P[i].decay_flag == 1 && P[i].decay_factor >= ln_decay_threshold)
        {
        local_decay_num++;
        create_daught_and_kick(i, local_decay_num, ti_step);
        }
     if(P[i].Type == 1 && P[i].decay_flag == 1)
        P[i].decay_factor += ti_step * All.Timebase_interval / hubble_a * ln_2 / DECAY_HALF_LIFE;
#endif 

#ifdef PETER_HALO_DECAY
     if(P[i].peter_decay_flag)
	   {
       random_kick_tab_update();
	   decay_p = 0.69314718056 / DECAY_HALF_LIFE 
	           * (P[i].Ti_endstep - P[i].Ti_begstep) * All.Timebase_interval / hubble_a;
		
       if(random_kick_tab[0] <= decay_p)
	     {
	     do_random_kick(P[i].Vel, DECAY_KICK_VEL);
         P[i].peter_decay_flag = 0;
	     peter_local_decay_num++;
	     }
	   }
#endif

	  P[i].Ti_begstep = P[i].Ti_endstep;
	  P[i].Ti_endstep = P[i].Ti_begstep + ti_step;


	  /* do the kick */

	  for(j = 0; j < 3; j++)
	    {
	      dv[j] = P[i].GravAccel[j] * dt_gravkick;
	      P[i].Vel[j] += dv[j];
#ifdef EXTRA_POTENTIAL            /* -dlcheng */
              P[i].Vel[j] += ex_pot_accel[j] * dt_gravkick;
#endif
	    }
#ifdef HALO_PART_RANDOM_KICK       /* -dlcheng */
          P[i].l_eff += update_l_eff(P[i].Pos, P[i].Vel, ti_step * All.Timebase_interval); //update l_eff after this step
#endif

	  if(P[i].Type == 0)	/* SPH stuff */
	    {
	      for(j = 0; j < 3; j++)
		{
		  dv[j] += SphP[i].HydroAccel[j] * dt_hydrokick;
		  P[i].Vel[j] += SphP[i].HydroAccel[j] * dt_hydrokick;

		  SphP[i].VelPred[j] =
		    P[i].Vel[j] - dt_gravkick2 * P[i].GravAccel[j] - dt_hydrokick2 * SphP[i].HydroAccel[j];
#ifdef PMGRID
		  SphP[i].VelPred[j] += P[i].GravPM[j] * dt_gravkickB;
#endif
		}

	      /* In case of cooling, we prevent that the entropy (and
	         hence temperature decreases by more than a factor 0.5 */

	      if(SphP[i].DtEntropy * dt_entr > -0.5 * SphP[i].Entropy)
		SphP[i].Entropy += SphP[i].DtEntropy * dt_entr;
	      else
		SphP[i].Entropy *= 0.5;

	      if(All.MinEgySpec)
		{
		  minentropy = All.MinEgySpec * GAMMA_MINUS1 / pow(SphP[i].Density * a3inv, GAMMA_MINUS1);
		  if(SphP[i].Entropy < minentropy)
		    {
		      SphP[i].Entropy = minentropy;
		      SphP[i].DtEntropy = 0;
		    }
		}

	      /* In case the timestep increases in the new step, we
	         make sure that we do not 'overcool' when deriving
	         predicted temperatures. The maximum timespan over
	         which prediction can occur is ti_step/2, i.e. from
	         the middle to the end of the current step */

	      dt_entr = ti_step / 2 * All.Timebase_interval;
	      if(SphP[i].Entropy + SphP[i].DtEntropy * dt_entr < 0.5 * SphP[i].Entropy)
		SphP[i].DtEntropy = -0.5 * SphP[i].Entropy / dt_entr;
	    }


	  /* if tree is not going to be reconstructed, kick parent nodes dynamically.
	   */
	  if(All.NumForcesSinceLastDomainDecomp < All.TotNumPart * All.TreeDomainUpdateFrequency)
	    {
	      no = Father[i];
	      while(no >= 0)
		{
		  for(j = 0; j < 3; j++)
		    Extnodes[no].vs[j] += dv[j] * P[i].Mass / Nodes[no].u.d.mass;

		  no = Nodes[no].u.d.father;
		}
	    }
	}
    }

/*! FINSH THE SYNCHRONIZATION AT THIS POINT -dlcheng */
#ifdef DECAY_DARK_MATTER                     /* -dlcheng */

  decay_tab_in[ThisTask] = local_decay_num;
  
  for(i = 0; i < NTask; i++)
    MPI_Allreduce(&decay_tab_in[i], &decay_tab_out[i], 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  for(i = 0; i < NTask; i++)
      new_part_num += decay_tab_out[i];

  if(ThisTask == 0)
    {
      printf("Using DL decay algorithm.\nDecayed particle number at this timestep: %d.\n", new_part_num);
      printf("Decay threshold: %g\n", decay_new_part_threshold);
      fflush(stdout);
    }

  new_part_num *= DECAY_DAUGHT_PART_NUM;

  for(i = 0; i < ThisTask; i++)
     id_start += decay_tab_out[i] * DECAY_DAUGHT_PART_NUM;     /* ID starts from 0, not 1 */
 
  id_start += All.TotNumPart;

  for(i = 0; i < local_decay_num * DECAY_DAUGHT_PART_NUM; i++)
    {
      P[NumPart + i].ID = id_start;
      id_start++;
    }

  Ntype[1] += new_part_num;
  All.TotNumPart += new_part_num;
  NtypeLocal[1] += local_decay_num * DECAY_DAUGHT_PART_NUM;
  NumPart += local_decay_num * DECAY_DAUGHT_PART_NUM;
  
  if(new_part_num > 0)
   {
    temp_flag = 1;
   }  /* if any particle decayed, set decay_tree_construct_flag = 1 for all CPUs */
   
  MPI_Allreduce(&temp_flag, &decay_tree_construct_flag, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
   
  if(decay_tree_construct_flag)
    All.NumForcesSinceLastDomainDecomp = 1 + All.TotNumPart * All.TreeDomainUpdateFrequency;
      /* ensure that the tree will be reconstructed */

  if(ThisTask == 0)
    {
      printf("Total number of particles: %ld.\n", (long int) All.TotNumPart);
      fflush(stdout);
    }
#endif

#ifdef PMGRID
  if(All.PM_Ti_endstep == All.Ti_Current)	/* need to do long-range kick */
    {
      ti_step = TIMEBASE;
      while(ti_step > (dt_displacement / All.Timebase_interval))
	ti_step >>= 1;

      if(ti_step > (All.PM_Ti_endstep - All.PM_Ti_begstep))	/* PM-timestep wants to increase */
	{
	  /* we only increase if an integer number of steps will bring us to the end */
	  if(((TIMEBASE - All.PM_Ti_endstep) % ti_step) > 0)
	    ti_step = All.PM_Ti_endstep - All.PM_Ti_begstep;	/* leave at old step */
	}

      if(All.Ti_Current == TIMEBASE)	/* we here finish the last timestep. */
	ti_step = 0;

      tstart = (All.PM_Ti_begstep + All.PM_Ti_endstep) / 2;
      tend = All.PM_Ti_endstep + ti_step / 2;

      if(All.ComovingIntegrationOn)
	dt_gravkick = get_gravkick_factor(tstart, tend);
      else
	dt_gravkick = (tend - tstart) * All.Timebase_interval;

      All.PM_Ti_begstep = All.PM_Ti_endstep;
      All.PM_Ti_endstep = All.PM_Ti_begstep + ti_step;

      if(All.ComovingIntegrationOn)
	dt_gravkickB = -get_gravkick_factor(All.PM_Ti_begstep, (All.PM_Ti_begstep + All.PM_Ti_endstep) / 2);
      else
	dt_gravkickB =
	  -((All.PM_Ti_begstep + All.PM_Ti_endstep) / 2 - All.PM_Ti_begstep) * All.Timebase_interval;

      for(i = 0; i < NumPart; i++)
	{
	  for(j = 0; j < 3; j++)	/* do the kick */
	    P[i].Vel[j] += P[i].GravPM[j] * dt_gravkick;

	  if(P[i].Type == 0)
	    {
	      if(All.ComovingIntegrationOn)
		{
		  dt_gravkickA = get_gravkick_factor(P[i].Ti_begstep, All.Ti_Current) -
		    get_gravkick_factor(P[i].Ti_begstep, (P[i].Ti_begstep + P[i].Ti_endstep) / 2);
		  dt_hydrokick = get_hydrokick_factor(P[i].Ti_begstep, All.Ti_Current) -
		    get_hydrokick_factor(P[i].Ti_begstep, (P[i].Ti_begstep + P[i].Ti_endstep) / 2);
		}
	      else
		dt_gravkickA = dt_hydrokick =
		  (All.Ti_Current - (P[i].Ti_begstep + P[i].Ti_endstep) / 2) * All.Timebase_interval;

	      for(j = 0; j < 3; j++)
		SphP[i].VelPred[j] = P[i].Vel[j]
		  + P[i].GravAccel[j] * dt_gravkickA
		  + SphP[i].HydroAccel[j] * dt_hydrokick + P[i].GravPM[j] * dt_gravkickB;
	    }
	}
    }
#endif

#ifdef HALO_PART_RANDOM_KICK  /* -dlcheng */
 int kick_num2;
 MPI_Allreduce(&kick_num, &kick_num2, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

 if(ThisTask == 0)
   printf("%d particles have been randomly kicked at this timestep.\n",kick_num2);
#endif

#ifdef PETER_HALO_DECAY   /* -dlcheng */
 MPI_Allreduce(&peter_local_decay_num, &peter_decay_num, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

 if(ThisTask == 0)
   printf("Using Peter et al. (2010) decay algorithm\n.Decayed particles number: %d.\n",peter_decay_num);
#endif

  t1 = second();
  All.CPU_TimeLine += timediff(t0, t1);
}




/*! This function normally (for flag==0) returns the maximum allowed timestep
 *  of a particle, expressed in terms of the integer mapping that is used to
 *  represent the total simulated timespan. The physical acceleration is
 *  returned in `aphys'. The latter is used in conjunction with the
 *  PSEUDOSYMMETRIC integration option, which also makes of the second
 *  function of get_timestep. When it is called with a finite timestep for
 *  flag, it returns the physical acceleration that would lead to this
 *  timestep, assuming timestep criterion 0.
 */
int get_timestep(int p,		/*!< particle index */
		 double *aphys,	/*!< acceleration (physical units) */
		 int flag	/*!< either 0 for normal operation, or finite timestep to get corresponding
				   aphys */ )
{
  double ax, ay, az, ac, csnd;
  double dt = 0, dt_courant = 0, dt_accel;
  int ti_step;

#ifdef RANDOM_KICK_TIMESTEP  /* -dlcheng */
  double dt_kick;
  double l_left; /* L_KICK - l_eff  */
  double w_kick; /* the weight at this position */

  w_kick = kick_weight(P[p].Pos);
  l_left = (L_KICK - P[p].l_eff) / w_kick; /* the real length  */
  dt_kick = find_timestep_random_kick(vel_t, accel_t, l_left);
#endif 

#ifdef CONDUCTION
  double dt_cond;
#endif

  if(flag == 0)
    {
      ax = fac1 * P[p].GravAccel[0];
      ay = fac1 * P[p].GravAccel[1];
      az = fac1 * P[p].GravAccel[2];

#ifdef PMGRID
      ax += fac1 * P[p].GravPM[0];
      ay += fac1 * P[p].GravPM[1];
      az += fac1 * P[p].GravPM[2];
#endif

#ifdef EXTRA_POTENTIAL  /* -dlcheng */
      ax += fac1 * ex_pot_accel[0];
      ay += fac1 * ex_pot_accel[1];
      az += fac1 * ex_pot_accel[2];
#endif

      if(P[p].Type == 0)
	{
	  ax += fac2 * SphP[p].HydroAccel[0];
	  ay += fac2 * SphP[p].HydroAccel[1];
	  az += fac2 * SphP[p].HydroAccel[2];
	}

      ac = sqrt(ax * ax + ay * ay + az * az);	/* this is now the physical acceleration */
      *aphys = ac;
    }
  else
    ac = *aphys;

  if(ac == 0)
    ac = 1.0e-30;

  switch (All.TypeOfTimestepCriterion)
    {
    case 0:
      if(flag > 0)
	{
	  dt = flag * All.Timebase_interval;
	  dt /= hubble_a;	/* convert dloga to physical timestep  */
	  ac = 2 * All.ErrTolIntAccuracy * atime * All.SofteningTable[P[p].Type] / (dt * dt);
	  *aphys = ac;
	  return flag;
	}

#if defined(DECAY_DARK_MATTER) && defined(UNEQUALSOFTENINGS)          /* -dlcheng */
      dt = dt_accel = sqrt(2 * All.ErrTolIntAccuracy * atime * All.SofteningTable[P[p].Type] / ac); /*! now determined by large mass particle only -dlcheng */
#else
      dt = dt_accel = sqrt(2 * All.ErrTolIntAccuracy * atime * All.SofteningTable[P[p].Type] / ac);
#endif

#ifdef ADAPTIVE_GRAVSOFT_FORGAS
      if(P[p].Type == 0)
	dt = dt_accel = sqrt(2 * All.ErrTolIntAccuracy * atime * SphP[p].Hsml / 2.8 / ac);
#endif
      break;
    default:
      endrun(888);
      break;
    }

  if(P[p].Type == 0)
    {
      csnd = sqrt(GAMMA * SphP[p].Pressure / SphP[p].Density);

      if(All.ComovingIntegrationOn)
	dt_courant = 2 * All.CourantFac * All.Time * SphP[p].Hsml / (fac3 * SphP[p].MaxSignalVel);
      else
	dt_courant = 2 * All.CourantFac * SphP[p].Hsml / SphP[p].MaxSignalVel;

      if(dt_courant < dt)
	dt = dt_courant;
    }

  /* convert the physical timestep to dloga if needed. Note: If comoving integration has not been selected,
     hubble_a=1.
   */
  dt *= hubble_a;

  if(dt >= All.MaxSizeTimestep)
    dt = All.MaxSizeTimestep;

  if(dt >= dt_displacement)
    dt = dt_displacement;

#ifdef RANDOM_KICK_TIMESTEP  /* -dlcheng */
  if(dt >= dt_kick)
    dt = dt_kick;
#endif 

  if(dt < All.MinSizeTimestep)
    {
#ifndef NOSTOP_WHEN_BELOW_MINTIMESTEP
      printf("warning: Timestep wants to be below the limit `MinSizeTimestep'\n");

      if(P[p].Type == 0)
	{
	  printf
	    ("Part-ID=%d  dt=%g dtc=%g ac=%g xyz=(%g|%g|%g)  hsml=%g  maxsignalvel=%g dt0=%g eps=%g\n",
	     (int) P[p].ID, dt, dt_courant * hubble_a, ac, P[p].Pos[0], P[p].Pos[1], P[p].Pos[2],
	     SphP[p].Hsml, SphP[p].MaxSignalVel,
	     sqrt(2 * All.ErrTolIntAccuracy * atime * All.SofteningTable[P[p].Type] / ac) * hubble_a,
	     All.SofteningTable[P[p].Type]);
	}
      else
	{
	  printf("Part-ID=%d  dt=%g ac=%g xyz=(%g|%g|%g)\n", (int) P[p].ID, dt, ac, P[p].Pos[0], P[p].Pos[1],
		 P[p].Pos[2]);
	}
      fflush(stdout);
      endrun(888);
#endif
      dt = All.MinSizeTimestep;
    }

  ti_step = dt / All.Timebase_interval;

  if(!(ti_step > 0 && ti_step < TIMEBASE))
    {
      printf("\nError: A timestep of size zero was assigned on the integer timeline!\n"
	     "We better stop.\n"
	     "Task=%d Part-ID=%d dt=%g tibase=%g ti_step=%d ac=%g xyz=(%g|%g|%g) tree=(%g|%g%g)\n\n",
	     ThisTask, (int) P[p].ID, dt, All.Timebase_interval, ti_step, ac,
	     P[p].Pos[0], P[p].Pos[1], P[p].Pos[2], P[p].GravAccel[0], P[p].GravAccel[1], P[p].GravAccel[2]);
#ifdef PMGRID
      printf("pm_force=(%g|%g|%g)\n", P[p].GravPM[0], P[p].GravPM[1], P[p].GravPM[2]);
#endif
      if(P[p].Type == 0)
	printf("hydro-frc=(%g|%g|%g)\n", SphP[p].HydroAccel[0], SphP[p].HydroAccel[1], SphP[p].HydroAccel[2]);

      fflush(stdout);
      endrun(818);
    }

  return ti_step;
}


/*! This function computes an upper limit ('dt_displacement') to the global
 *  timestep of the system based on the rms velocities of particles. For
 *  cosmological simulations, the criterion used is that the rms displacement
 *  should be at most a fraction MaxRMSDisplacementFac of the mean particle
 *  separation. Note that the latter is estimated using the assigned particle
 *  masses, separately for each particle type. If comoving integration is not
 *  used, the function imposes no constraint on the timestep.
 */
void find_dt_displacement_constraint(double hfac /*!<  should be  a^2*H(a)  */ )
{
  int i, j, type, *temp;
  int count[6];
  long long count_sum[6];
  double v[6], v_sum[6], mim[6], min_mass[6];
  double dt, dmean, asmth = 0;

  dt_displacement = All.MaxSizeTimestep;

  if(All.ComovingIntegrationOn)
    {
      for(type = 0; type < 6; type++)
	{
	  count[type] = 0;
	  v[type] = 0;
	  mim[type] = 1.0e30;
	}

      for(i = 0; i < NumPart; i++)
	{
	  v[P[i].Type] += P[i].Vel[0] * P[i].Vel[0] + P[i].Vel[1] * P[i].Vel[1] + P[i].Vel[2] * P[i].Vel[2];
	  if(mim[P[i].Type] > P[i].Mass)
	    mim[P[i].Type] = P[i].Mass;
	  count[P[i].Type]++;
	}

      MPI_Allreduce(v, v_sum, 6, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(mim, min_mass, 6, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

      temp = malloc(NTask * 6 * sizeof(int));
      MPI_Allgather(count, 6, MPI_INT, temp, 6, MPI_INT, MPI_COMM_WORLD);
      for(i = 0; i < 6; i++)
	{
	  count_sum[i] = 0;
	  for(j = 0; j < NTask; j++)
	    count_sum[i] += temp[j * 6 + i];
	}
      free(temp);

      for(type = 0; type < 6; type++)
	{
	  if(count_sum[type] > 0)
	    {
	      if(type == 0)
		dmean =
		  pow(min_mass[type] / (All.OmegaBaryon * 3 * All.Hubble * All.Hubble / (8 * M_PI * All.G)),
		      1.0 / 3);
	      else
		dmean =
		  pow(min_mass[type] /
		      ((All.Omega0 - All.OmegaBaryon) * 3 * All.Hubble * All.Hubble / (8 * M_PI * All.G)),
		      1.0 / 3);

	      dt = All.MaxRMSDisplacementFac * hfac * dmean / sqrt(v_sum[type] / count_sum[type]);

#ifdef PMGRID
	      asmth = All.Asmth[0];
#ifdef PLACEHIGHRESREGION
	      if(((1 << type) & (PLACEHIGHRESREGION)))
		asmth = All.Asmth[1];
#endif
	      if(asmth < dmean)
		dt = All.MaxRMSDisplacementFac * hfac * asmth / sqrt(v_sum[type] / count_sum[type]);
#endif

	      if(ThisTask == 0)
		printf("type=%d  dmean=%g asmth=%g minmass=%g a=%g  sqrt(<p^2>)=%g  dlogmax=%g\n",
		       type, dmean, asmth, min_mass[type], All.Time, sqrt(v_sum[type] / count_sum[type]), dt);

	      if(dt < dt_displacement)
		dt_displacement = dt;
	    }
	}

      if(ThisTask == 0)
	printf("displacement time constraint: %g  (%g)\n", dt_displacement, All.MaxSizeTimestep);
    }
}
