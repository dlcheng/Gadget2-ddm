#if defined(DECAY_DARK_MATTER) || defined(PETER_HALO_DECAY)

 #define MEM_INCREASE_FACTOR         24                          /*! the frequency of the splitting * daughter particle number */
 #define DECAY_DAUGHT_PART_NUM       1                           /*! the number of new particles generated after decay */
 #define DECAY_HALF_LIFE             (13.79 * All.HubbleParam)   /*! in unit of gadget time = 977792988.06 years */
                                                                 /* (13.79 * All.HubbleParam) for f=0.5
                                                                  * (26.80 * All.HubbleParam) for f=0.3
                                                                  * (90.72 * All.HubbleParam) for f=0.1
                                                                  */
 #define DECAY_KICK_VEL              200                         /*! in unit of km/s for peculiar velocity */
 #define INITIAL_DM_FRAC             0                           /*! the initial fraction of daughter particles */
#endif

#if defined(AD_EXTRA_POTENTIAL_GROWTH) || defined(AD_EXTRA_POTENTIAL_RELEASE)
 #define AD_TIME                     0.306601419451              /*! 300Myr in internal time unit, with h=1           */
#endif

#ifdef HALO_PART_RANDOM_KICK

 #define REF_WEIGHT_KICK             1.0                         /*! reference kick weight at 1kpc/h */
 #define SLOPE_KICK                  0.333333                    /*! the pow index of kick weight    */
 #define L_KICK                      15.0                        /*! define the max effective length to do a kick, in unit of kpc/h, 
                                                                     h=1 */
 #define VMAX_KICK                   200.0                       /*! the magnitude of kick velocity in interval unit km/s  */

#endif

#ifdef DECAY_DARK_MATTER
#undef  ADAPTIVE_GRAVSOFT_FORGAS  
#define UNEQUALSOFTENINGS
#endif
