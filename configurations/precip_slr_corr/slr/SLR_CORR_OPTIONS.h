C CPP options file for SLR_CORR package
C Use this file for selecting options within the SLR_CORR package

#ifdef ALLOW_SLR_CORR

C balance sea level rise (each time step by default)
#define ALLOW_SLR_CORR_BALANCE

C balance sea level rise over a running average (not by default)
#define ALLOW_SLR_CORR_SMOOTH_BALANCE


#endif /* ALLOW_SLR_CORR */
