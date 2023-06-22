# ABNORMAL END: S/R CALC_R_STAR

When using R_STAR coordinates, there is a potential for the model to crash if the layer thickness gets too small or too big. This is particularly prevalent when using sea ice next to the coast.

When this occurs, it's helpful to find where in the model the error occurred in order to adjust the model (e.g. potentially change the topography). The model will print a message where the error occurs, but only in the STDERR file for the tile. To find the tile(s) where the error happened, search for the error with `grep`:

```
grep CALC_R_STAR STDERR*
```
This will show which tiles the error occurred in:
```
STDERR.0112:STOP in CALC_R_STAR : too SMALL rStarFac[C,W,S] !
STDERR.0128:STOP in CALC_R_STAR : too SMALL rStarFac[C,W,S] !
```
Inside the STDERR files, you can see where the error occured, e.g.:
```
fail at i,j=  46   2 ; rStarFacW,H,eta =  0.097682  1.000000E+01 -5.381625E+00 -1.266474E+01
WARNING: r*FacW < hFacInf at       1 pts : bi,bj,Thid,Iter=   1   1   1    241260
STOP in CALC_R_STAR : too SMALL rStarFac[C,W,S] !
```
