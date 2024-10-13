new ;
closeall ;

/**************************************************/
/* Example that calls the procedure MIPROBIT.SRC  */
/**************************************************/
library myprocs pgraph ;

/****************/
/* 1. Constants */
/****************/
nobs = 5000 ;
npar = 3 ;

meanx = ones(1,npar) ;
sdx = 2 * rndu(1,npar) ;
bmat = (1 | 1 | 2 | 3) ;

/******************/
/* 2. Simulations */
/******************/
xobs = meanx + sdx.*rndn(nobs,npar) ;
xobs = ones(nobs,1)~xobs ;
eps = rndu(nobs,1) ;
eps = cdfni(eps) ;
yobs = xobs * bmat + eps ;
yobs = (yobs.>0) ;  

/*****************/
/* 3. Estimation */
/*****************/
namesb = "b0" | "b1" | "b2" | "b3" ;

{best, varest} = miprobit(yobs,xobs,namesb,1) ;

end ;