new ;

"--------------------------------------------" ;
" Program to check the procedure DCHOKERN.SRC" ;
"--------------------------------------------" ;

library pgraph myprocs ;

@ Simulating data: normal @
nobs = 10000 ;
meanx = 1 ;
sdx = 2 ;
a0 = 1 ;
a1 = 1.0 ;
a2 = -0.1 ;
sdeps = 1 ;

@ Simulations @
seed = 7387905 ;
x = meanx + sdx*rndns(nobs,1,seed) ;
eps = sdeps*rndns(nobs,1,seed) ;
y = a0 + a1*x + a2*x.*x - eps ;
y = (y.>0) ;


@ Values of x where we estimate Pr(Y=1|X) @
xval = pctiles(x,seqa(1,1,99)) ;

@ True P(X) @
ptrue = a0 + a1*xval + a2*xval.*xval ;
ptrue = cdfn(ptrue/sdeps) ;

@ Kernel m(x) @
pkern = dchokern(y,x,xval) ;

perror = 100*(pkern-ptrue)./ptrue ;

title("True and Kernel probabilities") ;
xlabel("x") ;
ylabel("P(x)") ;
xy(xval,ptrue~pkern) ;

/*
title("% estimation error") ;
xlabel("x") ;
ylabel("% error") ;
xy(xval,perror) ;
*/





