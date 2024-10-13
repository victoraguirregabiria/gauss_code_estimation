new ;
closeall ;
cls ;

//  ****************************************************************************
//
//  am_econometrica_2007_montecarlo.prg
//
//      GAUSS program that implements the Monte Carlo experiments reported in 
//      Section 4 of the paper: Aguirregabiria, V., and P. Mira (2007): 
//      "Sequential Estimation of Dynamic Discrete Games," Econometrica, 
//      Vol. 75, No. 1, pages 1–53.
//
//  by Victor Aguirregabiria
//
//  - Last version: May 2005
//  - Comments and explanations added on August 2011
//
//  ****************************************************************************


// ------------------------------------------------------------------------------
//                              MODEL
// ------------------------------------------------------------------------------
// 
//  MAIN FEATURES OF THE MODEL 
// 
//  - Dynamic game of firm entry-exit in a market
// 
//  -   i = Firm index, i belongs to (1,2,3,4,5)
//      t = Time index
//  
//  -   Decision variable:
//      a[it] = Indicator of the event 'firm i operates in the market at period t'
// 
//          a[it] = 0 ----> No active in market
//          a[it] = 1 ----> Active in market
// 
//  -   The game is dynamic because there is a sunk cost of entry
//
//  -   The state variables of this game are the indicators of whether
//      the firms were active or not in the market at previous period, 
//      and market size. These are payoff relevant state variables because
//      they determine whether a firm has to pay an entry cost to operate 
//      in the market.
// 
//  -   State variables
//      x[t] = ( s[t], a[1,t-1], a[2,t-1], a[3,t-1], a[4,t-1], a[5,t-1] )
//  
//      e0[it] and e1[it] = Firm i's private information.
// 
//  - Profit function
// 
//      If a[it] = 0 (firm i is not active at period t)
// 
//          Profit(0) = e0[it]
// 
//      If a[it] = 1 (firm i is not active at period t)
// 
//          Profit(1) = theta_fc_i - theta_ec * (1-a[i,t-1]) 
//                    + theta_rs * s[t]
//                    - theta_rn * ln(N[-it]) + 1) + e1[it]
//  where:
//          theta_fc_i, theta_ec, theta_rs, theta_rn are parameters
//
//          theta_fc_i  =   Firm i's fixed effect
//          theta_ec    =   Entry cost
//          N[-it]      =   Number of firms, other than i, operating in the market at period t. 
//
//  eps_i(0) and eps_i(1) are private information shocks that are
//  i.i.d. over time, markets, and players, and independent of each other,
//  with Extreme Value Type 1 distribution.
// 
//
// ------------------------------------------------------------------------------
//              MONTE CARLO EXPERIMENTS
// ------------------------------------------------------------------------------
//
//  - We implement 6 experiments.
//
//  - The following parameters are the same in the 6 experiments:
//
//      Number of local markets (M) =   400
//      Number of time periods  (T) =   1
//      Number of players (N)       =   5
//      Number of Monte Carlo simulations   =   1000
//
//      Discount factor     =   0.95
//      theta_fc_1          =   -1.9
//      theta_fc_2          =   -1.8
//      theta_fc_3          =   -1.7
//      theta_fc_4          =   -1.6
//      theta_fc_5          =   -1.5
//      theta_rs            =   1.0
//
//      Values for s[t]     =   (1,2,3,4,5)
//
//      Transition probability of s[t]  =   ( 0.8 ~ 0.2 ~ 0.0 ~ 0.0 ~ 0.0 )
//                                        | ( 0.2 ~ 0.6 ~ 0.2 ~ 0.0 ~ 0.0 )
//                                        | ( 0.0 ~ 0.2 ~ 0.6 ~ 0.2 ~ 0.0 )
//                                        | ( 0.0 ~ 0.0 ~ 0.2 ~ 0.6 ~ 0.2 )
//                                        | ( 0.0 ~ 0.0 ~ 0.0 ~ 0.2 ~ 0.8 ) 
//
//  - The following values of the parameters theta_rn and theta_ec
//    defines experiments 1 to 6:
//
//          Experiment 1:   theta_ec = 1.0      and     theta_rn = 0.0
//          Experiment 2:   theta_ec = 1.0      and     theta_rn = 1.0
//          Experiment 3:   theta_ec = 1.0      and     theta_rn = 2.0
//          Experiment 4:   theta_ec = 0.0      and     theta_rn = 1.0
//          Experiment 5:   theta_ec = 2.0      and     theta_rn = 1.0
//          Experiment 6:   theta_ec = 4.0      and     theta_rn = 1.0
//
//
//
// ------------------------------------------------------------------------------
//
//          STRUCTURE OF THIS PROGRAM
//
//              0. SELECTION OF THE MONTE CARLO EXPERIMENT TO IMPLEMENTS
//
//              1. VALUES OF PARAMETERS AND OTHER CONSTANTS
//
//              2. PROCEDURE FOR COMPUTATION OF MARKOV PERFECT EQUILIBRIUM
//                      - MPEPROB
//
//              3. COMPUTING A MARKOV PERFECT EQUILIBRIUM OF THIS DYNAMIC GAME
//
//              4. PROCEDURE TO SIMULATE DATA FROM THE COMPUTED EQUILIBRIUM
//                      - SIMDYGAM
//
//              5.  SIMULATING DATA (50,000 OBSERVATIONS) FROM THE EQUILIBRIUM 
//                  TO OBTAIN DESCRIPTIVE STATISTICS ON THE DYNAMICS OF MARKET STRUCTURE
//                  (TABLE 2 IN AGUIRREGABIRIA AND MIRA, 2007)
//
//              6.  PROCEDURES FOR THE ESTIMATION OF THE MODEL
//                      - FREQPROB
//                      - MILOGIT
//                      - CLOGIT
//                      - NPLDYGAM
//
//              7.  CHECKING FOR THE CONSISTENCY OF THE ESTIMATORS 
//                  (ESTIMATION WITH A VERY LARGE SAMPLE)
//
//              8.  MONTE CARLO EXPERIMENT
//
//              9.  SAVING RESULTS OF THE MONTE CARLO EXPERIMENT
//
//              10. STATISTICS FROM THE MONTE CARLO EXPERIMENTS
//                  TABLES 4 AND 5 IN AGUIRREGABIRIA AND MIRA (2007)
//
// --------------------------------------------------------------------------

// ------------------------------------------------------------------------------



// -----------------------------------------------------------
//  0. SELECTION OF THE MONTE CARLO EXPERIMENT TO IMPLEMENTS
// -----------------------------------------------------------
selexper = 1 ;  //  Value from 1 to 6 that represents the index of the Monte Carlo experiment
                //  to implement. A run of this program implements one experiment.

fileout  = "am_econometrica_2007_montecarlo_1.out" ;        // Name of output file
wdir = "c:\\mypapers\\dyngames_econometrica\\progau\\" ;    // Name of default directory

// Names of files for saving results of Monte Carlo experiment
nfile_bNP   = "bnp_exper_01" ; // Names of file with estimates when initial CCPs = frequency estimator
nfile_bSP   = "bsp_exper_01" ; // Names of file with estimates when initial CCPs = logit estimation
nfile_bR    = "br_exper_01"  ; // Names of file with estimates when initial CCPs = random U(0,1)
nfile_btrue = "btrue_exper_01"  ; // Names of file with estimates when initial CCPs = true CCPs


// -----------------------------------------------
//  1. VALUES OF PARAMETERS AND OTHER CONSTANTS
// -----------------------------------------------
nobs = 400 ;    // Number of markets (observations)
nrepli = 1000 ; // Number of Monte carlo replications
nplayer = 5 ;   //  Number of players
numexp = 6 ;    //  Total number of Monte Carlo experiments

theta_fc = zeros(numexp,nplayer) ;
theta_fc[.,1] = -1.9 * ones(numexp,1); // Vector with values of theta_fc_1 for each experiment
theta_fc[.,2] = -1.8 * ones(numexp,1); // Vector with values of theta_fc_2 for each experiment
theta_fc[.,3] = -1.7 * ones(numexp,1); // Vector with values of theta_fc_3 for each experiment
theta_fc[.,4] = -1.6 * ones(numexp,1); // Vector with values of theta_fc_4 for each experiment
theta_fc[.,5] = -1.5 * ones(numexp,1); // Vector with values of theta_fc_1 for each experiment

theta_rs = 1.0 * ones(numexp,1); // Vector with values of theta_rs for each experiment
disfact = 0.95 * ones(numexp,1); // Vector with values of discount factor for each experiment
sigmaeps = 1 * ones(numexp,1);  // Vector with values of std. dev. epsilon for each experiment

theta_rn = (0.0 | 1.0 | 2.0 | 1.0 | 1.0 | 1.0); // Vector with values of theta_rn for each experiment
theta_ec = (1.0 | 1.0 | 1.0 | 0.0 | 2.0 | 4.0); // Vector with values of theta_rn for each experiment

// Points of support and transition probability of state variable s[t], market size
sval = seqa(1,1,5) ;    // Support of market size
numsval = rows(sval) ;  // Number of possible market sizes
nstate = numsval * (2^nplayer) ;    // Number of points in the state space
ptrans = ( 0.8 ~ 0.2 ~ 0.0 ~ 0.0 ~ 0.0 )
       | ( 0.2 ~ 0.6 ~ 0.2 ~ 0.0 ~ 0.0 )
       | ( 0.0 ~ 0.2 ~ 0.6 ~ 0.2 ~ 0.0 )
       | ( 0.0 ~ 0.0 ~ 0.2 ~ 0.6 ~ 0.2 )
       | ( 0.0 ~ 0.0 ~ 0.0 ~ 0.2 ~ 0.8 ) ;

// Selecting the parameters for the experiment
theta_fc = theta_fc[selexper,.]' ;
theta_rs = theta_rs[selexper] ;
disfact = disfact[selexper] ;
sigmaeps = sigmaeps[selexper] ;
theta_ec = theta_ec[selexper] ;
theta_rn = theta_rn[selexper] ;

// Vector with true values of parameters
trueparam = (theta_fc | theta_rs | theta_rn | theta_ec | disfact | sigmaeps) ;

// Vector with names of parameters of profit function
namesb = ("FC_1" | "FC_2" | "FC_3" | "FC_4" | "FC_5" | "RS" | "RN" | "EC" ) ;

//  Openning output file and calling library of graphs 
buff = changedir(wdir) ;        // It sets default directory
output file = ^fileout reset ;  // It creates and opens output file 
format /mb1 /ros 16,4 ;         // It sets format of numerical results
library pgraph ;

// Seed for (pseudo) random number generation
rndseed 5333799  ;   

"*****************************************************************************************" ;
"*****************************************************************************************" ;
"               MONTE CARLO EXPERIMENT #";; selexper;
"*****************************************************************************************" ;
"*****************************************************************************************" ;

// -----------------------------------------------------------------------
//  2. PROCEDURE FOR COMPUTATION OF MARKOV PERFECT EQUILIBRIUM
//          - MPEPROB
// -----------------------------------------------------------------------

// ----------------------------------------------------------------------------
//
//  MPEPROB.SRC     
//      Procedure that computes players' equilibrium conditional choice probabilities (CCPs)
//      in a dynamic game of firms' market entry/exit with incomplete information.
//      Policy iteration algorithm (iterations in best response mapping).
//                      
//  by Victor Aguirregabiria
// 
// ----------------------------------------------------------------------------
// 
//  FORMAT:
//  { prob, psteady, mstate, dconv } = mpeprob(inip,maxiter)
// 
//  INPUTS:
//
//      inip    - (numx x nplayer) matrix of probabilities to initialize the algorithm
// 
//      maxiter -   Maximum number of policy iterations
// 
//  OUTPUTS:
//      prob   - (numx x nplayer) matrix with MPE probs of entry
//                The states in the rows of prob are ordered as follows.
//                Example: sval=(1|2) and 3 players:
//                         s[t]    a[1,t-1]    a[2,t-1]    a[3,t-1]
//                Row 1:     1           0           0           0
//                Row 2:     1           0           0           1
//                Row 3:     1           0           1           0
//                Row 4:     1           0           1           1
//                Row 5:     1           1           0           0
//                Row 6:     1           1           0           1
//                Row 7:     1           1           1           0
//                Row 8:     1           1           1           1
//                Row 9:     2           0           0           0
//                Row 10:    2           0           0           1
//                Row 11:    2           0           1           0
//                Row 12:    2           0           1           1
//                Row 13:    2           1           0           0
//                Row 14:    2           1           0           1
//                Row 15:    2           1           1           0
//                Row 16:    2           1           1           1
// 
//       psteady - (numx x 1) vector with steady-state distribution of {s[t],a[t-1]}
// 
//       mstate  - (numx x (nplayer+1)) matrix with values of state variables {s[t],a[t-1]}. 
//                 The states in the rows are ordered as in the matrix "prob" decribed above.
// 
//       dconv   - Indicator for convergence: 
//                   dconv = 1 ===> convergence achieved
//                   dconv = 0 ===> no convergence
// 
// ----------------------------------------------------------------------------

proc (4) = mpeprob(inip,maxiter);
  local nplayer, nums, numa, numx, aval, mstate,
        i, prob0, critconv, criter, dconv, iter, prob1,
        mi, ppi, ptranai0, ptranai1, hi,        
        profit1, profit0, ptrana, iptran0, iptran1,
        j, buff, v0, cbell, v1, maxv1, 
        psteady0, psteady1, ptran ;
        
  nums = rows(sval) ;  
  nplayer = cols(inip) ;
  numa = 2^nplayer ;
  numx = nums * numa ;

    "" ;
    "*****************************************************************************************" ;
    "   COMPUTING A MPE OF THE DYNAMIC GAME";
    "*****************************************************************************************" ;
    "" ;
    "----------------------------------------------------------------------------------------" ;
    "       Values of the structural parameters" ;
    "" ;
    "                       Fixed cost firm 1   =";; theta_fc[1];
    "                       Fixed cost firm 2   =";; theta_fc[2];
    "                       Fixed cost firm 3   =";; theta_fc[3];
    "                       Fixed cost firm 4   =";; theta_fc[4];
    "                       Fixed cost firm 5   =";; theta_fc[5];
    "       Parameter of market size (theta_rs) =";; theta_rs;
    "Parameter of competition effect (theta_rn) =";; theta_rn;
    "                     Entry cost (theta_ec) =";; theta_ec;
    "                       Discount factor     =";; disfact ;
    "                    Std. Dev. epsilons     =";; sigmaeps ;
    "" ;
    "----------------------------------------------------------------------------------------" ;
    "" ;
    "       BEST RESPONSE MAPPING ITERATIONS" ;
    "" ;

  // ---------------------------------------------------
  // a. Matrices with values of states of (s[t],a[t-1]) 
  // ---------------------------------------------------
  aval = zeros(numa,nplayer) ; 
  i=1 ;
  do while i<=nplayer ;
    aval[.,i] = ones(2^(i-1),1).*.(0|1).*.ones(2^(nplayer-i),1) ;               
    i=i+1 ;
  endo ; 
  mstate = zeros(numx,nplayer+1) ;
  mstate[.,1] = sval.*.ones(numa,1) ;
  mstate[.,2:nplayer+1] = ones(nums,1).*.aval ;

  // --------------------------------------- 
  // b. Initializing vector of probabilities 
  // ---------------------------------------
  prob0 = inip ;
  critconv = (1e-3)*(1/numx) ;
  criter = 1000 ;
  dconv = 1 ;

  // ---------------------- 
  // c. Iterative algorithm 
  // ---------------------- 
  iter=1 ;
  do while (criter>critconv).and(iter<=maxiter) ;
    "         Best response mapping iteration  =";; iter ;
    "         Convergence criterion =";; criter ;
    "" ;
    prob1 = prob0 ; 
    // ----------------------------------------------------
    // c.1. Matrix of transition probs Pr(a[t]|s[t],a[t-1]) 
    // ----------------------------------------------------
    ptrana = (prob1[.,1].^(aval[.,1]')).*((1-prob1[.,1]).^(1-aval[.,1]')) ;
    i=2 ;
    do while i<=nplayer ;
      ptrana = ptrana.*(prob1[.,i].^(aval[.,i]')).*((1-prob1[.,i]).^(1-aval[.,i]')) ;
      i=i+1 ;
    endo ;        
        
    i=1 ;    
    do while i<=nplayer ;
      // ----------------------------------------
      // c.2. Matrices Pr(a[t]|s[t],a[t-1],ai[t]) 
      // ----------------------------------------
      mi = aval[.,i]' ;
      ppi = prob1[.,i] ;
      iptran0 = (1-mi)./((ppi.^mi).*((1-ppi).^(1-mi))) ;
      iptran0 = ptrana .* iptran0 ;
      iptran1 = mi./((ppi.^mi).*((1-ppi).^(1-mi))) ;
      iptran1 = ptrana .* iptran1 ;
      clear mi ;
               
      // --------------------------------------
      // c.3. Computing hi = E[ln(N[-it]) + 1)]
      // --------------------------------------
      hi = aval ;
      hi[.,i] = ones(numa,1) ;
      hi= iptran1 * ln(sumc(hi')) ;
      
      // ----------------------------------------
      //  c.4. Matrices with expected profit firm i 
      // ---------------------------------------- 
      profit1 = theta_fc[i] 
              + theta_rs * mstate[.,1] 
              - theta_ec*(1-mstate[.,i+1])
              - theta_rn * hi ;              // Profit if firm is active
      profit0 = zeros(numx,1) ;           // Profit if firm is not active

                                
      // ----------------------------------------------------
      //  c.5. Transition probabilities for firm i              
      //          Pr(x[t+1],a[t] | x[t],a[t-1],ai[t]=0)       
      //          Pr(x[t+1],a[t] | x[t],a[t-1],ai[t]=1)       
      // ---------------------------------------------------- 
      iptran0 = (ptrans.*.ones(numa,numa)).*(ones(1,nums).*.iptran0);
      iptran1 = (ptrans.*.ones(numa,numa)).*(ones(1,nums).*.iptran1);
      
      // ---------------------------------------------
      //  c.6. Computing value function for firm i
      // ---------------------------------------------
      v0 = zeros(numx,1) ;
      cbell = 1000 ;
      do while cbell>critconv ;      
        v1 = (profit0 + disfact * iptran0 * v0)
           ~ (profit1 + disfact * iptran1 * v0) ;
        v1 = v1./sigmaeps ;
        maxv1 = maxc(v1') ;
        v1 = v1 - maxv1 ;
        v1 = sigmaeps * ( maxv1 + ln(exp(v1[.,1]) + exp(v1[.,2])) ) ;
        cbell = maxc( abs(v1-v0) ) ;
        v0 = v1 ;
      endo ;
      
      // --------------------------------------
      // c.7. Updating probabilities for firm i 
      // --------------------------------------
      v1 = (profit0 + disfact * iptran0 * v0)
         ~ (profit1 + disfact * iptran1 * v0) ;      
      v1 = v1./sigmaeps ;
      maxv1 = maxc(v1') ;
      v1 = v1 - maxv1 ;    
      prob1[.,i] = exp(v1[.,2])./(exp(v1[.,1])+exp(v1[.,2])) ;                   
      i=i+1 ;
    endo ;    

    criter = maxc(maxc(abs(prob1-prob0))) ;
    prob0 = prob1 ;
    iter = iter + 1 ;
  endo ;  
  clear iptran0, iptran1, v0, v1 ;

  if (criter>critconv) ;
    dconv = 0 ;
    psteady1 = 0 ;    
    "----------------------------------------------------------------------------------------" ;
    "         CONVERGENCE NOT ACHIEVED AFTER";; iter;; "BEST RESPONSE ITERATIONS";
    "----------------------------------------------------------------------------------------" ;
  else ;
    ptrana = (prob1[.,1].^(aval[.,1]')).*((1-prob1[.,1]).^(1-aval[.,1]')) ;
    i=2 ;
    do while i<=nplayer ;
      ptrana = ptrana.*(prob1[.,i].^(aval[.,i]')).*((1-prob1[.,i]).^(1-aval[.,i]')) ;
      i=i+1 ;
    endo ;
    ptrana = (ptrans.*.ones(numa,numa)).*(ones(1,nums).*.ptrana) ;   
    criter = 1000 ;
    psteady0 = (1/numx) * ones(numx,1) ;
    do while criter>critconv ;      
      psteady1 = ptrana' * psteady0 ;
      criter = maxc( abs(psteady1-psteady0) ) ;
      psteady0 = psteady1 ;
    endo ; 
    "----------------------------------------------------------------------------------------" ;
    "         CONVERGENCE ACHIEVED AFTER";; iter;; "BEST RESPONSE ITERATIONS";
    "----------------------------------------------------------------------------------------" ;
    "         EQUILIBRIUM PROBABILITIES" ;
    prob1 ;
    "----------------------------------------------------------------------------------------" ;
  endif ;  
    
  retp(prob1,psteady1,mstate,dconv) ;
endp ;


// -----------------------------------------------------------------------
//  3. COMPUTING A MARKOV PERFECT EQUILIBRIUM OF THIS DYNAMIC GAME
// -----------------------------------------------------------------------
maxiter = 200 ;    // Maximum number of Policy iterations
prob0 = 0.5*rndu(nstate,nplayer) ;
{pequil, psteady, vstate, dconv} = mpeprob(prob0,maxiter);


// -----------------------------------------------------------------------
//  4. PROCEDURE TO SIMULATE DATA FROM THE COMPUTED EQUILIBRIUM
//          - SIMDYGAM
// -----------------------------------------------------------------------

// -----------------------------------------------------------------------
//   SIMDYGAM.SRC    Simulates data of state and decision variables
//                   from the steady-state distribution of 
//                   a Markov Perfect Equilibrium in a dynamic
//                   game of firms' market entry/exit with 
//                   incomplete information
// 
//   by Victor Aguirregabiria
// 
// ----------------------------------------------------------------------------
// 
//   FORMAT:
//   { aobs, aobs_1, indsobs } = simdygam(nobs,pequil,psteady,mstate)
// 
//   INPUTS
//       nobs    - Number of simulations (markets)
// 
//       pchoice - (nstate x nplayer) matrix with MPE probs of entry
// 
//       psteady - (nstate x 1) vector with steady-state distribution of {s[t],a[t-1]}
// 
//       mstate  - (nstate x (nplayer+1)) matrix with values of state variables {s[t],a[t-1]}. 
//                The states in the rows of prob are ordered as follows.
//                Example: sval=(1|2) and 3 players:
//                         s[t]    a[1,t-1]    a[2,t-1]    a[3,t-1]
//                Row 1:     1           0           0           0
//                Row 2:     1           0           0           1
//                Row 3:     1           0           1           0
//                Row 4:     1           0           1           1
//                Row 5:     1           1           0           0
//                Row 6:     1           1           0           1
//                Row 7:     1           1           1           0
//                Row 8:     1           1           1           1
//                Row 9:     2           0           0           0
//                Row 10:    2           0           0           1
//                Row 11:    2           0           1           0
//                Row 12:    2           0           1           1
//                Row 13:    2           1           0           0
//                Row 14:    2           1           0           1
//                Row 15:    2           1           1           0
//                Row 16:    2           1           1           1
// 
// 
//   OUTPUTS:
//       aobs   - (nobs x nplayer) matrix with players' choices.
// 
//       aobs_1 - (nobs x nplayer) matrix with players' initial states.
// 
//       indsobs - (nobs x 1) vector with simulated values of s[t]
// 
// ----------------------------------------------------------------------------

proc (3) = simdygam(nobs,pchoice,psteady,mstate) ;
  local nplay, nums, numa, numx, i, pbuff0, pbuff1, uobs, xobs, aobs_1, aobs ;
        
  nplay = cols(pchoice) ;
  nums = rows(pchoice) ;
  numa = 2^nplay ;
  numx = nums/numa ;
 
  // ----------------------------------------------------------------------
  // a. Generating random draws from ergodic distribution of (s[t],a[t-1])
  // ----------------------------------------------------------------------
  pbuff1 = cumsumc(psteady) ;
  pbuff0 = cumsumc((0|psteady[1:nums-1])) ;
  uobs = rndu(nobs,1) ;
  uobs = ((uobs.>=(pbuff0')).*(uobs.<=(pbuff1'))) ;
  uobs = uobs * seqa(1,1,nums) ;
  xobs = mstate[uobs,1] ;
  aobs_1 = mstate[uobs,2:nplay+1] ;
  clear pbuff0, pbuff1 ;
  
  //  --------------------------------------------------------
  //  b. Generating random draws for a[t] (given s[t],a[t-1]) 
  //  --------------------------------------------------------
  pchoice = pchoice[uobs,.] ;
  uobs = rndu(nobs,nplay) ;
  aobs = (uobs.<=pchoice) ;
      
  retp(aobs, aobs_1, xobs) ;
endp ;

// --------------------------------------------------------------------------
//  5.  SIMULATING DATA (50,000 OBSERVATIONS) FROM THE EQUILIBRIUM 
//      TO OBTAIN DESCRIPTIVE STATISTICS ON THE DYNAMICS OF MARKET STRUCTURE
// --------------------------------------------------------------------------
nobsfordes = 50000 ;
{aobs,aobs_1,sobs} = simdygam(nobsfordes,pequil, psteady, vstate) ;

"" ;
"*****************************************************************************************" ;
"   DESCRIPTIVE STATISTICS FROM THE EQUILIBRIUM" ;
"   BASED ON";; nobsfordes "OBSERVATIONS";
"" ;
"   TABLE 2 OF THE PAPER AGUIRREGABIRIA AND MIRA (2007)" ;
"*****************************************************************************************" ;
"" ;
nf = sumc(aobs') ;      // Number of active firms in the market at t
nf_1 = sumc(aobs_1') ;  // Number of active firms in the market at t-1

//  Regression of (number of firms t) on (number of firms t-1)
__output = 0 ;  // For not printing the output of the regression below
{vnam,m,b,stb,vc,stderr,sigma,cx,rsq,resid,dwstat} = ols(0,nf,nf_1) ;
bareg_nf = b[2] ;  // Estimate of autorregressive parameter
entries = sumc((aobs.*(1-aobs_1))') ;   // Number of new entrants at t
exits = sumc(((1-aobs).*aobs_1)') ;     // Number of firm exits at t
excess = meanc(entries+exits-abs(entries-exits)) ; // Excess turnover
buff = corrx(entries~exits) ;    
corr_ent_exit = buff[1,2] ; // Correlation entries and exits
freq_active = meanc(aobs) ; // Frequencies of being active

"" ;
"----------------------------------------------------------------------------------------" ;
"       (1)    Average number of active firms   =";; meanc(nf) ;
"----------------------------------------------------------------------------------------" ;
"       (2)    Std. Dev. number of firms        =";; stdc(nf) ;
"----------------------------------------------------------------------------------------" ;
"       (3)    Regression N[t] on N[t-1]        =";; bareg_nf ;
"----------------------------------------------------------------------------------------" ;
"       (4)    Average number of entrants       =";; meanc(entries) ;
"----------------------------------------------------------------------------------------" ;
"       (5)    Average number of exits          =";; meanc(exits) ;
"----------------------------------------------------------------------------------------" ;
"       (6)    Excess turnover (in # of firms)  =";; excess ;
"----------------------------------------------------------------------------------------" ;
"       (7)    Correlation entries and exits    =";; corr_ent_exit ;
"----------------------------------------------------------------------------------------" ;
"       (8)    Frequencies of being active      =" ;; freq_active ;
"----------------------------------------------------------------------------------------" ;
"" ;

// --------------------------------------------------------------------------
//  6.  PROCEDURES FOR THE ESTIMATION OF THE MODEL
//          - FREQPROB
//          - MILOGIT
//          - CLOGIT
//          - NPLDYGAM
// --------------------------------------------------------------------------


// --------------------------------------------------------------------------
//  FREQPROB.SRC   Procedure that obtains a frequency estimation
//                 of Prob(Y|X) where Y is a vector of binary 
//                 variables and X is a vector of discrete variables
//
//  by Victor Aguirregabiria
//
// --------------------------------------------------------------------------
//
//  FORMAT:
//          freqp = freqprob(yobs,xobs,xval) 
//
//  INPUTS:
//
//      yobs    - (nobs x q) vector with sample observations 
//                of Y = Y1 ~ Y2 ~ ... ~ Yq
//            
//      xobs    - (nobs x k) matrix with sample observations of X
//
//      xval    - (numx x k) matrix with the values of X for which
//                we want to estimate Prob(Y|X).
//
//  OUTPUTS:
//
//      freqp   - (numx x q) vector with frequency estimates of
//                Pr(Y|X) for each value in xval.
//                Pr(Y1=1|X) ~ Pr(Y2=1|X) ~ ... ~ Pr(Yq=1|X) 
//
// --------------------------------------------------------------------------

proc (1) = freqprob(yobs,xobs,xval) ;
  local numx, numq, prob1, t, selx, denom, numer ; 
  numx = rows(xval) ;
  numq = cols(yobs) ;
  prob1 = zeros(numx,numq) ;
  t=1 ;
  do while t<=numx ;
    selx = prodc((xobs.==xval[t,.])') ;
    denom = sumc(selx) ;
    if (denom==0) ;
      prob1[t,.] = zeros(1,numq) ;
    else ;
      numer = sumc(selx.*yobs) ;
      prob1[t,.] = (numer')./denom ;
    endif ;
    t=t+1 ;
  endo ;
  retp(prob1) ;
endp ;


// -----------------------------------------------------------------
// MILOGIT - Estimation of a Logit Model by Maximum Likelihood
//           The optimization algorithm is a Newton's method.
//
// by Victor Aguirregabiria
//
// -----------------------------------------------------------------
//
// Format      {best,varest} = milogit(ydum,x,namesb)
//
// Input        ydum    - vector of observations of the dependent variable
//              x       - matrix of explanatory variables
//              namesb  - vector with names of parameters
//
// Output       best    - ML estimates
//              varest  - estimate of the covariance matrix
// -----------------------------------------------------------------

proc (1) = loglogit(ydum,x,b) ;
  local expxb, Fxb, llik, myzero ;
  myzero = 1E-12 ;
  expxb = exp(-x*b) ;
  Fxb = 1./(1+expxb) ;
  Fxb = Fxb + (myzero - Fxb).*(Fxb.<myzero)
            + (1-myzero - Fxb).*(Fxb.>1-myzero);
  llik = ydum'*ln(Fxb) + (1-ydum)'*ln(1-Fxb) ;
  retp(llik) ;
endp ;

proc (2) = milogit(ydum,x,namesb) ;
  local nobs, nparam, eps1, eps2, iter, criter1, 
        criter2, expxb0, Fxb0, phixb0, lamdab0, 
        dlogLb0, d2logLb0, b0, b1, lamda0, lamda1, 
        Avarb, sdb, tstat, llike, numy1, numy0, 
        logL0, LRI, pseudoR2, k ;

  format /mb1 /ros 16,6 ;

  nobs = rows(ydum) ;
  nparam = cols(x) ;
  eps1 = 1E-4 ;
  eps2 = 1E-2 ;
  b0 = zeros(nparam,1) ;  
  iter=1 ;
  criter1 = 1000 ;
  criter2 = 1000 ;

  do while (criter1>eps1).or(criter2>eps2) ;
/*
    "" ;
    "Iteration                = " iter ;
    "Log-Likelihood function  = " loglogit(ydum,x,b0) ;
    "Norm of b(k)-b(k-1)      = " criter1 ;
    "Norm of Gradient         = " criter2 ;
    "" ;
*/
    expxb0 = exp(-x*b0) ;
    Fxb0 = 1./(1+expxb0) ;
    dlogLb0 = x'*(ydum - Fxb0) ;
    d2logLb0 =  ( (Fxb0.*(1-Fxb0)).*x )'*x ;
    b1 = b0 + invpd(d2logLb0)*dlogLb0 ;
    criter1 = sqrt( (b1-b0)'*(b1-b0) ) ;
    criter2 = sqrt( dlogLb0'dlogLb0 ) ;
    b0 = b1 ;
    iter = iter + 1 ;
  endo ;

  expxb0 = exp(-x*b0) ;
  Fxb0 = 1./(1+expxb0) ;
  Avarb = - ( (Fxb0.*(1-Fxb0)).*x )'*x ;
  Avarb = inv(-Avarb) ;
  sdb    = sqrt(diag(Avarb)) ;
  tstat  = b0./sdb ;
  llike  = loglogit(ydum,x,b0) ;
  numy1  = sumc(ydum) ;
  numy0  = nobs - numy1 ;
  logL0  = numy1*ln(numy1) + numy0*ln(numy0) - nobs*ln(nobs) ;
  LRI    = 1 - llike/logL0 ;
  pseudoR2 = 1 - ( (ydum - Fxb0)'*(ydum - Fxb0) )/numy1 ;
/*
  "Number of Iterations     = " iter ;
  "Log-Likelihood function  = " llike ;
  "Akaike's AIC             = " 2*(nparam-llike) ;
  "Likelihood Ratio Index   = " LRI ;
  "Pseudo-R2                = " pseudoR2 ;
  "" ;
  "     ----------------------------------------------------------------";
  "         Parameter       Estimate        Standard        t-ratios";
  "                                         Errors" ;
  "     ----------------------------------------------------------------";
  k=1;
  do while k<=nparam;
    print $namesb[k];;b0[k];;sdb[k];;tstat[k];
    k=k+1 ;
  endo;
  "     ----------------------------------------------------------------";
*/
  retp(b0,avarb) ;
endp ;


// -------------------------------------------------------------------------------------
// CLOGIT -    Maximum Likelihood estimation of McFadden's Conditional Logit
//
//                  Optimization algorithm: Newton's method with analytical 
//                  gradient and hessian
//
// by Victor Aguirregabiria
//
// -------------------------------------------------------------------------------------
//
// Format      {best,varest} = clogit(ydum,x,restx,namesb)
//
// Input        ydum    - (nobs x 1) vector of observations of dependet variable
//                        Categorical variable with values: {1, 2, ..., nalt}
//
//              x       - (nobs x (k * nalt)) matrix of explanatory variables
//                        associated with unrestricted parameters.
//                        First k columns correspond to alternative 1, and so on
//
//              restx   - (nobs x nalt) vector of the sum of the explanatory
//                        variables whose parameters are restricted to be
//                        equal to 1.
//
//              namesb  - (k x 1) vector with names of parameters
//
//
//  Output      best    - (k x 1) vector with ML estimates.
//
//              varest  - (k x k) matrix with estimate of covariance matrix
//
// -------------------------------------------------------------------------------------

proc (2) = clogit(ydum,x,restx,namesb) ;
  local cconvb, myzero, nobs, nalt, npar, xysum, j,
        iter, criter, llike, b0, phat, sumpx, xxm, xbuff,
        d1llike, d2llike, b1, Avarb, sdb, tstat, 
        numyj, logL0, lrindex ;

  cconvb = 1e-6 ;
  myzero = 1e-16 ;
  nobs = rows(ydum) ;
  nalt = maxc(ydum) ;
  npar = cols(x)/nalt ;
  if npar/=rows(namesb) ;
    "ERROR: Dimensions of x";; npar;; "and of names(b0)";; rows(namesb) ;;
    "do not match " ;
    end ;
  endif;

  xysum = 0 ;
  j=1;
  do while j<=nalt ;
    xysum = xysum + sumc( (ydum.==j).*x[.,npar*(j-1)+1:npar*j] ) ;
    j=j+1 ;
  endo ;

  iter=1 ;
  criter = 1000 ;
  llike = -nobs ;
  b0 = zeros(npar,1) ;

  do while (criter>cconvb) ;
/*  
    "" ;
    "Iteration                = " iter ;
    "Log-Likelihood function  = " llike ;
    "Norm of b(k)-b(k-1)      = " criter ;
    "" ;
*/  
    // Computing probabilities
    phat = zeros(nobs,nalt) ;
    j=1 ;
    do while j<=nalt ;
      phat[.,j] = x[.,npar*(j-1)+1:npar*j]*b0 + restx[.,j] ;
      j=j+1 ;
    endo ;
    phat = phat - maxc(phat') ;
    phat = exp(phat)./sumc(exp(phat')) ;

    // Computing xmean
    sumpx = zeros(nobs,1) ;
    xxm = 0 ;
    llike = 0 ;
    j=1;
    do while j<=nalt ;
      xbuff = x[.,npar*(j-1)+1:npar*j] ; 
      sumpx = sumpx + phat[.,j] .*xbuff ;
      xxm = xxm + (phat[.,j].*xbuff)'*xbuff ;
      llike = llike
            + sumc( (ydum.==j)
                    .* ln( (phat[.,j].> myzero).*phat[.,j]
                         + (phat[.,j].<=myzero).*myzero    ) ) ;
      j=j+1 ;
    endo ;

    // Computing gradient
    d1llike = xysum - sumc(sumpx) ;

    // Computing hessian 
    d2llike = - (xxm - sumpx'*sumpx) ;
    
    // Gauss iteration
    b1 = b0 - inv(d2llike)*d1llike ;
    criter = sqrt( (b1-b0)'*(b1-b0) ) ;
    b0 = b1 ;
    iter = iter + 1 ;
  endo ;

  Avarb  = inv(-d2llike) ;
  sdb    = sqrt(diag(Avarb)) ;
  tstat  = b0./sdb ;
  
  numyj  = sumc(ydum.==(seqa(1,1,nalt)')) ;
  logL0  = sumc(numyj.*ln(numyj./nobs)) ;
  lrindex = 1 - llike/logL0 ;
/*
  "---------------------------------------------------------------------";
  "Number of Iterations     = " iter ;
  "Number of observations   = " nobs ;
  "Log-Likelihood function  = " llike ;
  "Likelihood Ratio Index   = " lrindex ;
  "---------------------------------------------------------------------";
  "       Parameter         Estimate        Standard        t-ratios";
  "                                         Errors" ;
  "---------------------------------------------------------------------";
  j=1;
  do while j<=npar;
    print $namesb[j];;b0[j];;sdb[j];;tstat[j];
    j=j+1 ;
  endo;
  "---------------------------------------------------------------------";
*/
  retp(b0,Avarb) ;
endp ;


//  ----------------------------------------------------------------------------
//  NPLDYGAM.SRC   Procedure that estimates the structural parameters 
//                 of dynamic game of firms' entry/exit
//                 using a Nested Pseudo-Likelihood (NPL) algorithm
//
//  by Victor Aguirregabiria
//
//  ----------------------------------------------------------------------------
//
//  FORMAT:
//  {best,varb} = npldygam(aobs,zobs,aobs_1,zval,ptranz,pchoice,bdisc,kiter)
//
//  INPUTS:
//      aobs    - (nobs x nplayer) matrix with observations of
//                firms' activity decisions (1=active ; 0=no active)
//      zobs    - (nobs x 1) vector with observations of market
//                exogenous characteristics.
//      aobs_1  - (nobs x nplayer) matrix with observations of
//                firms' initial states (1=incumbent; 0=potential entrant)
//      zval    - (numz x 1) vector with values of market characteristics
//      ptranz  - (numz x numz) matrix of transition probabilities
//                of market characteristics.
//      pchoice - ((numz*2^nplayer) x nplayer) matrix of players'
//                choice probabilities used to initialize the procedure.
//      bdisc   - Discount factor
//      kiter   - Number of NPL iterations
//
//  OUTPUTS:
//      best    - (kparam x kiter) matrix with estimates of parameters
//                (theta_fc_1 | theta_fc_2 | ...    | theta_fc_n | 
//                 theta_rs   | theta_rn   | theta_ec)
//                First column is the 1st-stage NPL, second column is
//                the 2nd-stage NPL, and so on.
//      varb    - (kparam x kaparam*kiter) matrix with estimated 
//                covariance matrices of the NPL estimators. 
//
//  ----------------------------------------------------------------------------

proc (2) = npldygam(aobs,zobs,aobs_1,zval,ptranz,pchoice,bdisc,kiter);
  local eulerc, myzero, nobs, nplayer, numa, numz, numx, 
        kparam, best, varb, mstate, twop, indzobs, indobs,
        i, namesb, aval, iter, ptrana, invi_bf, j, buff, 
        uobs0, uobs1, eobs0, eobs1, mi, ppi, ptranai0, 
        ptranai1, hi, umat0, umat1, sumu, sume, ww, 
        utilda0, utilda1, etilda0, etilda1, tetaest, 
        varest, u0, u1, e0, e1 ;        
        
  // ---------------
  // Some constants 
  // ---------------
  eulerc = 0.5772 ;
  myzero = 1e-16 ;
  nobs = rows(aobs) ;
  nplayer = cols(aobs) ;
  numa = 2^nplayer ;
  numz = rows(zval) ;  
  numx = numz*numa ;
  kparam = nplayer + 3 ;
  best = zeros(kparam,kiter) ;
  varb = zeros(kparam,kparam*kiter) ;  
  namesb = (0 $+ "theta_fc" $+ ftocv(seqa(1,1,nplayer),1,0))
         | "theta_rs" | "theta_rn" | "theta_ec" ; 
   
  // ----------------------------------------------
  // Matrix with values of states of (s[t],a[t-1]) 
  // ----------------------------------------------
  aval = zeros(numa,nplayer) ; 
  i=1 ;
  do while i<=nplayer ;
    aval[.,i] = ones(2^(i-1),1)
            .*.((0|1).*.ones(2^(nplayer-i),1)) ;               
    i=i+1 ;
  endo ; 
  mstate = zeros(numx,nplayer+1) ;
  mstate[.,1] = zval.*.ones(numa,1) ;
  mstate[.,2:nplayer+1] = ones(numz,1).*.aval ;
  
  // ------------------------------------------------
  // Matrix with observed indexes of state variables 
  // ------------------------------------------------
  indzobs = (zobs.==(zval'))*seqa(1,1,numz) ;  
  twop = 2.^(nplayer-seqa(1,1,nplayer)') ;
  indobs = sumc((aobs_1.*twop)') ;  
  indobs = (indzobs-1).*(2^nplayer) + indobs + 1 ;
   
  // ------------- 
  // NPL algorithm 
  // ------------- 
  aobs = 1 + reshape(aobs',nobs*nplayer,1) ;
  u0 = zeros(numx*nplayer,kparam) ;
  u1 = zeros(numx*nplayer,kparam) ;
  e0 = zeros(numx*nplayer,1) ;
  e1 = zeros(numx*nplayer,1) ;    
  iter=1 ;
  do while iter<=kiter ;
/*
    "" ;
    " -----------------------------------------------------" ;
    " POLICY ITERATION ESTIMATOR: STAGE =" ;; iter ;
    " -----------------------------------------------------" ;
    "" ;    
*/
    // -----------------------------------------------------------
    // (a) Matrix of transition probabilities Pr(a[t]|s[t],a[t-1]) 
    // -----------------------------------------------------------
    ptrana = (pchoice[.,1].^(aval[.,1]')).*((1-pchoice[.,1]).^(1-aval[.,1]')) ;
    i=2 ;
    do while i<=nplayer ;
      ptrana = ptrana.*(pchoice[.,i].^(aval[.,i]')).*((1-pchoice[.,i]).^(1-aval[.,i]')) ;
      i=i+1 ;
    endo ;          
    
    // -----------------------
    //  (b) Inverse of I-b*F
    // -----------------------   
    invi_bf = (ptranz.*.ones(numa,numa)).*(ones(1,numz).*.ptrana);
    invi_bf = crout(eye(numx)-bdisc*invi_bf) ;        
    
    // -----------------------------------------
    //  (c) Construction of explanatory variables 
    // -----------------------------------------
    uobs0 = zeros(nobs*nplayer,kparam) ;
    uobs1 = zeros(nobs*nplayer,kparam) ;
    eobs0 = zeros(nobs*nplayer,1) ;
    eobs1 = zeros(nobs*nplayer,1) ;   
    i=1 ;
    do while i<=nplayer ;
      // --------------------------------------------
      //  (c.1) Matrices Pr(a[t] | s[t],a[t-1], ai[t]) 
      // --------------------------------------------
      mi = aval[.,i]' ;
      ppi = pchoice[.,i] ;
      ppi= (ppi.>=myzero).*(ppi.<=(1-myzero)).*ppi
         + (ppi.<myzero).*myzero
         + (ppi.>(1-myzero)).*(1-myzero) ;      
      ptranai0 = (1-mi)./((ppi.^mi).*((1-ppi).^(1-mi))) ;
      ptranai0 = ptrana .* ptranai0 ;
      ptranai1 = mi./((ppi.^mi).*((1-ppi).^(1-mi))) ;
      ptranai1 = ptrana .* ptranai1 ;
      clear mi ;
               
      // ------------------------------------
      //  (c.2) Computing hi = E(ln(Sum(aj)+1)) 
      // ------------------------------------
      hi = aval ;
      hi[.,i] = ones(numa,1) ;
      hi = ptranai1 * ln(sumc(hi')) ;
      
      // ---------------------------
      //  (c.3) Creating U0 and U1 
      // ---------------------------
      umat0 = zeros(numx,nplayer+3) ;
      umat1 = eye(nplayer) ;
      umat1 = umat1[i,.] ;
      umat1 = (ones(numx,nplayer).*umat1)
            ~ mstate[.,1]~(-hi)~(mstate[.,i+1]-1) ;
      clear hi ;
      
      // ------------------------------
      //  (c.4) Creating sumu and sume 
      // ------------------------------
      sumu = (1-ppi).*umat0 + ppi.*umat1 ;
      sume = (1-ppi).*(eulerc-ln(1-ppi)) + ppi.*(eulerc-ln(ppi));
      clear ppi ;

      // -------------------
      //  (c.5) Creating ww 
      // -------------------
      ww = qrtsol(sumu~sume,lowmat(invi_bf)) ;
      ww = qrsol(ww,upmat1(invi_bf)) ;
      clear sumu, sume ;
            
      // ----------------------------------
      //  (c.6) Creating utilda and etilda 
      // ----------------------------------
      ptranai0 = (ptranz.*.ones(numa,numa))
               .*(ones(1,numz).*.ptranai0) ;
      utilda0 = umat0 + bdisc*(ptranai0*ww[.,1:kparam]) ;
      etilda0 = bdisc*(ptranai0*ww[.,kparam+1]) ;
      clear umat0, ptranai0 ;
      
      ptranai1 = (ptranz.*.ones(numa,numa))
               .*(ones(1,numz).*.ptranai1) ;
      utilda1 = umat1 + bdisc*(ptranai1*ww[.,1:kparam]) ;
      etilda1 = bdisc*(ptranai1*ww[.,kparam+1]) ;
      clear umat1, ptranai1 ;
            
      // -------------------------------------------
      //  (c.7) Creating observations uobs and eobs 
      // -------------------------------------------
      uobs0[(i-1)*nobs+1:i*nobs,.] = utilda0[indobs,.] ;
      uobs1[(i-1)*nobs+1:i*nobs,.] = utilda1[indobs,.] ;
      eobs0[(i-1)*nobs+1:i*nobs,.] = etilda0[indobs,.] ;
      eobs1[(i-1)*nobs+1:i*nobs,.] = etilda1[indobs,.] ;      
      u0[(i-1)*numx+1:i*numx,.] = utilda0 ;
      u1[(i-1)*numx+1:i*numx,.] = utilda1 ;
      e0[(i-1)*numx+1:i*numx,.] = etilda0 ;
      e1[(i-1)*numx+1:i*numx,.] = etilda1 ;
      clear utilda0, utilda1, etilda0, etilda1 ;            
      i=i+1 ;
    endo ;
      
    // ------------------------------------------
    //  (d) Pseudo Maximum Likelihood Estimation 
    // ------------------------------------------
    {tetaest,varest} = clogit(aobs,uobs0~uobs1,eobs0~eobs1,namesb) ;                  
    best[.,iter] = tetaest ;
    varb[.,(iter-1)*kparam+1:iter*kparam] = varest ;    
                   
    // ---------------------------- 
    //  (e) Updating probabilities 
    // ---------------------------- 
    i=1 ;
    do while i<=nplayer ;
      buff = (u1[(i-1)*numx+1:i*numx,.]-u0[(i-1)*numx+1:i*numx,.]) ;
      buff = buff*tetaest 
           + (e1[(i-1)*numx+1:i*numx,.]-e0[(i-1)*numx+1:i*numx,.]) ;
      pchoice[.,i] = exp(buff)./(1+exp(buff)) ;
      i=i+1 ;
    endo ;
    
    iter=iter+1 ;
  endo ;
    
  retp(best,varb) ;
endp ;

// --------------------------------------------------------------------------
//  7.  CHECKING FOR THE CONSISTENCY OF THE ESTIMATORS 
//      (ESTIMATION WITH A VERY LARGE SAMPLE)
// --------------------------------------------------------------------------

//  To check for consistency of the estimators (or for possible programming errors)
//  I estimated the model using each of the estimators using a large sample
//  of 400,000 markets. In all the experiments, and for each considered estimation
//  method, the estimates are equal to the true value up to the 4th decimal digit.
//
//  In this version of the program code, I have omitted this part to save memory
//  requirements and CPU time.
//
//  The user interested in checking for consistency can do it by using this program 
//  code with the following selections in PART 1:
//          nobs = 400000 ;    // Number of markets (observations)
//          nrepli = 1 ;       // Number of Monte carlo replications



// --------------------------------------------------------------------------
//  8.  MONTE CARLO EXPERIMENT
// --------------------------------------------------------------------------

"" ;
"*****************************************************************************************" ;
"       MONTE CARLO EXPERIMENT #";; selexper ;
"*****************************************************************************************" ;
"" ;
kparam = rows(trueparam)-2 ;    //  Number of parameters to estimate
npliter = 20 ;                  //  Maximum number of NPL iterations

//  Matrix that stores NPL fixed point estimates (for each replication and NPL iteration)
//  when we initialize the NPL algorithm with nonparametric frequency estimates of CCPs
bmatNP = zeros(kparam,nrepli*npliter) ; 

//  Matrix that stores NPL fixed point estimates (for each replication and NPL iteration)
//  when we initialize the NPL algorithm with logit estimates of CCPs
bmatSP = zeros(kparam,nrepli*npliter) ;

//  Matrix that stores NPL fixed point estimates (for each replication and NPL iteration)
//  when we initialize the NPL algorithm with U(0,1) random draws for the CCPs
bmatR = zeros(kparam,nrepli*npliter) ;

//  Matrix that stores NPL fixed point estimates (for each replication and NPL iteration)
//  when we initialize the NPL algorithm with the true CCPs
btruemat = zeros(kparam,nrepli) ;

redraws=0;  // We set the counter 'redraws' to zero. 
            // Note: When there are multicollinearity problems in a Monte Carlo sample
            // we ignore that sample and take a new one. We want to check for the 
            // number of times we have to make these re-draws and this is why we use
            // the counter 'redraws'

draw=1 ;    // Initializing the counter for the number of Monte Carlo replications
do while (draw<=nrepli) ;
  "     Replication =";; draw ;
  "         (a)     Simulations of x's and a's" ;
  flag=1 ; 
  do while (flag==1) ; 
    {aobs,aobs_1,sobs} = simdygam(nobs,pequil,psteady,vstate) ; 
    check = sumc(sumc(aobs~aobs_1).==zeros(2*nplayer,1)) ; 
    check = check + sumc(sumc(aobs~aobs_1).==(nobs.*ones(2*nplayer,1))) ; 
    if (check>0) ; 
      flag = 1 ; 
    elseif (check==0) ; 
      flag = 0 ; 
    endif ; 
    redraws=redraws+flag;   // Counts the number re-drawings
  endo ; 

  "         (b.1)   Estimation of initial CCPs (Non-Parametric)" ;
  prob0NP = freqprob(aobs,sobs~aobs_1,vstate) ; 

  "         (b.2)   NPL algorithm using frequency estimates as initial CCPs " ;  
  {best,varb} = npldygam(aobs,sobs,aobs_1,sval,ptrans,prob0NP,disfact,npliter);
  bmatNP[.,(draw-1)*npliter+1:draw*npliter] = best ;

  "         (c.1)   Estimation of initial CCPs (Semi-Parametric: Logit)" ;
  // Construct dependent (aobsSP) and explanatory variables (xobsSP)
    aobsSP = reshape(aobs,nobs*nplayer,1);
    alphai = ones(nobs,1).*.eye(nplayer); 
    xobsSP = sobs.*.ones(nplayer,1);
    nfirms_1 = sumc(aobs_1').*.ones(nplayer,1);
    aobs_1SP = reshape(aobs_1,nobs*nplayer,1);
    namesb = (0 $+ "alpha" $+ ftocv(seqa(1,1,nplayer),1,0))
           | "betaSP" | "thetaESP" | "deltaSP" ; 
    xobsSP = alphai~xobsSP~aobs_1SP~nfirms_1 ;
   // Logit estimation
    {best,varest} = milogit(aobsSP,xobsSP,namesb) ; 
   // Construct probabilities
        vstateSP = ones(rows(vstate),nplayer) ~ vstate 
                 ~ (sumc(vstate[.,2:nplayer+1]'));
        best=best[1:nplayer].*eye(nplayer)|ones(1,5).*.best[nplayer+1]
            |eye(nplayer).*best[nplayer+2]|ones(1,5).*.best[nplayer+3] ;
        prob0SP = 1 ./ (1+exp(-vstateSP*best)) ;

  "         (c.2)   NPL algorithm using Logit estimates as initial CCPs " ;       
  {best,varb} = npldygam(aobs,sobs,aobs_1,sval,ptrans,prob0SP,disfact,npliter);
  bmatSP[.,(draw-1)*npliter+1:draw*npliter] = best ;

  "         (d.1)   Estimation of initial CCPs (Completely Random)" ;
  prob0R = rndu(rows(vstate),nplayer) ; 

  "         (d.2)   NPL algorithm using U(0,1) random draws as initial CCPs " ;       
  {best,varb} = npldygam(aobs,sobs,aobs_1,sval,ptrans,prob0R,disfact,npliter);
  bmatR[.,(draw-1)*npliter+1:draw*npliter] = best ;

  "         (e)     NPL algorithm using true values as initial CCPs " ;       
  {best,varb} = npldygam(aobs,sobs,aobs_1,sval,ptrans,pequil,disfact,1);   
  btruemat[.,draw] = best ;

  draw=draw+1 ;

endo ;

"           Number of Re-drawings due to Multicollinearity =";;  redraws;

// --------------------------------------------------------------------------
//  9.  SAVING RESULTS OF THE MONTE CARLO EXPERIMENT
// --------------------------------------------------------------------------

save ^nfile_bNP = bmatNP ;
save ^nfile_bSP = bmatSP ;
save ^nfile_bR  = bmatR ;
save ^nfile_btrue = btruemat ;

// --------------------------------------------------------------------------
//  10. STATISTICS FROM THE MONTE CARLO EXPERIMENTS
//          TABLES 4 AND 5 IN AGUIRREGABIRIA AND MIRA (2007)
// --------------------------------------------------------------------------

// TABLE 4: EMPIRICAL MEANS AND EMPIRICAL STANDARD ERRORS OF THE ESTIMATES

bmatNP = reshape(bmatNP',nrepli,kparam*npliter) ;
bmatSP = reshape(bmatSP',nrepli,kparam*npliter) ;
bmatR  = reshape(bmatR' ,nrepli,kparam*npliter) ;
btruemat = btruemat' ;

// Empirical Means
mean_bmatNP = meanc(bmatNP) ;
mean_bmatNP = reshape(mean_bmatNP,npliter,kparam)' ;
mean_bmatSP = meanc(bmatSP) ;
mean_bmatSP = reshape(mean_bmatSP,npliter,kparam)' ;
mean_bmatR  = meanc(bmatR) ;
mean_bmatR = reshape(mean_bmatR,npliter,kparam)' ;
mean_bmattrue = meanc(btruemat) ;

// Empirical Medians
median_bmatNP = median(bmatNP) ;
median_bmatNP = reshape(median_bmatNP,npliter,kparam)' ;
median_bmatSP = median(bmatSP) ;
median_bmatSP = reshape(median_bmatSP,npliter,kparam)' ;
median_bmatR  = median(bmatR) ;
median_bmatR = reshape(median_bmatR,npliter,kparam)' ;
median_bmattrue = median(btruemat) ;

// Empirical Standard errors
se_bmatNP = stdc(bmatNP) ;
se_bmatNP = reshape(se_bmatNP,npliter,kparam)' ;
se_bmatSP = stdc(bmatSP) ;
se_bmatSP = reshape(se_bmatSP,npliter,kparam)' ;
se_bmatR  = stdc(bmatR) ;
se_bmatR = reshape(se_bmatR,npliter,kparam)' ;
se_bmattrue = stdc(btruemat) ;


"" ;
"*****************************************************************************************" ;
"   MONTE CARLO EXPERIMENT #";; selexper ;
"   EMPIRICAL MEANS AND STANDARD ERRORS"; 
"" ;
"   TABLE 4 OF THE PAPER AGUIRREGABIRIA AND MIRA (2007)" ;
"*****************************************************************************************" ;
"" ;
"----------------------------------------------------------------------------------------" ;
"                       theta_fc_1          theta_rs        theta_rn        theta_ec" ;
"----------------------------------------------------------------------------------------" ;
"TRUE VALUES       ";; trueparam[1];;  trueparam[6];;  trueparam[7];;  trueparam[8]; 
"----------------------------------------------------------------------------------------" ;
"" ;
"   MEAN 2step-True";;  mean_bmattrue[1];;  mean_bmattrue[6];;  mean_bmattrue[7];;  mean_bmattrue[8];
"" ;
" MEDIAN 2step-True";;  median_bmattrue[1];;median_bmattrue[6];;median_bmattrue[7];;median_bmattrue[8];
"" ;
"   S.E. 2step-True";;  se_bmattrue[1];;    se_bmattrue[6];;    se_bmattrue[7];;    se_bmattrue[8];
"" ;
"----------------------------------------------------------------------------------------" ;
"" ;
"   MEAN 2step-Freq";; mean_bmatNP[1,1];;   mean_bmatNP[6,1];;  mean_bmatNP[7,1];;  mean_bmatNP[8,1];
"" ;
" MEDIAN 2step-Freq";; median_bmatNP[1,1];; median_bmatNP[6,1];;median_bmatNP[7,1];;median_bmatNP[8,1];
"" ;
"   S.E. 2step-Freq";; se_bmatNP[1,1];;     se_bmatNP[6,1];;    se_bmatNP[7,1];;    se_bmatNP[8,1];
"" ;
"----------------------------------------------------------------------------------------" ;
"" ;
"   MEAN NPL-Freq  ";;  mean_bmatNP[1,npliter];;    mean_bmatNP[6,npliter];;    mean_bmatNP[7,npliter];;    mean_bmatNP[8,npliter];
"" ;
" MEDIAN NPL-Freq  ";;  median_bmatNP[1,npliter];;  median_bmatNP[6,npliter];;  median_bmatNP[7,npliter];;  median_bmatNP[8,npliter];
"" ;
"   S.E. NPL-Freq  ";;  se_bmatNP[1,npliter];;      se_bmatNP[6,npliter];;      se_bmatNP[7,npliter];;      se_bmatNP[8,npliter];
"" ;
"----------------------------------------------------------------------------------------" ;
"" ;
"  MEAN 2step-Logit";; mean_bmatSP[1,1];;  mean_bmatSP[6,1];;  mean_bmatSP[7,1];;  mean_bmatSP[8,1];
"" ;
"MEDIAN 2step-Logit";; median_bmatSP[1,1];;  median_bmatSP[6,1];;  median_bmatSP[7,1];;  median_bmatSP[8,1];
"" ;
"  S.E. 2step-Logit";; se_bmatSP[1,1];;    se_bmatSP[6,1];;    se_bmatSP[7,1];;    se_bmatSP[8,1];
"" ;
"----------------------------------------------------------------------------------------" ;
"" ;
"    MEAN NPL-Logit";; mean_bmatSP[1,npliter];;  mean_bmatSP[6,npliter];;  mean_bmatSP[7,npliter];;  mean_bmatSP[8,npliter];
"" ;
"  MEDIAN NPL-Logit";; median_bmatSP[1,npliter];;  median_bmatSP[6,npliter];;  median_bmatSP[7,npliter];;  median_bmatSP[8,npliter];
"" ;
"    S.E. NPL-Logit";; se_bmatSP[1,npliter];;    se_bmatSP[6,npliter];;    se_bmatSP[7,npliter];;    se_bmatSP[8,npliter];
"" ;
"----------------------------------------------------------------------------------------" ;
"" ;
" MEAN 2step-Random";; mean_bmatSP[1,1];;  mean_bmatSP[6,1];;  mean_bmatSP[7,1];;  mean_bmatSP[8,1];
"" ;
"MEDIAN 2step-Rando";; median_bmatSP[1,1];;  median_bmatSP[6,1];;  median_bmatSP[7,1];;  median_bmatSP[8,1];
"" ;
" S.E. 2step-Random";; se_bmatSP[1,1];;    se_bmatSP[6,1];;    se_bmatSP[7,1];;    se_bmatSP[8,1];
"" ;
"----------------------------------------------------------------------------------------" ;
"" ;
"   MEAN NPL-Random";; mean_bmatR[1,npliter];;  mean_bmatR[6,npliter];;  mean_bmatR[7,npliter];;  mean_bmatR[8,npliter];
"" ;
" MEDIAN NPL-Random";; median_bmatR[1,npliter];;  median_bmatR[6,npliter];;  median_bmatR[7,npliter];;  median_bmatR[8,npliter];
"" ;
"   S.E. NPL-Random";; se_bmatR[1,npliter];;    se_bmatR[6,npliter];;    se_bmatR[7,npliter];;    se_bmatR[8,npliter];
"" ;
"----------------------------------------------------------------------------------------" ;


// TABLE 5: SQUARE-ROOT MEAN SQUARE ERRORS OF DIFFERENT ESTIMATORS
//          RATIOS OVER THE SQUARE-ROOT MSE OF THE 2-STEP PML USING THE TRUE CCPs

// Empirical Squere-root MSE
trueparam = trueparam[1:kparam] ;

srmse_true = sqrt((mean_bmattrue-trueparam).^2 + se_bmattrue.^2) ;
srmse_NP   = sqrt((mean_bmatNP-trueparam).^2   + se_bmatNP.^2) ;
srmse_SP   = sqrt((mean_bmatSP-trueparam).^2   + se_bmatSP.^2) ;
srmse_R    = sqrt((mean_bmatR-trueparam).^2    + se_bmatR.^2) ;

// Ratios
srmse_NP   = srmse_NP./srmse_true ;
srmse_SP   = srmse_SP./srmse_true ;
srmse_R    = srmse_R./srmse_true ;

"" ;
"*****************************************************************************************" ;
"   MONTE CARLO EXPERIMENT #";; selexper ;
"   SQUARE-ROOT MEAN SQUARE ERRORS"; 
"   RATIOS OVER THE SQUARE-ROOT MSE OF THE 2-STEP PML USING THE TRUE CCPs" ;
"" ;
"   TABLE 5 OF THE PAPER AGUIRREGABIRIA AND MIRA (2007)" ;
"*****************************************************************************************" ;
"" ;
"----------------------------------------------------------------------------------------" ;
"                       theta_fc_1          theta_rs        theta_rn        theta_ec" ;
"----------------------------------------------------------------------------------------" ;
"SQ-MSE 2-step-TRUE";; srmse_true[1];;  srmse_true[6];;  srmse_true[7];;  srmse_true[8]; 
"----------------------------------------------------------------------------------------" ;
"" ;
" RATIO: 2step-Freq";; srmse_NP[1,1];;   srmse_NP[6,1];;  srmse_NP[7,1];;  srmse_NP[8,1];
"" ;
"----------------------------------------------------------------------------------------" ;
"" ;
"   RATIO: NPL-Freq";; srmse_NP[1,npliter];;   srmse_NP[6,npliter];;  srmse_NP[7,npliter];;  srmse_NP[8,npliter];
"" ;
"----------------------------------------------------------------------------------------" ;
"" ;
"RATIO: 2step-Logit";; srmse_SP[1,1];;   srmse_SP[6,1];;  srmse_SP[7,1];;  srmse_SP[8,1];
"" ;
"----------------------------------------------------------------------------------------" ;
"" ;
"  RATIO: NPL-Logit";; srmse_SP[1,npliter];;   srmse_SP[6,npliter];;  srmse_SP[7,npliter];;  srmse_SP[8,npliter];
"" ;
"----------------------------------------------------------------------------------------" ;
"" ;
"RATIO: 2step-Rando";; srmse_R[1,1];;   srmse_R[6,1];;  srmse_R[7,1];;  srmse_R[8,1];
"" ;
"----------------------------------------------------------------------------------------" ;
"" ;
" RATIO: NPL-Random";; srmse_R[1,npliter];;   srmse_R[6,npliter];;  srmse_R[7,npliter];;  srmse_R[8,npliter];
"" ;
"----------------------------------------------------------------------------------------" ;

output off ;
end ;

