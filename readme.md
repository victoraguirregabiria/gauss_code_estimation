MONTE CARLO EXPERIMENTS IN AGUIRREGABIRIA & MIRA (ECONOMETRICA, 2007)

 am_econometrica_2007_montecarlo.prg      Program

am_econometrica_2007_montecarlo_1.pdf  Results 1

am_econometrica_2007_montecarlo_2.pdf  Results 2

am_econometrica_2007_montecarlo_3.pdf  Results 3

am_econometrica_2007_montecarlo_4.pdf  Results 4

am_econometrica_2007_montecarlo_5.pdf  Results 5

am_econometrica_2007_montecarlo_6.pdf  Results 6

---------------------------------------------------------------

ESTIMATION OF (STANDARD) DISCRETE CHOICE MODELS

milogit.src             Procedure for the Maximum Likelihood estimation of a Logit model.

milogit example    Program that runs an example calling the procedure milogit.src.

miprobit.src           Procedure for the Maximum Likelihood estimation of a Probit model.

miprobit example Program that runs an example calling the procedure miprobit.src.

multilog.src           Procedure for the Maximum Likelihood estimation of a Multinomial Logit.

multilog example Program that runs an example calling the procedure multilog.src.

clogit.src               Procedure for the Maximum Likelihood estimation of McFadden's Conditional Logit.

clogit example       Program that runs an example calling the procedure clogit.src.

---------------------------------------------------------------
ESTIMATION OF SINGLE-AGENT DISCRETE-CHOICE DYNAMIC-PROGRAMMING MODELS

npl_sing.src           Procedure that estimates the structural parameters of a discrete-choice single-agent dynamic programming model using the Nested Pseudo Likelihood (NPL) algorithm in Aguirregabiria and Mira (Econometrica, 2002). It calls the following procedures:

·        clogit.src          Procedure for the Maximum Likelihood estimation of McFadden's Conditional Logit.

npl_sing.e             Program that estimates the bus replacement model in Rust (Econometrica, 1987). This program calls the library nplprocs.lcg and the procedures npl_sing.src, multilog.src , discthre.src and the GAUSS dataset:

·        bus1234.dat      Rust’s bus replacement data set (bus engine groups 1, 2, 3 and 4).

·        bus1234.dht      

---------------------------------------------------------------
ESTIMATION OF STATIC DISCRETE GAMES

nplgame.src
Procedure that estimates the structural parameters of an entry model of incomplete information using a Nested Pseudo-Likelihood (NPL) algorithm.

equiprob.src
Procedure that computes players' choice probabilities of Bayesian Nash Equilibrium in a static game of firms' market entry.

probmat.src
Procedure for the mapping from probabilities to probabilities in a model of entry with incomplete information.

simgame.src
Procedure that simulates observations from a model of entry.

dynprob.src
Procedure for the mapping from probabilities to probabilities in a static model of entry with incomplete information.

kpiegame.src
Procedure that estimates the structural parameters of an entry model of incomplete information using 1-stage of the Nested Pseudo-Likelihood (NPL) algorithm.

psumbern.src
Procedure that obtains the probability distribution of a sum of independent but heterogeneous Bernoullis.

---------------------------------------------------------------
ESTIMATION OF DYNAMIC DISCRETE GAMES

mpeprob.src
Procedure that computes players' choice probabilities of Markov Perfect Equilibrium in a dynamic game of firms' market entry/exit with incomplete information.

mlegame.src
Procedure that estimates the structural parameters of a model of entry of incomplete information using Maximum Likelihood.

npldygam.src
Procedure that estimates the structural parameters of dynamic game of firms' entry/exit using a Nested Pseudo-Likelihood (NPL) algorithm.

simdygam.src
Compute the steady-state distribution of state and decision variables in a Markov Perfect Equilibrium of a dynamic game of firms' market entry/exit with incomplete information.  

---------------------------------------------------------------
NONPARAMETRIC METHODS

kernel1.src            Procedure for the kernel estimation of a univariate density function.

nadaraya_cv.src    Procedure for Nadaraya-Watson estimation of a nonparametric regression function. Cross-Validation for bandwidth choice.

nadaraya_cv e       Program that runs an example calling the procedure nadaraya_cv.src.

isonpreg.src          Procedure that obtains a nonparametric isotonic (monotonic) regression using the max-min estimator first proposed by Brunk (AS, 1958).

mono_si.src          Procedure that estimates a nonparametric regression function using the SI (Smoothing- Isotonising) estimator in Mammen (AS, 1991).

mono_is.src          Procedure that estimates a nonparametric regression function using the IS (Isotonising-Smoothing) estimator in Mammen (AS, 1991).

freqprob.src          Procedure that obtains a frequency estimation of Prob(Y|X) where Y is discrete and X is a vector of discrete variables.

dchokern.src         Procedure for the kernel estimation of Prob(Y|X) where Y is a binary variable and X is a vector of continuous variables.

dchokern.e            Program that runs an example calling the procedure dchokern.src.

---------------------------------------------------------------
OTHERS

iniprob.src
Obtains initial reduced form estimates of conditional choice probabilities to be used as starting value for the K-stage PIE.

nplld.src
Procedure that estimates a dynamic labor demand model with lump-sum and linear hiring and firing costs, using a nested fixed point algorithm.

disckpie.src
It discretizes a vector of variables according to a criterion selected by the user.

ccpld.src
Procedure that estimates a dynamic labor demand model with lump-sum and linear hiring and firing costs, using a nested fixed point algorithm.

disckld.src
Discretizes decision and state variables and obtains matrices of transition probabilities in a labor demand model with linear hiring and firing costs.

solvld.src
Procedure that solves a labor demand model with with two types of labor. The algorithm exploits a "smooth" Bellman's equation.

simld.src
Procedure that simulates data of state and decision variables in a labour demand model with linear hiring and firing costs.

dissolld.src
Procedure that discretizes decisions and state variables and obtains matrices of transition probabilities in a dynamic labor demand model with two types of labor.
