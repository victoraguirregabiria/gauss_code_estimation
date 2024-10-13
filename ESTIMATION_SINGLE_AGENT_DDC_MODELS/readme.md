ESTIMATION OF SINGLE-AGENT DISCRETE-CHOICE DYNAMIC-PROGRAMMING MODELS
---------------------------------------------------------------------
npl_sing.src           Procedure that estimates the structural parameters of a discrete-choice single-agent dynamic programming model using the Nested Pseudo Likelihood (NPL) algorithm in Aguirregabiria and Mira (Econometrica, 2002). It calls the following procedures:
      clogit.src          Procedure for the Maximum Likelihood estimation of McFadden's Conditional Logit.

npl_sing.e             Program that estimates the bus replacement model in Rust (Econometrica, 1987). This program calls the library nplprocs.lcg and the procedures npl_sing.src, multilog.src , discthre.src and the GAUSS dataset:
·        bus1234.dat      Rust’s bus replacement data set (bus engine groups 1, 2, 3 and 4).
·        bus1234.dht      
