This is the code used to produce the results of the paper: [Statistical mechanics of the maximum-average submatrix problem](https://arxiv.org/abs/2303.05237)

# Requirements

The data and plots were generated using a [Jupyter](https://jupyter.org) notebook and the packages

    LinearAlgebra
    IJulia v1.24.0
    Plots v1.38.4
    FastGaussQuadrature v0.5.0
    Roots v2.0.8
    LogExpFunctions v0.3.21
    CSV v0.10.9
    DataFrames v1.4.4
    LinearRegression v0.2.1

# The dataset

```m_less_0.1.csv``` and ```m_great_0.1.csv``` contain all data used in the analysis and in the generation of the figures.
Each line contains the following fields

    type                # finte M / small M     
    beta                # inverse temperature     
    m                   # magnetization 
    h                   # magnetic field
    q0                  # overlap 0  
    q1                  # overlap 1  
    p1                  # Parisi parameter 
    err                 # max between the xerr and rerr, see below  
    h_p1                # magnetic field at p = 1      
    q0_p1               # overlap 0 at p = 1        
    q1_p1               # overlap 1 at p = 1        
    p1_p1               # Parisi parameter fixed to 1      
    err_p1              # max between the xerr and rerr at p = 1, see below         
    complexity_p1       # complexity at p=1              
    average             # submatrix-average      
    replicon            # type-I 1-RSB stability (stable if positive)         
    entropy             # total entropy      
    average_p1          # submatrix-average at p=1          
    replicon_p1         # type-I 1-RSB stability at p=1 (stable if positive)          
    entropy_p           # total entropy at p=1          

The data can be checked by rerunning the SP solver initialized at each lines values of the parameters.
```data_processed.csv``` contains the data cleaned and processed, with the following fields

    type                # M denotes finite m, S small m, see sections below      
    beta                # inverse temperature     
    m                   # magnetization 
    h                   # magnetic field
    q0                  # overlap 0  
    q1                  # overlap 1  
    p1                  # Parisi parameter 
    err                 # max between the xerr and rerr, see below  
    complexity_p1       # complexity at p=1
    average             # submatrix-average 
    replicon            # type-I 1-RSB stability (stable if positive) 
    entropy             # total entropy 
    computed_m          # cross check that the obtained value of h produces the correct magnetisation

# Generating figures 

The figures are generated from ```data_processed.csv``` in the ```analysis.ipynb``` Jupyter notebook.

# Generating data

The data is generated from in the ```analysis.ipynb``` Jupyter notebook, we provide good starting conditions for the solver for many values of $m$.

## Saddle point equation solver

To solve the 1-RSB saddle point equations, use the function 

    function RSB1_SP_fixM(
        beta,                   # inverse temperature
        m_target                # magnetization
        ;                       # -- from here optional parameteres with their default values
        init_q0 = 0.2,          # initialization for q0
        init_q1 = 0.5,          # initialization for q1
        init_p1 = 0.5,          # initialization for p
        init_h = 0.,            # initialization for h
        damping = 0.5,          # damping for the iteration of the q0 and q1 equations
        maxsteps = 1000,          # maximum number of overall iterations
        tol = 1e-6,             # the iterations stop when the L2 distance between 
                                # updated and old parameters, as well as the residuals of 
                                # all equations, are less than tol
        updatestep_p1 = 5,      # number of iterations in between each solution of the equation for p
        damping_p1 = 0.25,      # damping for the update of p
        fix_p1 = false,         # if false, solve Σ(p) = 0, else solve at fixed p = init_p1   
        updatestep_h = 1,       # number of iterations in between each solution of the the equation for h
        damping_h = 0.1         # damping for the update of h
    )

It returns a tuple with the following informations

    beta, 
    h at convergence, 
    m_target, 
    q0 at convergence, 
    q1 at convergence, 
    p1 at convergence, 
    [norm of residuals, L2 distance between last and second-to-last order parameters], 
    1 if converged before maxsteps, 0 otherwise, 
    steps used to reach convergence

## Finite magnetization

To generate the data for finite $m$, we first find an inverse temperature $\beta_0$ already in the RSB region $q_0 \neq q_1$ and solve the SP equations. This may require a large amount of steps, as well as some fine-tuning of all optional parameters (even thought the defaults should work just fine).
Then, fix $\beta_1 = \beta_0 + \beta_{\rm step}$ for some small step $\beta_{\rm step}$ (for example 0.1), and solve the SP equations at $\beta_1$ using the solution at $\beta_0$ as initialization.
Repeat up (or down) to the desired inverse temperature.

## Small magnetization

To simulate data at small $m$ and $Σ(p) = 0$, use the initialization 

    beta = intensivebeta * sqrt(log(1/m)/m)
    init_h = -1.5 * sqrt(log(1/m)*m)
    init_q0 = 1.001 * m^2
    init_q1 = 0.999 * m
    init_p1 = 2 / intensivebeta

while for $p=1$ use 

    beta = intensivebeta * sqrt(log(1/m)/m)
    init_h = - (intensivebeta^2 + 2)/(2 * intensivebeta) * sqrt(log(1/m)*m)
    init_q0 = 1.001 * m^2
    init_q1 = 0.999 * m
    init_p1 = 1.

Then, after finding a fist solution at intensive beta = 2.1, use the same stepping procedure as that described in the previous section.
It is particularly important to simulate the first SP equations at any $m$ at intensivebeta ~ 2.1
and then descend in intensive inverse temperature very slowly into the region of positive complexity. 
Indeed, the SP equations are very delicate in that region, and tend to default quickly to the RS solution.

These scalings allows to see the rich phase diagram for 1. < intensivebeta < 4.

## Observables 

The following observables can be computed 

    RSB1_free_entropy(beta, h, m, q0, q1, p1)
    RSB1_energy(beta, h, m, q0, q1, p1)
    RSB1_maxavg(beta, h, m, q0, q1, p1)
    RSB1_entropy(beta, h, m, q0, q1, p1)
    RSB1_complexity(beta, h, m, q0, q1, p1)
    RSB1_replicon(beta, h, m, q0, q1, p1)

# RS

The code provides also a RS solver. 
It works very similarly to the 1-RSB solver, so refer to the sections above for the "documentation".
