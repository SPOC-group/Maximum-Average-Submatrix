This is the code used to produce the results of the paper:

# Generating data

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

A good set of starting inverse temperatures is 

    [m, beta]
    [0.1, 30.]
    [0.2, 22.]
    [0.3, 18.]
    [0.4, 15.]
    [0.5, 15.]
    [0.6, 15.]
    [0.7, 20.] 
    [0.8, 22.]
    [0.9, 35.]

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

These scalings allows to see the rich phase diagram for 1. < intensivebeta < 4.


Then, after finding a fist solution at intensive beta = 2.5, use the same stepping procedure as that described in the previous section.
It is particularly important to simulate the first SP equations at any $m$ at intensivebeta ~ 2.5
and then descend in intensive inverse temperature very slowly into the region of positive complexity. 
Indeed, the SP equations are very delicate in that region, and tend to default quickly to the RS solution.


## Saving data

The data was then appended into a simple dataframe generated with 

    DataFrame( type=String[], beta=Float64[], h=Float64[], m=Float64[], q0=Float64[], q1=Float64[], p1=Float64[], err=Float64[])

in which err = max(rerr, xerr), and saved as a CSV file using

    CSV.write(database_file, database)

We do not provide a unified script to perform this operation as the data was generated in multiple different occasions