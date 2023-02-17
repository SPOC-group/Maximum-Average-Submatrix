include("saddle-point-solver.jl")

### TODO check RS stability and 1RSB stability. They should match.

# test 1: 0-th to 4-th even moments of a gaussian
begin
    println("--------------")
    println("Answer: [1, 1, 3]")
    logIntGauss.([x->0., x->log(x^2), x->log(x^4)]) .|> exp |> println
    intGauss.([x->1, x->x^2, x->x^4]) |> println
    println("--------------")
end

# test 2: check that 1RSB equations give residual of order 1e-5 on precomputed solutions
begin
    println("--------------")
    println("Answer: all ~ 1e-5")
    beta, h, m, q0, q1, p1 = 17.25, -0.5461060715852927, 0.35, 0.2820352755178014, 0.3274186640689256, 0.34478193694074744
    RSB1eqM(beta, h, m, q0, q1, p1) |> abs |> println
    RSB1eqQ0(beta, h, m, q0, q1, p1) |> abs |> println
    RSB1eqQ1(beta, h, m, q0, q1, p1) |> abs |> println
    RSB1eqP1(beta, h, m, q0, q1, p1) |> abs |> println

    beta, h, m, q0, q1, p1 = 18.691525423728812, -0.28983199171338925, 0.5, 0.4480748731055634, 0.48117688708959494, 0.3682023999488861
    RSB1eqM(beta, h, m, q0, q1, p1) |> abs |> println
    RSB1eqQ0(beta, h, m, q0, q1, p1) |> abs |> println
    RSB1eqQ1(beta, h, m, q0, q1, p1) |> abs |> println
    RSB1eqP1(beta, h, m, q0, q1, p1) |> abs |> println
    println("--------------")
end

# test 3: check  1RSB solvers on a precomputed solution
begin
    println("--------------")
    println("Answer: all ~ 1e-5")
    beta, h, m, q0, q1, p1 = 17.25, -0.5461060715852927, 0.35, 0.2820352755178014, 0.3274186640689256, 0.34478193694074744
    RSB1solveForH(beta, -0.2, m, q0, q1, p1) - h |> abs |> println
    RSB1solveForP1(beta, h, m, q0, q1, 0.7) - p1 |> abs |> println

    beta, h, m, q0, q1, p1 = 18.691525423728812, -0.28983199171338925, 0.5, 0.4480748731055634, 0.48117688708959494, 0.3682023999488861
    RSB1solveForH(beta, -0.5, m, q0, q1, p1) - h |> abs |> println
    RSB1solveForP1(beta, h, m, q0, q1, 0.7) - p1 |> abs |> println
    println("--------------")
end

# test 4: check that  1RSB solver one finds a precomputed solution
begin
    println("--------------")
    println("Answer: [17.25, -0.55, 0.35, 0.28, 0.33, 0.34]")
    beta, h, m, q0, q1, p1 = 17.25, -0.5461060715852927, 0.35, 0.2820352755178014, 0.3274186640689256, 0.34478193694074744

    RSB1_SP_fixM(beta, m;
        init_q0 = 0.25,
        init_q1 = 0.35,
        init_p1 = 0.35,
        init_h = -0.5,
        damping = 0.5, 
        maxsteps = 500, 
        tol = 1e-5, 
        updatestep_p1 = 5, 
        updatestep_h = 1
    ) |> println
    println("--------------")
end

# test 5: check that  1RSB solver one finds a precomputed solution at p1=1 fixed
begin
    println("--------------")
    println("Answer: [23., -0.52, 0.035, 0.0088, 0.028, 1.]")
    # beta, h, m, q0, q1, p1 = 19.6, -0.5089834385201364, 0.035, 0.0029700560065282506, 0.017881130578959736, 1.
    beta, h, m, q0, q1, p1 = 23., -0.5231493708442925, 0.035, 0.008876092238851033, 0.0282599747742614, 1.

    RSB1_SP_fixM(beta, m;
        init_q0 = m^2,
        init_q1 = m,
        init_p1 = p1,
        init_h = -0.1,
        damping = 0.5, 
        maxsteps = 500, 
        tol = 1e-5, 
        updatestep_p1 = 5, fix_p1 = true,
        updatestep_h = 1
    )  |> println
    println("--------------")
end

# test 6: check that RS equations give residual of order 1e-5 on precomputed solutions
begin
    println("--------------")
    println("Answer: all ~ 1e-11")
    beta, h, m, q = 19.11764705882353, -0.5766777743725129, 0.05, 0.017749085445828854
    RSeqM(beta, h, m, q) |> abs |> println
    RSeqQ(beta, h, m, q) |> abs |> println

    # beta, h, m, q0, q1, p1 = 18.691525423728812, -0.28983199171338925, 0.5, 0.4480748731055634, 0.48117688708959494, 0.3682023999488861
    # RSB1eqM(beta, h, m, q0, q1, p1) |> abs |> println
    # RSB1eqQ0(beta, h, m, q0, q1, p1) |> abs |> println
    # RSB1eqQ1(beta, h, m, q0, q1, p1) |> abs |> println
    # RSB1eqP1(beta, h, m, q0, q1, p1) |> abs |> println
    println("--------------")
end

# test 7: check that  RS solver one finds a precomputed solution
begin
    println("--------------")
    println("Answer: [19.12, -0.58, 0.05, 0.018]")
    beta, h, m, q = 19.11764705882353, -0.5766777743725129, 0.05, 0.017749085445828854

    RS_SP_fixM(beta, m;
        init_q = 0.025,
        damping = 0.5, 
        maxsteps = 500, 
        tol = 1e-5, 
        updatestep_h = 1
    ) |> println
    println("--------------")
end