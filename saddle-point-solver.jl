using FastGaussQuadrature, LinearAlgebra, Roots
using LogExpFunctions: log1pexp, logistic, logsumexp

# Gaussian pdf and cdf
gaussPDF(z; μ=0., σ2=1.) = exp(- (z - μ)^2 / (2 * σ2)) / sqrt(2 * pi * σ2)
gaussCDF(z; μ=0., σ2=1.) = 0.5 * (1. + erf( (z - μ) / sqrt(2. * σ2) ))

# Regularised sqrt
mysqrt(x) = sqrt(abs(x))

# precomputed Gauss-Hermite quadrature weights and points
const xGauss = gausshermite(201)[1]
const wGauss = gausshermite(201)[2]

# given f, this computes log int gauss(x) exp(f(x)) dx in a stable way using quadrature + logsumexp
function logIntGauss(f)
    # I have to compute (modulo conversion factor 
    # due to quadrature being defined for exp(-x^2))
    
    # log sum(wGauss * exp(f(xGauss)))
    # = log sum( exp( f(xGauss) + log(wGauss) ) )
    # = logsumexp( f(xGauss) + log(wGauss) )

    integrand(x, w) = f(sqrt(2) * x) + log(w)
    zs = integrand.(xGauss, wGauss)

    # now I compute log(sum(exp(zs)))
    logsumexp(zs) - 0.5*log(pi)
end

# given f, this computes int gauss(x) f(x) dx
function intGauss(f)
    dot(wGauss, f.(sqrt(2) .* xGauss)) / sqrt(pi)
end

################# untested RS 

# RS saddle point equations

function RSeqM(beta, h, m, q)
    betaH(u) = beta * h + beta^2/2 * (m-q) + beta * u * mysqrt(q)

    # use log(logistic(x)) = x - log1pexp(x)
    # and log((1+exp(x))^y) = y * log1pexp(x)
    # and A / B = exp(log(A) - log(B))
    residual = exp(logIntGauss( u -> betaH(u) - log1pexp(betaH(u))))

    residual - m
end

function RSeqQ(beta, h, m, q)
    betaH(u) = beta * h + beta^2/2 * (m-q) + beta * u * mysqrt(q)

    # use log(logistic(x)) = x - log1pexp(x)
    # and log((1+exp(x))^y) = y * log1pexp(x)
    # and A / B = exp(log(A) - log(B))
    residual = exp(logIntGauss( u -> 2 * betaH(u) - 2 * log1pexp(betaH(u))))

    residual - q
end

# scalar solvers for h and p1 at fixed values of other parameters

function RSsolveForH(beta, h_init, m, q)
    find_zero(h -> RSeqM(beta, h, m, q), h_init, Order1())
end

# RS_SP_fixM iterator
function RS_SP_fixM(beta, m_target;
        init_q = (m_target^2 + m_target)/2,
        init_h = 0.,
        damping = 0.5, 
        maxsteps = 50, 
        tol = 1e-5, 
        updatestep_h = 1, damping_h = 0.1
    )

    current_h = init_h
    current_q = init_q

    eqq = RSeqQ(beta, current_h, m_target, current_q)

    rerr = 200.
    xerr = 200.

    for i in 1:maxsteps

        # update q0 and q1
        # recall that eqQ0 returns integral - q0, so that the update of the fp equation is integral = q0 + eqQ0
        new_q = current_q + (1-damping) * eqq

        # update h
        new_h = if i % updatestep_h == 0
            try
                damping_h * current_h + (1 - damping_h) * RSsolveForH(beta, current_h, m_target, new_q)
            catch
                current_h + 0.1 * (rand()-0.5)
            end
        else
            current_h
        end

        

        # check for tolerance 
        eqm  = RSeqM(beta, new_h, m_target, new_q)
        eqq = RSeqQ(beta, new_h, m_target, new_q)

        rerr = norm([eqm, eqq])
        xerr = norm([ new_h - current_h, new_q - current_q])

        current_h = new_h
        current_q = new_q

        if rerr < tol && xerr < tol
            return (beta, current_h, m_target, current_q, [rerr, xerr], 1., i)
        end
    end

    return (beta, current_h, m_target, current_q, [rerr, xerr], 0., maxsteps)
end

# RS stability 
function RS_replicon(beta, h, m, q)
    betaH(u) = beta * h + beta^2/2 * (m-q) + beta * u * mysqrt(q)
    1/beta^2 - exp(logIntGauss( u -> 2 * betaH(u) - 4 * log1pexp(betaH(u))))
end

################# 1RSB 

# 1RSB saddle point equations

function RSB1eqM(beta, h, m, q0, q1, p1)
    betaH(u, v) = beta * h + beta^2/2 * (m-q1) + beta * v * mysqrt(q1-q0) + beta * u * mysqrt(q0)

    # use log(logistic(x)) = x - log1pexp(x)
    # and log((1+exp(x))^y) = y * log1pexp(x)
    # and A / B = exp(log(A) - log(B))
    logIntegrand(u) = 
        logIntGauss( v -> (p1 - 1) * log1pexp(betaH(u, v)) + betaH(u, v) ) - 
        logIntGauss( v -> p1 * log1pexp(betaH(u, v)) )

    residual = exp(logIntGauss( u -> logIntegrand(u) ))
    residual - m
end

function RSB1eqQ0(beta, h, m, q0, q1, p1)
    betaH(u, v) = beta * h + beta^2/2 * (m-q1) + beta * v * mysqrt(q1-q0) + beta * u * mysqrt(q0)

    # use log(logistic(x)) = x - log1pexp(x)
    # and log((1+exp(x))^y) = y * log1pexp(x)
    # and A / B = exp(log(A) - log(B))
    logIntegrand(u) = 
        2 * logIntGauss( v -> (p1 - 1) * log1pexp(betaH(u, v)) + betaH(u, v) ) - 
        2 * logIntGauss( v -> p1 * log1pexp(betaH(u, v)) )

    residual = exp(logIntGauss( u -> logIntegrand(u) ))
    residual - q0
end

function RSB1eqQ1(beta, h, m, q0, q1, p1)
    betaH(u, v) = beta * h + beta^2/2 * (m-q1) + beta * v * mysqrt(q1-q0) + beta * u * mysqrt(q0)

    # use log(logistic(x)) = x - log1pexp(x)
    # and log((1+exp(x))^y) = y * log1pexp(x)
    # and A / B = exp(log(A) - log(B))
    logIntegrand(u) = 
        logIntGauss( v -> (p1 - 2) * log1pexp(betaH(u, v)) + 2 * betaH(u, v) ) - 
        logIntGauss( v -> p1 * log1pexp(betaH(u, v)) )

    residual = exp(logIntGauss( u -> logIntegrand(u) ))
    residual - q1
end

function RSB1eqP1(beta, h, m, q0, q1, p1)
    betaH(u, v) = beta * h + beta^2/2 * (m-q1) + beta * v * mysqrt(q1-q0) + beta * u * mysqrt(q0)

    # use log((1+exp(x))^y) = y * log1pexp(x)
    # and A / B = exp(log(A) - log(B))
    logIntegrand1(u) = 
        logIntGauss( v -> p1 * log1pexp(betaH(u, v)) + log(log1pexp(betaH(u, v))) ) - 
        logIntGauss( v -> p1 * log1pexp(betaH(u, v)) )

    integral1 = exp(logIntGauss( u -> logIntegrand1(u) ))

    logIntegrand2(u) = logIntGauss( v -> p1 * log1pexp(betaH(u, v)) )

    integral2 = intGauss( u -> logIntegrand2(u) )

    p1^2 * beta^2 / 4 * (q0^2 - q1^2) + p1 * integral1 - integral2
end

# scalar solvers for h and p1 at fixed values of other parameters

function RSB1solveForH(beta, h_init, m, q0, q1, p1)
    find_zero(h -> RSB1eqM(beta, h, m, q0, q1, p1), h_init, Order1())
end

function RSB1solveForP1(beta, h, m, q0, q1, p1_init)
    find_zero(p1 -> RSB1eqP1(beta, h, m, q0, q1, p1), p1_init, Order1())
end

# SP iteration

function RSB1_SP_fixM(beta, m_target;
        init_q0 = 0.2,
        init_q1 = 0.5,
        init_p1 = 0.5,
        init_h = 0.,
        damping = 0.5, 
        maxsteps = 1000, 
        tol = 1e-6, 
        updatestep_p1 = 5, damping_p1 = 0.25, fix_p1 = false,
        updatestep_h = 1, damping_h = 0.1
    )

    current_h = init_h
    current_q0 = init_q0
    current_q1 = init_q1
    current_p1 = init_p1

    # new_h = 0.
    # new_q0 = 0.
    # new_q1 = 0.
    # new_p1 = 0.

    eqq0 = RSB1eqQ0(beta, current_h, m_target, current_q0, current_q1, current_p1)
    eqq1 = RSB1eqQ1(beta, current_h, m_target, current_q0, current_q1, current_p1)

    rerr = 200.
    xerr = 200.

    for i in 1:maxsteps

        # update q0 and q1
        # recall that eqQ0 returns integral - q0, so that the update of the fp equation is integral = q0 + eqQ0
        new_q0 = current_q0 + (1-damping) * eqq0
        new_q1 = current_q1 + (1-damping) * eqq1

        # update h
        new_h = if i % updatestep_h == 0
            try
                damping_h * current_h + (1 - damping_h) * RSB1solveForH(beta, current_h, m_target, new_q0, new_q1, current_p1)
            catch
                current_h + 0.1 * (rand()-0.5)
            end
        else
            current_h
        end

        # update p1
        new_p1 = if fix_p1
            current_p1
        elseif i % updatestep_p1 == 0
            try
                damping_p1 * current_p1 + (1 - damping_p1) * RSB1solveForP1(beta, new_h, m_target, new_q0, new_q1, current_p1)
            catch
                0.5#clamp(current_p1 + 0.1 * (rand()-0.5), 0.1, 0.9)
            end
        else
            current_p1
        end
          
 

        # check for tolerance 
        eqm  = RSB1eqM(beta, new_h, m_target, new_q0, new_q1, new_p1)
        eqq0 = RSB1eqQ0(beta, new_h, m_target, new_q0, new_q1, new_p1)
        eqq1 = RSB1eqQ1(beta, new_h, m_target, new_q0, new_q1, new_p1)
        eqp1 = fix_p1 ? 0. : RSB1eqP1(beta, new_h, m_target, new_q0, new_q1, new_p1)

        rerr = norm([eqm, eqq0, eqq1, eqp1])
        xerr = norm([ new_h - current_h, new_q0 - current_q0, new_q1 - current_q1, new_p1 - current_p1 ])

        current_h = new_h
        current_q0 = new_q0
        current_q1 = new_q1
        current_p1 = new_p1

        if rerr < tol && xerr < tol
            return (beta, current_h, m_target, current_q0, current_q1, current_p1, [rerr, xerr], 1., i)
        end
    end

    return (beta, current_h, m_target, current_q0, current_q1, current_p1, [rerr, xerr], 0., maxsteps)
end

function RSB1_residual(beta, h, m, q0, q1, p1)
    eqm  = RSB1eqM(beta, h, m, q0, q1, p1)
    eqq0 = RSB1eqQ0(beta, h, m, q0, q1, p1)
    eqq1 = RSB1eqQ1(beta, h, m, q0, q1, p1)
    eqp1 = RSB1eqP1(beta, h, m, q0, q1, p1)

    [eqm, eqq0, eqq1, eqp1, norm([eqm, eqq0, eqq1, eqp1])]
end

# observables 

function RSB1_free_entropy(beta, h, m, q0, q1, p1)
    betaH(u, v) = beta * h + beta^2/2 * (m-q1) + beta * v * mysqrt(q1-q0) + beta * u * mysqrt(q0)

    logIntegrand(u) = 
        logIntGauss( v -> p1 * log1pexp(betaH(u, v)) )

    integral = intGauss( u -> logIntegrand(u) )

    - beta^2 / 4 * (m^2 - p1 * q0^2 + (p1-1) * q1^2) + integral / p1
end

function RSB1_energy(beta, h, m, q0, q1, p1)
    - beta / 2 * (m^2 - p1 * q0^2 + (p1-1) * q1^2)
end

function RSB1_maxavg(beta, h, m, q0, q1, p1)
    beta / m^2 * (m^2 - p1 * q0^2 + (p1-1) * q1^2)
end

function RSB1_entropy(beta, h, m, q0, q1, p1)
    energy = RSB1_energy(beta, h, m, q0, q1, p1)
    free_entropy = RSB1_free_entropy(beta, h, m, q0, q1, p1)
    free_entropy + beta * energy - beta * h * m
end

function RSB1_complexity(beta, h, m, q0, q1, p1)
    - RSB1eqP1(beta, h, m, q0, q1, p1)
end

# stability 

function RSB1_replicon(beta, h, m, q0, q1, p1)
    betaH(u, v) = beta * h + beta^2/2 * (m-q1) + beta * v * mysqrt(q1-q0) + beta * u * mysqrt(q0)

    # use log(logistic(x)) = x - log1pexp(x)
    # and log(1-logistic(x)) =  - log1pexp(x)
    # and log((1+exp(x))^y) = y * log1pexp(x)
    # and A / B = exp(log(A) - log(B))
    logIntegrand(u) = 
        logIntGauss( v -> (p1 - 4) * log1pexp(betaH(u, v)) + 2 * betaH(u, v) ) - 
        logIntGauss( v -> p1 * log1pexp(betaH(u, v)) )

    integral = exp(logIntGauss( u -> logIntegrand(u) ))
    1/beta^2 - integral
end

####### 