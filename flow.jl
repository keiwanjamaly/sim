using Logging: global_logger
using TerminalLoggers: TerminalLogger
global_logger(TerminalLogger())
using DifferentialEquations
using Sundials
using BenchmarkTools
using LinearAlgebra

mutable struct FlowConfig
    sigmaMax::Float64
    grid::Array{Float64}
    dx::Float64
    Lambda::Float64
    mu::Float64
    T::Float64
    beta::Float64
    NFlavor::Float64
    function FlowConfig(sigmaMax::Float64, nGrid::Int64, Lambda::Float64, mu::Float64, T::Float64, NFlavor::Float64)
        self = new()

        self.sigmaMax = 6.0
        self.grid = LinRange(0, sigmaMax, nGrid)
        self.dx = sigmaMax / nGrid
        self.Lambda = Lambda
        self.mu = mu
        self.T = T
        self.beta = 1 / T
        self.NFlavor = NFlavor

        return self
    end
end

function calK(config::FlowConfig, t::Float64)
    config.Lambda * exp(-t)
end

function calT(config::FlowConfig, k::Float64)
    -log(k / config.Lambda)
end

function e_f(k::Float64, sigma::Float64)
    sqrt(k^2 + sigma^2)
end

function e_b(k::Float64, ux::Float64)
    sqrt(k^2 + ux)
end

function n_b(x)
    1 / (exp(x) + 1)
end

function differentiateU(u::Array{Float64}, dx::Float64)
    u_x = Array{Float64}(undef, length(uInitial) + 1)
    u_x[1] = u[2] / dx
    for i = 1:(length(u)-1)
        u_x[i+1] = (u[i+1] - u[i]) / dx
    end
    u_x[end] = u_x[end-1]
    u_x
end

function initialCondition(x::Float64, config::FlowConfig)
    # intermediate = 1 / sqrt(1 + (1 / config.Lambda)^2)
    # (x / pi) * (atanh(intermediate) - intermediate)
    intermediate = (2 + config.Lambda^2 - 2 * sqrt(1 + config.Lambda^2)) / (2 * pi * sqrt(1 + config.Lambda^2))
    x * intermediate
end

function S(k::Float64, sigma::Float64, config::FlowConfig)
    e = e_f(k, sigma)
    beta = config.beta
    minus = beta * (e - config.mu) / 2
    plus = beta * (e + config.mu) / 2

    # sigma * k^3 / (4 * e^3 * pi) * (e * beta * (sech(minus)^2 + sech(plus)^2) - 2 * (tanh(minus) + tanh(plus)))
    sigma * k^4 / (8 * e^3 * pi) * (e * beta * (sech(minus)^2 + sech(plus)^2) - 2 * (tanh(minus) + tanh(plus)))
end

function Q(k::Float64, ux::Float64, config::FlowConfig)
    e = e_b(k, ux)
    beta = config.beta
    N = config.NFlavor

    # -k^3 / (2 * pi * e * N) * (1 + 2 * n_b(beta * e))
    -k^4 / (4 * pi * e * N) * (1 + 2 * n_b(beta * e))
end


function f(du::Array{Float64}, u::Array{Float64}, p::FlowConfig, t::Float64)
    k = calK(config, t)
    for i = 1:length(u)
        du[i] = S(k, p.grid[i], p)
    end

    ux = differentiateU(u, config.dx)
    QArr = Array{Float64}(undef, length(u) + 1)
    for i = 1:length(ux)
        QArr[i] = Q(k, ux[i], config)
    end

    for i = 1:length(u)
        du[i] += (QArr[i+1] - QArr[i]) / config.dx
    end

    println(t)
end


# u0 = 1 / 2
# tspan = (0.0, 1.0)
# prob = ODEProblem(f, u0, tspan)
# print("integrating\n")d
# sol = solve(prob, Tsit5(), reltol=1e-2, abstol=1e-2)

# print("done :)\n")

sigmaMax = 6.0
nGrid = 1000
Lambda = 1e3
kir = 1e-4
mu = 0.1
T = 0.1
NFlavor = 2.0
grid = LinRange(0, sigmaMax, nGrid)

config = FlowConfig(sigmaMax, nGrid, Lambda, mu, T, NFlavor)

uInitial = map(x -> initialCondition(x, config), grid)

tspan = (0.0, calT(config, kir))
vCenter = Array{Float64}(undef, nGrid)
vSide = Array{Float64}(undef, nGrid - 1)
jac_bandity = Tridiagonal(vSide, vCenter, vSide)
f_p = ODEFunction(f; jac_prototype=jac_bandity)
prob = ODEProblem(f_p, uInitial, tspan, config)

# sol = solve(prob, TRBDF2(), reltol=1e-10, abstol=1e-10)

