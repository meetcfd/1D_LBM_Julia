using Plots
# LBM - 1D, diffusion equation D1Q2

m = 100 #lattices
n = 5 # time_steps

x = collect(range(0,m)) # spatial domain

ϕ = zeros(m+1) 

α = 0.25 # thermal diffusivity
ω = 1. / (α + .5)

ϕ_wall = 1. # boundary condition at high temperature wall

f1 = 0.5 .* ϕ
f2 = 0.5 .* ϕ

function collision!(non_dim_dep_var, relax_fac, f1, f2)
    feq = 0.5 .* non_dim_dep_var
    f1 .= (1. - relax_fac) .* f1 .+ relax_fac .* feq
    f2 .= (1. - relax_fac) .* f2 .+ relax_fac .* feq 
end

function streaming!(f1, f2)
   f1[begin + 1:end] .= f1[begin:end - 1]
   f2[begin:end - 1] .= f2[begin + 1:end]
end

function boundary_conditions!(wall_temperature, f1, f2)
    f1[begin] = wall_temperature - f2[begin]
    f2[end] = f2[end - 1]
end

function time_step!(f1, f2, ϕ, ω, ϕ_wall)
    collision!(ϕ, ω, f1, f2)
    streaming!(f1, f2)
    boundary_conditions!(ϕ_wall, f1, f2)
    ϕ .= f1 + f2
end

P = plot(layout=(2,1));

for t in 1:n
    time_step!(f1, f2, ϕ, ω, ϕ_wall)
    if mod(t,50) == 0
        flux = 2 * ω .* (f1 - f2)
        plot!(P, x, [ϕ,flux], layout=(2,1), label=t)
    end
end

display(P)
