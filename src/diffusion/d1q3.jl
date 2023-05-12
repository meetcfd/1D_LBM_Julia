using Plots
# LBM - 1D, diffusion equation D1Q3

m = 100 #lattices
n = 200 # time_steps

x = collect(range(0,m)) # spatial domain
ϕ = zeros(m+1) 
ϕ_wall = 1. # boundary condition at high temperature wall

α = 0.25 # thermal diffusivity
ω = 1. / (3. * α + .5)

w0 = 4. / 6.
w = 1. / 6.

f0 = w0 .* ϕ
f1 = w .* ϕ
f2 = w .* ϕ

function collision!(non_dim_dep_var, relax_fac, f0, f1, f2, w0, w)
    f0eq = w0 .* non_dim_dep_var
    feq = w .* non_dim_dep_var
    f0 .= (1. - relax_fac) .* f0 .+ relax_fac .* f0eq
    f1 .= (1. - relax_fac) .* f1 .+ relax_fac .* feq
    f2 .= (1. - relax_fac) .* f2 .+ relax_fac .* feq 
    #return f0, f1, f2
end

function streaming!(f1, f2)
   f1[begin + 1:end] .= f1[begin:end - 1]
   f2[begin:end - 1] .= f2[begin + 1:end]
   #return f1, f2
end

function boundary_conditions!(wall_temperature, f0, f1, f2)
    f0[end] = f0[end - 1]
    f1[begin] = wall_temperature - f2[begin] - f0[begin]
    f2[end] = f2[end - 1]
end

function time_step!(f0, f1, f2, ϕ, ω, ϕ_wall)
    collision!(ϕ, ω, f0, f1, f2, w0, w)
    streaming!(f1, f2)
    boundary_conditions!(ϕ_wall, f0, f1, f2)
    # f0[end] = f0[end - 1]
    # f1[begin] = ϕ_wall - f2[begin] - f0[begin]     
    # f2[end] = f2[end - 1]
    ϕ .= f0 + f1 + f2
    #return f0, f1, f2, ϕ
end

P = plot(layout=(2,1));

for t in 1:n
    time_step!(f0, f1, f2, ϕ, ω, ϕ_wall)
    if mod(t,50) == 0
        flux = 3. * ω .* (f1 - f2)
        plot!(P, x, [ϕ,flux], layout=(2,1), label=t)
    end
end

display(P)