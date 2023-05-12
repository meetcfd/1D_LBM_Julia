using Plots

m = 100 # lattices
n = 200 # timesteps

x = collect(range(0,1,m+1)) #indep variable
ϕ = exp.(-100 .* (x .- 0.25).^2) # dep variable

w0 = 4/6
w = 1/6
c² = 1/3

u = 0.1 # advec velocity
α = 0. # diff coeff
ϕ_wall = 0.

ω = 1 / (α + .5)

f0 = w0 .* ϕ
f1 = w .* ϕ
f2 = w .* ϕ 

feq0 = zeros(m+1)
feq1 = zeros(m+1)
feq2 = zeros(m+1)

function feq!(eq_dis, w, ϕ, u, c², k)
    if k == 0
        eq_dis .= w .* ϕ 
    elseif k == 1
        eq_dis .= w .* ϕ .* (1 + u / c²)
    elseif k == 2
        eq_dis .= w .* ϕ .* (1 - u / c²)
    end
end

function collision!(f0, f1, f2, feq0, feq1, feq2, w, w0, ϕ, u, c², ω)
    feq!(feq0, w0, ϕ, u, c², 0)
    feq!(feq1, w, ϕ, u, c², 1)
    feq!(feq2, w, ϕ, u, c², 2)
    f0 .= (1 - ω) .* f0 + ω .* feq0
    f1 .= (1 - ω) .* f1 + ω .* feq1
    f2 .= (1 - ω) .* f2 + ω .* feq2 
end

function streaming!(f1, f2)
    f1[begin + 1:end] .= f1[begin:end-1]
    f2[begin:end - 1] .= f2[begin + 1: end]
end

function boundary_conditions!(f0, f1, f2, ϕ_wall)
    f0[end] = f0[end - 1]
    f1[begin] = ϕ_wall - f0[begin] - f2[begin]
    f2[end] = f2[end - 1] 
end

function time_step!(ϕ, f0, f1, f2, feq0, feq1, feq2, w, w0, u, c², ω, ϕ_wall)
    collision!(f0, f1, f2, feq0, feq1, feq2, w, w0, ϕ, u, c², ω)
    streaming!(f1, f2)
    boundary_conditions!(f0, f1, f2, ϕ_wall)
    ϕ .= f0 .+ f1 .+ f2
end

P = plot();

for t in 1:n
    time_step!(ϕ, f0, f1, f2, feq0, feq1, feq2, w, w0, u, c², ω, ϕ_wall)
    if mod(t,50) == 0
        plot!(P, x, ϕ, layout=(2,1), label=t)
    end
end

display(P)
