using Plots

# 1. Smooth Algebraic Primitives
@inline smooth_max(x::T, ϵ::T) where {T <: AbstractFloat} = T(0.5) * (x + sqrt(x * x + ϵ * ϵ))
@inline smooth_min(x::T, ϵ::T) where {T <: AbstractFloat} = T(0.5) * (x - sqrt(x * x + ϵ * ϵ))

# 2. Z0HR Scheme Core Function
function z0hr_stability(
    Ri::T;
    α_θ::T = 1.0,
    B_um::T = 16.0,
    B_uh::T = 16.0,
    Ri_c::T = 0.25,
    ϵ::T    = 1e-3
) where {T <: AbstractFloat}

    # Smooth coordinates
    Ri_pos = smooth_max(Ri, ϵ)
    Ri_neg = smooth_min(Ri, ϵ)

    # Smooth stable factor
    g_raw         = one(T) - (Ri_pos / Ri_c)
    stable_factor = smooth_max(g_raw, ϵ)

    # Unstable radicands (strictly > 1.0)
    unstab_m = one(T) - (B_um * Ri_neg)
    unstab_h = one(T) - (B_uh * Ri_neg)

    # Functions
    S_m = (stable_factor^2) * sqrt(unstab_m)
    S_h = (one(T) / α_θ) * (stable_factor^2) * (unstab_h^T(0.75))

    return S_m, S_h
end

# 3. Domain Sampling (Convective to Strongly Stable)
Ri_grid = range(-1.0, 0.4, length=1000)

# Evaluate over the grid
results = [z0hr_stability(ri; B_um=16.0, B_uh=16.0) for ri in Ri_grid]
S_m_vals = getindex.(results, 1)
S_h_vals = getindex.(results, 2)

# 4. Generate Plot
p = plot(
    Ri_grid, S_m_vals,
    label="S_m (Momentum, B_u = 16)",
    lw=2.5,
    color=:blue,
    xlabel="Gradient Richardson Number (Ri)",
    ylabel="Stability Function Value",
    title="Z0HR Scheme Diagnostics (B_u,m = B_u,h = 16.0)",
    legend=:topright,
    grid=:true,
    framestyle=:box
)

plot!(
    p, Ri_grid, S_h_vals,
    label="S_h (Heat, B_u = 16)",
    lw=2.5,
    color=:red,
    linestyle=:dash
)

# Reference guidelines
vline!(p, [0.0], color=:gray, linestyle=:dot, label="Neutral (Ri = 0)")
vline!(p, [0.25], color=:black, linestyle=:dashdot, label="Critical Threshold (Ri_c = 0.25)")
hline!(p, [0.0], color=:black, lw=0.5, label="")

# Display and Save
display(p)
savefig(p, "z0hr_bu16_diagnostics.png")