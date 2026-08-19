using Plots

# 1. Smooth Primitives & Z0HR Scheme
@inline smooth_max(x::T, ϵ::T) where {T<:AbstractFloat} = T(0.5) * (x + sqrt(x * x + ϵ * ϵ))
@inline smooth_min(x::T, ϵ::T) where {T<:AbstractFloat} = T(0.5) * (x - sqrt(x * x + ϵ * ϵ))

function z0hr_stability(Ri::T; α_θ::T=1.0, B_um::T=16.0, B_uh::T=16.0, Ri_c::T=0.25, ϵ::T=1e-3) where {T<:AbstractFloat}
    Ri_pos = smooth_max(Ri, ϵ)
    Ri_neg = smooth_min(Ri, ϵ)

    g_raw = one(T) - (Ri_pos / Ri_c)
    stable_factor = smooth_max(g_raw, ϵ)

    unstab_m = one(T) - (B_um * Ri_neg)
    unstab_h = one(T) - (B_uh * Ri_neg)

    S_m = (stable_factor^2) * sqrt(unstab_m)
    S_h = (one(T) / α_θ) * (stable_factor^2) * (unstab_h^T(0.75))

    return S_m, S_h
end

# 2. Classic Piecewise Businger-Dyer Formulation (C^0 Kinks)
function classic_businger_dyer(Ri::T; α_θ::T=1.0, B_um::T=16.0, B_uh::T=16.0, Ri_c::T=0.25) where {T<:AbstractFloat}
    if Ri < zero(T)
        S_m = sqrt(one(T) - B_um * Ri)
        S_h = (one(T) / α_θ) * (one(T) - B_uh * Ri)^T(0.75)
    elseif Ri <= Ri_c
        stable_factor = one(T) - Ri / Ri_c
        S_m = stable_factor^2
        S_h = (one(T) / α_θ) * (stable_factor^2)
    else
        S_m = zero(T)
        S_h = zero(T)
    end
    return S_m, S_h
end

# 3. Grid Sampling
Ri_grid = range(-0.8, 0.35, length=1000)

z0hr_res = [z0hr_stability(ri; B_um=16.0, B_uh=16.0) for ri in Ri_grid]
classic_res = [classic_businger_dyer(ri; B_um=16.0, B_uh=16.0) for ri in Ri_grid]

Sm_z0hr = getindex.(z0hr_res, 1)
Sh_z0hr = getindex.(z0hr_res, 2)
Sm_classic = getindex.(classic_res, 1)
Sh_classic = getindex.(classic_res, 2)

# 4. Plot Generation
p = plot(
    xlabel="Gradient Richardson Number (Ri)",
    ylabel="Stability Function Value",
    title="Z0HR Regularization vs. Classic Piecewise Businger-Dyer (B_u = 16)",
    legend=:topright,
    grid=true,
    framestyle=:box
)

# Plot Classic Piecewise (Dashed / Reference)
plot!(p, Ri_grid, Sm_classic, label="Classic S_m (Piecewise)", lw=2, color=:black, linestyle=:dash)
plot!(p, Ri_grid, Sh_classic, label="Classic S_h (Piecewise)", lw=2, color=:gray, linestyle=:dashdot)

# Plot Z0HR (Solid / Continuous)
plot!(p, Ri_grid, Sm_z0hr, label="Z0HR S_m (Regularized)", lw=2.5, color=:blue)
plot!(p, Ri_grid, Sh_z0hr, label="Z0HR S_h (Regularized)", lw=2.5, color=:red)

# Threshold Markers
vline!(p, [0.0], color=:gray, linestyle=:dot, label="Neutral (Ri=0)")
vline!(p, [0.25], color=:darkgreen, linestyle=:dot, label="Cutoff (Ri_c=0.25)")

display(p)
savefig(p, "z0hr_vs_classic_bu16.png")