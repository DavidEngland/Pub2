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
Ri_grid = range(-0.6, 0.35, length=2000)

z0hr_res = [z0hr_stability(ri; B_um=16.0, B_uh=16.0) for ri in Ri_grid]
classic_res = [classic_businger_dyer(ri; B_um=16.0, B_uh=16.0) for ri in Ri_grid]

Sm_z0hr = getindex.(z0hr_res, 1)
Sh_z0hr = getindex.(z0hr_res, 2)
Sm_classic = getindex.(classic_res, 1)
Sh_classic = getindex.(classic_res, 2)

# 4. Main Panel: Full Domain Overview
p1 = plot(
    xlabel="Gradient Richardson Number (Ri)",
    ylabel="Stability Function Value",
    title="Overview: Z0HR vs Classic (B_u = 16)",
    legend=:topright,
    grid=true,
    framestyle=:box
)

# Plot Classic underneath (Thick, semi-transparent black/gray)
plot!(p1, Ri_grid, Sm_classic, label="Classic S_m", lw=5, alpha=0.35, color=:black)
plot!(p1, Ri_grid, Sh_classic, label="Classic S_h", lw=5, alpha=0.35, color=:gray)

# Plot Z0HR on top (Crisp, thinner dashed lines)
plot!(p1, Ri_grid, Sm_z0hr, label="Z0HR S_m", lw=2, linestyle=:dash, color=:blue)
plot!(p1, Ri_grid, Sh_z0hr, label="Z0HR S_h", lw=2, linestyle=:dash, color=:red)

# 5. Subpanel: Zoomed Inset at Neutral State (Ri ≈ 0)
p2 = plot(
    xlims=(-0.04, 0.04), ylims=(0.7, 1.3),
    title="Zoom: Neutral State (Ri ≈ 0)",
    xlabel="Ri", ylabel="S",
    grid=true, legend=false, framestyle=:box
)
plot!(p2, Ri_grid, Sm_classic, lw=4, alpha=0.35, color=:black)
plot!(p2, Ri_grid, Sm_z0hr, lw=2, linestyle=:dash, color=:blue)
plot!(p2, Ri_grid, Sh_classic, lw=4, alpha=0.35, color=:gray)
plot!(p2, Ri_grid, Sh_z0hr, lw=2, linestyle=:dash, color=:red)

# 6. Subpanel: Zoomed Inset at Critical Threshold (Ri ≈ 0.25)
p3 = plot(
    xlims=(0.20, 0.28), ylims=(-0.01, 0.08),
    title="Zoom: Cutoff (Ri ≈ 0.25)",
    xlabel="Ri", ylabel="S",
    grid=true, legend=false, framestyle=:box
)
plot!(p3, Ri_grid, Sm_classic, lw=4, alpha=0.35, color=:black)
plot!(p3, Ri_grid, Sm_z0hr, lw=2, linestyle=:dash, color=:blue)
plot!(p3, Ri_grid, Sh_classic, lw=4, alpha=0.35, color=:gray)
plot!(p3, Ri_grid, Sh_z0hr, lw=2, linestyle=:dash, color=:red)

# Combine into a 2-row multi-panel layout
final_plot = plot(p1, plot(p2, p3, layout=(1, 2)), layout=(2, 1), size=(800, 700))

display(final_plot)
savefig(final_plot, "z0hr_vs_classic_bu16_multi.png")