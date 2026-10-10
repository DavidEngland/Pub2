module Z0HR

export smooth_max, smooth_min, z0hr_split
export z0hr_stability, compute_lagrange_coeffs, eval_lagrange_kernel, solve_zeta_z0hr
export intrinsic_derivatives

"""
    smooth_max(x, ϵ)

Infinitely differentiable C^∞ smooth maximum operator:
smooth_max(x; ϵ) = (x + √(x² + ϵ²)) / 2
"""
@inline smooth_max(x::T, ϵ::T) where {T<:AbstractFloat} = T(0.5) * (x + sqrt(x * x + ϵ * ϵ))

"""
    smooth_min(x, ϵ)

Infinitely differentiable C^∞ smooth minimum operator:
smooth_min(x; ϵ) = (x - √(x² + ϵ²)) / 2
"""
@inline smooth_min(x::T, ϵ::T) where {T<:AbstractFloat} = T(0.5) * (x - sqrt(x * x + ϵ * ϵ))

"""
    z0hr_split(Ri; ϵ=1e-3)

Zero-Offset Hyperbolic Regularization (Z0HR) coordinate split.
Decomposes `Ri` into smooth positive (stable) `Ri_pos` and negative (unstable) `Ri_neg`
components preserving origin invariance (Ri⁺(0) = 0, Ri⁻(0) = 0).
"""
@inline function z0hr_split(Ri::T; ϵ::T=T(1e-3)) where {T<:AbstractFloat}
    sqrt_disc_ri = sqrt(Ri * Ri + ϵ * ϵ)
    Ri_pos = T(0.5) * (Ri + sqrt_disc_ri - ϵ)
    Ri_neg = T(0.5) * (Ri - sqrt_disc_ri + ϵ)
    w_stable = T(0.5) * (one(T) + Ri / sqrt_disc_ri)
    return Ri_pos, Ri_neg, w_stable
end

"""
    z0hr_stability(Ri; α_θ=1.0, B_um=4.0/Ri_c, B_uh=8.0/(3Ri_c), Ri_c=0.25, ϵ=1e-3)

Branch-free, C^∞-smooth evaluation of momentum (`S_m`) and heat/scalar (`S_h`)
stability functions across all stability regimes without if/else branching.
"""
@inline function z0hr_stability(
    Ri::T;
    α_θ::T  = one(T),
    Ri_c::T = T(0.25),
    B_um::T = T(4.0) / Ri_c,
    B_uh::T = T(8.0) / (T(3.0) * Ri_c),
    ϵ::T    = T(1e-3)
) where {T<:AbstractFloat}

    # 1. Smooth Hyperbolic Coordinates
    sqrt_disc_ri = sqrt(Ri * Ri + ϵ * ϵ)
    Ri_pos       = T(0.5) * (Ri + sqrt_disc_ri - ϵ)
    Ri_neg       = T(0.5) * (Ri - sqrt_disc_ri + ϵ)

    # 2. Smooth Stable Extinction Factor (Algebraic Smooth Max)
    g_raw         = one(T) - (Ri_pos / Ri_c)
    sqrt_disc_g   = sqrt(g_raw * g_raw + ϵ * ϵ)
    stable_factor = T(0.5) * (g_raw + sqrt_disc_g)

    # 3. Unstable Terms (Naturally > 1.0 everywhere for B_u > 0)
    unstab_m = one(T) - (B_um * Ri_neg)
    unstab_h = one(T) - (B_uh * Ri_neg)

    # 4. Unified C^∞ Stability Functions
    S_m = (stable_factor^2) * sqrt(unstab_m)
    S_h = (one(T) / α_θ) * (stable_factor^2) * (unstab_h^T(0.75))

    return S_m, S_h
end

"""
    compute_lagrange_coeffs(Pr_t0, β_m, β_h)

Derives forward Taylor expansion coefficients c_k for Ri(ζ) = c₁ζ + c₂ζ² + c₃ζ³ + O(ζ⁴)
and universal Lagrange reversion coefficients a_k for ζ(Ri) = a₁Ri + a₂Ri² + a₃Ri³ + O(Ri⁴).
"""
@inline function compute_lagrange_coeffs(Pr_t0::T, β_m::T, β_h::T) where {T<:AbstractFloat}
    c1 = Pr_t0
    c2 = β_h - T(2.0) * β_m * Pr_t0
    c3 = T(3.0) * β_m * β_m * Pr_t0 - T(2.0) * β_m * β_h

    c1_inv = one(T) / c1
    a1 = c1_inv
    a2 = -c2 * (c1_inv^3)
    a3 = (T(2.0) * c2 * c2 - c1 * c3) * (c1_inv^5)

    return (a1, a2, a3)
end

"""
    eval_lagrange_kernel(Ri, (a1, a2, a3))

Evaluates 3rd-order near-neutral Lagrange series reversion via Horner's scheme.
"""
@inline function eval_lagrange_kernel(Ri::T, a::NTuple{3, T}) where {T<:AbstractFloat}
    a1, a2, a3 = a
    return Ri * (a1 + Ri * (a2 + Ri * a3))
end

"""
    solve_zeta_z0hr(Ri, a_coeffs; Ri_c=0.19, Pr_t0=1.0, ϵ=1e-3)

Branch-free solver evaluating similarity parameter ζ(Ri) across all stability regimes
using smooth Z0HR weighting to blend Lagrange reversion and stable quadratic closures.
"""
@inline function solve_zeta_z0hr(
    Ri::T,
    a_coeffs::NTuple{3, T};
    Ri_c::T  = T(0.19),
    Pr_t0::T = one(T),
    ϵ::T     = T(1e-3)
) where {T<:AbstractFloat}

    Ri_pos, Ri_neg, w_stable = z0hr_split(Ri; ϵ=ϵ)
    zeta_neutral = eval_lagrange_kernel(Ri, a_coeffs)

    denom_stable = max(one(T) - Ri_pos / Ri_c, ϵ)
    zeta_stable = (Ri_pos / Pr_t0) / denom_stable

    zeta = (one(T) - w_stable) * zeta_neutral + w_stable * zeta_stable
    return zeta
end

"""
    intrinsic_derivatives(ζ; α_θ=1.0, β_m=4.7, β_h=2.35)

Computes analytical intrinsic similarity derivatives R⁽ⁿ⁾(ζ) = dⁿR/dζⁿ for n ∈ {1,2,3,4}
populating intrinsic similarity vector R_ζ in Partial Bell transformations.
Detects Prandtl-matched manifold (α_θ β_m = 2β_h) for optimized power-law evaluation.
"""
function intrinsic_derivatives(
    ζ::T;
    α_θ::T = one(T),
    β_m::T = T(4.7),
    β_h::T = T(2.35)
) where {T<:AbstractFloat}

    u = one(T) + β_m * ζ
    u2 = u * u
    u3 = u2 * u
    u4 = u3 * u
    u5 = u4 * u
    u6 = u5 * u

    # Check if parameters lie on the Prandtl-matched manifold (α_θ β_m = 2β_h)
    if isapprox(α_θ * β_m, T(2.0) * β_h; atol=T(1e-6))
        # Pure constant-numerator power law collapse
        R_prime = α_θ / u3
        R_dprime = -T(3.0) * α_θ * β_m / u4
        R_triple = T(12.0) * α_θ * (β_m^2) / u5
        R_quad   = -T(60.0) * α_θ * (β_m^3) / u6
    else
        # General linear derivative theorem
        term_lin = β_m * (T(2.0) * β_h - α_θ * β_m) * ζ

        R_prime  = (α_θ + (T(2.0) * β_h - α_θ * β_m) * ζ) / u3
        R_dprime = (T(2.0) * (β_h - T(2.0) * α_θ * β_m) - T(2.0) * term_lin) / u4
        R_triple = (T(6.0) * β_m * (T(3.0) * α_θ * β_m - T(2.0) * β_h) + T(6.0) * β_m * term_lin) / u5
        R_quad   = (-T(24.0) * (β_m^2) * (T(4.0) * α_θ * β_m - T(3.0) * β_h) - T(24.0) * (β_m^2) * term_lin) / u6
    end

    return (R_prime, R_dprime, R_triple, R_quad)
end

end # module
