module PartialBellMapping

using LinearAlgebra

export build_bell_matrix, invert_bell_matrix, det_bell_matrix
export curvature_decomposition, fold_knee_scaling, subgrid_knee_detection

"""
    build_bell_matrix(g1, g2, g3, g4)

Constructs the 4x4 lower-triangular Partial Bell Transformation Matrix B
mapping intrinsic similarity derivatives R_ζ to physical spatial derivatives Ri_z:
Ri_z = B * R_ζ
where g_j = dʲζ / dzʲ.
"""
function build_bell_matrix(g1::T, g2::T, g3::T, g4::T) where {T<:AbstractFloat}
    B = zeros(T, 4, 4)
    # Row 1 (m=1)
    B[1, 1] = g1
    # Row 2 (m=2)
    B[2, 1] = g2
    B[2, 2] = g1 * g1
    # Row 3 (m=3)
    B[3, 1] = g3
    B[3, 2] = T(3.0) * g1 * g2
    B[3, 3] = g1^3
    # Row 4 (m=4)
    B[4, 1] = g4
    B[4, 2] = T(4.0) * g1 * g3 + T(3.0) * (g2^2)
    B[4, 3] = T(6.0) * (g1^2) * g2
    B[4, 4] = g1^4
    return B
end

"""
    det_bell_matrix(g1)

Computes the analytical determinant of the Partial Bell Matrix: det(B) = g₁¹⁰.
"""
@inline det_bell_matrix(g1::T) where {T<:AbstractFloat} = g1^10

"""
    invert_bell_matrix(g1, g2, g3, g4)

Computes the exact closed-form analytical lower-triangular inverse Partial Bell Matrix B⁻¹:
R_ζ = B⁻¹ * Ri_z
without dynamic matrix LU factorization or floating-point pivoting errors.
"""
function invert_bell_matrix(g1::T, g2::T, g3::T, g4::T) where {T<:AbstractFloat}
    invB = zeros(T, 4, 4)

    g1_sq = g1 * g1
    g1_3  = g1_sq * g1
    g1_4  = g1_3 * g1
    g1_5  = g1_4 * g1
    g1_6  = g1_5 * g1
    g1_7  = g1_6 * g1

    g2_sq = g2 * g2
    g2_3  = g2_sq * g2

    # Row 1
    invB[1, 1] = one(T) / g1

    # Row 2
    invB[2, 1] = -g2 / g1_3
    invB[2, 2] = one(T) / g1_sq

    # Row 3
    invB[3, 1] = (T(3.0) * g2_sq - g1 * g3) / g1_5
    invB[3, 2] = -T(3.0) * g2 / g1_4
    invB[3, 3] = one(T) / g1_3

    # Row 4
    invB[4, 1] = (-g1_sq * g4 + T(10.0) * g1 * g2 * g3 - T(15.0) * g2_3) / g1_7
    invB[4, 2] = (T(15.0) * g2_sq - T(4.0) * g1 * g3) / g1_6
    invB[4, 3] = -T(6.0) * g2 / g1_5
    invB[4, 4] = one(T) / g1_4

    return invB
end

"""
    curvature_decomposition(R_prime, R_dprime, g1, g2)

Decomposes total physical profile curvature into three GSPT layers:
- C_const: Intrinsic constitutive similarity curvature R''(ζ) * g₁²
- C_coord: Geometric coordinate/flux curvature R'(ζ) * g₂
- R_coord: Geometric diagnostic ratio (C_coord / C_const)
"""
function curvature_decomposition(R_prime::T, R_dprime::T, g1::T, g2::T) where {T<:AbstractFloat}
    C_const = R_dprime * (g1^2)
    C_coord = R_prime * g2
    Ri_zz_total = C_const + C_coord
    R_coord = iszero(C_const) ? T(Inf) : C_coord / C_const

    return (C_const = C_const, C_coord = C_coord, Total_Curvature = Ri_zz_total, Diagnostic_Ratio = R_coord)
end

"""
    fold_knee_scaling(R_prime, g2, p)

Evaluates the double-factorial sequence scaling B_{2p, p}|_{g₁=0} = (2p-1)!! g₂ᵖ
at a coordinate fold knee (g₁ = ζ_z = 0, g₂ ≠ 0).
"""
function fold_knee_scaling(R_prime::T, g2::T, p::Int) where {T<:AbstractFloat}
    double_factorial = 1
    for k in 1:2:2p-1
        double_factorial *= k
    end
    return T(double_factorial) * (g2^p) * R_prime
end

"""
    subgrid_knee_detection(z_grid, curvature)

Locates sub-grid coordinate fold knee heights (z_knee) using linear zero-crossing
interpolation across profile curvature extrema.
"""
function subgrid_knee_detection(z_grid::Vector{T}, curvature::Vector{T}) where {T<:AbstractFloat}
    n = length(z_grid)
    k = argmin(abs.(curvature[2:n-1])) + 1
    i1, i2 = (k > 2 && curvature[k-1] * curvature[k] <= zero(T)) ? (k - 1, k) : (k, k + 1)

    c1, c2 = curvature[i1], curvature[i2]
    denom = c2 - c1
    if iszero(denom)
        return z_grid[k]
    else
        return z_grid[i1] - c1 * (z_grid[i2] - z_grid[i1]) / denom
    end
end

end # module
