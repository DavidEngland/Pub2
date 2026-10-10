"""
Z0HR_Framework.py - Unified C^∞ Hyperbolic Regularization & Partial Bell Derivative Mapping
Atmospheric Boundary Layer Hyperdiffusion Framework
"""

import numpy as np

__all__ = [
    'smooth_max', 'smooth_min', 'z0hr_split', 'z0hr_stability',
    'compute_lagrange_coeffs', 'eval_lagrange_kernel', 'solve_zeta_z0hr',
    'intrinsic_derivatives', 'build_bell_matrix', 'det_bell_matrix',
    'invert_bell_matrix', 'curvature_decomposition', 'subgrid_knee_detection'
]

def smooth_max(x, eps=1e-3):
    """C^∞ smooth maximum operator: smooth_max(x; ϵ) = (x + √(x² + ϵ²)) / 2"""
    return 0.5 * (x + np.sqrt(x**2 + eps**2))

def smooth_min(x, eps=1e-3):
    """C^∞ smooth minimum operator: smooth_min(x; ϵ) = (x - √(x² + ϵ²)) / 2"""
    return 0.5 * (x - np.sqrt(x**2 + eps**2))

def z0hr_split(Ri, eps=1e-3):
    """
    Zero-Offset Hyperbolic Regularization (Z0HR) coordinate split.
    Decomposes Ri into smooth Ri+ and Ri- projections preserving origin invariance: Ri⁺(0) = 0, Ri⁻(0) = 0.
    """
    sqrt_disc_ri = np.sqrt(Ri**2 + eps**2)
    Ri_pos = 0.5 * (Ri + sqrt_disc_ri - eps)
    Ri_neg = 0.5 * (Ri - sqrt_disc_ri + eps)
    w_stable = 0.5 * (1.0 + Ri / sqrt_disc_ri)
    return Ri_pos, Ri_neg, w_stable

def z0hr_stability(Ri, alpha_theta=1.0, Ri_c=0.25, B_um=16.0, B_uh=16.0, eps=1e-3):
    """
    Branch-free C^∞-smooth evaluation of momentum (S_m) and heat/scalar (S_h)
    stability functions across all stability regimes without conditional branching.
    """
    Ri_pos, Ri_neg, _ = z0hr_split(Ri, eps=eps)
    g_raw = 1.0 - (Ri_pos / Ri_c)
    stable_factor = smooth_max(g_raw, eps)
    unstab_m = 1.0 - (B_um * Ri_neg)
    unstab_h = 1.0 - (B_uh * Ri_neg)
    S_m = (stable_factor**2) * np.sqrt(unstab_m)
    S_h = (1.0 / alpha_theta) * (stable_factor**2) * (unstab_h**0.75)
    return S_m, S_h

def compute_lagrange_coeffs(Pr_t0, beta_m, beta_h):
    """
    Derives forward Taylor expansion coefficients c_k for Ri(ζ) = c₁ζ + c₂ζ² + c₃ζ³
    and universal Lagrange reversion coefficients a_k for ζ(Ri) = a₁Ri + a₂Ri² + a₃Ri³.
    """
    c1 = Pr_t0
    c2 = beta_h - 2.0 * beta_m * Pr_t0
    c3 = 3.0 * beta_m**2 * Pr_t0 - 2.0 * beta_m * beta_h
    c1_inv = 1.0 / c1
    a1 = c1_inv
    a2 = -c2 * (c1_inv**3)
    a3 = (2.0 * c2**2 - c1 * c3) * (c1_inv**5)
    return a1, a2, a3

def eval_lagrange_kernel(Ri, a_coeffs):
    """Evaluates 3rd-order near-neutral Lagrange series reversion via Horner's scheme."""
    a1, a2, a3 = a_coeffs
    return Ri * (a1 + Ri * (a2 + Ri * a3))

def solve_zeta_z0hr(Ri, a_coeffs, Ri_c=0.19, Pr_t0=1.0, eps=1e-3):
    """
    Branch-free solver evaluating similarity parameter ζ(Ri) across all stability regimes
    using smooth Z0HR weighting to blend Lagrange reversion and stable quadratic closures.
    """
    Ri_pos, Ri_neg, w_stable = z0hr_split(Ri, eps=eps)
    zeta_neutral = eval_lagrange_kernel(Ri, a_coeffs)
    denom_stable = np.maximum(1.0 - Ri_pos / Ri_c, eps)
    zeta_stable = (Ri_pos / Pr_t0) / denom_stable
    return (1.0 - w_stable) * zeta_neutral + w_stable * zeta_stable

def intrinsic_derivatives(zeta, alpha_theta=1.0, beta_m=5.0, beta_h=2.5):
    """
    Computes analytical intrinsic similarity derivatives R⁽ⁿ⁾(ζ) = dⁿR/dζⁿ for n ∈ {1,2,3,4}
    populating intrinsic similarity vector R_ζ in Partial Bell transformations.
    """
    u = 1.0 + beta_m * zeta
    u3, u4, u5, u6 = u**3, u**4, u**5, u**6
    if np.isclose(alpha_theta * beta_m, 2.0 * beta_h):
        R_prime = alpha_theta / u3
        R_dprime = -3.0 * alpha_theta * beta_m / u4
        R_triple = 12.0 * alpha_theta * (beta_m**2) / u5
        R_quad = -60.0 * alpha_theta * (beta_m**3) / u6
    else:
        term_lin = beta_m * (2.0 * beta_h - alpha_theta * beta_m) * zeta
        R_prime = (alpha_theta + (2.0 * beta_h - alpha_theta * beta_m) * zeta) / u3
        R_dprime = (2.0 * (beta_h - 2.0 * alpha_theta * beta_m) - 2.0 * term_lin) / u4
        R_triple = (6.0 * beta_m * (3.0 * alpha_theta * beta_m - 2.0 * beta_h) + 6.0 * beta_m * term_lin) / u5
        R_quad = (-24.0 * (beta_m**2) * (4.0 * alpha_theta * beta_m - 3.0 * beta_h) - 24.0 * (beta_m**2) * term_lin) / u6
    return R_prime, R_dprime, R_triple, R_quad

def build_bell_matrix(g1, g2, g3, g4):
    """
    Constructs the 4x4 lower-triangular Partial Bell Transformation Matrix B
    mapping intrinsic similarity derivatives R_ζ to physical spatial derivatives Ri_z: Ri_z = B * R_ζ
    """
    B = np.zeros((4, 4))
    B[0, 0] = g1
    B[1, 0] = g2
    B[1, 1] = g1**2
    B[2, 0] = g3
    B[2, 1] = 3.0 * g1 * g2
    B[2, 2] = g1**3
    B[3, 0] = g4
    B[3, 1] = 4.0 * g1 * g3 + 3.0 * (g2**2)
    B[3, 2] = 6.0 * (g1**2) * g2
    B[3, 3] = g1**4
    return B

def det_bell_matrix(g1):
    """Computes the analytical determinant of the Partial Bell Matrix: det(B) = g₁¹⁰."""
    return g1**10

def invert_bell_matrix(g1, g2, g3, g4):
    """
    Computes the exact closed-form analytical lower-triangular inverse Partial Bell Matrix B⁻¹:
    R_ζ = B⁻¹ * Ri_z
    """
    invB = np.zeros((4, 4))
    g1_sq, g1_3, g1_4, g1_5, g1_6, g1_7 = g1**2, g1**3, g1**4, g1**5, g1**6, g1**7
    g2_sq, g2_3 = g2**2, g2**3

    invB[0, 0] = 1.0 / g1
    invB[1, 0] = -g2 / g1_3
    invB[1, 1] = 1.0 / g1_sq
    invB[2, 0] = (3.0 * g2_sq - g1 * g3) / g1_5
    invB[2, 1] = -3.0 * g2 / g1_4
    invB[2, 2] = 1.0 / g1_3
    invB[3, 0] = (-g1_sq * g4 + 10.0 * g1 * g2 * g3 - 15.0 * g2_3) / g1_7
    invB[3, 1] = (15.0 * g2_sq - 4.0 * g1 * g3) / g1_6
    invB[3, 2] = -6.0 * g2 / g1_5
    invB[3, 3] = 1.0 / g1_4
    return invB

def curvature_decomposition(R_prime, R_dprime, g1, g2):
    """
    Decomposes total physical profile curvature into three GSPT layers:
    C_const (constitutive similarity), C_coord (coordinate geometry), and diagnostic ratio R_coord.
    """
    C_const = R_dprime * (g1**2)
    C_coord = R_prime * g2
    Total_Curvature = C_const + C_coord
    Diagnostic_Ratio = np.inf if C_const == 0 else C_coord / C_const
    return C_const, C_coord, Total_Curvature, Diagnostic_Ratio

def subgrid_knee_detection(z_grid, curvature):
    """Locates sub-grid coordinate fold knee heights (z_knee) using linear zero-crossing interpolation."""
    curvature = np.array(curvature)
    n = len(z_grid)
    k = np.argmin(np.abs(curvature[1:n-1])) + 1
    i1, i2 = (k-1, k) if (k > 1 and curvature[k-1] * curvature[k] <= 0) else (k, k+1)
    c1, c2 = curvature[i1], curvature[i2]
    denom = c2 - c1
    return z_grid[k] if denom == 0 else z_grid[i1] - c1 * (z_grid[i2] - z_grid[i1]) / denom
