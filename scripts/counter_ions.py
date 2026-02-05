import numpy as np
import argparse

# ----------------------------
# Spherical <-> Cartesian
# ----------------------------
def sph_to_cart(theta, phi):
    """theta in [0, pi], phi in [0, 2pi) -> (N,3) unit vectors."""
    st, ct = np.sin(theta), np.cos(theta)
    sp, cp = np.sin(phi), np.cos(phi)
    x = cp * st
    y = sp * st
    z = ct
    return np.stack([x, y, z], axis=-1)

def cart_to_sph(xyz):
    """(N,3) -> (theta, phi) canonical angles."""
    x, y, z = xyz.T
    # numerical safety: re-normalize
    r = np.linalg.norm(xyz, axis=1, keepdims=True)
    xyz = xyz / np.clip(r, 1e-15, None)
    x, y, z = xyz.T
    theta = np.arccos(np.clip(z, -1.0, 1.0))  # [0, pi]
    phi = np.mod(np.arctan2(y, x), 2*np.pi)   # [0, 2pi)
    return theta, phi


# ----------------------------
# Fibonacci-sphere initializer
# ----------------------------
def fibonacci_sphere_angles(N, seed=None):
    """
    Even-ish spread on S^2 (unit). Returns theta, phi.
    Uses the 'equal area' latitude distribution + golden angle longitudes.
    """
    rng = np.random.default_rng(seed)
    k = np.arange(N, dtype=float)
    z = (2.0*k + 1.0)/N - 1.0  # in (-1,1)
    theta = np.arccos(z)       # [0, pi]
    golden_angle = np.pi * (3.0 - np.sqrt(5.0))
    phi = np.mod(k * golden_angle, 2*np.pi)
    # tiny jitter to avoid exact symmetries/singularities
    phi = np.mod(phi + 1e-6 * rng.standard_normal(N), 2*np.pi)
    return theta, phi


# ----------------------------
# Energy and gradient on angles
# ----------------------------
def energy_and_grad(theta, phi):
    """
    Coulomb energy E = sum_{i<j} 1/||ri - rj|| on unit sphere
    and gradient wrt (theta_i, phi_i) computed via Euclidean forces
    projected onto the spherical basis (∂r/∂θ, ∂r/∂φ).
    Returns: E, g_theta, g_phi
    """
    N = theta.shape[0]
    R = sph_to_cart(theta, phi)  # (N,3)

    # Pairwise differences and distances
    diff = R[:, None, :] - R[None, :, :]                     # (N,N,3)
    r2 = np.einsum('ijk,ijk->ij', diff, diff)                # (N,N)
    np.fill_diagonal(r2, 1.0)                                # avoid /0; won't be used
    r = np.sqrt(r2)
    inv_r = 1.0 / r
    # Energy: sum over i<j
    E = np.sum(np.triu(inv_r, 1))

    # Forces: Fi = dE/dri = sum_{j!=i} (ri - rj)/|ri - rj|^3
    inv_r3 = inv_r**3
    np.fill_diagonal(inv_r3, 0.0)
    F = np.einsum('ij,ijk->ik', inv_r3, diff)  # (N,3)

    # Local spherical basis derivatives
    st, ct = np.sin(theta), np.cos(theta)
    sp, cp = np.sin(phi), np.cos(phi)

    # ∂r/∂θ and ∂r/∂φ (not unit; that's intended for correct chain rule)
    e_theta = np.stack([cp*ct, sp*ct, -st], axis=-1)          # (N,3)
    e_phi   = np.stack([-sp*st, cp*st, np.zeros_like(st)], axis=-1)

    # Gradients via chain rule: dE/dθ_i = Fi · ∂r_i/∂θ_i, similarly for φ
    g_theta = np.einsum('ij,ij->i', F, e_theta)
    g_phi   = np.einsum('ij,ij->i', F, e_phi)

    return E, g_theta, g_phi


def energy_and_forces_xyz(theta, phi, eps2=0.0):
    """Return energy and Euclidean forces Fi = dE/dri for softened 1/sqrt(r^2+eps2)."""
    R = sph_to_cart(theta, phi)                     # (N,3)
    diff = R[:, None, :] - R[None, :, :]            # (N,N,3)
    r2 = np.einsum('ijk,ijk->ij', diff, diff)       # (N,N)
    np.fill_diagonal(r2, 1.0)                       # won't be used
    r2_soft = r2 + eps2
    inv_r = 1.0 / np.sqrt(r2_soft)
    E = np.sum(np.triu(inv_r, 1))                   # energy
    inv_r3 = inv_r**3
    np.fill_diagonal(inv_r3, 0.0)
    F = np.einsum('ij,ijk->ik', inv_r3, diff)       # (N,3)
    return E, R, F

def geodesic_step(R, V):
    """
    Move each point along tangent vector V by geodesic amount |V|.
    R, V: (N,3); R assumed unit; V tangent (orthogonal to R).
    Returns new R'.
    """
    # step angle per particle (arc length on S^2)
    w = np.linalg.norm(V, axis=1)                   # (N,)
    # avoid division by zero
    mask = w > 0
    Vhat = np.zeros_like(V)
    Vhat[mask] = V[mask] / w[mask, None]
    # Rodrigues on the sphere: cos(w) R + sin(w) Vhat
    R_new = np.copy(R)
    cw = np.cos(w)
    sw = np.sin(w)
    R_new = cw[:, None]*R + sw[:, None]*Vhat
    # re-normalize (numerical safety)
    R_new /= np.linalg.norm(R_new, axis=1, keepdims=True)
    return R_new

def minimize_coulomb_on_sphere(
    Q, radius, center,
    steps=2000,
    eta0=0.2,
    decay=800.0,
    eta_min=1e-4,
    step_max=0.25,
    backtrack_max=20,
    armijo_c=1e-4,
    eps0=1e-3,
    eps_final=0.0,
    seed=None,
    verbose=False,
    verbose_every=100,
    grad_tol=1e-3,     # <<< NEW: stop when max_i ||Ft_i|| <= grad_tol
    patience=3         # <<< NEW: require this many consecutive passes
):
    theta, phi = fibonacci_sphere_angles(Q, seed=seed)
    eta = float(eta0)

    def softening_at(k):
        t = min(1.0, k / max(1, steps-1))
        eps = (1-t)*eps0 + t*eps_final
        return eps*eps

    E, R, F = energy_and_forces_xyz(theta, phi, eps2=softening_at(0))
    calm = 0  # counts consecutive iterations under tol

    for k in range(steps):
        # project to tangent plane (F = -∇E)
        proj = np.einsum('ij,ij->i', F, R)[:, None]
        Ft = F - proj * R
        Ft_norm = np.linalg.norm(Ft, axis=1)
        max_norm = np.max(Ft_norm)
        rms_norm = np.sqrt(np.mean(Ft_norm**2))

        # --- early stopping check BEFORE attempting a step ---
        if max_norm <= grad_tol:
            calm += 1
            if calm >= patience:
                if verbose:
                    print(f"Early stop at iter {k}: max||Ft||={max_norm:.3e} (patience {patience})")
                break
        else:
            calm = 0

        # step size schedule + trust region
        step_scale = max(eta_min, eta0 / (1.0 + k/decay))
        w = step_scale * Ft_norm
        clip = np.minimum(1.0, step_max / (w + 1e-15))
        dir_t = Ft / (Ft_norm[:, None] + 1e-15)
        V = dir_t * (step_scale * clip)[:, None]   # downhill since F = -∇E

        # Armijo backtracking
        eps2 = softening_at(k)
        E0 = E
        alpha = 1.0
        accepted = False

        for _ in range(backtrack_max):
            R_try = geodesic_step(R, alpha*V)
            th_try, ph_try = cart_to_sph(R_try)
            E_try, _, F_try = energy_and_forces_xyz(th_try, ph_try, eps2=eps2)
            # directional decrease prediction: dpred = (-∇E)·(αV) with -∇E = F
            dpred = np.sum(F * (alpha * V))
            if E_try <= E0 - armijo_c * dpred:
                accepted = True
                break
            alpha *= 0.5

        if not accepted:
            # nothing helpful found; treat as converged if forces are tiny
            if verbose:
                print(f"Backtrack failed at iter {k}; max||Ft||={max_norm:.3e}")
            break

        # accept
        R, E, F = R_try, E_try, F_try
        theta, phi = cart_to_sph(R)

        if verbose and (k % verbose_every == 0 or k == steps-1):
            print(f"iter {k:5d} | E={E:.8f} | eta={step_scale:.5f} | ||Ft||_max={max_norm:.3e} | ||Ft||_rms={rms_norm:.3e}")

    # optional SciPy polish unchanged ...
    # (keep your L-BFGS-B block)
    xyz = (sph_to_cart(theta, phi)*radius)+center

    return xyz


if __name__ == "__main__":
    parser = argparse.ArgumentParser()#description="Process some arguments.")
    parser.add_argument("-Q", "--charge", type=int, help="Total charge in sphere")
    parser.add_argument("-r", "--radius", type=float, help="Radius of the sphere")
    parser.add_argument("-c", "--center", nargs=3, type=float, help="Center of the sphere (x, y, z)")
    parser.add_argument("-p", "--plot", action="store_true", default=False, help="plot 3D sphere")
    args = parser.parse_args()
    Q = args.charge
    r = args.radius
    c = args.center
    theta, phi, xyz, E = minimize_coulomb_on_sphere(
        Q, r, c,
        steps=2000,
        grad_tol=1e-6,
        patience=3,
        verbose=False,
        seed=42
    )
    print(f"Final energy: {E:.10f}")
    print(xyz)
    if args.plot == True:
        plot_points_on_sphere(xyz, title=f"N={Q}, E={E:.6f}")   # should remove the E from plotting function
