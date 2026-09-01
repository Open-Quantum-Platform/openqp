# Analytic nuclear Hessian for MRSF-TDDFT

## Provenance and scope

This development begins from Hiroya Nakata's original MRSF response and
analytic-gradient formulation and implementation, as reported with Seunghoon
Lee, Emma Eunji Kim, Sangyoub Lee, and Cheol Ho Choi in J. Chem. Phys. **150**,
184111 (2019).  The response space and spin-pairing conventions are those of
Lee, Filatov, Lee, and Choi, J. Chem. Phys. **149**, 104101 (2018), including
its Supporting Information.  The closed-shell TDHF/TDDFT Hessian already in
OpenQP supplies verified nuclear-integral derivatives and response solvers; it
does not define the MRSF equations.

The initial analytic method has the following validity range:

- a restricted open-shell triplet reference with exactly two singly occupied
  spatial orbitals;
- singlet or triplet collinear MRSF response in the Tamm--Dancoff
  approximation;
- MRSF-TDHF and local or global-hybrid LDA/GGA MRSF-TDDFT;
- one MPI rank until distributed response contractions have an independent
  verification;
- no range separation, meta-GGA dependence, unrestricted MRSF reference, or
  quintet response in the initial implementation.

The four missing type-IV configurations of the founding two-SOMO method remain
outside the response space.  A large response norm or a state belonging to a
near-degenerate manifold is reported as an ill-conditioned state-specific
Hessian, not silently interpreted as an isolated-state result.

## Reference and response space

Let (C), (O=\{O_1,O_2\}), and (V) denote the closed, open, and virtual
orbital spaces.  The reference density is

\[
 \rho_0^{\mathrm{MR}}=\frac12\left(\rho_0^{M_S=+1}
                                  +\rho_0^{M_S=-1}\right).
\]

The complete parent-derived response topology is

\[
 \mathrm{CO}_{+},\ \mathrm{OV}_{+},\ \mathrm{CV}_{+},\
 \mathrm{OO}_{\pm},\ \mathrm{CV}_{-},\ \mathrm{OV}_{-},\
 \mathrm{CO}_{-}.
\]

The OO parent phases are aligned before spin adaptation.  With (L) and (R)
denoting the two diagonal OO parents,

\[
 X_{\mathrm{OO}}^{T}=\frac{L+R}{\sqrt2},\qquad
 X_{\mathrm{OO}}^{S}=\frac{L-R}{\sqrt2}.
\]

The off-diagonal OO configurations (G=(O_2,O_1)) and
(D=(O_1,O_2)) occur only in the singlet response.  Spin-pairing coupling is
restricted to CO--CO, CO--OV, and OV--OV contractions.  Its nuclear derivative
is part of the Hessian; the three user scale factors are independent of nuclear
coordinates.

## Stationary MRSF Lagrangian

For spin multiplicity (k\in\{S,T\}), the symmetric MRSF eigenproblem is

\[
 A^{(k)}X^{(k)}=\Omega_k X^{(k)},\qquad
 X^{(k)\mathsf T}X^{(k)}=1,
\]

where the physical spin-pairing convention of Nakata *et al.* is
(A^{(T)}=A^{(0)}+C_{\mathrm{SPC}}) and
(A^{(S)}=A^{(0)}-C_{\mathrm{SPC}}).  An implementation may store a response
channel with the opposite sign, but that internal channel must then be named
and converted explicitly; it does not change this physical definition.

OpenQP defines the special-channel action
(K_{\mathrm{raw}}[X]\equiv
\mathcal M[f_{1:6}[d_{1:6}(X)]]=-C_{\mathrm{SPC}}[X]).  Channel 7, conventionally
called `ball`, contributes to (A^{(0)}) and is not an SPC channel.  Hence the
raw special-channel multiplier is +1 for a singlet and -1 for a triplet.  All
first- and second-order derivatives retain this conversion.

The stationary Lagrangian is

\[
 \mathcal L_k = X^{\mathsf T}A^{(k)}X
 -\Omega_k(X^{\mathsf T}X-1)
 + Z^{\mathsf T}f_{\mathrm{ROHF}}
 - W:(C^{\mathsf T}SC-I).
\]

Stationarity with respect to orbital rotations gives the Guest--Saunders ROHF
equation

\[
 \bar J Z^{(k)}=-\bar R^{(k)},
\]

with the CO, CV, and OV rotations retained.  The relaxed difference density is
(P^{(k)}=T^{(k)}+Z^{(k)}), and (W^{(k)}) is the corresponding
energy-weighted density.

## First nuclear response

For a Cartesian coordinate (R_\eta), differentiation of the MRSF
eigenproblem gives

\[
 (A^{(k)}-\Omega_k I)X^{(k),\eta}
 =-Q_k\left(A^{(k),\eta}-\Omega_k^{\eta}I\right)X^{(k)},
\]

where (Q_k=I-X^{(k)}X^{(k)\mathsf T}) for an isolated root and
(X^{(k)\mathsf T}X^{(k),\eta}=0).  For a degenerate or nearly degenerate
manifold, (Q_k) is replaced by the complement of the complete state
projector, and the reported result is the projector response unless a diabatic
state definition is supplied.

The differentiated MRSF operator contains all of the following terms:

1. derivatives of the alpha and beta ROHF Fock matrices;
2. AO overlap and orbital-response terms;
3. derivatives of exact-exchange and exchange-correlation response kernels;
4. the seven CO/OV/CV/OO parent-derived density responses;
5. derivatives of every CO--CO, CO--OV, and OV--OV spin-pairing contraction.

The differentiated Z-vector equation is

\[
 \bar J Z^{(k),\eta}
 =-\left(\bar R^{(k),\eta}+\bar J^{\eta}Z^{(k)}\right).
\]

The derivatives (T^{(k),\eta}), (P^{(k),\eta}), and
(W^{(k),\eta}) are then formed from the orbital response,
(X^{(k),\eta}), and (Z^{(k),\eta}).  The two-particle difference density
(\Gamma^{(k),\eta}) includes the derivative of the ordinary MRSF term and of
the spin-pairing term.  Omitting the latter is equivalent to differentiating
MRSF(0), not MRSF-TDDFT.

## Cartesian Hessian

The excitation-energy gradient is written in the AO basis as

\[
 \Omega_k^\xi=h^\xi:P^{(k)}-S^\xi:W^{(k)}
                 +g^\xi:\Gamma^{(k)}+G_{\mathrm{xc}}^\xi.
\]

Its derivative with respect to (R_\eta) is

\[
\begin{aligned}
 \Omega_k^{\xi\eta}={}&h^{\xi\eta}:P^{(k)}
 +h^\xi:P^{(k),\eta}
 -S^{\xi\eta}:W^{(k)}-S^\xi:W^{(k),\eta}\\
 &+g^{\xi\eta}:\Gamma^{(k)}+g^\xi:\Gamma^{(k),\eta}
 +G_{\mathrm{xc}}^{\xi\eta}.
\end{aligned}
\]

The total state Hessian is

\[
 H_{k,\xi\eta}=E_{\mathrm{ROHF}}^{\xi\eta}
                  +\Omega_k^{\xi\eta}.
\]

Both Cartesian index orders are evaluated before the final symmetric average.
The unsymmetrized difference is retained as a response residual.  Translational
projection is applied only after the complete reference and excitation terms
have been assembled.

For LDA/GGA functionals, (G_{\mathrm{xc}}^{\xi\eta}) contains the fixed-grid
functional derivatives, basis-function motion, atom-centred grid motion, and
first and second derivatives of normalized partition weights.  Global exact
exchange is treated by the same two-electron derivative equations as
MRSF-TDHF.  A functional is not enabled until each of these terms agrees with
finite differences of the analytic MRSF gradient.

## Independent verification criteria

Molecular calculations begin only after the following finite model criteria
are satisfied:

1. an isolated second-quantized audit recovers the complete spin-adapted
   seven-part response topology, with no duplicate CSFs and with the stated
   L/R phases; this audit is not the production representation;
2. an independently constructed Slater--Condon Hamiltonian and the MRSF
   matrix-vector product agree for every unit-vector column to
   (10^{-10}\ E_h);
3. singlet and triplet OO functions are normalized, mutually orthogonal, and
   eigenfunctions of \(\hat S^2\) with eigenvalues 0 and 2;
4. analytic first and second eigenvalue derivatives agree with central finite
   differences at (h), (h/2), and (h/4), with the expected convergence;
5. the MRSF(0) limit is recovered when all spin-pairing scales are zero, and
   the SF limit is recovered when the second parent-derived terms are removed;
6. the state-specific solver has overlap at least 0.99 for an isolated root;
   a degenerate projector agrees to (10^{-8});
7. molecular Hessians agree with finite differences of the analytic OpenQP
   MRSF gradient at three displacement sizes, and the residual translational
   and rotational components decrease with the response tolerances and DFT
   grid density.

No skipped or expected-failure result counts as satisfying a criterion.
