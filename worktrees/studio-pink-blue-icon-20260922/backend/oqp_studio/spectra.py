"""Line-shape broadening for computed spectra.

A calculation gives sticks: a position and an intensity per transition. A
measured spectrum is those sticks convolved with a line shape. Lorentzian is
the default here because homogeneous (lifetime) broadening produces a
Lorentzian, so it reproduces the wings of a real band far better than a
Gaussian does; Gaussian (inhomogeneous broadening) and pseudo-Voigt (a mix of
the two) are offered for comparison.
"""

from __future__ import annotations

import math

import numpy as np

SHAPES = ("lorentzian", "gaussian", "voigt")

# Molar absorptivity prefactor for a line shape normalized in cm^-1
# (Hirata/Casida convention): eps(nu) = 2.174e8 * f * g(nu).
UV_PREFACTOR = 2.174e8


def _lorentzian(x: np.ndarray, center: float, fwhm: float) -> np.ndarray:
    half = fwhm / 2.0
    return (half / math.pi) / ((x - center) ** 2 + half * half)


def _gaussian(x: np.ndarray, center: float, fwhm: float) -> np.ndarray:
    sigma = fwhm / (2.0 * math.sqrt(2.0 * math.log(2.0)))
    return np.exp(-((x - center) ** 2) / (2.0 * sigma * sigma)) / (
        sigma * math.sqrt(2.0 * math.pi)
    )


def profile(x: np.ndarray, center: float, fwhm: float, shape: str,
            eta: float = 0.5) -> np.ndarray:
    """One normalized line shape, in whatever unit x carries."""
    if shape == "gaussian":
        return _gaussian(x, center, fwhm)
    if shape == "voigt":
        # Pseudo-Voigt: the usual cheap stand-in for the true convolution.
        mix = min(max(eta, 0.0), 1.0)
        return mix * _lorentzian(x, center, fwhm) + (1 - mix) * _gaussian(x, center, fwhm)
    return _lorentzian(x, center, fwhm)


def broaden(centers: list[float], weights: list[float], *,
            lo: float, hi: float, fwhm: float, shape: str = "lorentzian",
            points: int = 1200, eta: float = 0.5,
            scale: float = 1.0) -> tuple[list[float], list[float]]:
    """Convolve sticks with `shape` over a regular grid from lo to hi."""
    x = np.linspace(lo, hi, points)
    y = np.zeros_like(x)
    for center, weight in zip(centers, weights):
        if weight == 0.0:
            continue
        y += weight * profile(x, center, fwhm, shape, eta)
    return [float(v) for v in x], [float(v) for v in y * scale]


EV_PER_HARTREE = 27.211386245988
NM_EV = 1239.841984
CM_PER_EV = 8065.543937


def vibrational_spectrum(frequencies: list[float], intensities: list[float], *,
                         fwhm: float = 20.0, shape: str = "lorentzian",
                         eta: float = 0.5, points: int = 1400) -> dict:
    """IR or Raman band shape on a cm^-1 axis."""
    real = [(f, i) for f, i in zip(frequencies, intensities) if f > 0]
    if not real:
        return {"x": [], "y": [], "sticks": []}
    lo = max(0.0, min(f for f, _ in real) - 12 * fwhm)
    hi = max(f for f, _ in real) + 12 * fwhm
    x, y = broaden([f for f, _ in real], [i for _, i in real],
                   lo=lo, hi=hi, fwhm=fwhm, shape=shape, points=points, eta=eta)
    return {
        "x": x, "y": y,
        "sticks": [{"position": f, "intensity": i} for f, i in real],
    }


def electronic_spectrum(energies_ev: list[float], strengths: list[float], *,
                        fwhm_ev: float = 0.3, shape: str = "lorentzian",
                        eta: float = 0.5, points: int = 1400) -> dict:
    """UV/Vis band shape, computed in cm^-1 and reported against eV and nm.

    Broadening a spectrum is only meaningful on a linear-in-energy axis, so
    the convolution runs in wavenumbers and the nm axis is derived from it.
    """
    pairs = [(e, f) for e, f in zip(energies_ev, strengths) if e > 0]
    if not pairs:
        return {"x": [], "y": [], "sticks": [], "x_nm": []}
    fwhm_cm = fwhm_ev * CM_PER_EV
    centers = [e * CM_PER_EV for e, _ in pairs]
    # A Lorentzian formally extends all the way to zero energy.  That tail is
    # harmless on an energy plot, but its wavelength conversion diverges at
    # zero and can make an ESA plot span millions of nanometres.  Retain the
    # relevant low-energy wing while keeping the wavelength axis finite.
    lo = max(min(centers) * 0.5, min(centers) - 8 * fwhm_cm)
    hi = max(centers) + 8 * fwhm_cm
    x_cm, y = broaden(centers, [f for _, f in pairs], lo=lo, hi=hi,
                      fwhm=fwhm_cm, shape=shape, points=points, eta=eta,
                      scale=UV_PREFACTOR)
    return {
        "x": [v / CM_PER_EV for v in x_cm],          # eV
        "x_nm": [NM_EV / (v / CM_PER_EV) for v in x_cm],
        "y": y,                                       # L mol^-1 cm^-1
        "sticks": [
            {"position": e, "position_nm": NM_EV / e, "intensity": f}
            for e, f in pairs
        ],
    }


def energy_spectrum(energies_ev: list[float], strengths: list[float], *,
                    fwhm_ev: float = 0.3, shape: str = "lorentzian",
                    eta: float = 0.5, points: int = 1400) -> dict:
    """Broaden sticks on a direct electron-energy (eV) axis.

    EKT IP/EA roots are electron removal/addition energies, not optical
    transitions.  Keeping this axis in eV avoids the physically misleading
    wavelength conversion used for UV/Vis spectra.
    """
    pairs = [(e, f) for e, f in zip(energies_ev, strengths) if e > 0]
    if not pairs:
        return {"x": [], "y": [], "sticks": []}
    lo = max(0.0, min(e for e, _ in pairs) - 8 * fwhm_ev)
    hi = max(e for e, _ in pairs) + 8 * fwhm_ev
    x, y = broaden([e for e, _ in pairs], [f for _, f in pairs],
                   lo=lo, hi=hi, fwhm=fwhm_ev, shape=shape,
                   points=points, eta=eta)
    return {
        "x": x,
        "y": y,
        "sticks": [{"position": e, "intensity": f} for e, f in pairs],
    }
