"""LAMP-tau: ab initio spin-phonon relaxation analysis for single-ion magnets.

The package is organised in two parts:

``mag_anal``
    Analysis of ab initio results: effective spin Hamiltonian, g tensors,
    spin-phonon coupling and relaxation times.
``SPC_Tool``
    Scripts that generate displaced structures, run the DMET+CASSCF-SO /
    CASSCF calculations and extract the crystal-field parameters.

Importing the package itself is intentionally lightweight; the heavy classes
are resolved lazily so that ``import mag_anal`` does not require PySCF.
"""

__version__ = "1.0.0"

__all__ = ["CAL_RELAX", "PSH", "__version__"]


def __getattr__(name):
    if name == "CAL_RELAX":
        from .Anal_Tau import CAL_RELAX

        return CAL_RELAX
    if name == "PSH":
        from .psh import PSH

        return PSH
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
