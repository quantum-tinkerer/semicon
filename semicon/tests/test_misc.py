import numpy as np

from semicon.kp_models.symbols import sigma_x, sigma_y, sigma_z, Jx, Jy, Jz
from semicon.misc import spin_matrices


sigma_x = np.array(sigma_x.tolist(), dtype=complex)
sigma_y = np.array(sigma_y.tolist(), dtype=complex)
sigma_z = np.array(sigma_z.tolist(), dtype=complex)

Jx = np.array(Jx.tolist(), dtype=complex)
Jy = np.array(Jy.tolist(), dtype=complex)
Jz = np.array(Jz.tolist(), dtype=complex)


def test_spin_matrices_explicit():
    S = spin_matrices(1/2)
    assert np.allclose(S, 0.5 * np.array([sigma_x, sigma_y, sigma_z]))

    S = spin_matrices(3/2)
    assert np.allclose(S, np.array([Jx, Jy, Jz]))


def test_spin_matrices_commutation():
    max_s = 10
    for s in np.arange(1/2, max_s + 1/2, 1/2):
        Sx, Sy, Sz = spin_matrices(1/2)

        np.allclose((Sx @ Sy - Sy @ Sx), 1j * Sz)
        np.allclose((Sy @ Sz - Sz @ Sy), 1j * Sx)
        np.allclose((Sz @ Sx - Sx @ Sz), 1j * Sy)
