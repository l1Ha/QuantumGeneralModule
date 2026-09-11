#!/usr/bin/env python3
"""
Unit tests for pygenmod Python companion package.
"""
import sys
import os
import unittest
import numpy as np

# Add parent directory to sys.path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from pygenmod import (
    to_au, from_au, AU2EV, FS2AU,
    PulseConfig, PULSE_GAUSSIAN, PULSE_SIN2,
    pulse_envelope, pulse_electric_field, pulse_stark_shift,
    dvr_sinc_init, fgh_solve_bound_states,
    get_atom_config, keldysh_parameter, hhg_cutoff_energy,
    soft_core_coulomb_potential, adk_ionization_rate,
    landau_zener_probability, calculate_channel_populations,
    plot_wavefunctions, plot_pulses
)


class TestPyGenMod(unittest.TestCase):

    def test_constants_and_units(self):
        # Roundtrip fs
        val = 150.0
        val_au = to_au(val, 'fs')
        val_back = from_au(val_au, 'fs')
        self.assertAlmostEqual(val, val_back, places=9)

        # Roundtrip eV
        val_ev = 13.605693
        val_au = to_au(val_ev, 'eV')
        val_back = from_au(val_au, 'eV')
        self.assertAlmostEqual(val_ev, val_back, places=9)

    def test_pulse_generation(self):
        cfg = PulseConfig(
            shape_type=PULSE_GAUSSIAN,
            field_peak=0.05,
            freq_central=0.057,
            duration=50.0 * FS2AU,
            t_center=0.0
        )
        t = np.array([0.0, 25.0 * FS2AU, -25.0 * FS2AU])
        env = pulse_envelope(t, cfg)
        self.assertAlmostEqual(env[0], 1.0, places=7)
        self.assertAlmostEqual(env[1], 0.5, places=7)
        self.assertAlmostEqual(env[2], 0.5, places=7)

        # Electric field at t=0
        ef = pulse_electric_field(np.array([0.0]), cfg)
        self.assertAlmostEqual(ef[0], 0.05, places=7)

    def test_dvr_fgh_harmonic_oscillator(self):
        # HO: V(x) = 0.5 * m * omega^2 * x^2, m=1, omega=1 -> E_v = v + 0.5
        dvr = dvr_sinc_init(x_min=-8.0, x_max=8.0, n_points=128, mass=1.0)
        v_pot = 0.5 * (dvr.x**2)
        eig_vals, wavefuncs = fgh_solve_bound_states(dvr, v_pot)

        # Ground state ~ 0.5
        self.assertAlmostEqual(eig_vals[0], 0.5, places=3)
        # First excited state ~ 1.5
        self.assertAlmostEqual(eig_vals[1], 1.5, places=3)
        # Second excited state ~ 2.5
        self.assertAlmostEqual(eig_vals[2], 2.5, places=3)

        # Normalization check
        norm_v0 = np.sum(np.abs(wavefuncs[:, 0])**2) * dvr.dx
        self.assertAlmostEqual(norm_v0, 1.0, places=5)

    def test_coulomb_atomic(self):
        atom = get_atom_config("Ar")
        self.assertAlmostEqual(atom.ip_au, 0.57915, places=4)

        # Keldysh gamma
        gamma = keldysh_parameter(0.057, 0.05, 0.5)
        self.assertAlmostEqual(gamma, 1.14, places=2)

        # Cutoff law
        e_cut = hhg_cutoff_energy(0.5, 0.05, 0.057)
        self.assertGreater(e_cut, 1.10)

        # Soft-core potential
        v = soft_core_coulomb_potential(np.array([0.0]), soft_a=1.0)
        self.assertAlmostEqual(v[0], -1.0, places=7)

        # ADK rate positive
        w_adk = adk_ionization_rate(0.05, 0.5)
        self.assertGreater(w_adk, 0.0)

    def test_hhg_and_multistate(self):
        # Landau Zener
        p_lz = landau_zener_probability(0.1, 1.0, 1.0)
        self.assertAlmostEqual(p_lz, np.exp(-2.0 * np.pi * 0.01), places=6)

        # Populations
        psi1 = np.ones(10) / np.sqrt(10.0)
        psi2 = np.zeros(10)
        pop1, pop2, ratio = calculate_channel_populations(psi1, psi2, dx=1.0)
        self.assertAlmostEqual(pop1, 1.0, places=6)
        self.assertAlmostEqual(pop2, 0.0, places=6)
        self.assertAlmostEqual(ratio, 0.0, places=6)

    def test_visualizer_smoke(self):
        # Quick check that visualization runs without error
        dvr = dvr_sinc_init(x_min=-5.0, x_max=5.0, n_points=64, mass=1.0)
        v_pot = 0.5 * (dvr.x**2)
        eig_vals, wavefuncs = fgh_solve_bound_states(dvr, v_pot)

        plot_wavefunctions(dvr.x, v_pot, eig_vals, wavefuncs, n_states=3,
                           title="HO Test Wavefunctions", filename="test_ho.png")
        self.assertTrue(os.path.exists("test_ho.png"))
        os.remove("test_ho.png")


if __name__ == '__main__':
    unittest.main(verbosity=2)
