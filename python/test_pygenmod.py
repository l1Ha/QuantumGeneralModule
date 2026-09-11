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
    rovibrational_state_index, rovibrational_state_unindex,
    rot_matrix_cos_theta, calc_franck_condon_factors,
    calc_rotational_constants_bv, build_rovibrational_hamiltonian,
    build_rovibrational_dipole_matrix,
    van_der_waals_mean_length, square_well_scattering_length_exact,
    calc_scattering_length_numerov, plot_scattering_length_wavefunction,
    calc_differential_cross_section, calc_differential_cross_section_identical,
    plot_differential_cross_sections,
    calc_multichannel_close_coupling, plot_multichannel_smatrix, plot_feshbach_resonance,
    calc_scattering_wavefunction_ti, plot_scattering_wavefunction,
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

    def test_rovibrational(self):
        # Index roundtrip
        idx = rovibrational_state_index(v=2, j=3, j_max=5)
        v, j = rovibrational_state_unindex(idx, j_max=5)
        self.assertEqual(v, 2)
        self.assertEqual(j, 3)

        # Transition dipole selection rule <0|cos|1> = 1/sqrt(3)
        me = rot_matrix_cos_theta(0, 1, 0)
        self.assertAlmostEqual(me, 1.0 / np.sqrt(3.0), places=6)
        self.assertEqual(rot_matrix_cos_theta(0, 2, 0), 0.0)

        # Franck-Condon factor
        chi = np.zeros((10, 2))
        chi[2, 0] = 1.0
        chi[2, 1] = 1.0
        fc = calc_franck_condon_factors(chi, chi, dx=1.0)
        self.assertAlmostEqual(fc[0, 0], 1.0, places=6)

    def test_scattering(self):
        # 1. Van der Waals mean length
        a_bar = van_der_waals_mean_length(mass=1.0, c6_au=50.0)
        self.assertAlmostEqual(a_bar, 1.51153, places=4)

        # 2. Square well analytical scattering length
        as_exact = square_well_scattering_length_exact(r_well=2.0, v0_well=1.0, mass=1.0)
        self.assertAlmostEqual(as_exact, 2.22898, places=4)

        # 3. Numerov scattering length
        r = np.linspace(0.01, 10.0, 500)
        v = np.where(r <= 2.0, -1.0, 0.0)
        as_num, u = calc_scattering_length_numerov(r, v, mass=1.0)
        self.assertAlmostEqual(as_num, as_exact, delta=0.03)

        # 4. Smoke test for plot_scattering_length_wavefunction
        plot_scattering_length_wavefunction(r, v, u, as_num, filename="test_scat.png")
        self.assertTrue(os.path.exists("test_scat.png"))
        os.remove("test_scat.png")

    def test_differential_cross_sections(self):
        # Phase shifts for l=0, 1, 2
        delta = np.array([0.45, 0.15, 0.05])
        energy = 0.05
        mass = 1.0
        theta = np.linspace(0.0, np.pi, 181)

        # Distinguishable
        ds_dist = calc_differential_cross_section(energy, mass, delta, theta)
        self.assertEqual(len(ds_dist), 181)
        self.assertTrue(np.all(ds_dist >= 0.0))

        # Quantum statistics: Boson vs Fermion at theta = pi/2 (index 90)
        ds_boson = calc_differential_cross_section_identical(energy, mass, delta, theta, 'boson')
        ds_fermion = calc_differential_cross_section_identical(energy, mass, delta, theta, 'fermion')
        ds_unpol = calc_differential_cross_section_identical(energy, mass, delta, theta, 'fermion_unpolarized')

        # At theta = pi/2:
        # 1) Fermion with only odd l has P_odd(0) = 0 -> strictly 0
        self.assertAlmostEqual(ds_fermion[90], 0.0, places=7)

        # 2) Boson has constructive interference: ds_boson(pi/2) = 4 * ds_even(pi/2)
        # Check that ds_boson(pi/2) > 0 and ds_boson >= 0 everywhere
        self.assertGreater(ds_boson[90], 0.0)
        self.assertTrue(np.all(ds_boson >= 0.0))

        # 3) Unpolarized spin-1/2: 0.25 * ds_boson + 0.75 * ds_fermion
        self.assertAlmostEqual(ds_unpol[90], 0.25 * ds_boson[90] + 0.75 * ds_fermion[90], places=7)

        # Smoke test for plotting
        ds_map = {
            "Distinguishable": ds_dist,
            "Bosons": ds_boson,
            "Polarized Fermions": ds_fermion
        }
        plot_differential_cross_sections(theta, ds_map, filename="test_dcs.png")
        self.assertTrue(os.path.exists("test_dcs.png"))
        os.remove("test_dcs.png")

    def test_multichannel_close_coupling(self):
        # 2-channel system: Channel 1 open (thresh=0), Channel 2 closed (thresh=0.4)
        np_pts = 200
        r_grid = np.linspace(0.8, 8.0, np_pts)
        v_mat = np.zeros((2, 2, np_pts))
        for i, r in enumerate(r_grid):
            v_mat[0, 0, i] = -0.5 * np.exp(-(r - 2.0)**2)
            v_mat[1, 1, i] = -1.0 * np.exp(-(r - 2.2)**2)
            v_mat[0, 1, i] = 0.1 * np.exp(-(r - 2.1)**2)
            v_mat[1, 0, i] = v_mat[0, 1, i]

        thresh = np.array([0.0, 0.4])
        # Energy E = 0.15 -> Ch 1 open, Ch 2 closed
        res = calc_multichannel_close_coupling(r_grid, v_mat, mass=1.0, total_energy=0.15, thresholds=thresh)
        self.assertEqual(res['n_open'], 1)
        self.assertEqual(res['n_closed'], 1)
        # Unitarity of open subspace: |S_11|^2 = 1.0
        self.assertAlmostEqual(res['prob_matrix'][0, 0], 1.0, places=5)

        # 2 open channels test: E = 0.5
        res_open2 = calc_multichannel_close_coupling(r_grid, v_mat, mass=1.0, total_energy=0.5, thresholds=thresh)
        self.assertEqual(res_open2['n_open'], 2)
        self.assertEqual(res_open2['n_closed'], 0)
        # Unitarity sum = 1.0
        self.assertAlmostEqual(np.sum(res_open2['prob_matrix'][:, 0]), 1.0, places=4)
        self.assertAlmostEqual(np.sum(res_open2['prob_matrix'][:, 1]), 1.0, places=4)

        # Smoke test for plotting S-matrix heatmap and Feshbach resonance
        plot_multichannel_smatrix(res_open2['prob_matrix'], ["Ch 1", "Ch 2"], filename="test_smat.png")
        self.assertTrue(os.path.exists("test_smat.png"))
        os.remove("test_smat.png")

        plot_feshbach_resonance(np.array([0.05, 0.10, 0.15]), np.array([1.2, 5.8, -2.1]), filename="test_fb.png")
        self.assertTrue(os.path.exists("test_fb.png"))
        os.remove("test_fb.png")

    def test_continuous_scattering_wavefunction(self):
        r_grid = np.linspace(0.01, 25.0, 1000)
        v_pot = -0.8 * np.exp(-(r_grid - 2.0)**2)
        mass = 1.0
        energy = 0.20
        k = np.sqrt(2.0 * mass * energy)
        target_amp = np.sqrt(2.0 * mass / (np.pi * k))

        # Energy normalized
        u_wf_e, delta = calc_scattering_wavefunction_ti(r_grid, v_pot, mass=mass, energy=energy, l=0, norm_type='energy')
        dr = r_grid[1] - r_grid[0]
        du = (u_wf_e[-1] - u_wf_e[-2]) / dr
        amp_num = np.sqrt(u_wf_e[-1]**2 + (du / k)**2)
        self.assertAlmostEqual(amp_num, target_amp, delta=0.03 * target_amp)
        self.assertAlmostEqual(u_wf_e[0], 0.0, places=4)

        # Unit amplitude
        u_wf_1, _ = calc_scattering_wavefunction_ti(r_grid, v_pot, mass=mass, energy=energy, l=0, norm_type='unit')
        amp_1 = np.sqrt(u_wf_1[-1]**2 + (((u_wf_1[-1] - u_wf_1[-2]) / dr) / k)**2)
        self.assertAlmostEqual(amp_1, 1.0, delta=0.03)

        # Plot test
        plot_scattering_wavefunction(r_grid, v_pot, u_wf_e, energy=energy, phase_shift=delta, filename="test_scat_wf.png")
        self.assertTrue(os.path.exists("test_scat_wf.png"))
        os.remove("test_scat_wf.png")

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
