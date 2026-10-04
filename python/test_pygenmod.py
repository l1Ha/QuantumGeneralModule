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
    to_au, from_au, FS2AU,
    PulseConfig, PULSE_GAUSSIAN,
    pulse_envelope, pulse_electric_field,
    dvr_sinc_init, fgh_solve_bound_states,
    get_atom_config, keldysh_parameter, hhg_cutoff_energy,
    soft_core_coulomb_potential, adk_ionization_rate,
    landau_zener_probability, calculate_channel_populations,
    rovibrational_state_index, rovibrational_state_unindex,
    rot_matrix_cos_theta, calc_franck_condon_factors,
    van_der_waals_mean_length, square_well_scattering_length_exact,
    calc_scattering_length_numerov, plot_scattering_length_wavefunction,
    calc_differential_cross_section, calc_differential_cross_section_identical,
    plot_differential_cross_sections,
    calc_multichannel_close_coupling, plot_multichannel_smatrix, plot_feshbach_resonance,
    calc_scattering_wavefunction_ti, plot_scattering_wavefunction,
    plot_wavefunctions,
    BASIS_UNCOUPLED, BASIS_F_COUPLED, BASIS_TOTAL_SPIN, BASIS_FIELD_DRESSED,
    AU2GHZ, clebsch_gordan_half, get_cold_atom_preset,
    calc_breit_rabi_energies, build_field_collision_channels,
    calc_basis_transform_matrix,
    build_spin_exchange_matrix, fit_feshbach_resonance_parameters,
    plot_breit_rabi_diagram, plot_magnetic_feshbach_resonance
)



class TestPublicAPI(unittest.TestCase):
    """包根 __all__ 声明的每个名字都必须真实可导入（公开 API 表面）。"""

    def test_all_names_importable(self):
        import pygenmod
        missing = [name for name in pygenmod.__all__ if not hasattr(pygenmod, name)]
        self.assertEqual(missing, [])
        # 第 17/18 章参考实现的关鍵入口必须出现在 __all__ 中
        for name in ("run_qct_ensemble", "SOPHamiltonian", "tt_round",
                     "smolyak_build", "mctdh_propagate"):
            self.assertIn(name, pygenmod.__all__)


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

    def test_field_scattering_module(self):
        # 1. Clebsch-Gordan half integer
        cg1 = clebsch_gordan_half(1, 1, 1, -1, 2, 0)
        self.assertAlmostEqual(cg1, 1.0 / np.sqrt(2.0), places=12)

        # 2. Breit-Rabi Rb87 zero field hyperfine splitting
        rb87 = get_cold_atom_preset("87Rb")
        evals, _ = calc_breit_rabi_energies(rb87, 0.0)
        hfs_split_ghz = (evals[3] - evals[0]) * AU2GHZ
        self.assertAlmostEqual(hfs_split_ghz, 6.83468261, delta=1e-4)

        # 3. Channels and Unitary transformation across 4 bases
        ch_unc = build_field_collision_channels(rb87, rb87, BASIS_UNCOUPLED, two_Mtot=2, l_max=0)
        ch_f = build_field_collision_channels(rb87, rb87, BASIS_F_COUPLED, two_Mtot=2, l_max=0)
        ch_spin = build_field_collision_channels(rb87, rb87, BASIS_TOTAL_SPIN, two_Mtot=2, l_max=0)
        n_ch = len(ch_unc)
        self.assertGreater(n_ch, 0)

        # U(f <- unc)
        u_f_unc = calc_basis_transform_matrix(rb87, rb87, ch_unc, ch_f, BASIS_UNCOUPLED, BASIS_F_COUPLED, 50.0)
        np.testing.assert_allclose(u_f_unc @ u_f_unc.T, np.eye(n_ch), atol=1e-12)

        # U(spin <- unc)
        u_spin_unc = calc_basis_transform_matrix(rb87, rb87, ch_unc, ch_spin, BASIS_UNCOUPLED, BASIS_TOTAL_SPIN, 50.0)
        np.testing.assert_allclose(u_spin_unc @ u_spin_unc.T, np.eye(n_ch), atol=1e-12)

        # U(dress <- unc)
        u_dress_unc = calc_basis_transform_matrix(rb87, rb87, ch_unc, ch_unc, BASIS_UNCOUPLED, BASIS_FIELD_DRESSED, 50.0)
        np.testing.assert_allclose(u_dress_unc @ u_dress_unc.T, np.eye(n_ch), atol=1e-12)

        # Chain rule: U(f <- spin) == U(f <- unc) @ U(unc <- spin)
        u_f_spin = calc_basis_transform_matrix(rb87, rb87, ch_spin, ch_f, BASIS_TOTAL_SPIN, BASIS_F_COUPLED, 50.0)
        np.testing.assert_allclose(u_f_spin, u_f_unc @ u_spin_unc.T, atol=1e-12)

        # 4. Spin exchange operator s1 . s2 is diagonal in Total Spin basis
        p_exc = build_spin_exchange_matrix(ch_spin, BASIS_TOTAL_SPIN)
        np.testing.assert_allclose(p_exc, np.diag(np.diag(p_exc)), atol=1e-14)

        # 5. Feshbach fit
        # 避开 b0 = 80 G 共振极点：解析式在 B = b0 处发散，采样网格需错开该点。
        b_test = np.linspace(50.0, 110.0, 61)
        b_test = b_test[np.abs(b_test - 80.0) > 1.0e-9]
        a_test = 100.0 * (1.0 - 5.0 / (b_test - 80.0))
        b0_fit, delta_b_fit, a_bg_fit = fit_feshbach_resonance_parameters(b_test, a_test)
        self.assertAlmostEqual(b0_fit, 80.0, delta=1.0)
        self.assertAlmostEqual(delta_b_fit, 5.0, delta=1.0)
        self.assertAlmostEqual(a_bg_fit, 100.0, delta=5.0)

        # 6. Plotting smoke test
        plot_breit_rabi_diagram(rb87, b_max_gauss=50.0, num_b=30, save_path="test_br.png")
        self.assertTrue(os.path.exists("test_br.png"))
        os.remove("test_br.png")

        plot_magnetic_feshbach_resonance(b_test, a_test, b0_fit, delta_b_fit, a_bg_fit, save_path="test_mfb.png")
        self.assertTrue(os.path.exists("test_mfb.png"))
        os.remove("test_mfb.png")

        # 7. Heteronuclear 40K + 87Rb multi-basis test
        k40 = get_cold_atom_preset("40K")
        ch_u_het = build_field_collision_channels(k40, rb87, BASIS_UNCOUPLED, two_Mtot=-7, l_max=0)
        ch_f_het = build_field_collision_channels(k40, rb87, BASIS_F_COUPLED, two_Mtot=-7, l_max=0)
        ch_s_het = build_field_collision_channels(k40, rb87, BASIS_TOTAL_SPIN, two_Mtot=-7, l_max=0)
        self.assertEqual(len(ch_u_het), 12)
        self.assertEqual(len(ch_f_het), 12)
        self.assertEqual(len(ch_s_het), 12)

        u_f_het = calc_basis_transform_matrix(k40, rb87, ch_u_het, ch_f_het, BASIS_UNCOUPLED, BASIS_F_COUPLED, 540.0)
        u_s_het = calc_basis_transform_matrix(k40, rb87, ch_u_het, ch_s_het, BASIS_UNCOUPLED, BASIS_TOTAL_SPIN, 540.0)
        np.testing.assert_allclose(u_f_het @ u_f_het.T, np.eye(12), atol=1e-12)
        np.testing.assert_allclose(u_s_het @ u_s_het.T, np.eye(12), atol=1e-12)


# ---------------------------------------------------------------------------
# Chapter 17: QCT dynamics (pygenmod.qct)
# ---------------------------------------------------------------------------

class TestQCT(unittest.TestCase):
    """QCT 轨迹引擎：守恒律、统计一致性、EBK 指认（书中第 17 章）。"""

    def test_ebk_action_quantization(self):
        # EBK 对 Morse 振子是精确的：j=0 时应逐位复现 G(v)
        import pygenmod.qct as Q
        par = Q.LEPSParameters()
        for v in (0, 1, 2, 3):
            e_ebk = Q.ebk_internal_energy(v, 0, par)
            g = Q.morse_vibrational_energy(v, par)
            self.assertAlmostEqual(e_ebk + par.d_e, g, delta=1e-9 * par.d_e)

    def test_trajectory_conservation(self):
        # Velocity-Verlet: 总能量漂移 < 1e-4 Ha，总角动量守恒到 1e-10
        import pygenmod.qct as Q
        from pygenmod.constants import EV2AU
        par = Q.LEPSParameters()
        cfg = Q.QCTConfig(e_coll=1.5 * EV2AU, b_max=3.0, r_start=9.0, r_end=11.0,
                          dt=3.0, max_steps=8000, v_initial=0, j_initial=0, seed=7)
        rng = np.random.default_rng(cfg.seed)
        sampler = Q.InitialSampler.build(cfg, par)
        for _ in range(5):
            tr = Q.qct_init_trajectory(cfg, rng, par, sampler)
            e0, j0 = tr.e_total, np.linalg.norm(tr.j_total)
            Q.propagate_trajectory(tr, cfg, par)
            e1 = Q._total_energy(tr.q, tr.p, cfg.mass, par)
            j1 = np.linalg.norm(Q._total_angular_momentum(tr.q, tr.p, cfg.mass))
            self.assertLess(abs(e1 - e0), 1.0e-4)
            self.assertLess(abs(j1 - j0), 1.0e-10 * max(j0, 1.0))
            self.assertTrue(tr.converged)

    def test_internal_energy_bound(self):
        # 严格不等式：|dE_int| <= E_coll（能量只在内部与平动之间流动）
        import pygenmod.qct as Q
        from pygenmod.constants import EV2AU
        par = Q.LEPSParameters()
        cfg = Q.QCTConfig(e_coll=1.5 * EV2AU, b_max=3.0, r_start=9.0, r_end=11.0,
                          dt=3.0, n_traj=80, max_steps=8000, seed=42)
        res = Q.run_qct_ensemble(cfg, par)
        e0 = -par.d_e + Q.morse_vibrational_energy(0, par)
        for t in res.trajectories:
            if not t.converged:
                continue
            pair, atom = Q._product_jacobi(t.q)
            e_int, _, _, _ = Q._product_internal_state(
                t.q, t.p, pair, atom, cfg.mass, par, t.r_dir_initial)
            self.assertLessEqual(abs(e_int - e0), cfg.e_coll + 1.0e-6)

    def test_ensemble_statistics(self):
        # MC 估计量、分层估计量与 Wilson 区间的一致性；种子可复现
        import pygenmod.qct as Q
        from pygenmod.constants import EV2AU
        par = Q.LEPSParameters()
        cfg = Q.QCTConfig(e_coll=1.5 * EV2AU, b_max=3.0, r_start=9.0, r_end=11.0,
                          dt=3.0, n_traj=200, max_steps=8000, seed=42)
        res = Q.run_qct_ensemble(cfg, par)
        sig, err = Q.qct_cross_section(res.n_reactive, res.n_traj, cfg.b_max)
        self.assertAlmostEqual(sig, res.cross_section, delta=1e-12)
        self.assertAlmostEqual(err, res.stat_error, delta=1e-12)
        self.assertLessEqual(res.wilson_lo, res.cross_section + 1e-12)
        self.assertGreaterEqual(res.wilson_hi, res.cross_section - 1e-12)
        strat = Q.stratified_cross_section(res.trajectories, cfg)
        # 分层估计量应在 3 倍统计误差内与 MC 估计量一致
        self.assertLess(abs(strat - res.cross_section),
                        3.0 * res.stat_error + 1.0e-9)
        # 产物态截面之和不超过总截面
        s_sum = sum(res.state_cross_section.values())
        self.assertLessEqual(s_sum, res.cross_section + 1e-9)
        # 确定性
        res2 = Q.run_qct_ensemble(cfg, par)
        self.assertEqual(res.n_reactive, res2.n_reactive)

    def test_zpe_constraint_and_thermal_rate(self):
        # 被动 ZPE 约束剔除不可指认轨迹；MB 积分用解析 sigma 验证
        import pygenmod.qct as Q
        from pygenmod.constants import EV2AU
        par = Q.LEPSParameters()
        cfg = Q.QCTConfig(e_coll=1.5 * EV2AU, b_max=3.0, r_start=9.0, r_end=11.0,
                          dt=3.0, n_traj=100, max_steps=8000, seed=42,
                          enforce_zpe=True)
        res = Q.run_qct_ensemble(cfg, par)
        self.assertTrue(all(v >= 0 for (v, _) in res.state_cross_section))
        # 解析检验：sigma(E) = 常数时 k(T) = pref * sigma * T^2（网格须覆盖 E ~ kT 峰）
        mu = 0.5 * 1837.15
        temp = 500.0 * 3.1668e-6   # 500 K in Hartree
        energies = np.linspace(0.05 * temp, 40.0 * temp, 4000)
        sigma0 = 1.3
        k = Q.qct_thermal_rate(energies, np.full_like(energies, sigma0), temp, mu)
        pref = np.sqrt(8.0 / (np.pi * mu * temp ** 3))
        self.assertAlmostEqual(k / (pref * sigma0 * temp ** 2), 1.0, delta=5e-3)

    def test_parallel_ensemble_execution(self):
        # 多进程并行采样（针对服务器多核硬件加速）：轨迹总数严格守恒
        import pygenmod.qct as Q
        from pygenmod.constants import EV2AU
        par = Q.LEPSParameters()
        cfg_seq = Q.QCTConfig(e_coll=1.5 * EV2AU, b_max=3.0, n_traj=160,
                              max_steps=6000, seed=42, n_workers=1)
        res_seq = Q.run_qct_ensemble(cfg_seq, par)
        cfg_par = Q.QCTConfig(e_coll=1.5 * EV2AU, b_max=3.0, n_traj=160,
                              max_steps=6000, seed=42, n_workers=2)
        res_par = Q.run_qct_ensemble(cfg_par, par)
        self.assertEqual(res_par.n_traj, 160)
        self.assertLess(abs(res_par.cross_section - res_seq.cross_section),
                        3.0 * res_seq.stat_error + 0.5)


# ---------------------------------------------------------------------------
# Chapter 18: high-dimensional methods (SOP / TT / Smolyak / MCTDH)
# ---------------------------------------------------------------------------

def _random_sop(rng, d, n, m, hermitian=False):
    import pygenmod.sop_hamiltonian as S
    terms = []
    for _ in range(m):
        facs = []
        for _ in range(d):
            a = rng.normal(size=(n, n))
            if hermitian:
                a = a + a.T
            facs.append(a)
        terms.append(S.SOPTerm(coeff=float(rng.normal()), factors=facs))
    return S.SOPHamiltonian(dims=(n,) * d, terms=terms)


def _sop_dense_kron(sop):
    n_tot = int(np.prod(sop.dims))
    h = np.zeros((n_tot, n_tot), dtype=complex)
    for t in sop.terms:
        op = np.array([[t.coeff]])
        for f in t.factors:
            op = np.kron(op, f)
        h += op
    return h


class TestSOPHamiltonian(unittest.TestCase):
    """SOP 矩阵自由作用与 POTFIT 分解（书中 18.1 节）。"""

    def test_apply_matches_dense(self):
        rng = np.random.default_rng(0)
        sop = _random_sop(rng, d=3, n=4, m=5)
        v = rng.normal(size=(4, 4, 4))
        hv = sop.apply(v)
        href = (_sop_dense_kron(sop) @ v.ravel()).reshape(4, 4, 4)
        self.assertLess(np.abs(hv - href).max() / np.abs(href).max(), 1e-12)

    def test_potfit_reconstruction(self):
        import pygenmod.sop_hamiltonian as S
        rng = np.random.default_rng(1)
        vg = rng.normal(size=(4, 5, 6))
        sop = S.sop_from_potfit(vg, eps=1e-10)
        # 重构张量 = 所有项的外积和
        out = np.zeros(vg.shape, dtype=complex)
        grids = np.meshgrid(*[np.arange(n) for n in vg.shape], indexing="ij")
        for term in sop.terms:
            acc = np.ones(vg.shape, dtype=complex)
            for k, f in enumerate(term.factors):
                acc = acc * np.asarray(f, dtype=complex)[grids[k]]
            out += term.coeff * acc
        self.assertLess(np.abs(out - vg).max() / np.abs(vg).max(), 1e-8)


class TestTensorTrain(unittest.TestCase):
    """TT-SVD、舍入、内积与 TT 算符（书中 18.3 节）。"""

    def test_roundtrip_and_dot(self):
        import pygenmod.tensor_train as T
        rng = np.random.default_rng(2)
        a = rng.normal(size=(4, 5, 6, 3))
        tt = T.tt_from_dense(a, eps=1e-12)
        self.assertLess(np.abs(T.tt_to_dense(tt) - a).max() / np.abs(a).max(), 1e-10)
        ttr = T.tt_round(tt, eps=1e-10)
        self.assertLess(np.abs(T.tt_to_dense(ttr) - a).max() / np.abs(a).max(), 1e-8)
        b = rng.normal(size=a.shape)
        dot_tt = T.tt_dot(T.tt_from_dense(a, eps=1e-14), T.tt_from_dense(b, eps=1e-14))
        self.assertLess(abs(dot_tt - np.sum(a * b.conj())), 1e-10)

    def test_operator_from_sop(self):
        import pygenmod.tensor_train as T
        rng = np.random.default_rng(3)
        sop = _random_sop(rng, d=3, n=4, m=5)
        cores = T.sop_to_tt_operator(sop)
        hd = _sop_dense_kron(sop)
        self.assertLess(np.abs(T.tt_operator_to_dense(cores) - hd).max()
                        / np.abs(hd).max(), 1e-12)
        v = rng.normal(size=(4, 4, 4))
        hv = sop.apply(v)
        w = T.tt_operator_apply(cores, T.tt_from_dense(v, eps=1e-14), eps=1e-12)
        self.assertLess(np.abs(T.tt_to_dense(w) - hv).max() / np.abs(hv).max(), 1e-10)


class TestSparseGrid(unittest.TestCase):
    """Smolyak 稀疏网格（书中 18.4 节）。"""

    def test_coefficient_equals_difference_form(self):
        import pygenmod.sparse_grid as G
        for (d, n) in ((2, 5), (3, 6)):
            g1 = G.smolyak_build(d, n)
            g2 = G._difference_form_grid(d, n)
            self.assertEqual(len(g1.points), len(g2.points))
            d1 = {tuple(np.round(p, 10)): w for p, w in zip(g1.points, g1.weights)}
            d2 = {tuple(np.round(p, 10)): w for p, w in zip(g2.points, g2.weights)}
            self.assertEqual(set(d1), set(d2))
            for key in d1:
                self.assertAlmostEqual(d1[key], d2[key], delta=1e-10)
            # 常数必须精确：权重和 = 2^D
            self.assertAlmostEqual(g1.weights.sum(), 2.0 ** d, delta=1e-12)

    def test_integration_convergence(self):
        import pygenmod.sparse_grid as G
        cs = np.array([2.0, 0.7])
        exact = np.prod([2.0 * np.sin(c) / c for c in cs])
        errs = []
        for n_level in (4, 6, 8):
            g = G.smolyak_build(2, n_level)
            vals = np.prod(np.cos(g.points * cs), axis=-1)
            errs.append(abs(g.integrate(vals) - exact))
        self.assertLess(errs[-1], 1e-8)
        self.assertLess(errs[-1], errs[0])   # 收敛


class TestMCTDHCore(unittest.TestCase):
    """MCTDH 核心：全空间精确性、规范与守恒、SPF 收敛（书中 18.2 节）。"""

    @classmethod
    def setUpClass(cls):
        import pygenmod.sop_hamiltonian as S
        rng = np.random.default_rng(1)
        d, n = 3, 6
        def herm(nn):
            a = rng.normal(size=(nn, nn)) + 1j * rng.normal(size=(nn, nn))
            return a + a.conj().T
        h1, h2, h3 = herm(n), herm(n), herm(n)
        cls.sop = S.SOPHamiltonian(dims=(n,) * d, terms=[
            S.SOPTerm(1.0, [h1, np.eye(n), np.eye(n)]),
            S.SOPTerm(1.0, [np.eye(n), h2, np.eye(n)]),
            S.SOPTerm(1.0, [np.eye(n), np.eye(n), h3]),
            S.SOPTerm(0.3, [h1, h2, np.eye(n)]),
            S.SOPTerm(0.3, [np.eye(n), h2, h3]),
        ])
        cls.hd = 0.5 * (cls.sop.to_dense() + cls.sop.to_dense().conj().T)
        cls.w_eig, cls.U = np.linalg.eigh(cls.hd)
        cls.a0 = rng.normal(size=(n,) * 3) + 1j * rng.normal(size=(n,) * 3)
        cls.a0 /= np.linalg.norm(cls.a0)
        cls.r = np.random.default_rng(42)
        cls.qs = []
        for k in range(3):
            q, _ = np.linalg.qr(cls.r.normal(size=(n, n)) + 1j * cls.r.normal(size=(n, n)))
            cls.qs.append(q)

    @staticmethod
    def _project(a, spfs):
        # primitive -> coefficients: A[i] = <phi_i|psi>，需要 conj(spf)
        b = a
        for k in range(len(spfs)):
            b = np.moveaxis(np.tensordot(spfs[k].conj(), b, axes=([1], [k])), 0, k)
        return b

    def _psi(self, st):
        psi = st.A
        for k in range(len(st.spf)):
            psi = np.moveaxis(np.tensordot(st.spf[k], psi, axes=([0], [k])), 0, k)
        return psi

    def _energy(self, st):
        # <Psi|H|Psi> = <A|H_A A>，H_A 通过 SPF 基矩阵 eta 作用
        import pygenmod.mctdh_core as M
        d = st.A.ndim
        e = 0.0 + 0.0j
        for term in self.sop.terms:
            c = st.A
            for k in range(d):
                eta = M._spf_apply(st.spf[k], term.factors[k])
                c = M._apply_along_axis(eta, c, k)
            e += term.coeff * np.sum(st.A.conj() * c)
        return float(np.real(e))

    def test_full_space_exact(self):
        # 单位基与随机旋转基都须复现精确演化（rotated 基检验 SPF 方程）；
        # 精确参考必须从同一初态出发（rotated 基下的投影态）。
        import pygenmod.mctdh_core as M
        t_end, dt = 0.5, 0.002
        n_steps = int(t_end / dt)
        aex = self.U.conj().T @ self.a0.reshape(-1)
        aex = np.einsum("ij,j->i", self.U,
                        np.exp(-1j * self.w_eig * t_end) * aex).reshape(self.a0.shape)
        st = M.MCTDHState(A=self.a0.copy(),
                          spf=[np.eye(self.a0.shape[0], dtype=complex)] * 3)
        M.mctdh_propagate(st, self.sop, M.MCTDHConfig(dt=dt, n_steps=n_steps))
        fid_id = abs(np.sum(aex.conj() * st.A)) ** 2 / np.linalg.norm(st.A) ** 2
        spfs = [self.qs[k].T.copy() for k in range(3)]
        a_rot = self._project(self.a0, spfs)
        st2 = M.MCTDHState(A=a_rot.copy(), spf=[s.copy() for s in spfs])
        psi0 = self._psi(st2)
        arex = self.U.conj().T @ psi0.reshape(-1)
        arex = np.einsum("ij,j->i", self.U,
                         np.exp(-1j * self.w_eig * t_end) * arex).reshape(self.a0.shape)
        M.mctdh_propagate(st2, self.sop, M.MCTDHConfig(dt=dt, n_steps=n_steps))
        psi = self._psi(st2)
        fid_rot = abs(np.sum(arex.conj() * psi)) ** 2 / (np.linalg.norm(arex) ** 2 *
                                                         np.linalg.norm(psi) ** 2)
        self.assertGreater(fid_id, 1 - 1e-7)
        self.assertGreater(fid_rot, 1 - 1e-6)

    def test_reduced_conservation_and_convergence(self):
        # 归约流形上：范数、正交归一与能量守恒；保真度随 SPF 数单调提高
        import pygenmod.mctdh_core as M
        t_end, dt = 0.5, 0.002
        n_steps = int(t_end / dt)
        fids = []
        for n_spf in (3, 4):
            spfs = [self.qs[k][:, :n_spf].T.copy() for k in range(3)]
            a_red = self._project(self.a0, spfs)
            a_red /= np.linalg.norm(a_red)
            # 匹配的精确参考：从同一投影初态出发
            psi0 = self._psi(M.MCTDHState(A=a_red, spf=spfs))
            c0 = self.U.conj().T @ psi0.reshape(-1)
            arex = np.einsum("ij,j->i", self.U,
                             np.exp(-1j * self.w_eig * t_end) * c0).reshape(self.a0.shape)
            st = M.MCTDHState(A=a_red.copy(), spf=[s.copy() for s in spfs])
            e0 = self._energy(st)
            for _ in range(n_steps // 2):
                M.mctdh_propagate_step(st, self.sop, M.MCTDHConfig(dt=dt, n_steps=1))
            e_mid = self._energy(st)
            for _ in range(n_steps // 2):
                M.mctdh_propagate_step(st, self.sop, M.MCTDHConfig(dt=dt, n_steps=1))
            psi = self._psi(st)
            e1 = self._energy(st)
            orth = max(np.abs(s @ s.conj().T - np.eye(s.shape[0])).max() for s in st.spf)
            self.assertLess(orth, 1e-10)
            self.assertLess(abs(st.norm() - 1.0), 1e-8)
            self.assertLess(abs(e_mid - e0) / max(abs(e0), 1.0), 5e-3)
            self.assertLess(abs(e1 - e0) / max(abs(e0), 1.0), 5e-3)
            fid = abs(np.sum(arex.conj() * psi)) ** 2 / (np.linalg.norm(arex) ** 2 *
                                                         np.linalg.norm(psi) ** 2)
            fids.append(fid)
        self.assertGreater(fids[1], fids[0])

    def test_separable_product_tracking(self):
        # 可分离哈密顿量 + 乘积初态：MCTDH 应高精度跟随精确解
        import pygenmod.sop_hamiltonian as S
        import pygenmod.mctdh_core as M
        rng = np.random.default_rng(5)
        nn = 6
        def herm(nn):
            a = rng.normal(size=(nn, nn)) + 1j * rng.normal(size=(nn, nn))
            return a + a.conj().T
        g1, g2 = herm(nn), herm(nn)
        sop = S.SOPHamiltonian(dims=(nn, nn), terms=[
            S.SOPTerm(1.0, [g1, np.eye(nn)]), S.SOPTerm(1.0, [np.eye(nn), g2])])
        w1, u1 = np.linalg.eigh(g1)
        w2, u2 = np.linalg.eigh(g2)
        t_end, dt = 1.0, 0.001
        v1 = u1 @ (np.exp(-1j * w1 * t_end) * (u1.conj().T @ np.eye(nn)[:, 0]))
        v2 = u2 @ (np.exp(-1j * w2 * t_end) * (u2.conj().T @ np.eye(nn)[:, 0]))
        psi_ex = np.outer(v1, v2)
        spfs = []
        rr = np.random.default_rng(50)
        for _ in range(2):
            vecs = np.column_stack([np.eye(nn)[:, 0],
                                    rr.normal(size=nn) + 1j * rr.normal(size=nn)])
            q, _ = np.linalg.qr(vecs)
            spfs.append(q.T.copy())
        a = np.zeros((2, 2), dtype=complex)
        a[0, 0] = 1.0
        st = M.MCTDHState(A=a, spf=spfs)
        M.mctdh_propagate(st, sop, M.MCTDHConfig(dt=dt, n_steps=int(t_end / dt)))
        psi = self._psi(st)
        fid = abs(np.sum(psi_ex.conj() * psi)) ** 2 / np.linalg.norm(psi) ** 2
        self.assertGreater(fid, 0.999)



if __name__ == '__main__':
    unittest.main(verbosity=2)
