"""
2D Bicubic Tensor-Product Spline (PES interpolation) — Python reference.

Mirrors `mod_interpolation::spline_2d_*` (src/mod_interpolation.f90):
natural-boundary 1D cubic splines along both axes build the node second
derivatives zxx = d2z/dx2, zyy = d2z/dy2 and the mixed derivative
zxxyy = d4z/(dx2 dy2); inside a cell the Hermite-type basis (T linear hats,
U cubic spline bases) reproduces node values and second derivatives exactly
and is C1 continuous.
"""
import numpy as np


def _spline_1d_y2(x, y):
    """Natural cubic spline second derivatives at the nodes (tridiagonal solve)."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    n = x.size
    if n < 3:
        return np.zeros(n)
    h = np.diff(x)
    # natural boundary: y2[0] = y2[n-1] = 0
    a = np.zeros(n)
    b = np.zeros(n)
    c = np.zeros(n)
    d = np.zeros(n)
    b[1:-1] = 2.0 * (h[:-1] + h[1:])
    a[1:-1] = h[:-1]
    c[1:-1] = h[1:]
    d[1:-1] = 6.0 * ((y[2:] - y[1:-1]) / h[1:] - (y[1:-1] - y[:-2]) / h[:-1])
    # Thomas algorithm (b[0] and b[n-1] unused with natural BCs)
    cp = np.zeros(n)
    dp = np.zeros(n)
    for i in range(1, n - 1):
        m = b[i] - a[i] * cp[i - 1]
        cp[i] = c[i] / m
        dp[i] = (d[i] - a[i] * dp[i - 1]) / m
    y2 = np.zeros(n)
    for i in range(n - 2, 0, -1):
        y2[i] = dp[i] - cp[i] * y2[i + 1]
    return y2


class Spline2D:
    """Bicubic tensor-product spline on a regular grid (natural boundaries)."""

    def __init__(self, x, y, z):
        self.x = np.asarray(x, dtype=float)
        self.y = np.asarray(y, dtype=float)
        z = np.asarray(z, dtype=float)
        nx, ny = self.x.size, self.y.size
        if z.shape != (nx, ny):
            raise ValueError("z must have shape (len(x), len(y))")
        self.z = z
        # (1) per row (fixed x): spline along y -> zyy
        self.zyy = np.zeros((nx, ny))
        for i in range(nx):
            self.zyy[i, :] = _spline_1d_y2(self.y, z[i, :])
        # (2) per column (fixed y): spline along x -> zxx
        self.zxx = np.zeros((nx, ny))
        for j in range(ny):
            self.zxx[:, j] = _spline_1d_y2(self.x, z[:, j])
        # (3) per column: spline of zyy along x -> mixed derivative zxxyy
        self.zxxyy = np.zeros((nx, ny))
        for j in range(ny):
            self.zxxyy[:, j] = _spline_1d_y2(self.x, self.zyy[:, j])

    @staticmethod
    def _locate(g, xv):
        n = g.size
        if xv <= g[0]:
            return 0
        if xv >= g[-1]:
            return n - 2
        return int(np.searchsorted(g, xv, side="right") - 1)

    @staticmethod
    def _basis(gi, gip1, xv):
        h = gip1 - gi
        d0 = gip1 - xv
        d1 = xv - gi
        t0, t1 = d0 / h, d1 / h
        u0 = (d0 ** 3 / h - h * d0) / 6.0
        u1 = (d1 ** 3 / h - h * d1) / 6.0
        dt0, dt1 = -1.0 / h, 1.0 / h
        du0 = (h - 3.0 * d0 ** 2 / h) / 6.0
        du1 = (3.0 * d1 ** 2 / h - h) / 6.0
        return t0, t1, u0, u1, dt0, dt1, du0, du1

    def eval(self, xv, yv):
        i = self._locate(self.x, xv)
        j = self._locate(self.y, yv)
        bx = self._basis(self.x[i], self.x[i + 1], xv)
        by = self._basis(self.y[j], self.y[j + 1], yv)
        t0x, t1x, u0x, u1x = bx[:4]
        t0y, t1y, u0y, u1y = by[:4]
        return (
            self.z[i, j] * t0x * t0y + self.z[i + 1, j] * t1x * t0y
            + self.z[i, j + 1] * t0x * t1y + self.z[i + 1, j + 1] * t1x * t1y
            + self.zxx[i, j] * u0x * t0y + self.zxx[i + 1, j] * u1x * t0y
            + self.zxx[i, j + 1] * u0x * t1y + self.zxx[i + 1, j + 1] * u1x * t1y
            + self.zyy[i, j] * t0x * u0y + self.zyy[i + 1, j] * t1x * u0y
            + self.zyy[i, j + 1] * t0x * u1y + self.zyy[i + 1, j + 1] * t1x * u1y
            + self.zxxyy[i, j] * u0x * u0y + self.zxxyy[i + 1, j] * u1x * u0y
            + self.zxxyy[i, j + 1] * u0x * u1y + self.zxxyy[i + 1, j + 1] * u1x * u1y
        )

    def grad(self, xv, yv):
        """Analytic gradient (dS/dx, dS/dy) of the bicubic spline."""
        i = self._locate(self.x, xv)
        j = self._locate(self.y, yv)
        bx = self._basis(self.x[i], self.x[i + 1], xv)
        by = self._basis(self.y[j], self.y[j + 1], yv)
        t0x, t1x, u0x, u1x, dt0x, dt1x, du0x, du1x = bx
        t0y, t1y, u0y, u1y, dt0y, dt1y, du0y, du1y = by

        dzdx = (
            self.z[i, j] * dt0x * t0y + self.z[i + 1, j] * dt1x * t0y
            + self.z[i, j + 1] * dt0x * t1y + self.z[i + 1, j + 1] * dt1x * t1y
            + self.zxx[i, j] * du0x * t0y + self.zxx[i + 1, j] * du1x * t0y
            + self.zxx[i, j + 1] * du0x * t1y + self.zxx[i + 1, j + 1] * du1x * t1y
            + self.zyy[i, j] * dt0x * u0y + self.zyy[i + 1, j] * dt1x * u0y
            + self.zyy[i, j + 1] * dt0x * u1y + self.zyy[i + 1, j + 1] * dt1x * u1y
            + self.zxxyy[i, j] * du0x * u0y + self.zxxyy[i + 1, j] * du1x * u0y
            + self.zxxyy[i, j + 1] * du0x * u1y + self.zxxyy[i + 1, j + 1] * du1x * u1y
        )
        dzdy = (
            self.z[i, j] * t0x * dt0y + self.z[i + 1, j] * t1x * dt0y
            + self.z[i, j + 1] * t0x * dt1y + self.z[i + 1, j + 1] * t1x * dt1y
            + self.zxx[i, j] * u0x * dt0y + self.zxx[i + 1, j] * u1x * dt0y
            + self.zxx[i, j + 1] * u0x * dt1y + self.zxx[i + 1, j + 1] * u1x * dt1y
            + self.zyy[i, j] * t0x * du0y + self.zyy[i + 1, j] * t1x * du0y
            + self.zyy[i, j + 1] * t0x * du1y + self.zyy[i + 1, j + 1] * t1x * du1y
            + self.zxxyy[i, j] * u0x * du0y + self.zxxyy[i + 1, j] * u1x * du0y
            + self.zxxyy[i, j + 1] * u0x * du1y + self.zxxyy[i + 1, j + 1] * u1x * du1y
        )
        return dzdx, dzdy


class Spline3D:
    """Tricubic tensor-product spline on a regular 3D grid (natural boundaries).

    Mirrors `spline_3d_t` (src/mod_interpolation.f90): the 1D natural spline
    kernel is applied axis by axis — along x (-> fxx), then along y
    (-> fyy, fxy), then along z (-> fzz, fxz, fyz, fxyz) — so that the
    in-cell 2x2x2 Hermite expansion reproduces all eight node tensor
    components exactly (f, fxx, fyy, fzz, fxy, fxz, fyz, fxyz).
    """

    def __init__(self, x, y, z, f):
        self.x = np.asarray(x, dtype=float)
        self.y = np.asarray(y, dtype=float)
        self.z = np.asarray(z, dtype=float)
        f = np.asarray(f, dtype=float)
        nx, ny, nz = self.x.size, self.y.size, self.z.size
        if f.shape != (nx, ny, nz):
            raise ValueError("f must have shape (len(x), len(y), len(z))")
        self.f = f

        fxx = np.zeros((nx, ny, nz))
        for j in range(ny):
            for k in range(nz):
                fxx[:, j, k] = _spline_1d_y2(self.x, f[:, j, k])
        fyy = np.zeros((nx, ny, nz))
        fxy = np.zeros((nx, ny, nz))
        for i in range(nx):
            for k in range(nz):
                fyy[i, :, k] = _spline_1d_y2(self.y, f[i, :, k])
                fxy[i, :, k] = _spline_1d_y2(self.y, fxx[i, :, k])
        fzz = np.zeros((nx, ny, nz))
        fxz = np.zeros((nx, ny, nz))
        fyz = np.zeros((nx, ny, nz))
        fxyz = np.zeros((nx, ny, nz))
        for i in range(nx):
            for j in range(ny):
                fzz[i, j, :] = _spline_1d_y2(self.z, f[i, j, :])
                fxz[i, j, :] = _spline_1d_y2(self.z, fxx[i, j, :])
                fyz[i, j, :] = _spline_1d_y2(self.z, fyy[i, j, :])
                fxyz[i, j, :] = _spline_1d_y2(self.z, fxy[i, j, :])
        self.fxx, self.fyy, self.fzz = fxx, fyy, fzz
        self.fxy, self.fxz, self.fyz, self.fxyz = fxy, fxz, fyz, fxyz

    def _cell(self, xv, yv, zv):
        i = Spline2D._locate(self.x, xv)
        j = Spline2D._locate(self.y, yv)
        k = Spline2D._locate(self.z, zv)
        bx = Spline2D._basis(self.x[i], self.x[i + 1], xv)
        by = Spline2D._basis(self.y[j], self.y[j + 1], yv)
        bz = Spline2D._basis(self.z[k], self.z[k + 1], zv)
        return i, j, k, bx, by, bz

    def eval(self, xv, yv, zv):
        i, j, k, bx, by, bz = self._cell(xv, yv, zv)
        val = 0.0
        for p in (0, 1):
            for q in (0, 1):
                for r in (0, 1):
                    ii, jj, kk = i + p, j + q, k + r
                    tx, ux = bx[0 + p], bx[2 + p]
                    ty, uy = by[0 + q], by[2 + q]
                    tz, uz = bz[0 + r], bz[2 + r]
                    val += (
                        self.f[ii, jj, kk] * tx * ty * tz
                        + self.fxx[ii, jj, kk] * ux * ty * tz
                        + self.fyy[ii, jj, kk] * tx * uy * tz
                        + self.fzz[ii, jj, kk] * tx * ty * uz
                        + self.fxy[ii, jj, kk] * ux * uy * tz
                        + self.fxz[ii, jj, kk] * ux * ty * uz
                        + self.fyz[ii, jj, kk] * tx * uy * uz
                        + self.fxyz[ii, jj, kk] * ux * uy * uz
                    )
        return val

    def grad(self, xv, yv, zv):
        """Analytic gradient (df/dx, df/dy, df/dz)."""
        i, j, k, bx, by, bz = self._cell(xv, yv, zv)
        gx = gy = gz = 0.0
        for p in (0, 1):
            for q in (0, 1):
                for r in (0, 1):
                    ii, jj, kk = i + p, j + q, k + r
                    tx, ux, dtx, dux = bx[0 + p], bx[2 + p], bx[4 + p], bx[6 + p]
                    ty, uy, dty, duy = by[0 + q], by[2 + q], by[4 + q], by[6 + q]
                    tz, uz, dtz, duz = bz[0 + r], bz[2 + r], bz[4 + r], bz[6 + r]
                    gx += (
                        self.f[ii, jj, kk] * dtx * ty * tz
                        + self.fxx[ii, jj, kk] * dux * ty * tz
                        + self.fyy[ii, jj, kk] * dtx * uy * tz
                        + self.fzz[ii, jj, kk] * dtx * ty * uz
                        + self.fxy[ii, jj, kk] * dux * uy * tz
                        + self.fxz[ii, jj, kk] * dux * ty * uz
                        + self.fyz[ii, jj, kk] * dtx * uy * uz
                        + self.fxyz[ii, jj, kk] * dux * uy * uz
                    )
                    gy += (
                        self.f[ii, jj, kk] * tx * dty * tz
                        + self.fxx[ii, jj, kk] * ux * dty * tz
                        + self.fyy[ii, jj, kk] * tx * duy * tz
                        + self.fzz[ii, jj, kk] * tx * dty * uz
                        + self.fxy[ii, jj, kk] * ux * duy * tz
                        + self.fxz[ii, jj, kk] * ux * dty * uz
                        + self.fyz[ii, jj, kk] * tx * duy * uz
                        + self.fxyz[ii, jj, kk] * ux * duy * uz
                    )
                    gz += (
                        self.f[ii, jj, kk] * tx * ty * dtz
                        + self.fxx[ii, jj, kk] * ux * ty * dtz
                        + self.fyy[ii, jj, kk] * tx * uy * dtz
                        + self.fzz[ii, jj, kk] * tx * ty * duz
                        + self.fxy[ii, jj, kk] * ux * uy * dtz
                        + self.fxz[ii, jj, kk] * ux * ty * duz
                        + self.fyz[ii, jj, kk] * tx * uy * duz
                        + self.fxyz[ii, jj, kk] * ux * uy * duz
                    )
        return gx, gy, gz
