#!/usr/bin/env python3
"""Lab 1 system identification from the suspension simulator logs.

Reads the CSV the simulator writes (Time(s), Input, Output_Displacement, 100 ms
sampling) and works out the gain, the damping ratio, the natural frequency and
from those the spring and damping coefficients of

    x'' = -b x' - k x + F   ->   G(s) = 1 / (s^2 + b s + k)

Usage:
    analyse_lab1.py step.csv --sine-dir .
    analyse_lab1.py step.csv --sine 0.5=sine_0p5.csv --sine 1.0=sine_1p0.csv ...
    analyse_lab1.py --selftest

Every figure it writes is already annotated the way the report template asks for.
"""

import argparse
import csv
import json
import math
import os

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

matplotlib.rcParams.update({
    "font.family": "sans-serif",
    "font.size": 12,
    "axes.labelsize": 13,
    "axes.titlesize": 14,
    "figure.dpi": 160,
    "savefig.bbox": "tight",
})

ANNOT = dict(color="#b3261e", fontsize=11)


def _norm(name):
    """Header text reduced to letters and digits, so spacing and brackets stop mattering."""
    return "".join(c for c in name.lower() if c.isalnum())


def _column(row, *candidates):
    """The key in row whose header matches one of the candidate spellings."""
    table = {_norm(k): k for k in row if k}
    for cand in candidates:
        if _norm(cand) in table:
            return table[_norm(cand)]
    return None


def read_log(path):
    """Return time, input and output columns from a simulator CSV.

    The exact header spelling the simulator writes is matched loosely, because
    Time(s) and Time (s), or Output_Displacement and Output Displacement, are the
    same column and guessing wrong at the bench costs a session.
    """
    t, u, y = [], [], []
    cols = None
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh):
            if cols is None:
                ct = _column(row, "Time(s)", "Time", "Time s", "t")
                cu = _column(row, "Input", "Input(V)", "u")
                cy = _column(row, "Output_Displacement", "Output Displacement",
                             "Output", "Displacement", "x")
                missing = [n for n, c in (("time", ct), ("input", cu),
                                          ("output displacement", cy)) if c is None]
                if missing:
                    raise SystemExit(
                        f"{path}: cannot find the {', '.join(missing)} column.\n"
                        f"  headers in the file: {[k for k in row if k]}\n"
                        f"  add the spelling to _column() and rerun")
                cols = (ct, cu, cy)
            t.append(float(row[cols[0]]))
            u.append(float(row[cols[1]]))
            y.append(float(row[cols[2]]))
    if len(t) < 3:
        raise SystemExit(f"{path}: only {len(t)} samples, nothing to identify")
    return np.array(t), np.array(u), np.array(y)


def step_index(u):
    """First sample where the input leaves its starting value."""
    change = np.flatnonzero(np.abs(u - u[0]) > 1e-9)
    if change.size == 0:
        raise SystemExit("no step found in the input column")
    return int(change[0])


def turning_points(y):
    """Indices of local maxima and minima, ignoring flat runs."""
    d = np.diff(y)
    sign = np.sign(d)
    nz = sign != 0
    idx = np.flatnonzero(nz)
    if idx.size < 2:
        return np.array([], dtype=int), np.array([], dtype=int)
    s = sign[idx]
    flips = np.flatnonzero(s[:-1] != s[1:])
    peaks = [idx[i] + 1 for i in flips if s[i] > 0]
    troughs = [idx[i] + 1 for i in flips if s[i] < 0]
    return np.array(peaks, dtype=int), np.array(troughs, dtype=int)


def analyse_step(t, u, y, settle_fraction=0.1):
    """Identify the plant from one step log."""
    i0 = step_index(u)
    y0 = float(np.mean(y[:i0])) if i0 > 0 else float(y[0])
    u0 = float(u[0])
    du = float(np.mean(u[-max(1, int(len(u) * settle_fraction)):]) - u0)
    tail = max(1, int(len(y) * settle_fraction))
    y_inf = float(np.mean(y[-tail:]))
    dy = y_inf - y0

    res = {
        "step_time_s": float(t[i0]),
        "input_step": du,
        "output_baseline": y0,
        "output_final": y_inf,
        "output_change": dy,
        "gain_A": dy / du,
    }

    resp = y[i0:] - y0
    t_resp = t[i0:] - t[i0]
    peaks, troughs = turning_points(resp)
    peaks = peaks[resp[peaks] > dy] if peaks.size else peaks

    if peaks.size:
        res["damping_case"] = "underdamped"
        overshoot = (resp[peaks[0]] - dy) / dy
        res["first_peak_time_s"] = float(t_resp[peaks[0]])
        res["overshoot"] = float(overshoot)
        ln = math.log(overshoot)
        zeta = -ln / math.sqrt(math.pi ** 2 + ln ** 2)

        # damped period, averaged over whatever peaks the log actually contains
        if peaks.size >= 2:
            period = float(np.mean(np.diff(t_resp[peaks])))
        elif troughs.size:
            period = float(2 * (t_resp[troughs[0]] - t_resp[peaks[0]]))
        else:
            raise SystemExit("only one turning point, log a longer response")
        wd = 2 * math.pi / period
        wn = wd / math.sqrt(1 - zeta ** 2)
        res.update({
            "damped_period_s": period,
            "wd_rad_s": wd,
            "zeta": float(zeta),
            "wn_rad_s": float(wn),
            "peak_times_s": [float(v) for v in t_resp[peaks]],
        })

        # log decrement as a second, independent reading of zeta
        if peaks.size >= 2:
            ratios = np.log((resp[peaks[:-1]] - dy) / (resp[peaks[1:]] - dy))
            delta = float(np.mean(ratios))
            res["zeta_log_decrement"] = float(delta / math.sqrt(4 * math.pi ** 2 + delta ** 2))
    else:
        res["damping_case"] = "overdamped or critically damped"
        res.update(fit_overdamped(t_resp, resp, dy))

    zeta, wn = res["zeta"], res["wn_rad_s"]
    res["k"] = float(wn ** 2)
    res["b"] = float(2 * zeta * wn)
    res["poles"] = poles(zeta, wn)
    res["corner_frequency_rad_s"] = float(wn)

    # k arrives twice by different routes. They agree only if the input column is
    # force in the same units the model assumes, so the gap between them is the
    # scaling check, not a rounding artefact.
    res["k_from_wn"] = float(wn ** 2)
    res["k_from_gain"] = float(1.0 / res["gain_A"]) if res["gain_A"] else float("nan")
    if math.isfinite(res["k_from_gain"]) and res["k_from_wn"]:
        res["k_mismatch_percent"] = float(
            100 * (res["k_from_gain"] - res["k_from_wn"]) / res["k_from_wn"])
    return res


def fit_overdamped(t, y, dy):
    """Grid search for zeta >= 1 and wn when the step shows no overshoot.

    A no-overshoot response gives no peak to read, so the parameters come from a
    least squares sweep over the analytic step response instead.
    """
    wn_grid = np.geomspace(0.05, 50.0, 400)
    zeta_grid = np.linspace(1.0, 6.0, 300)
    best = None
    for wn in wn_grid:
        for zeta in zeta_grid:
            err = float(np.sum((step_model(t, dy, zeta, wn) - y) ** 2))
            if best is None or err < best[0]:
                best = (err, zeta, wn)
    _, zeta, wn = best
    # 63.2 percent point, quoted as the reduced first order model the brief allows
    reach = np.flatnonzero(y >= 0.632 * dy)
    out = {"zeta": float(zeta), "wn_rad_s": float(wn), "fit_method": "least squares sweep"}

    # A best fit sitting on the edge of the sweep means the true value is outside it,
    # so the number reported is the boundary and not an identification.
    edge = []
    if zeta <= zeta_grid[1]:
        edge.append(f"zeta pinned at the lower limit {zeta_grid[0]:.2f}")
    if zeta >= zeta_grid[-2]:
        edge.append(f"zeta pinned at the upper limit {zeta_grid[-1]:.2f}")
    if wn <= wn_grid[1]:
        edge.append(f"wn pinned at the lower limit {wn_grid[0]:.3f}")
    if wn >= wn_grid[-2]:
        edge.append(f"wn pinned at the upper limit {wn_grid[-1]:.1f}")
    if edge:
        out["grid_warning"] = "; ".join(edge)
        print("WARNING: overdamped fit hit the edge of the search grid: "
              + out["grid_warning"] + "\n  widen the ranges in fit_overdamped and rerun")
    if reach.size:
        tau = float(t[reach[0]])
        out["tau_63_s"] = tau
        out["reduced_order_pole"] = -1.0 / tau
    return out


def step_model(t, dy, zeta, wn):
    """Analytic unit-ish step response of A wn^2 / (s^2 + 2 zeta wn s + wn^2)."""
    if zeta < 1:
        wd = wn * math.sqrt(1 - zeta ** 2)
        phi = math.atan2(math.sqrt(1 - zeta ** 2), zeta)
        env = np.exp(-zeta * wn * t) / math.sqrt(1 - zeta ** 2)
        return dy * (1 - env * np.sin(wd * t + phi))
    if abs(zeta - 1) < 1e-9:
        return dy * (1 - np.exp(-wn * t) * (1 + wn * t))
    r = math.sqrt(zeta ** 2 - 1)
    s1, s2 = -wn * (zeta - r), -wn * (zeta + r)
    return dy * (1 + (s2 * np.exp(s1 * t) - s1 * np.exp(s2 * t)) / (s1 - s2))


def poles(zeta, wn):
    if zeta < 1:
        wd = wn * math.sqrt(1 - zeta ** 2)
        return [[-zeta * wn, wd], [-zeta * wn, -wd]]
    r = math.sqrt(zeta ** 2 - 1)
    return [[-wn * (zeta - r), 0.0], [-wn * (zeta + r), 0.0]]


def fit_sinusoid(t, y, w):
    """Least squares amplitude and phase of y at a known frequency."""
    basis = np.column_stack([np.cos(w * t), np.sin(w * t), np.ones_like(t)])
    coef, *_ = np.linalg.lstsq(basis, y, rcond=None)
    a, b, _ = coef
    return math.hypot(a, b), math.atan2(-b, a)


def sine_dir_logs(directory):
    """Sine logs in a folder, sorted, so the frequencies need not be typed out."""
    if not directory:
        return []
    names = [n for n in sorted(os.listdir(directory))
             if n.lower().startswith("sine_") and n.lower().endswith(".csv")]
    if not names:
        raise SystemExit(f"{directory}: no files named sine_<frequency>.csv")
    return [os.path.join(directory, n) for n in names]


def frequency_from_name(path):
    """Drive frequency from a name like sine_1p6.csv, where 1p6 means 1.6 rad/s."""
    stem = os.path.splitext(os.path.basename(path))[0]
    text = stem.split("_", 1)[1].replace("p", ".")
    try:
        return float(text)
    except ValueError:
        raise SystemExit(f"{path}: cannot read a frequency out of the file name")


def analyse_sine(path, w, skip_fraction=0.5):
    """One frequency point: magnitude ratio in dB and phase shift in degrees."""
    t, u, y = read_log(path)
    start = int(len(t) * skip_fraction)
    t, u, y = t[start:], u[start:], y[start:]
    amp_u, ph_u = fit_sinusoid(t, u, w)
    amp_y, ph_y = fit_sinusoid(t, y, w)
    if amp_u < 1e-9:
        raise SystemExit(f"{path}: input amplitude is zero at w = {w}")
    # Wrap into (-180, 180] first. A second order lowpass only ever lags, so a
    # reading well above zero is a wrapped lag, while a few degrees above zero is
    # just noise at a frequency far below the corner and must stay near zero.
    phase = math.degrees(ph_y - ph_u)
    phase = (phase + 180) % 360 - 180
    if phase > 10:
        phase -= 360
    return {
        "file": os.path.basename(path),
        "w_rad_s": w,
        "input_amplitude": amp_u,
        "output_amplitude": amp_y,
        "magnitude": amp_y / amp_u,
        "magnitude_dB": 20 * math.log10(amp_y / amp_u),
        "phase_deg": phase,
        "samples_used": int(len(t)),
    }


def model_response(w, gain, zeta, wn):
    s = 1j * w
    g = gain * wn ** 2 / (s ** 2 + 2 * zeta * wn * s + wn ** 2)
    return 20 * np.log10(np.abs(g)), np.degrees(np.unwrap(np.angle(g)))


def frequency_summary(points, res):
    """Answer the four Stage 1 questions a second time, from the measured points.

    The brief is built around comparing these against the step test answers, so
    every quantity here has a counterpart in analyse_step.
    """
    if len(points) < 3:
        return {"note": "fewer than three frequency points, no independent reading"}
    pts = sorted(points, key=lambda p: p["w_rad_s"])
    w = np.array([p["w_rad_s"] for p in pts])
    mag = np.array([p["magnitude"] for p in pts])
    mag_db = np.array([p["magnitude_dB"] for p in pts])
    ph = np.array([p["phase_deg"] for p in pts])
    wn_step = res["wn_rad_s"]
    out = {}

    # Question 1 again: the gain is the low frequency asymptote, read off the lowest
    # point rather than an average, because the curve is already rising towards
    # resonance and averaging over a rising curve only adds bias.
    a_freq = float(mag[0])
    ratio = float(w[0] / wn_step)
    out["gain_A_from_asymptote"] = a_freq
    out["gain_lowest_w_over_wn"] = ratio
    out["gain_mismatch_percent"] = float(100 * (a_freq - res["gain_A"]) / res["gain_A"])
    # a lightly damped plant has not settled to A yet at w = w[0], so quote how much
    # of the mismatch the model itself predicts before blaming the measurement
    denom = math.hypot(1 - ratio ** 2, 2 * res["zeta"] * ratio)
    out["gain_expected_bias_percent"] = float(100 * (1 / denom - 1)) if denom else 0.0
    if ratio > 0.3:
        out["gain_note"] = (f"lowest frequency is {ratio:.2f} wn, too close to the corner "
                            "for a clean asymptote, so this reads high")

    # Questions 2 and 3 again: only a complex pair lifts the magnitude above the asymptote
    i_peak = int(np.argmax(mag_db))
    interior = 0 < i_peak < len(w) - 1
    mr = float(mag[i_peak] / a_freq)
    out.update({"peak_w_rad_s": float(w[i_peak]), "Mr": mr, "peak_is_interior": interior})
    if interior and mr > 1.0:
        q = 1.0 / (2 * mr)
        disc = 1 - 4 * q ** 2
        if disc >= 0:
            z_sq = (1 - math.sqrt(disc)) / 2
            out["zeta_from_Mr"] = float(math.sqrt(z_sq))
            if 1 - 2 * z_sq > 0:
                out["wn_from_peak"] = float(w[i_peak] / math.sqrt(1 - 2 * z_sq))
    else:
        out["resonance_note"] = ("no interior magnitude peak, so the poles are real or "
                                 "the resonance sits outside the frequencies tested")

    # Question 4 again: the phase crosses -90 degrees at wn whatever zeta is
    for i in range(len(w) - 1):
        if (ph[i] + 90) * (ph[i + 1] + 90) <= 0 and ph[i] != ph[i + 1]:
            f = (-90 - ph[i]) / (ph[i + 1] - ph[i])
            lo, hi = math.log10(w[i]), math.log10(w[i + 1])
            out["wn_from_phase"] = float(10 ** (lo + f * (hi - lo)))
            break
    else:
        out["phase_note"] = "the phase never reaches -90 degrees, test a higher frequency"

    # Two poles and no zeros roll off at -40 dB/decade, but only well above wn. Just
    # past resonance the curve falls far steeper than that, so a slope measured from
    # points at 2 wn and 4 wn reads nearer -60 and is not evidence against the model.
    hf = w >= 2.0 * wn_step
    out["highest_w_over_wn"] = float(w[-1] / wn_step)
    if int(hf.sum()) >= 2:
        out["hf_slope_dB_per_decade"] = float(np.polyfit(np.log10(w[hf]), mag_db[hf], 1)[0])
        out["hf_points_used"] = int(hf.sum())
        if out["highest_w_over_wn"] < 4:
            out["slope_note"] = (f"highest frequency is only {out['highest_w_over_wn']:.1f} wn, "
                                 "so the -40 dB/decade asymptote is not reached yet and this "
                                 "slope reads steeper. Do not call that a failed model")
    else:
        out["slope_note"] = "fewer than two points above 2 wn, roll off not measurable"

    # validation in the frequency domain, against a model built from the step alone
    mmag, mph = model_response(w, res["gain_A"], res["zeta"], wn_step)
    mph = np.where(mph > 10, mph - 360, mph)
    out["model_mag_rms_dB"] = float(np.sqrt(np.mean((mag_db - mmag) ** 2)))
    out["model_phase_rms_deg"] = float(np.sqrt(np.mean((ph - mph) ** 2)))
    return out


def plot_step_input(t, u, outdir):
    fig, ax = plt.subplots(figsize=(7, 3.2))
    ax.plot(t, u, color="#1f3b73", lw=2)
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Input [V]")
    ax.set_title("Step input applied to the suspension")
    ax.grid(alpha=0.3)
    save(fig, outdir, "input_step.png")


def plot_step_response(t, u, y, res, outdir):
    i0 = step_index(u)
    fig, ax = plt.subplots(figsize=(7.5, 4.4))
    ax.plot(t, y, color="#1f3b73", lw=2, label="measured displacement")
    ax.axhline(res["output_final"], color="#b3261e", ls="--", lw=1.2)
    ax.annotate(f"final value = {res['output_final']:.3f} m",
                xy=(t[-1], res["output_final"]), xytext=(-10, 8),
                textcoords="offset points", ha="right", **ANNOT)
    ax.axvline(t[i0], color="#4a4a4a", ls=":", lw=1.2)
    ax.annotate("step applied", xy=(t[i0], y.min()), xytext=(6, 10),
                textcoords="offset points", color="#4a4a4a", fontsize=11)

    if res["damping_case"] == "underdamped":
        tp = res["step_time_s"] + res["first_peak_time_s"]
        peak = res["output_final"] + res["overshoot"] * res["output_change"]
        ax.annotate("", xy=(tp, peak), xytext=(tp, res["output_final"]),
                    arrowprops=dict(arrowstyle="<->", color="#b3261e", lw=1.4))
        ax.annotate(f"overshoot = {100 * res['overshoot']:.1f}%",
                    xy=(tp, (peak + res['output_final']) / 2), xytext=(10, 0),
                    textcoords="offset points", va="center", **ANNOT)
        if len(res.get("peak_times_s", [])) >= 2:
            t1 = res["step_time_s"] + res["peak_times_s"][0]
            t2 = res["step_time_s"] + res["peak_times_s"][1]
            lvl = res["output_final"] + 1.15 * res["overshoot"] * res["output_change"]
            ax.annotate("", xy=(t1, lvl), xytext=(t2, lvl),
                        arrowprops=dict(arrowstyle="<->", color="#b3261e", lw=1.4))
            ax.annotate(f"period T = {res['damped_period_s']:.2f} s",
                        xy=(t2, lvl), xytext=(8, -4),
                        textcoords="offset points", ha="left", **ANNOT)

    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Displacement [m]")
    ax.margins(y=0.16)
    ax.set_title("Measured step response")
    ax.grid(alpha=0.3)
    ax.legend(loc="lower right")
    save(fig, outdir, "step_response.png")


def plot_bode(points, res, outdir, name="bode_plot.png", title="Measured frequency response"):
    w = np.array([p["w_rad_s"] for p in points])
    mag = np.array([p["magnitude_dB"] for p in points])
    ph = np.array([p["phase_deg"] for p in points])
    grid = np.geomspace(min(w.min(), res["wn_rad_s"]) / 5, max(w.max(), res["wn_rad_s"]) * 5, 400)
    mmag, mph = model_response(grid, res["gain_A"], res["zeta"], res["wn_rad_s"])

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(7.5, 6), sharex=True)
    ax1.semilogx(grid, mmag, color="#4a4a4a", lw=1.6, label="identified model")
    ax1.semilogx(w, mag, "o", color="#1f3b73", ms=7, label="measured")
    ax1.axvline(res["wn_rad_s"], color="#b3261e", ls="--", lw=1.2)
    ax1.annotate(f"$\\omega_n$ = {res['wn_rad_s']:.2f} rad/s",
                 xy=(res["wn_rad_s"], mag.max()), xytext=(8, -4),
                 textcoords="offset points", **ANNOT)
    ax1.set_ylabel("Magnitude [dB]")
    ax1.grid(which="both", alpha=0.3)
    ax1.legend(loc="lower left")
    ax1.set_title(title)

    ax2.semilogx(grid, mph, color="#4a4a4a", lw=1.6)
    ax2.semilogx(w, ph, "o", color="#1f3b73", ms=7)
    ax2.axvline(res["wn_rad_s"], color="#b3261e", ls="--", lw=1.2)
    ax2.axhline(-90, color="#b3261e", ls=":", lw=1.2)
    ax2.annotate("-90 deg at the natural frequency", xy=(w.min(), -90),
                 xytext=(4, 6), textcoords="offset points", **ANNOT)
    ax2.set_xlabel("Frequency [rad/s]")
    ax2.set_ylabel("Phase [deg]")
    ax2.grid(which="both", alpha=0.3)
    save(fig, outdir, name)


def plot_validation(t, u, y, res, outdir):
    i0 = step_index(u)
    t_resp = t[i0:] - t[i0]
    model = res["output_baseline"] + step_model(t_resp, res["output_change"],
                                                res["zeta"], res["wn_rad_s"])
    fig, ax = plt.subplots(figsize=(7.5, 4.4))
    ax.plot(t[i0:], y[i0:], color="#1f3b73", lw=2, label="measured")
    ax.plot(t[i0:], model, color="#b3261e", lw=1.8, ls="--", label="identified model")
    rms = float(np.sqrt(np.mean((y[i0:] - model) ** 2)))
    ax.annotate(f"RMS error = {rms:.4f} m", xy=(0.98, 0.06), xycoords="axes fraction",
                ha="right", **ANNOT)
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Displacement [m]")
    ax.set_title("Simulated against measured step response")
    ax.grid(alpha=0.3)
    ax.legend(loc="lower right")
    save(fig, outdir, "validation_step.png")
    return rms


def save(fig, outdir, name):
    path = os.path.join(outdir, name)
    fig.savefig(path)
    plt.close(fig)
    print(f"wrote {path}")


def summarise(res, points, rms, freq=None, repeats=None):
    zeta, wn = res["zeta"], res["wn_rad_s"]
    print("\nIdentified model")
    print(f"  gain A          {res['gain_A']:.4f} m/V   ->  k = 1/A = {1 / res['gain_A']:.4f}")
    print(f"  damping case    {res['damping_case']}")
    print(f"  zeta            {zeta:.4f}")
    if "zeta_log_decrement" in res:
        print(f"  zeta (log dec)  {res['zeta_log_decrement']:.4f}")
    print(f"  wn              {wn:.4f} rad/s")
    print(f"  k = wn^2        {res['k']:.4f}")
    print(f"  b = 2 zeta wn   {res['b']:.4f}")
    p = res["poles"]
    print(f"  poles           {p[0][0]:.4f} + {p[0][1]:.4f}j, {p[1][0]:.4f} + {p[1][1]:.4f}j")
    print(f"  G(s) = 1 / (s^2 + {res['b']:.4f} s + {res['k']:.4f})")
    if rms is not None:
        print(f"  step RMS error  {rms:.5f} m")

    if "k_mismatch_percent" in res:
        print("\nGain consistency check")
        print(f"  k from the final value   {res['k_from_gain']:.4f}")
        print(f"  k from wn^2              {res['k_from_wn']:.4f}")
        print(f"  difference               {res['k_mismatch_percent']:+.1f}%")
        if abs(res["k_mismatch_percent"]) > 10:
            print("  the two disagree by more than 10%, so the input column is probably"
                  "\n  not force in the units the model assumes. Say so in the report.")

    if repeats and len(repeats) > 1:
        print("\nRepeatability across step runs")
        print("  file            A         zeta      wn        k         b")
        for name, r in repeats:
            print(f"  {name:14s}  {r['gain_A']:8.4f}  {r['zeta']:8.4f}  "
                  f"{r['wn_rad_s']:8.4f}  {r['k']:8.4f}  {r['b']:8.4f}")
        for key, label in (("gain_A", "A"), ("zeta", "zeta"), ("wn_rad_s", "wn"),
                           ("k", "k"), ("b", "b")):
            vals = np.array([r[key] for _, r in repeats])
            spread = 100 * (vals.max() - vals.min()) / abs(vals.mean())
            print(f"  {label:5s} mean {vals.mean():8.4f}   spread {spread:5.1f}% of mean")

    if points:
        print("\nFrequency response points")
        print("  w [rad/s]   |G| [dB]   phase [deg]   file")
        for pt in points:
            print(f"  {pt['w_rad_s']:9.3f}  {pt['magnitude_dB']:9.2f}  {pt['phase_deg']:11.1f}   {pt['file']}")

    if freq:
        print("\nRead a second time off the frequency response")
        if "gain_A_from_asymptote" in freq:
            print(f"  gain A                   {freq['gain_A_from_asymptote']:.4f} m/V "
                  f"({freq['gain_mismatch_percent']:+.1f}% against the step)")
        if "zeta_from_Mr" in freq:
            print(f"  zeta from the peak Mr    {freq['zeta_from_Mr']:.4f} "
                  f"(step gave {zeta:.4f})")
        if "wn_from_peak" in freq:
            print(f"  wn from the peak         {freq['wn_from_peak']:.4f} rad/s")
        if "wn_from_phase" in freq:
            print(f"  wn from the -90 crossing {freq['wn_from_phase']:.4f} rad/s "
                  f"({100 * (freq['wn_from_phase'] - wn) / wn:+.1f}% against the step)")
        if "hf_slope_dB_per_decade" in freq:
            print(f"  high frequency slope     {freq['hf_slope_dB_per_decade']:.1f} dB/decade "
                  f"over {freq['hf_points_used']} points, expect -40 for two poles")
        print(f"  model vs measured        {freq['model_mag_rms_dB']:.2f} dB RMS, "
              f"{freq['model_phase_rms_deg']:.1f} deg RMS")
        for k in ("gain_note", "resonance_note", "phase_note", "slope_note", "note"):
            if k in freq:
                print(f"  note: {freq[k]}")


def selftest(outdir):
    """Run the whole pipeline on data generated from a known plant.

    This checks the maths, nothing else. The numbers here are synthetic and must
    never end up in the report.
    """
    k_true, b_true = 4.0, 0.8
    wn = math.sqrt(k_true)
    zeta = b_true / (2 * wn)
    gain = 1 / k_true
    t = np.arange(0, 60, 0.1)
    u = np.where(t >= 5, 1.0, 0.0)
    y = np.zeros_like(t)
    mask = t >= 5
    y[mask] = 5.0 + step_model(t[mask] - 5, gain * 1.0, zeta, wn)
    y[~mask] = 5.0

    res = analyse_step(t, u, y)
    print(f"selftest: true k = {k_true}, b = {b_true}")
    print(f"selftest: found k = {res['k']:.4f}, b = {res['b']:.4f}")
    assert abs(res["k"] - k_true) < 0.15, res["k"]
    assert abs(res["b"] - b_true) < 0.15, res["b"]

    # the no-overshoot branch takes a different route to zeta and wn, check it too
    k_od, b_od = 1.0, 3.0
    wn_od = math.sqrt(k_od)
    zeta_od = b_od / (2 * wn_od)
    y_od = np.full_like(t, 5.0)
    y_od[mask] = 5.0 + step_model(t[mask] - 5, 1 / k_od, zeta_od, wn_od)
    res_od = analyse_step(t, u, y_od)
    print(f"selftest overdamped: true k = {k_od}, b = {b_od}, "
          f"found k = {res_od['k']:.4f}, b = {res_od['b']:.4f}")
    assert res_od["damping_case"].startswith("overdamped"), res_od["damping_case"]
    assert abs(res_od["k"] - k_od) < 0.1, res_od["k"]
    assert abs(res_od["b"] - b_od) < 0.2, res_od["b"]

    points = []
    for w in [0.4, 0.5, 1.0, 1.6, 2.0, 2.5, 4.0, 8.0, 16.0]:
        ts = np.arange(0, 120, 0.1)
        us = np.sin(w * ts)
        mag, ph = model_response(np.array([w]), gain, zeta, wn)
        ys = 5.0 + 10 ** (mag[0] / 20) * np.sin(w * ts + math.radians(ph[0]))
        path = os.path.join(outdir, f"selftest_sine_{w}.csv")
        with open(path, "w", newline="") as fh:
            wr = csv.writer(fh)
            wr.writerow(["Time(s)", "Input", "Output_Displacement"])
            wr.writerows(zip(ts, us, ys))
        pt = analyse_sine(path, w)
        os.remove(path)
        assert abs(pt["magnitude_dB"] - mag[0]) < 0.1, pt
        assert abs(pt["phase_deg"] - ph[0]) < 1.0, pt
        points.append(pt)

    # the frequency reading must land on the same plant the step test found
    freq = frequency_summary(points, res)
    assert abs(freq["gain_A_from_asymptote"] - gain) / gain < 0.10, freq
    assert abs(freq["zeta_from_Mr"] - zeta) < 0.05, freq
    assert abs(freq["wn_from_phase"] - wn) < 0.05, freq
    assert abs(freq["hf_slope_dB_per_decade"] + 40) < 5, freq

    # header spelling must not matter, this is the one thing that cannot be retried
    for header in (["Time (s)", "Input", "Output Displacement"],
                   ["time", "INPUT", "output_displacement"]):
        path = os.path.join(outdir, "selftest_headers.csv")
        with open(path, "w", newline="") as fh:
            wr = csv.writer(fh)
            wr.writerow(header)
            wr.writerows(zip(t, u, y))
        tt, uu, yy = read_log(path)
        os.remove(path)
        assert len(tt) == len(t) and abs(analyse_step(tt, uu, yy)["k"] - k_true) < 0.15

    plot_step_input(t, u, outdir)
    plot_step_response(t, u, y, res, outdir)
    plot_bode(points, res, outdir)
    rms = plot_validation(t, u, y, res, outdir)
    summarise(res, points, rms, freq, [("selftest", res)])
    print("\nselftest passed, all figures written from synthetic data")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("step_csv", nargs="*",
                    help="step test logs from the simulator. Pass all the repeat runs "
                         "to get a spread; the first one is the source of the figures")
    ap.add_argument("--sine", action="append", default=[], metavar="W=FILE",
                    help="sine test log and the frequency in rad/s it was run at")
    ap.add_argument("--sine-dir", metavar="DIR",
                    help="folder of sine logs named sine_<frequency>.csv, "
                         "where 1p6 means 1.6 rad/s")
    ap.add_argument("--outdir", default="report/figures", help="where the figures go")
    ap.add_argument("--skip", type=float, default=0.5,
                    help="fraction of each sine log to discard as transient")
    ap.add_argument("--selftest", action="store_true",
                    help="check the analysis against a plant with known k and b")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    if args.selftest:
        selftest(args.outdir)
        return
    if not args.step_csv:
        ap.error("give a step log, or --selftest")

    repeats = [(os.path.basename(p), analyse_step(*read_log(p))) for p in args.step_csv]
    t, u, y = read_log(args.step_csv[0])
    res = repeats[0][1]
    points = []
    for spec in args.sine:
        w, _, path = spec.partition("=")
        points.append(analyse_sine(path, float(w), args.skip))
    for path in sine_dir_logs(args.sine_dir):
        points.append(analyse_sine(path, frequency_from_name(path), args.skip))
    points.sort(key=lambda p: p["w_rad_s"])

    plot_step_input(t, u, args.outdir)
    plot_step_response(t, u, y, res, args.outdir)
    rms = plot_validation(t, u, y, res, args.outdir)
    freq = {}
    if points:
        plot_bode(points, res, args.outdir)
        freq = frequency_summary(points, res)

    summarise(res, points, rms, freq, repeats)
    out = {"step": res, "frequency_points": points, "step_rms_error": rms,
           "frequency_summary": freq,
           "step_repeats": {name: r for name, r in repeats}}
    with open(os.path.join(args.outdir, "results.json"), "w") as fh:
        json.dump(out, fh, indent=2)
    print(f"wrote {os.path.join(args.outdir, 'results.json')}")


if __name__ == "__main__":
    main()
