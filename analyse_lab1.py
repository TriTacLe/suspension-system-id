#!/usr/bin/env python3
"""Lab 1 system identification from the suspension simulator logs.

Reads the CSV the simulator writes (Time(s), Input, Output_Displacement, 100 ms
sampling) and works out the gain, the damping ratio, the natural frequency and
from those the spring and damping coefficients of

    x'' = -b x' - k x + F   ->   G(s) = 1 / (s^2 + b s + k)

Usage:
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


def read_log(path):
    """Return time, input and output columns from a simulator CSV."""
    t, u, y = [], [], []
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh):
            keys = {k.strip().lower(): k for k in row}
            t.append(float(row[keys["time(s)"]]))
            u.append(float(row[keys["input"]]))
            y.append(float(row[keys["output_displacement"]]))
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
    return res


def fit_overdamped(t, y, dy):
    """Grid search for zeta >= 1 and wn when the step shows no overshoot.

    A no-overshoot response gives no peak to read, so the parameters come from a
    least squares sweep over the analytic step response instead.
    """
    best = None
    for wn in np.geomspace(0.05, 50.0, 400):
        for zeta in np.linspace(1.0, 6.0, 300):
            err = float(np.sum((step_model(t, dy, zeta, wn) - y) ** 2))
            if best is None or err < best[0]:
                best = (err, zeta, wn)
    _, zeta, wn = best
    # 63.2 percent point, quoted as the reduced first order model the brief allows
    reach = np.flatnonzero(y >= 0.632 * dy)
    out = {"zeta": float(zeta), "wn_rad_s": float(wn), "fit_method": "least squares sweep"}
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


def analyse_sine(path, w, skip_fraction=0.5):
    """One frequency point: magnitude ratio in dB and phase shift in degrees."""
    t, u, y = read_log(path)
    start = int(len(t) * skip_fraction)
    t, u, y = t[start:], u[start:], y[start:]
    amp_u, ph_u = fit_sinusoid(t, u, w)
    amp_y, ph_y = fit_sinusoid(t, y, w)
    if amp_u < 1e-9:
        raise SystemExit(f"{path}: input amplitude is zero at w = {w}")
    phase = math.degrees(ph_y - ph_u)
    while phase > 0:
        phase -= 360
    while phase < -360:
        phase += 360
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


def summarise(res, points, rms):
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
    if points:
        print("\nFrequency response points")
        print("  w [rad/s]   |G| [dB]   phase [deg]   file")
        for pt in points:
            print(f"  {pt['w_rad_s']:9.3f}  {pt['magnitude_dB']:9.2f}  {pt['phase_deg']:11.1f}   {pt['file']}")


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
    for w in [0.5, 1.0, 1.6, 2.0, 2.5, 4.0]:
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

    plot_step_input(t, u, outdir)
    plot_step_response(t, u, y, res, outdir)
    plot_bode(points, res, outdir)
    rms = plot_validation(t, u, y, res, outdir)
    summarise(res, points, rms)
    print("\nselftest passed, all figures written from synthetic data")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("step_csv", nargs="?", help="step test log from the simulator")
    ap.add_argument("--sine", action="append", default=[], metavar="W=FILE",
                    help="sine test log and the frequency in rad/s it was run at")
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

    t, u, y = read_log(args.step_csv)
    res = analyse_step(t, u, y)
    points = []
    for spec in args.sine:
        w, _, path = spec.partition("=")
        points.append(analyse_sine(path, float(w), args.skip))
    points.sort(key=lambda p: p["w_rad_s"])

    plot_step_input(t, u, args.outdir)
    plot_step_response(t, u, y, res, args.outdir)
    rms = plot_validation(t, u, y, res, args.outdir)
    if points:
        plot_bode(points, res, args.outdir)

    summarise(res, points, rms)
    out = {"step": res, "frequency_points": points, "step_rms_error": rms}
    with open(os.path.join(args.outdir, "results.json"), "w") as fh:
        json.dump(out, fh, indent=2)
    print(f"wrote {os.path.join(args.outdir, 'results.json')}")


if __name__ == "__main__":
    main()
