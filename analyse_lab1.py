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

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt

matplotlib.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.size": 12,
        "axes.labelsize": 13,
        "axes.titlesize": 14,
        "figure.dpi": 160,
        "savefig.bbox": "tight",
    }
)

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
                cy = _column(
                    row,
                    "Output_Displacement",
                    "Output Displacement",
                    "Output",
                    "Displacement",
                    "x",
                )
                missing = [
                    n
                    for n, c in (
                        ("time", ct),
                        ("input", cu),
                        ("output displacement", cy),
                    )
                    if c is None
                ]
                if missing:
                    raise SystemExit(
                        f"{path}: cannot find the {', '.join(missing)} column.\n"
                        f"  headers in the file: {[k for k in row if k]}\n"
                        f"  add the spelling to _column() and rerun"
                    )
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


def identify_underdamped(t, y, dy, resolution):
    """Damping ratio and natural frequency from the overshoot and the peak timing.

    t and y are measured from the step, dy is the total change in the output.
    Returns None when no peak clears the resolution, which is how the caller finds
    out the response has no usable overshoot and needs the overdamped fit instead.
    """
    peaks, troughs = turning_points(y)

    # The output is quantised, so a single sample one step above the final value
    # is indistinguishable from a late peak and drags the averaged period badly.
    # Require a peak to clear the final value by several resolution steps first.
    peaks = peaks[y[peaks] > dy + 3 * resolution]
    if not peaks.size:
        return None

    overshoot = (y[peaks[0]] - dy) / dy
    ln = math.log(overshoot)
    zeta = -ln / math.sqrt(math.pi**2 + ln**2)

    # Averaging over several peaks is better when the log has them. When only the
    # first overshoot clears the resolution, the time to it is half a damped period
    # and that is the whole of the timing information available.
    if peaks.size >= 2:
        period = float(np.mean(np.diff(t[peaks])))
        source = f"mean spacing of {peaks.size} peaks"
    else:
        period = float(2 * t[peaks[0]])
        source = "twice the time to the only peak clearing the resolution"
    wd = 2 * math.pi / period
    wn = wd / math.sqrt(1 - zeta**2)

    res = {
        "first_peak_time_s": float(t[peaks[0]]),
        "overshoot": float(overshoot),
        "damped_period_s": period,
        "damped_period_source": source,
        "wd_rad_s": wd,
        "zeta": float(zeta),
        "wn_rad_s": float(wn),
        "peak_times_s": [float(v) for v in t[peaks]],
    }

    # Log decrement as a second, independent reading of zeta. Extrema alternate
    # max, min, max at t = n*pi/wd, so the first trough is a usable second extremum
    # even when the second peak has already decayed under the resolution. Its
    # timing is not usable though: quantisation flattens the bottom over nearly a
    # second, so the period stays on the sharp first peak.
    after = troughs[troughs > peaks[0]]
    after = after[np.abs(y[after] - dy) > 3 * resolution]
    if after.size:
        extrema = np.array([peaks[0], after[0]])
        dev = np.abs(y[extrema] - dy)
        ln_r = math.log(dev[1] / dev[0])
        zeta_r = -ln_r / math.sqrt(math.pi**2 + ln_r**2)
        res["zeta_log_decrement"] = float(zeta_r)
        res["extrema_deviations_m"] = [float(v) for v in dev]
        res["extrema_times_s"] = [float(v) for v in t[extrema]]
    return res


def analyse_step(t, u, y, settle_fraction=0.1):
    """Identify the plant from one step log."""
    i0 = step_index(u)
    y0 = float(np.mean(y[:i0])) if i0 > 0 else float(y[0])
    u0 = float(u[0])
    du = float(np.mean(u[-max(1, int(len(u) * settle_fraction)) :]) - u0)
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
    levels = np.unique(y)
    resolution = float(np.min(np.diff(levels))) if levels.size > 1 else 0.0
    res["output_resolution_m"] = resolution

    under = identify_underdamped(t_resp, resp, dy, resolution)
    if under:
        res["damping_case"] = "underdamped"
        res.update(under)
    else:
        res["damping_case"] = "overdamped or critically damped"
        res.update(fit_overdamped(t_resp, resp, dy))

    # What a two real pole fit can manage on the same record, as the numerical
    # counterpart to reading the pole type off the shape of the curve.
    res["two_exponential_fit"] = fit_two_exponentials(t_resp, resp)
    res["two_exponential_fit"]["model_rms_m"] = float(
        np.sqrt(np.mean((step_model(t_resp, dy, res["zeta"], res["wn_rad_s"]) - resp) ** 2))
    )

    zeta, wn = res["zeta"], res["wn_rad_s"]
    res["k"] = float(wn**2)
    res["b"] = float(2 * zeta * wn)
    res["poles"] = poles(zeta, wn)
    res["corner_frequency_rad_s"] = float(wn)

    # k arrives twice by different routes. They agree only if the input column is
    # force in the same units the model assumes, so the gap between them is the
    # scaling check, not a rounding artefact.
    res["k_from_wn"] = float(wn**2)
    res["k_from_gain"] = float(1.0 / res["gain_A"]) if res["gain_A"] else float("nan")
    if math.isfinite(res["k_from_gain"]) and res["k_from_wn"]:
        res["k_mismatch_percent"] = float(
            100 * (res["k_from_gain"] - res["k_from_wn"]) / res["k_from_wn"]
        )
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
    out = {
        "zeta": float(zeta),
        "wn_rad_s": float(wn),
        "fit_method": "least squares sweep",
    }

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
        print(
            "WARNING: overdamped fit hit the edge of the search grid: "
            + out["grid_warning"]
            + "\n  widen the ranges in fit_overdamped and rerun"
        )
    if reach.size:
        tau = float(t[reach[0]])
        out["tau_63_s"] = tau
        out["reduced_order_pole"] = -1.0 / tau
    return out


def fit_two_exponentials(t, y):
    """Best fit of c + A exp(-bt) + B exp(-ct) with b and c real.

    This is the two real pole step response. Fitting it is the direct way to
    answer whether the poles are real: if the record needs a complex pair, no
    choice of real b and c reaches it. A and B and the constant are linear once
    b and c are fixed, so only the two rates are searched.
    """
    rates = np.geomspace(0.02, 20.0, 240)
    best = None
    for i, b in enumerate(rates):
        for c in rates[i + 1 :]:
            basis = np.column_stack([np.ones_like(t), np.exp(-b * t), np.exp(-c * t)])
            coeff, *_ = np.linalg.lstsq(basis, y, rcond=None)
            fit = basis @ coeff
            err = float(np.sqrt(np.mean((fit - y) ** 2)))
            if best is None or err < best[0]:
                best = (err, b, c, fit)
    err, b, c, fit = best
    return {
        "rms_m": err,
        "rate_1": float(b),
        "rate_2": float(c),
        "peak_above_final_m": float(np.max(fit) - fit[-1]),
    }


def step_model(t, dy, zeta, wn):
    """Analytic unit-ish step response of A wn^2 / (s^2 + 2 zeta wn s + wn^2)."""
    if zeta < 1:
        wd = wn * math.sqrt(1 - zeta**2)
        phi = math.atan2(math.sqrt(1 - zeta**2), zeta)
        env = np.exp(-zeta * wn * t) / math.sqrt(1 - zeta**2)
        return dy * (1 - env * np.sin(wd * t + phi))
    if abs(zeta - 1) < 1e-9:
        return dy * (1 - np.exp(-wn * t) * (1 + wn * t))
    r = math.sqrt(zeta**2 - 1)
    s1, s2 = -wn * (zeta - r), -wn * (zeta + r)
    return dy * (1 + (s2 * np.exp(s1 * t) - s1 * np.exp(s2 * t)) / (s1 - s2))


def poles(zeta, wn):
    if zeta < 1:
        wd = wn * math.sqrt(1 - zeta**2)
        return [[-zeta * wn, wd], [-zeta * wn, -wd]]
    r = math.sqrt(zeta**2 - 1)
    return [[-wn * (zeta - r), 0.0], [-wn * (zeta + r), 0.0]]


def fit_sinusoid(t, y, w):
    """Least squares amplitude and phase of y at a known frequency."""
    basis = np.column_stack([np.cos(w * t), np.sin(w * t), np.ones_like(t)])
    coef, *_ = np.linalg.lstsq(basis, y, rcond=None)
    a, b, _ = coef
    return math.hypot(a, b), math.atan2(-b, a)


def sine_dir_logs(directory):
    """Sine logs in a folder, split into the sweep and the amplitude repeats.

    A name with a tag after the frequency, like sine_1p0_a3.csv, is the same
    frequency driven at another amplitude. It belongs to the linearity check
    rather than to the sweep, so it comes back separately instead of turning into
    a second point at a frequency already covered.
    """
    if not directory:
        return [], []
    names, tagged = [], []
    for n in sorted(os.listdir(directory)):
        low = n.lower()
        if not (low.startswith("sine_") and low.endswith(".csv")):
            continue
        (tagged if "_" in low[5:-4] else names).append(n)
    if not names:
        raise SystemExit(f"{directory}: no files named sine_<frequency>.csv")
    return (
        [os.path.join(directory, n) for n in names],
        [os.path.join(directory, n) for n in tagged],
    )


def frequency_from_name(path):
    """Drive frequency and amplitude from a name like sine_1p6.csv or sine_1p0_a3.csv.

    1p6 means 1.6 rad/s. A trailing _a3 means the run was driven at amplitude 3,
    which is the linearity repeat rather than a frequency point.
    """
    stem = os.path.splitext(os.path.basename(path))[0]
    text = stem.split("_", 1)[1]
    amp = 1.0
    if "_a" in text:
        text, _, tail = text.partition("_a")
        try:
            amp = float(tail.replace("p", "."))
        except ValueError:
            raise SystemExit(f"{path}: cannot read a drive amplitude out of the name")
    try:
        return float(text.replace("p", ".")), amp
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


def linearity_check(base, amp_step, points, amp_points):
    """What changes when the same test is driven at a different amplitude.

    A linear plant scales the output by whatever the input was scaled by and
    leaves every identified parameter alone, so the ratios below are the test and
    the parameter pairs are what it costs to be wrong about it.
    """
    out = {}
    if amp_step:
        out["step"] = {
            "input_ratio": amp_step["input_step"] / base["input_step"],
            "output_ratio": amp_step["output_change"] / base["output_change"],
            "gain_A": [base["gain_A"], amp_step["gain_A"]],
            "zeta": [base["zeta"], amp_step["zeta"]],
            "wn_rad_s": [base["wn_rad_s"], amp_step["wn_rad_s"]],
            "overshoot": [base["overshoot"], amp_step["overshoot"]],
            "first_peak_time_s": [
                base["first_peak_time_s"],
                amp_step["first_peak_time_s"],
            ],
        }
    sines = []
    for p in amp_points:
        near = [q for q in points if abs(q["w_rad_s"] - p["w_rad_s"]) < 1e-9]
        if not near:
            continue
        q = near[0]
        sines.append(
            {
                "w_rad_s": p["w_rad_s"],
                "input_ratio": p["input_amplitude"] / q["input_amplitude"],
                "output_ratio": p["output_amplitude"] / q["output_amplitude"],
                "magnitude": [q["magnitude"], p["magnitude"]],
                "phase_deg": [q["phase_deg"], p["phase_deg"]],
            }
        )
    if sines:
        out["sine"] = sines
    return out


def model_response(w, gain, zeta, wn):
    s = 1j * w
    g = gain * wn**2 / (s**2 + 2 * zeta * wn * s + wn**2)
    return 20 * np.log10(np.abs(g)), np.degrees(np.unwrap(np.angle(g)))


def freq_residual(p, w, mag_db, phase_deg):
    """Magnitude and phase residuals of a delayed second order model, stacked."""
    gain, zeta, wn, delay = p
    mag, ph = model_response(w, gain, zeta, wn)
    return np.concatenate([mag - mag_db, ph - np.degrees(w * delay) - phase_deg])


def fit_frequency_model(w, mag_db, phase_deg, guess, fixed_delay=None):
    """Least squares gain, zeta, wn and delay fitted to the frequency points alone.

    Magnitude in dB and phase in degrees carry equal weight, since the two
    residuals scatter by about the same amount on this data. Gauss-Newton with a
    numerical Jacobian; the model is smooth in all four parameters and a step
    estimate starts close enough to converge. The guess is the only thing this
    borrows from the step test, so the answer is independent of it.

    With fixed_delay set the delay is held there and only the three plant
    parameters move, which shows what the fit has to do to the plant to cover a
    delay it is not allowed to name.
    """
    p = np.array(guess, float)
    free = 3 if fixed_delay is not None else len(p)
    if fixed_delay is not None:
        p[3] = fixed_delay
    for _ in range(200):
        r = freq_residual(p, w, mag_db, phase_deg)
        cols = []
        for i in range(free):
            h = 1e-6 * max(abs(p[i]), 1e-3)
            q = p.copy()
            q[i] += h
            cols.append((freq_residual(q, w, mag_db, phase_deg) - r) / h)
        delta, *_ = np.linalg.lstsq(np.column_stack(cols), -r, rcond=None)
        p[:free] += delta
        if np.max(np.abs(delta / np.maximum(np.abs(p[:free]), 1e-9))) < 1e-12:
            break
    gain, zeta, wn, delay = p
    mag, ph = model_response(w, gain, zeta, wn)
    return {
        "gain_A": float(gain),
        "zeta": float(zeta),
        "wn_rad_s": float(wn),
        "delay_s": float(delay),
        "k": float(wn**2),
        "b": float(2 * zeta * wn),
        "alpha": float(gain * wn**2),
        "mag_rms_dB": float(np.sqrt(np.mean((mag - mag_db) ** 2))),
        "phase_rms_deg": float(
            np.sqrt(np.mean((ph - np.degrees(w * delay) - phase_deg) ** 2))
        ),
    }


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
    denom = math.hypot(1 - ratio**2, 2 * res["zeta"] * ratio)
    out["gain_expected_bias_percent"] = float(100 * (1 / denom - 1)) if denom else 0.0
    if ratio > 0.3:
        out["gain_note"] = (
            f"lowest frequency is {ratio:.2f} wn, too close to the corner "
            "for a clean asymptote, so this reads high"
        )

    # The lowest point is not the DC gain: the curve has already begun to rise towards
    # resonance there, by the amount just computed. Dividing that lift back out gives
    # what the run would have read at zero frequency, and Mr has to be measured against
    # that rather than against the lowest point, or the same bias lands in zeta.
    a_dc = float(a_freq * denom) if denom else a_freq
    out["gain_A_at_dc"] = a_dc
    out["gain_dc_mismatch_percent"] = float(
        100 * (a_dc - res["gain_A"]) / res["gain_A"]
    )

    # Questions 2 and 3 again: only a complex pair lifts the magnitude above the asymptote
    i_peak = int(np.argmax(mag_db))
    interior = 0 < i_peak < len(w) - 1
    mr = float(mag[i_peak] / a_dc)
    out.update(
        {"peak_w_rad_s": float(w[i_peak]), "Mr": mr, "peak_is_interior": interior}
    )
    if interior and mr > 1.0:
        q = 1.0 / (2 * mr)
        disc = 1 - 4 * q**2
        if disc >= 0:
            z_sq = (1 - math.sqrt(disc)) / 2
            out["zeta_from_Mr"] = float(math.sqrt(z_sq))
            if 1 - 2 * z_sq > 0:
                out["wn_from_peak"] = float(w[i_peak] / math.sqrt(1 - 2 * z_sq))
    else:
        out["resonance_note"] = (
            "no interior magnitude peak, so the poles are real or "
            "the resonance sits outside the frequencies tested"
        )

    # Question 4 again: the phase crosses -90 degrees at wn whatever zeta is
    for i in range(len(w) - 1):
        if (ph[i] + 90) * (ph[i + 1] + 90) <= 0 and ph[i] != ph[i + 1]:
            f = (-90 - ph[i]) / (ph[i + 1] - ph[i])
            lo, hi = math.log10(w[i]), math.log10(w[i + 1])
            out["wn_from_phase"] = float(10 ** (lo + f * (hi - lo)))
            break
    else:
        out["phase_note"] = (
            "the phase never reaches -90 degrees, test a higher frequency"
        )

    # Two poles and no zeros roll off at -40 dB/decade, but only well above wn. Just
    # past resonance the curve falls far steeper than that, so a slope measured from
    # points at 2 wn and 4 wn reads nearer -60 and is not evidence against the model.
    hf = w >= 2.0 * wn_step
    out["highest_w_over_wn"] = float(w[-1] / wn_step)
    if int(hf.sum()) >= 2:
        out["hf_slope_dB_per_decade"] = float(
            np.polyfit(np.log10(w[hf]), mag_db[hf], 1)[0]
        )
        out["hf_points_used"] = int(hf.sum())
        if out["highest_w_over_wn"] < 4:
            out["slope_note"] = (
                f"highest frequency is only {out['highest_w_over_wn']:.1f} wn, "
                "so the -40 dB/decade asymptote is not reached yet and this "
                "slope reads steeper. Do not call that a failed model"
            )
    else:
        out["slope_note"] = "fewer than two points above 2 wn, roll off not measurable"

    # validation in the frequency domain, against a model built from the step alone
    mmag, mph = model_response(w, res["gain_A"], res["zeta"], wn_step)
    mph = np.where(mph > 10, mph - 360, mph)
    out["model_mag_rms_dB"] = float(np.sqrt(np.mean((mag_db - mmag) ** 2)))
    out["model_phase_rms_deg"] = float(np.sqrt(np.mean((ph - mph) ** 2)))

    # the same four parameters read out of these points alone. Starting from the step
    # answer only sets where the search begins; the result owes it nothing else.
    if len(w) >= 4:
        guess = [res["gain_A"], res["zeta"], wn_step, 0.1]
        out["fit"] = fit_frequency_model(w, mag_db, ph, guess)
        # The gain the full fit returns is pulled up by the points past the corner,
        # where the second order model is already leaving the measurement. Repeating
        # it below 2 rad/s is the reading the report quotes for A and alpha.
        low = w <= 2.0
        if low.sum() >= 4:
            out["fit_below_2"] = fit_frequency_model(
                w[low], mag_db[low], ph[low], guess
            )
            out["fit_below_2"]["points_used"] = int(low.sum())
    return out


def refit_at_delays(points, res, delays):
    """Refit the plant to the frequency points with the delay held at each value.

    Held at zero the fit has to bend zeta and wn to cover the missing lag, so the
    spread across these refits is a direct measure of what leaving the delay out
    costs the plant parameters.
    """
    pts = sorted(points, key=lambda p: p["w_rad_s"])
    w = np.array([p["w_rad_s"] for p in pts])
    mag_db = np.array([p["magnitude_dB"] for p in pts])
    ph = np.array([p["phase_deg"] for p in pts])
    guess = [res["gain_A"], res["zeta"], res["wn_rad_s"], 0.0]
    return [fit_frequency_model(w, mag_db, ph, guess, fixed_delay=d) for d in delays]


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
    ax.annotate(
        f"final value = {res['output_final']:.3f} m",
        xy=(t[-1], res["output_final"]),
        xytext=(-8, -16),
        textcoords="offset points",
        ha="right",
        **ANNOT,
    )
    # the rise itself is the gain, so mark it the way the brief's example does
    ax.axhline(res["output_baseline"], color="#b3261e", ls=":", lw=1.0)
    xg = t[i0] - 0.45 * (t[i0] - t[0])
    ax.annotate(
        "",
        xy=(xg, res["output_final"]),
        xytext=(xg, res["output_baseline"]),
        arrowprops=dict(arrowstyle="<->", color="#b3261e", lw=1.4),
    )
    ax.annotate(
        f"$A$ = {res['gain_A']:.3f} m/V",
        xy=(xg, res["output_final"]),
        xytext=(0, 14),
        textcoords="offset points",
        ha="center",
        va="bottom",
        **ANNOT,
    )
    ax.axvline(t[i0], color="#4a4a4a", ls=":", lw=1.2)
    ax.annotate(
        "step applied",
        xy=(t[i0], y.min()),
        xytext=(6, 10),
        textcoords="offset points",
        color="#4a4a4a",
        fontsize=11,
    )

    if res["damping_case"] == "underdamped":
        tp = res["step_time_s"] + res["first_peak_time_s"]
        peak = res["output_final"] + res["overshoot"] * res["output_change"]
        rise = res["output_change"]
        ax.annotate(
            "",
            xy=(tp, peak),
            xytext=(tp, res["output_final"]),
            arrowprops=dict(arrowstyle="<->", color="#b3261e", lw=1.4),
        )
        ax.annotate(
            f"overshoot = {100 * res['overshoot']:.1f}%",
            xy=(tp, peak),
            xytext=(12, 12),
            textcoords="offset points",
            va="bottom",
            ha="left",
            **ANNOT,
        )

        # time to the first peak, which is where the damped period comes from
        lvl = res["output_final"] - 0.28 * rise
        ax.annotate(
            "",
            xy=(res["step_time_s"], lvl),
            xytext=(tp, lvl),
            arrowprops=dict(arrowstyle="<->", color="#b3261e", lw=1.4),
        )
        # label at the right end of the bar, since a centred one lands on the rising edge
        ax.annotate(
            f"$T_p$ = {res['first_peak_time_s']:.2f} s",
            xy=(tp, lvl),
            xytext=(10, -5),
            textcoords="offset points",
            ha="left",
            **ANNOT,
        )

        # the undershoot is only millimetres deep, so it gets a pointer rather than a bar
        if "extrema_deviations_m" in res:
            tu = res["step_time_s"] + res["extrema_times_s"][1]
            depth = res["extrema_deviations_m"][1]
            ax.annotate(
                f"undershoot = {1000 * depth:.0f} mm,\n"
                f"decay ratio gives $\\zeta$ = {res['zeta_log_decrement']:.3f}",
                xy=(tu, res["output_final"] - depth),
                xytext=(18, -96),
                textcoords="offset points",
                arrowprops=dict(arrowstyle="->", color="#b3261e", lw=1.2),
                **ANNOT,
            )

        if len(res.get("peak_times_s", [])) >= 2:
            t1 = res["step_time_s"] + res["peak_times_s"][0]
            t2 = res["step_time_s"] + res["peak_times_s"][1]
            top = res["output_final"] + 1.35 * res["overshoot"] * rise
            ax.annotate(
                "",
                xy=(t1, top),
                xytext=(t2, top),
                arrowprops=dict(arrowstyle="<->", color="#b3261e", lw=1.4),
            )
            ax.annotate(
                f"$T_d$ = {res['damped_period_s']:.2f} s",
                xy=((t1 + t2) / 2, top),
                xytext=(0, 6),
                textcoords="offset points",
                ha="center",
                **ANNOT,
            )

    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Displacement [m]")
    ax.margins(y=0.16)
    ax.set_title("Measured step response")
    ax.grid(alpha=0.3)
    ax.legend(loc="lower right")
    save(fig, outdir, "step_response.png")


def fit_delay(points, res):
    """Least squares transport delay from the phase the model does not account for.

    A second order system cannot lag past 180 degrees, so measurements beyond that
    need something else. An extra pole adds lag that saturates; a delay adds lag that
    grows without bound and leaves the magnitude alone. Fitting -wT to the residual
    tells the two apart and gives the delay if that is what it is.
    """
    w = np.array([p["w_rad_s"] for p in points])
    ph = np.array([p["phase_deg"] for p in points])
    _, model_ph = model_response(w, res["gain_A"], res["zeta"], res["wn_rad_s"])
    excess = np.radians(ph - model_ph)
    delay = float(-np.sum(excess * w) / np.sum(w * w))
    return {
        "delay_s": delay,
        "phase_rms_no_delay_deg": float(np.sqrt(np.mean((ph - model_ph) ** 2))),
        "phase_rms_with_delay_deg": float(
            np.sqrt(np.mean((ph - model_ph + np.degrees(w * delay)) ** 2))
        ),
    }


def plot_bode(
    points,
    res,
    outdir,
    name="bode_plot.png",
    title="Measured frequency response",
    delay=None,
):
    w = np.array([p["w_rad_s"] for p in points])
    mag = np.array([p["magnitude_dB"] for p in points])
    ph = np.array([p["phase_deg"] for p in points])
    grid = np.geomspace(
        min(w.min(), res["wn_rad_s"]) / 5, max(w.max(), res["wn_rad_s"]) * 5, 400
    )
    mmag, mph = model_response(grid, res["gain_A"], res["zeta"], res["wn_rad_s"])

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(7.5, 6.4), sharex=True)
    ax1.semilogx(grid, mmag, color="#4a4a4a", lw=1.6, label="identified model")
    ax1.semilogx(w, mag, "o", color="#1f3b73", ms=7, label="measured")
    ax1.axvline(res["wn_rad_s"], color="#b3261e", ls="--", lw=1.2)
    ax1.annotate(
        f"$\\omega_n$ = {res['wn_rad_s']:.3f} rad/s",
        xy=(res["wn_rad_s"], mmag.min() + 0.15 * (mmag.max() - mmag.min())),
        xytext=(8, 0),
        textcoords="offset points",
        **ANNOT,
    )

    # the resonant peak is a second reading of zeta, so mark where it was measured
    ipk = int(np.argmax(mag))
    ax1.annotate(
        f"resonant peak, {10 ** (mag[ipk] / 20):.3f} at {w[ipk]:.2f} rad/s",
        xy=(w[ipk], mag[ipk]),
        xytext=(0, -52),
        textcoords="offset points",
        ha="center",
        arrowprops=dict(arrowstyle="->", color="#b3261e", lw=1.2),
        **ANNOT,
    )
    # the low frequency asymptote is the gain, which is where the frequency test reads
    # it a second time, so mark the level the model settles onto
    a_db = 20 * math.log10(res["gain_A"])
    ax1.axhline(a_db, color="#b3261e", ls=":", lw=1.0)
    ax1.annotate(
        f"$A$ = {res['gain_A']:.3f} m/V",
        xy=(grid.max(), a_db),
        xytext=(-4, -14),
        textcoords="offset points",
        ha="right",
        **ANNOT,
    )
    ax1.set_ylabel("Magnitude [dB]")
    ax1.grid(which="both", alpha=0.3)
    ax1.legend(loc="lower left")
    ax1.set_title(title)

    ax2.semilogx(grid, mph, color="#4a4a4a", lw=1.6, label="identified model")
    if delay:
        ax2.semilogx(
            grid,
            mph - np.degrees(grid * delay["delay_s"]),
            color="#1a7f4b",
            lw=1.6,
            ls="--",
            label=f"with {delay['delay_s'] * 1000:.0f} ms delay",
        )
    ax2.semilogx(w, ph, "o", color="#1f3b73", ms=7, label="measured")
    ax2.axvline(res["wn_rad_s"], color="#b3261e", ls="--", lw=1.2)
    ax2.axhline(-90, color="#b3261e", ls=":", lw=1.2)
    # both reference labels start at the left edge, ahead of where either curve drops
    ax2.annotate(
        "$-90^\\circ$ at $\\omega_n$",
        xy=(grid.min(), -90),
        xytext=(4, 6),
        textcoords="offset points",
        **ANNOT,
    )
    ax2.axhline(-180, color="#8a8a8a", ls=":", lw=1.0)
    ax2.annotate(
        "$-180^\\circ$, the limit for two poles",
        xy=(grid.min(), -180),
        xytext=(4, 6),
        textcoords="offset points",
        color="#5a5a5a",
        fontsize=11,
    )
    ax2.set_xlabel("Frequency [rad/s]")
    ax2.set_ylabel("Phase [deg]")
    ax2.grid(which="both", alpha=0.3)
    ax2.legend(loc="lower left", fontsize=10)
    save(fig, outdir, name)


def step_response_model(t_resp, res, delay=0.0):
    return res["output_baseline"] + step_model(
        np.clip(t_resp - delay, 0, None),
        res["output_change"],
        res["zeta"],
        res["wn_rad_s"],
    )


def fit_step_delay(t, u, y, res, span=0.3, grid=0.001):
    """The pure time shift of the model that fits the step record best.

    Same question the phase residual answers in the frequency domain, asked of the
    step record instead, so the two answers can be compared. The shift is read on a
    1 ms grid, far finer than the sampling, but the step onset is only known to the
    sample it was logged in, so the answer carries a full sample period of doubt.
    """
    i0 = step_index(u)
    t_resp = t[i0:] - t[i0]
    meas = y[i0:]
    best = None
    for delay in np.arange(0.0, span + grid, grid):
        r = meas - step_response_model(t_resp, res, delay)
        rms = float(np.sqrt(np.mean(r**2)))
        if best is None or rms < best[1]:
            best = (float(delay), rms, float(np.max(np.abs(r))))
    return {"delay_s": best[0], "rms_m": best[1], "peak_residual_m": best[2]}


def plot_validation(t, u, y, res, outdir, delay=None):
    i0 = step_index(u)
    t_resp = t[i0:] - t[i0]
    model = step_response_model(t_resp, res)
    fig, (ax, axr) = plt.subplots(
        2, 1, figsize=(7.5, 6.2), gridspec_kw={"height_ratios": [2.1, 1]}
    )
    # both panels count from the step, so a feature in the residual sits under the
    # part of the response that produced it
    ax.plot(t_resp, y[i0:], color="#1f3b73", lw=2, label="measured")
    ax.plot(t_resp, model, color="#b3261e", lw=1.8, ls="--", label="identified model")
    rms = float(np.sqrt(np.mean((y[i0:] - model) ** 2)))
    # clear of the legend, which sits in the lower right corner
    ax.annotate(
        f"RMS error = {rms:.4f} m",
        xy=(0.98, 0.30),
        xycoords="axes fraction",
        ha="right",
        **ANNOT,
    )
    ax.set_xlabel("")
    ax.set_ylabel("Displacement [m]")
    ax.set_title("Simulated against measured step response")
    ax.grid(alpha=0.3)
    ax.legend(loc="lower right")

    # the two curves above sit on top of each other at this scale, so the residual
    # gets its own axes over the window where it is larger than the logger can resolve
    resid = 1000 * (y[i0:] - model)
    step_mm = 1000 * res["output_resolution_m"]
    window = 1.1 * res["damped_period_s"]
    ax.axvspan(0, window, color="#1f3b73", alpha=0.06, zorder=0)
    axr.axhspan(
        -step_mm,
        step_mm,
        color="#8a8a8a",
        alpha=0.25,
        label=f"$\\pm{step_mm:.0f}$ mm logger resolution",
    )
    axr.axhline(0, color="#8a8a8a", lw=0.8)
    axr.plot(t_resp, resid, color="#1f3b73", lw=2, label="measured $-$ model")
    i_peak = int(np.argmax(np.abs(resid)))
    axr.annotate(
        f"{resid[i_peak]:.1f} mm at {t_resp[i_peak]:.1f} s",
        xy=(t_resp[i_peak], resid[i_peak]),
        xytext=(10, 4),
        textcoords="offset points",
        **ANNOT,
    )
    if delay is not None:
        shifted = 1000 * (y[i0:] - step_response_model(t_resp, res, delay))
        axr.plot(
            t_resp,
            shifted,
            color="#1a7f4b",
            lw=1.6,
            ls="--",
            label=f"model shifted {1000 * delay:.0f} ms later",
        )
    axr.set_xlim(0, window)
    axr.set_xlabel("Time since the step [s]")
    axr.set_ylabel("Residual [mm]")
    axr.grid(alpha=0.3)
    axr.legend(loc="lower right", fontsize=10)
    save(fig, outdir, "validation_step.png")

    # the sample after the last one out of the band, or the end of the record when the
    # residual never settles inside it
    outside = np.flatnonzero(np.abs(resid) > step_mm)
    recovery = (
        float(t_resp[min(outside[-1] + 1, len(t_resp) - 1)]) if outside.size else 0.0
    )
    return {
        "rms_m": rms,
        "peak_residual_mm": float(resid[i_peak]),
        "peak_residual_time_s": float(t_resp[i_peak]),
        "residual_within_resolution_s": recovery,
    }


def plot_frequency_validation(points, res, outdir, delay=None):
    """Residual of the identified model against every measured frequency point.

    The Bode plot shows agreement. This shows what is left over, which is where the
    magnitude runs out of signal and the phase runs past what two poles can do.
    """
    pts = sorted(points, key=lambda p: p["w_rad_s"])
    w = np.array([p["w_rad_s"] for p in pts])
    mag = np.array([p["magnitude_dB"] for p in pts])
    ph = np.array([p["phase_deg"] for p in pts])
    steps = np.array([p["output_amplitude"] for p in pts]) / res["output_resolution_m"]
    mmag, mph = model_response(w, res["gain_A"], res["zeta"], res["wn_rad_s"])
    dmag, dph = mag - mmag, ph - mph
    mag_rms = float(np.sqrt(np.mean(dmag**2)))

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(7.5, 6.4), sharex=True)
    ax1.axhspan(
        -mag_rms, mag_rms, color="#8a8a8a", alpha=0.25, label=f"${mag_rms:.2f}$ dB RMS"
    )
    ax1.axhline(0, color="#8a8a8a", lw=0.8)
    ax1.semilogx(
        w, dmag, "o-", color="#1f3b73", ms=7, lw=1.4, label="measured $-$ model"
    )
    # the top of the sweep is where the output stops being resolvable, not where the
    # model starts failing, so say how many resolution steps each of those points has
    for x, d, n in zip(w, dmag, steps):
        if n < 30:
            ax1.annotate(
                f"{n:.0f} steps",
                xy=(x, d),
                xytext=(-9, 2),
                textcoords="offset points",
                ha="right",
                va="center",
                **ANNOT,
            )
    ax1.margins(y=0.2)
    ax1.set_ylabel("Magnitude error [dB]")
    ax1.grid(which="both", alpha=0.3)
    ax1.legend(loc="upper left", fontsize=10)
    ax1.set_title("Identified model against the measured frequency points")

    ax2.axhline(0, color="#8a8a8a", lw=0.8)
    ax2.semilogx(
        w, dph, "o-", color="#1f3b73", ms=7, lw=1.4, label="measured $-$ model"
    )
    out = {
        "mag_rms_dB": mag_rms,
        "phase_rms_deg": float(np.sqrt(np.mean(dph**2))),
        "mag_error_dB": [float(v) for v in dmag],
        "phase_error_deg": [float(v) for v in dph],
        "output_amplitude_steps": [float(v) for v in steps],
    }
    if delay is not None:
        grid = np.geomspace(w.min(), w.max(), 200)
        ax2.semilogx(
            grid,
            -np.degrees(grid * delay),
            color="#b3261e",
            ls=":",
            lw=1.4,
            label=f"$-\\omega T$ at $T = {1000 * delay:.0f}$ ms",
        )
        left = dph + np.degrees(w * delay)
        ax2.semilogx(
            w,
            left,
            "s--",
            color="#1a7f4b",
            ms=6,
            lw=1.4,
            label="with the delay included",
        )
        out["phase_rms_with_delay_deg"] = float(np.sqrt(np.mean(left**2)))
    ax2.set_xlabel("Frequency [rad/s]")
    ax2.set_ylabel("Phase error [deg]")
    ax2.grid(which="both", alpha=0.3)
    ax2.legend(loc="lower left", fontsize=10)
    save(fig, outdir, "validation_freq.png")
    return out


def save(fig, outdir, name):
    path = os.path.join(outdir, name)
    fig.savefig(path)
    plt.close(fig)
    print(f"wrote {path}")


def summarise(res, points, val, freq=None, repeats=None, lin=None):
    zeta, wn = res["zeta"], res["wn_rad_s"]
    rms = val["rms_m"] if val else None
    print("\nIdentified model")
    print(
        f"  gain A          {res['gain_A']:.4f} m/V   ->  k = 1/A = {1 / res['gain_A']:.4f}"
    )
    print(f"  damping case    {res['damping_case']}")
    if "two_exponential_fit" in res:
        te = res["two_exponential_fit"]
        print(
            f"  two real poles  best fit {te['rms_m'] * 1000:.1f} mm RMS against "
            f"{te['model_rms_m'] * 1000:.1f} mm for the complex pair, and it clears "
            f"the final value by only {te['peak_above_final_m'] * 1000:.2f} mm"
        )
    print(f"  zeta            {zeta:.4f}")
    if "zeta_log_decrement" in res:
        d = res["extrema_deviations_m"]
        print(
            f"  zeta (log dec)  {res['zeta_log_decrement']:.4f}"
            f"  from extrema {1000 * d[0]:.0f} and {1000 * d[1]:.0f} mm off final"
        )
    print(f"  wn              {wn:.4f} rad/s")
    print(f"  k = wn^2        {res['k']:.4f}")
    print(f"  b = 2 zeta wn   {res['b']:.4f}")
    p = res["poles"]
    print(
        f"  poles           {p[0][0]:.4f} + {p[0][1]:.4f}j, {p[1][0]:.4f} + {p[1][1]:.4f}j"
    )
    print(f"  G(s) = 1 / (s^2 + {res['b']:.4f} s + {res['k']:.4f})")
    if rms is not None:
        print(f"  step RMS error  {rms:.5f} m")
        print(
            f"  worst residual  {val['peak_residual_mm']:.1f} mm at "
            f"{val['peak_residual_time_s']:.1f} s, inside the logger resolution "
            f"from {val['residual_within_resolution_s']:.1f} s on"
        )
    if val and "step_delay" in val:
        d = val["step_delay"]
        print(
            f"  best time shift {1000 * d['delay_s']:.0f} ms, which takes the step "
            f"RMS to {d['rms_m'] * 1000:.2f} mm and the worst point to "
            f"{d['peak_residual_m'] * 1000:.2f} mm"
        )

    if "k_mismatch_percent" in res:
        print("\nGain consistency check")
        print(f"  k from the final value   {res['k_from_gain']:.4f}")
        print(f"  k from wn^2              {res['k_from_wn']:.4f}")
        print(f"  difference               {res['k_mismatch_percent']:+.1f}%")
        if abs(res["k_mismatch_percent"]) > 10:
            print(
                "  the two disagree by more than 10%, so the input column is probably"
                "\n  not force in the units the model assumes. Say so in the report."
            )

    if repeats and len(repeats) > 1:
        print("\nRepeatability across step runs")
        print("  file            A         zeta      wn        k         b")
        for name, r in repeats:
            print(
                f"  {name:14s}  {r['gain_A']:8.4f}  {r['zeta']:8.4f}  "
                f"{r['wn_rad_s']:8.4f}  {r['k']:8.4f}  {r['b']:8.4f}"
            )
        for key, label in (
            ("gain_A", "A"),
            ("zeta", "zeta"),
            ("wn_rad_s", "wn"),
            ("k", "k"),
            ("b", "b"),
        ):
            vals = np.array([r[key] for _, r in repeats])
            spread = 100 * (vals.max() - vals.min()) / abs(vals.mean())
            print(
                f"  {label:5s} mean {vals.mean():8.4f}   spread {spread:5.1f}% of mean"
            )

    if points:
        print("\nFrequency response points")
        print("  w [rad/s]   |G| [dB]   phase [deg]   file")
        for pt in points:
            print(
                f"  {pt['w_rad_s']:9.3f}  {pt['magnitude_dB']:9.2f}  {pt['phase_deg']:11.1f}   {pt['file']}"
            )

    if freq:
        print("\nRead a second time off the frequency response")
        if "gain_A_from_asymptote" in freq:
            print(
                f"  gain A at the lowest w   {freq['gain_A_from_asymptote']:.4f} m/V "
                f"({freq['gain_mismatch_percent']:+.1f}% against the step)"
            )
        if "gain_A_at_dc" in freq:
            print(
                f"  gain A with lift removed {freq['gain_A_at_dc']:.4f} m/V "
                f"({freq['gain_dc_mismatch_percent']:+.2f}% against the step)"
            )
        if "zeta_from_Mr" in freq:
            print(
                f"  zeta from the peak Mr    {freq['zeta_from_Mr']:.4f} "
                f"(step gave {zeta:.4f})"
            )
        if "wn_from_peak" in freq:
            print(f"  wn from the peak         {freq['wn_from_peak']:.4f} rad/s")
        if "wn_from_phase" in freq:
            print(
                f"  wn from the -90 crossing {freq['wn_from_phase']:.4f} rad/s "
                f"({100 * (freq['wn_from_phase'] - wn) / wn:+.1f}% against the step)"
            )
        if "hf_slope_dB_per_decade" in freq:
            print(
                f"  high frequency slope     {freq['hf_slope_dB_per_decade']:.1f} dB/decade "
                f"over {freq['hf_points_used']} points, expect -40 for two poles"
            )
        print(
            f"  model vs measured        {freq['model_mag_rms_dB']:.2f} dB RMS, "
            f"{freq['model_phase_rms_deg']:.1f} deg RMS"
        )
        if "fit" in freq:
            f = freq["fit"]
            print("  fitted to these points alone, using no step data:")
            print(
                f"    A    {f['gain_A']:.4f} m/V   ({100 * (f['gain_A'] - res['gain_A']) / res['gain_A']:+.1f}% against the step)"
            )
            print(
                f"    zeta {f['zeta']:.4f}         ({100 * (f['zeta'] - zeta) / zeta:+.1f}%)"
            )
            print(
                f"    wn   {f['wn_rad_s']:.4f} rad/s ({100 * (f['wn_rad_s'] - wn) / wn:+.1f}%)"
            )
            print(f"    T    {1000 * f['delay_s']:.1f} ms")
            print(f"    alpha {f['alpha']:.4f}")
        if "fit_below_2" in freq:
            f = freq["fit_below_2"]
            print(
                f"  the same fit over the {f['points_used']} points at or below 2 rad/s:"
            )
            print(
                f"    A    {f['gain_A']:.4f} m/V   ({100 * (f['gain_A'] - res['gain_A']) / res['gain_A']:+.1f}% against the step)"
            )
            print(
                f"    zeta {f['zeta']:.4f}         wn {f['wn_rad_s']:.4f} rad/s   "
                f"alpha {f['alpha']:.4f}"
            )
        if "delay_s" in freq:
            print(
                f"  fitted transport delay   {freq['delay_s'] * 1000:.1f} ms, "
                f"phase {freq['phase_rms_no_delay_deg']:.1f} -> "
                f"{freq['phase_rms_with_delay_deg']:.1f} deg RMS"
            )
        if freq.get("fixed_delay_refits"):
            print("  refitted with the delay held fixed:")
            print("    T [ms]    zeta        wn        mag RMS    phase RMS")
            for f in freq["fixed_delay_refits"]:
                print(
                    f"    {1000 * f['delay_s']:6.0f}    {f['zeta']:.4f}    "
                    f"{f['wn_rad_s']:.4f}    {f['mag_rms_dB']:6.2f} dB   "
                    f"{f['phase_rms_deg']:6.2f} deg"
                )
        for k in ("gain_note", "resonance_note", "phase_note", "slope_note", "note"):
            if k in freq:
                print(f"  note: {freq[k]}")

    if lin:
        print("\nLinearity, against runs at another drive amplitude")
        if "step" in lin:
            s = lin["step"]
            print(
                f"  step   input x{s['input_ratio']:.2f}  output x{s['output_ratio']:.4f}  "
                f"A {s['gain_A'][0]:.4f} -> {s['gain_A'][1]:.4f}  "
                f"zeta {s['zeta'][0]:.4f} -> {s['zeta'][1]:.4f}  "
                f"wn {s['wn_rad_s'][0]:.4f} -> {s['wn_rad_s'][1]:.4f}"
            )
            print(
                f"         first peak at {s['first_peak_time_s'][0]:.1f} s -> "
                f"{s['first_peak_time_s'][1]:.1f} s"
            )
        for s in lin.get("sine", []):
            print(
                f"  sine   {s['w_rad_s']:.2f} rad/s  input x{s['input_ratio']:.2f}  "
                f"output x{s['output_ratio']:.4f}  "
                f"|G| {s['magnitude'][0]:.4f} -> {s['magnitude'][1]:.4f}  "
                f"phase {s['phase_deg'][0]:.2f} -> {s['phase_deg'][1]:.2f} deg"
            )


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
    print(
        f"selftest overdamped: true k = {k_od}, b = {b_od}, "
        f"found k = {res_od['k']:.4f}, b = {res_od['b']:.4f}"
    )
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
    for header in (
        ["Time (s)", "Input", "Output Displacement"],
        ["time", "INPUT", "output_displacement"],
    ):
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
    delay = fit_delay(points, res)
    plot_bode(points, res, outdir, delay=delay)
    plot_frequency_validation(points, res, outdir, delay=delay["delay_s"])
    freq.update(delay)
    step_delay = fit_step_delay(t, u, y, res)
    val = plot_validation(t, u, y, res, outdir, delay=step_delay["delay_s"])
    val["step_delay"] = step_delay
    assert step_delay["delay_s"] < 0.02, step_delay
    summarise(res, points, val, freq, [("selftest", res)])
    print("\nselftest passed, all figures written from synthetic data")


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument(
        "step_csv",
        nargs="*",
        help="step test logs from the simulator. Pass all the repeat runs "
        "to get a spread; the first one is the source of the figures",
    )
    ap.add_argument(
        "--sine",
        action="append",
        default=[],
        metavar="W=FILE",
        help="sine test log and the frequency in rad/s it was run at",
    )
    ap.add_argument(
        "--sine-dir",
        metavar="DIR",
        help="folder of sine logs named sine_<frequency>.csv, "
        "where 1p6 means 1.6 rad/s",
    )
    ap.add_argument(
        "--amp-step",
        metavar="FILE",
        help="step log taken at a different input amplitude, for the linearity check",
    )
    ap.add_argument("--outdir", default="report/figures", help="where the figures go")
    ap.add_argument(
        "--skip",
        type=float,
        default=0.5,
        help="fraction of each sine log to discard as transient",
    )
    ap.add_argument(
        "--selftest",
        action="store_true",
        help="check the analysis against a plant with known k and b",
    )
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
    sweep, tagged = sine_dir_logs(args.sine_dir)
    for path in sweep:
        points.append(analyse_sine(path, frequency_from_name(path)[0], args.skip))
    points.sort(key=lambda p: p["w_rad_s"])
    amp_points = []
    for path in tagged:
        w, amp = frequency_from_name(path)
        print(
            f"{os.path.basename(path)} is the {w:g} rad/s run at drive amplitude "
            f"{amp:g}, held out of the sweep and used for the linearity check"
        )
        amp_points.append(analyse_sine(path, w, args.skip))

    plot_step_input(t, u, args.outdir)
    plot_step_response(t, u, y, res, args.outdir)
    step_delay = fit_step_delay(t, u, y, res)
    val = plot_validation(t, u, y, res, args.outdir, delay=step_delay["delay_s"])
    val["step_delay"] = step_delay
    freq, delay = {}, {}
    if points:
        delay = fit_delay(points, res)
        plot_bode(points, res, args.outdir, delay=delay)
        freq = frequency_summary(points, res)
        freq.update(delay)
        freq["validation"] = plot_frequency_validation(
            points, res, args.outdir, delay=delay["delay_s"]
        )
        freq["fixed_delay_refits"] = refit_at_delays(
            points, res, [0.0, step_delay["delay_s"], delay["delay_s"]]
        )

    amp_step = analyse_step(*read_log(args.amp_step)) if args.amp_step else None
    lin = linearity_check(res, amp_step, points, amp_points)

    summarise(res, points, val, freq, repeats, lin)
    out = {
        "step": res,
        "frequency_points": points,
        "step_rms_error": val["rms_m"],
        "step_validation": val,
        "frequency_summary": freq,
        "step_repeats": {name: r for name, r in repeats},
        "linearity": lin,
    }
    with open(os.path.join(args.outdir, "results.json"), "w") as fh:
        json.dump(out, fh, indent=2)
    print(f"wrote {os.path.join(args.outdir, 'results.json')}")


if __name__ == "__main__":
    main()
