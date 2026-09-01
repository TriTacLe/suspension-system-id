# Suspension system identification

Recovers the transfer function of an unknown spring, mass and damper suspension from two
bench measurements: one step response and a set of single frequency sinusoids. Written for
the F-34B suspension simulator used in EEE3094S Control Systems Engineering at UCT, where
the plant parameters are generated per student and cannot be looked up.

The plant is

```
m x'' = -b x' - k x + F,      m = 1      ->      G(s) = 1 / (s^2 + b s + k)
```

so identifying it means measuring the steady state gain, the damping ratio and the natural
frequency, then converting those to k and b.

## What the analysis does

From the step log it reads the gain off the final value, the damping ratio off the first
overshoot with the logarithmic decrement as a cross-check, and the damped period off the
mean spacing of every peak in the record. If the response shows no overshoot it sweeps
zeta and the natural frequency for the least squares best fit to the analytic overdamped
step response instead, and reports the reduced first order model from the 63.2 percent
point.

From the sine logs it reads the same quantities a second time and independently: the gain
from the low frequency end, the damping ratio from the height of the resonant peak, and the
natural frequency from where the phase crosses -90 degrees. Amplitude and phase at each
drive frequency come from a least squares fit against a sinusoid at the known frequency,
which works on records far too coarsely sampled for peak picking.

It then validates the model in both domains, as RMS error against the step record and as
RMS magnitude and phase error against the measured frequency points, and prints the spread
across repeat step runs so the uncertainty of the method can be told apart from the error
of the model.

## Running it

Needs Python with numpy and matplotlib.

```
python -m venv .venv
.venv/bin/pip install numpy matplotlib
```

Check it still works before trusting it on real data. The self-test generates responses
from plants with known k and b and asserts that the analysis recovers them, in both the
underdamped and overdamped cases, along with the frequency readings and two alternative CSV
header spellings.

```
.venv/bin/python analyse_lab1.py --selftest --outdir /tmp/check
```

Then, with the logs in one folder:

```
.venv/bin/python analyse_lab1.py step.csv step2.csv step3.csv --sine-dir .
```

Sine logs are named for the frequency they were driven at, with `p` standing in for the
decimal point, so `sine_1p6.csv` is 1.6 rad/s. The frequency is read out of the file name.
Figures come from the first step log; the rest give the repeatability spread.

The run writes `results.json` and four figures into `report/figures/`: the applied input,
the annotated step response, the Bode plot with the identified model overlaid, and the
simulated response against the measured one.

## Report

`report/` holds a LaTeX skeleton that consumes those figures. Values that only exist once
the bench data is captured are wrapped in `\gap{}` and print in red, so anything still red
is unfinished. Build with `pdflatex main.tex` twice, from inside `report/`.

## Bench notes

Two things are easy to get wrong on the day and expensive to discover afterwards. The
simulator appends every run to one fixed CSV file name, so every run has to be copied,
renamed and cleared immediately, or runs stack up in one file that still looks valid. And
the frequency plan has to reach far enough above the corner for the magnitude to reach its
asymptote, otherwise the roll-off measures near -60 dB per decade rather than -40 and reads
as though a second order model had failed.
