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
mean spacing of every peak in the record. Only peaks clearing the final value by several
times the logger resolution count, since a sample one step above the settled value is
quantisation noise and counting it wrecks the averaged period.

It also fits the two real pole form, a constant plus two decaying exponentials, to the same
record. That fit cannot produce an overshoot at all, which is the numerical version of
reading the pole type off the shape of the curve rather than judging it by eye.

From the sine logs it reads the same quantities a second time and independently: the gain
from the low frequency end, the damping ratio from the height of the resonant peak, and the
natural frequency from where the phase crosses -90 degrees. Amplitude and phase at each
drive frequency come from a least squares fit against a sinusoid at the known frequency,
which works on records far too coarsely sampled for peak picking.

A second order model cannot lag past -180 degrees. Where the measured phase runs past it,
the excess is fitted as a transport delay, which grows without bound and leaves the
magnitude alone, rather than as a third pole, which would bend both.

It then validates the model in both domains, as RMS error against the step record and as
RMS magnitude and phase error against the measured frequency points, and prints the spread
across repeat step runs so the uncertainty of the method can be told apart from the error
of the model. Runs at a second input amplitude check that the plant is linear at all.

## Running it

`simulate_lab1.m` does the whole identification and draws the figures. It needs the Control
System Toolbox, for `tf` and `step`.

```
matlab -batch simulate_lab1
```

No arguments. The log folder is set by the `datadir` line at the top of the script, which is
resolved relative to the file, so it runs from anywhere. Point it at your own session
folder.

The arithmetic was checked once against a second implementation written from the same
equations without sharing any code, which agreed to four decimal places on every parameter.
That one was in Python and has been removed. Two things it established are worth recording,
since the code that showed them is gone: the two implementations could disagree about a bug
but not about the plant, and fitting the delay to synthetic data from a plant with no delay
returns 0 ms, so a delay fitted to real logs is in the measurements rather than invented by
the fit.

## Log naming

Sine logs are named for the frequency they were driven at, with `p` standing in for the
decimal point, so `sine_1p6.csv` is 1.6 rad/s. The frequency is read out of the file name.
A name carrying a tag after the frequency, like `sine_1p0_a3.csv`, is that frequency driven
at another amplitude: it is held out of the sweep and used for the linearity check instead.
Figures come from the first step log; the rest give the repeatability spread.

## Report

`report/` holds a LaTeX skeleton that consumes the figures. Values that only exist once the
bench data is captured are wrapped in `\gap{}` and print in red, so anything still red is
unfinished. Build with `pdflatex main.tex` twice, from inside `report/`.

## Bench notes

Two things are easy to get wrong on the day and expensive to discover afterwards. The
simulator appends every run to one fixed CSV file name, so each run has to be copied,
renamed and cleared before the next is started, or several runs stack up in one file that
still looks valid. And the sine frequencies have to reach far enough above the corner for
the magnitude to settle onto its asymptote, otherwise the measured roll-off reads far too
shallow and looks as though a second order model had failed.
