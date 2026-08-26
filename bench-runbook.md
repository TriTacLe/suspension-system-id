---
type: note
status: active
project: uct
course: EEE3094S
tags: [uct, control-systems, lab]
---

# Lab 1 bench runbook

What to do at the Control Lab machine so the session is capture only. Everything after the
capture already exists: `analyse_lab1.py` does the identification and draws the figures, and
`report/` is the write-up with the measured numbers marked in red.

Due 4 September, one PDF to Gradescope through the LTI link in the Lab 1 folder on Amathuba.

## Before you sit down

The plant is keyed to your student number and only exists on the lab machines, so nothing
here can be done in advance. Bring a USB stick for the CSV files.

## Capture, in order

1. Sign on with LXXTRI004 so your own transfer function loads. Wait for **Ready**.
2. Click **Laboratory Testing**. Do not use **Start Step Test**, it steps the height of the
   aircraft instead of the suspension input.
3. Step test. Enter a step in **Offset Input**, press **Apply Input**, then **Start Logging**
   before the response begins and **Stop Logging** well after it has gone flat. A record that
   stops while the curve is still moving gives a wrong final value, and the final value is
   the one number everything else scales off.
4. Copy `C:\F34BSuspensionTestData\F34BSuspensionTestData.csv` to the stick as `step.csv`.
   The simulator overwrites the same file on every run, so copy it before the next test.
5. Repeat the step at least twice. 100 ms sampling can miss the true peak, and two records
   let you say how much the overshoot reading moves.
6. Sine tests. Run seven frequencies, from the natural frequency you can already estimate at
   the bench: peaks in the step response are one damped period apart, so
   `wn ~ 2*pi/T`. Test at 0.2, 0.5, 0.8, 1.0, 1.25, 2 and 4 times that. Enter the frequency
   in rad/s and an amplitude in **Sine Wave Geni**, **Apply Input**, log for at least ten
   full cycles, then copy the CSV as `sine_<frequency>.csv`.
7. Watch the magnitude as you go. If the output amplitude peaks between two of your points,
   add a frequency there. That peak is the sharpest evidence of the damping ratio you will
   get.

## After the session

Run the analysis from `assignments/lab-1/`:

```
../../.venv/bin/python analyse_lab1.py step.csv \
    --sine 0.4=sine_0p4.csv --sine 1.0=sine_1p0.csv ...
```

It needs numpy and matplotlib. The virtualenv at `../../.venv` in the course folder already
has them; anywhere else, `python -m venv .venv && .venv/bin/pip install numpy matplotlib`.

It prints the gain, damping ratio, natural frequency, k, b and the poles, writes
`results.json` next to the figures, and puts the four annotated figures straight into
`report/figures/`. Check `--selftest` still passes if you change the script; it identifies a
plant with known coefficients and fails loudly if the maths breaks.

Then fill in the report. Every value the bench produces is marked `\gap{...}` and prints in
red, so anything still red is unfinished. Build with `pdflatex main.tex` twice, the second
pass resolves the cross-references.
