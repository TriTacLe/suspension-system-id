# Lab 1: System Identification

What the lab asks for, condensed from the brief. The brief and the simulator manual are the
authority; this file exists so the repo is self-explanatory without them.

## The task

Identify the transfer function of a simulated spring, mass and damper suspension using two
experiments. The plant parameters are randomised per student and cannot be looked up, so
every number has to come off a measurement.

The system is

```
m x'' = -b x' - k x + F
```

with `x` the displacement of the mass and `F` an applied force proportional to the input
voltage. Everything is scaled to the mass, so `m = 1` and initial conditions are zero.

The randomised plant can land anywhere on the damping range. If it comes out overdamped the
response will look first order, and proposing a reduced order model is the right answer
rather than a failure. The brief says outright that it may not be possible to place both
poles or recover both coefficients from tests this rough, and that the marks are for the
approach, not the accuracy of the final model.

## Stage 1, step test

Apply a step and record the response over a sensible interval.

The simulator's `Start Step Test` steps the height of the plane, not the suspension, which
is not the response wanted here. Use `Laboratory Testing` and apply a unit step with the
offset input.

Four questions to answer off the step response, and the answers set up Stage 2:

1. What is the gain of the system?
2. Are the poles real or complex?
3. Where is the dominant pole or complex pole pair?
4. Given 3, what corner frequency would you expect?

## Stage 2, frequency response test

Drive the system with single sinusoids across a range chosen from the Stage 1 answers, and
read the magnitude and phase change between input and output at each one. The simulator does
not sweep, so the Bode plot is assembled point by point.

Use the frequency response to confirm all four Stage 1 answers. Depending on the damping it
may also expose a non-dominant pole that the step response could not show, which extends the
answer to question 3.

At least five frequencies, more of them clustered near the estimated corner. Plot the points
as they come in so the gaps are visible while there is still time to fill them.

## Report

Written to the provided template, which sets the minimum content rather than a form to fill
in. The stated marking emphasis is on communication: the marker sees the report, not the
work behind it.

What the brief asks for specifically:

- Calculations typed with an equation editor, one step per line, each step introduced in
  words, with the formulas and assumptions named.
- Every symbol defined in the text, including the obvious ones.
- No unexplained numbers. Where a value came from, usually a plot, has to be stated.
- Figures exported properly rather than screenshotted, axes labelled and legible, captions
  that describe the feature that matters.
- Every figure referenced in the text. An unreferenced figure should be deleted instead.
- Anything important in a plot marked in three places: annotated on the figure, named in the
  caption, and used as evidence in the text.

## Source files

`Lab 1.pdf` (brief and report template) and `Lab 1 User manual.pdf` (simulator) are on
Amathuba under the course's lab material. They are not committed here.
