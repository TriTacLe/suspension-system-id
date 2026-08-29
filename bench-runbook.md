---
type: note
status: active
project: uct
course: EEE3094S
tags: [uct, control-systems, lab]
---

# Lab 1 lab day: everything to do, in order

Follow this top to bottom at the machine. No electronics, no wiring, no oscilloscope. Lab 1
is a Windows program on a Control Lab PC. You click buttons, watch a graph, and copy CSV
files onto a USB stick. Everything after that is already written and waiting in this folder.

Due 4 September. One PDF to Gradescope through the LTI link in the Lab 1 folder on Amathuba.

## What the experiment actually is

The program simulates the suspension of a jet: a spring and a shock absorber carrying a
mass. Two numbers describe it, the spring coefficient `k` and the damping coefficient `b`,
and they are generated from your student number, so nobody else has your plant and you
cannot look the answer up.

You cannot measure `k` and `b` directly. You can only push the suspension and watch how it
moves. So you push it in two ways:

- Push once and hold, a step. How far it settles tells you the spring. How much it bounces
  past that point, and how fast the bouncing dies away, tells you the damping.
- Push it up and down at a steady rate, a sine. Do that at several rates and you find the
  rate where the suspension responds most, which is a second, independent way to get the
  same numbers.

That is the whole lab. Ten logs, one hour.

## Before you go

- Book your slot under the Groups tab in the Lab 1 folder on Amathuba if you have not.
- Bring a USB stick. The CSV files only exist on that machine.
- Bring something to write on. Four numbers get written down while the graph is still on
  screen, and you cannot recover them later.
- No preparation is possible on your laptop. The plant does not exist until you sign on.

## Stage 0: sign on

1. Open the suspension simulator.
2. Top box, **Select your student number**, click the dropdown, pick **LXXTRI004**.
3. Wait until the panel reads **Ready**. Your plant is now loaded. If you skip this you are
   testing somebody else's system and the marks go with it.
4. Click **Laboratory Testing** in the **Test Functions** box.

Do not click **Start Step Test**. That button drops the aircraft and steps its height, not
the suspension input. It looks like the right button and it is not.

A **Laboratory Tests** panel opens with an **Offset Input** box, a **Sine Wave Geni** box
holding Rads/sec and Amplitude, a green **Apply Input** and a red **Stop lab. Testing**.

## Stage 1: the step test

Do this three times. The whole run is about five minutes.

1. Click **Show Graph** so you can see the response, and **Reset Graph** so it starts clean.
2. Type `1` into **Offset Input**.
3. Click **Start Logging**. Logging first, always. A log that starts after the step has
   already happened is useless, because the analysis needs the flat part before the step to
   know where the suspension was resting.
4. Click **Apply Input**.
5. Watch the graph. The displacement rises, probably overshoots, and settles.
6. Wait until the trace has been visibly flat for a good while. Longer than feels
   necessary. The value it settles at is the single most important number in the lab, since
   the spring coefficient comes straight from it. A log cut while the curve is still
   drifting gives a wrong spring coefficient and everything downstream inherits it.
7. Click **Stop Logging**.
8. Copy `C:\F34BSuspensionTestData\F34BSuspensionTestData.csv` to your stick and rename it
   `step.csv`. Copy it now. The program overwrites that same file on the next run.
9. Repeat twice more, saving `step2.csv` and `step3.csv`. Sampling is coarse, and the
   repeats tell you how much the reading moves between runs, which is worth a sentence in
   the report.

### Write these down before you touch anything else

While the graph is still on screen:

| What | How to read it | Your value |
| --- | --- | --- |
| Resting level | Where the trace sat before the step | |
| Final level | Where the trace settles after the step | |
| Highest point | The top of the first hump, if there is one | |
| Hump spacing | Time from the first hump top to the second hump top, in seconds | |

If there is no hump at all, if the trace just rises and flattens, write "no overshoot" and
carry on. That is a valid result, it means the suspension is heavily damped, and the
analysis script handles it by a different route.

### Work out one number at the bench

Call the hump spacing T seconds. Then

```
W = 6.28 / T
```

W is roughly the frequency the suspension naturally wants to oscillate at, in radians per
second, and it is what tells you which sine frequencies are worth testing. For example a
hump spacing of 3.2 s gives W = 6.28 / 3.2 = 2.0.

If there was no overshoot, use W = 1 and expect to add frequencies as you go.

W is the damped natural frequency, not the undamped one. They differ by the factor
sqrt(1 - zeta squared), which is a few percent for a lightly damped suspension. That is
close enough to choose test frequencies with, which is all W is for here. Do not write W
into the report as the natural frequency: the script computes both properly from the log.

## Stage 2: the sine tests

Eight runs, roughly forty five minutes. Test at these frequencies, worked out from your W:

| Point | Frequency to enter | Example, W = 2.0 |
| --- | --- | --- |
| 1 | 0.2 x W | 0.4 |
| 2 | 0.5 x W | 1.0 |
| 3 | 0.8 x W | 1.6 |
| 4 | 1.0 x W | 2.0 |
| 5 | 1.25 x W | 2.5 |
| 6 | 2 x W | 4.0 |
| 7 | 4 x W | 8.0 |
| 8 | 6 x W | 12.0 |

Two points well below W and three well above are what prove the suspension behaves like a
second order system rather than a simpler one. The three around W are where the interesting
behaviour lives.

Point 8 is the one that earns the roll-off claim. Two poles and no zeros give -40 dB per
decade, but only well above W. Measured between 2W and 4W the slope comes out nearer -60,
because just past resonance the curve falls much steeper than its asymptote. Without a
point out at 6W the report would show about -60 and read as though the second order model
had failed, when it has not. If you are short of time, point 8 matters more than point 5.

Logging is fixed at 100 ms, so the highest frequency worth entering is about 15 rad/s. Above
that there are too few samples per cycle to fit anything. If 6 x W lands above 15, enter 15
and note in the report that the sampler set the ceiling.

For each frequency:

1. Click **Reset Graph**.
2. Type the frequency into **Rads/sec** in the **Sine Wave Geni** box.
3. Type `1` into **Amplitude**.
4. Click **Start Logging**.
5. Click **Apply Input**.
6. Let it run until you have counted at least ten full up and down cycles. At low
   frequencies this takes a while. Do not cut it short: the first few cycles are settling
   and get discarded by the analysis, so a short log leaves almost nothing usable.
7. Click **Stop Logging**.
8. Copy the CSV to your stick, named after the frequency with `p` in place of the decimal
   point: 1.6 rad/s becomes `sine_1p6.csv`, 0.4 becomes `sine_0p4.csv`. The analysis script
   reads the frequency straight out of the file name, so this naming matters.

### Watch one thing while they run

Note roughly how big the output swing is at each frequency. It should grow as you approach
W and shrink after. If it looks biggest somewhere between two of your frequencies, run one
extra test in that gap. That peak is the strongest evidence of the damping you will get all
day, and an extra log costs three minutes.

Also note whether the output looks shifted in time against the input. It should lag more
and more as frequency rises. At W it should be about a quarter cycle behind.

## Before you leave

- Eleven CSV files on the stick: `step.csv`, `step2.csv`, `step3.csv`, and eight `sine_*.csv`.
- Open one of them in Excel on the lab machine and check it has three columns, a time, an
  input and an output displacement, with numbers under them. An empty or single column file
  means logging was not running and the test has to be redone.
- Write down the exact column headings as spelled in the file. The script matches them
  loosely, so spacing and brackets do not matter, but if it still cannot find a column it
  prints the headings it did find and they go straight into `_column()`.
- The four handwritten numbers from Stage 1.
- Click **Stop lab. Testing**.

## Things that go wrong

**Every CSV looks the same.** The program writes to one fixed file name and overwrites it.
Copy and rename after every single run.

**The file is empty or has only headers.** Logging was not started before the input was
applied. Redo that run.

**The response never settles inside the log.** Log for longer. If it still will not settle,
say so in the report and use the last part of the record, but say what you did.

**No overshoot at all.** Fine. Skip the hump spacing, note "no overshoot", use W = 1 and
spread the sine frequencies wider, say 0.1 to 10 rad/s.

**Wrong student number.** Sign off, sign on again as LXXTRI004, redo everything. Data from
another plant is worth nothing.

## After the lab

From `assignments/lab-1/`, with the CSVs in one folder:

```
../../.venv/bin/python analyse_lab1.py step.csv step2.csv step3.csv --sine-dir .
```

That is the whole analysis. It reads every step log you pass and every `sine_*.csv` in the
folder, prints the gain, damping ratio, natural frequency, spring and damping coefficients,
poles and transfer function, writes `results.json`, and puts four finished figures into
`report/figures/`. The figures come from the first step log; the others are there to give
the spread across repeat runs, which the report quotes as the repeatability of the method.

It also reads the four Stage 1 answers a second time straight off the sine data: the gain
from the low frequency end, the damping ratio from the height of the resonant peak, and the
natural frequency from where the phase crosses -90 degrees. Comparing those against the step
answers is what the brief is built around, so both sets go in the report.

It needs numpy and matplotlib. The virtualenv at `../../.venv` in the course folder has
them already. Anywhere else: `python -m venv .venv && .venv/bin/pip install numpy matplotlib`.

Check it still passes its own tests before trusting it on your data:

```
../../.venv/bin/python analyse_lab1.py --selftest --outdir /tmp/lab1check
```

Then:

1. Open `report/`. Every value the lab produces is marked in red. Fill each one from what
   the script printed. Anything still red is unfinished.
2. Several red gaps are sentences, not numbers, asking what your own data showed. Those are
   where the marks are. The brief says outright that you are marked on approach and on how
   clearly you explain it, not on how accurate the model is.
3. Build with `pdflatex main.tex` twice from inside `report/`. The second pass fills in the
   cross-references.
4. Check the PDF has no red left in it.
5. Submit `main.pdf` to Gradescope through the LTI link in the Lab 1 folder on Amathuba.

## One thing worth doing for marks

Before you run the analysis, write down what you expect: is the system bouncy or sluggish,
roughly what spring and damping coefficients do you expect from the shape you saw, where do
you expect the sine response to peak. Then compare against what the script gives. A
prediction that turns out wrong, explained honestly, reads far better than numbers with no
reasoning attached, and reasoning is exactly where marks are usually lost in this course.
