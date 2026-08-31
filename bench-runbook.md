---
type: note
status: active
project: uct
course: EEE3094S
tags: [uct, control-systems, lab]
---

# Lab 1: everything to do, in order

Follow this top to bottom. Part 1 is the lab day. Part 2 is the report. No electronics, no
wiring, no oscilloscope. Lab 1 is a Windows program on a Control Lab PC. You click buttons,
watch a graph, and copy CSV files onto a USB stick.

Due 4 September. One PDF to Gradescope through the LTI link in the Lab 1 folder on Amathuba.

## What the experiment actually is

The program simulates the suspension of a jet: a spring and a shock absorber carrying a
mass. Two numbers describe it, the spring coefficient `k` and the damping coefficient `b`,
and they are generated from your student number, so nobody else has your plant and you
cannot look the answer up.

The brief gives you the equation of the system:

```
m x'' = -b x' - k x + F,   with m = 1,   so   G(s) = 1 / (s^2 + b s + k)
```

You cannot measure `k` and `b` directly. You can only push the suspension and watch how it
moves. So you push it in two ways:

- Push once and hold, a step. How far it settles tells you the spring. How much it bounces
  past that point, and how fast the bouncing dies away, tells you the damping.
- Push it up and down at a steady rate, a sine. Do that at several rates and you find the
  rate where the suspension responds most, which is a second, independent way to get the
  same numbers.

Ten to eleven logs, about an hour.

## How this is marked, read this before anything else

The brief says it twice, in bold: **you are marked primarily on your approach, not the
accuracy and detail of your eventual model**, and **your markers see your report, not your
actual work**. A rough model explained clearly beats a precise model dumped on the page.

That has one practical consequence for the lab day: write down what you see and what you
think it means while you are at the bench. You cannot reconstruct your reasoning three days
later, and the reasoning is what is being marked.

# Part 1: the lab day

## Before you go

- Book your slot under the Groups tab in the Lab 1 folder on Amathuba if you have not.
- Bring a USB stick. The CSV files only exist on that machine.
- Bring something to write on. Several numbers get written down while the graph is still on
  screen, and you cannot recover them later.
- No preparation is possible on your laptop. The plant does not exist until you sign on.

## Stage 0: sign on

1. Open the suspension simulator.
2. Top box, **Select your student number**, click the dropdown, pick **LXXTRI004**.
3. Wait until the panel reads **Ready**. Your plant is now loaded. If you skip this you are
   testing somebody else's system and the marks go with it.
4. Click **Laboratory Testing** in the **Test Functions** box.

Do not click **Start Step Test**. The brief warns about this directly: that button steps the
height of the plane, not the suspension. It looks like the right button and it is not.

A **Laboratory Tests** panel opens with an **Offset Input** box, a **Sine Wave Geni** box
holding Rads/sec and Amplitude, a green **Apply Input** and a red **Stop lab. Testing**.

### If LXXTRI004 is not in the dropdown

It was missing as of 29 August, the list running from LTSNKU002 straight to LZZALB001.
Nothing in this runbook works without it, because no selection means no plant. Do not sit
through the slot hoping it appears.

Check the list first, before anything else. If the number is not there:

1. Find the lab staff member on duty and ask them to add LXXTRI004 to the simulator. It is
   a list on that machine, so there is a chance it can be done while you wait.
2. If they cannot, photograph the dropdown showing the gap between LTSNKU002 and
   LZZALB001, with the date visible if you can manage it. That photo is the evidence for
   the extension request.
3. Mail Dr Shield from the lab, not later. Say the slot has been used and the number is
   still missing, attach the photo, and ask her to confirm the extension.
4. Rebook the next available slot before you leave the building.

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
| Shape in words | Bouncy, or slow rise with no bounce | |

If there is no hump at all, if the trace just rises and flattens, write "no overshoot" and
carry on. That is a valid result, it means the suspension is heavily damped, and the
analysis script handles it by a different route.

### The four questions the brief asks at the end of Stage 1

These are graded deliverables, not warm-up. Answer them roughly at the bench, in pen, then
properly at home with the script. Your bench answers decide which frequencies you test, so
getting them roughly right matters more than getting them exactly right.

**1. What is the gain of the system?**

Gain is how far the output moves for a one unit input, once everything has settled.

```
gain A = (final level - resting level) / size of the step
```

You typed `1` into Offset Input, so the step size is 1 and the gain is just the change in
level. And since the plant is `1/(s^2 + bs + k)`, the gain at zero frequency is `1/k`, so
`k = 1 / A`. That is your spring coefficient already.

**2. Does the system have real or complex poles?**

Look at the shape. Nothing else.

- It bounced past the final value at least once, so complex poles. Underdamped.
- It rose and flattened with no bounce at all, so real poles. Overdamped or critically
  damped.

Write down which one you saw and the evidence for it, because that sentence is the start of
your report.

**3. Where do you estimate the dominant pole or complex pole pair to be?**

If it bounced (complex poles):

```
T  = hump spacing in seconds
wd = 6.28 / T                       damped natural frequency
OS = (highest point - final level) / (final level - resting level)
```

`OS` is the overshoot as a fraction. Roughly, a 10 percent overshoot means a damping ratio
near 0.6, 25 percent near 0.4, 50 percent near 0.2. The exact relation is
`zeta = -ln(OS) / sqrt(pi^2 + ln(OS)^2)`, and the script does it properly. Then

```
wn    = wd / sqrt(1 - zeta^2)
poles = -zeta*wn  +/-  j*wd
```

The real part is how fast the bouncing dies away, the imaginary part is how fast it bounces.

If it did not bounce (real poles): you probably cannot see both poles from a step. Say so.
Take the time to reach 63 percent of the final value, call it `tau`, and report a reduced
first order model `A / (tau s + 1)` with the dominant pole at `-1/tau`. The brief explicitly
invites this: "you may wish to propose a reduced-order model". Proposing it and explaining
why is worth marks. Pretending you found two poles you could not see loses them.

**4. From 3, what corner frequency do you expect?**

The corner frequency is where the Bode magnitude stops being flat and starts falling.

- Complex poles: it is `wn`, and there will be a peak near it.
- Real dominant pole at `-1/tau`: it is `1/tau`.

Call this number W. It is the centre of your sine testing.

If there was no overshoot and you have no idea, use W = 1 and expect to add frequencies as
you go.

Note W is not exactly the natural frequency when you get it from the hump spacing, since
`wd` and `wn` differ by `sqrt(1 - zeta^2)`, a few percent for a lightly damped suspension.
That is close enough to choose test frequencies with, which is all W is for here. Do not
write the bench W into the report as the natural frequency: the script computes both
properly from the log.

## Stage 2: the sine tests

Eight runs, roughly forty five minutes. The brief says you need at least five frequencies
and should put extra ones near the corner frequency. Eight is what makes the report solid.

Work them out from your W:

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

### Plot as you go, the brief asks for this

Keep a rough table on paper and sketch it as you work:

| Frequency | Output swing size | Behind the input by |
| --- | --- | --- |
| | | |

The swing should grow as you approach W and shrink after. If it looks biggest somewhere
between two of your frequencies, run one extra test in that gap. That peak is the strongest
evidence of the damping you will get all day, and an extra log costs three minutes.

The lag should grow steadily as frequency rises. At W it should be about a quarter cycle
behind. Far above W it should approach half a cycle.

### Confirming the four answers, which is what Stage 2 is for

The brief asks you to use the frequency response to confirm your Stage 1 answers. The script
does all of this from your logs, but know what it is doing:

- Gain, from the flat low frequency end of the magnitude plot. Compare against the step gain.
- Complex or real poles, from whether there is a resonant peak. A peak means complex poles.
  No peak, just a steady bend, means real poles.
- Pole positions, from the height of the peak, which gives `zeta`, and the frequency where
  the phase passes -90 degrees, which gives `wn` exactly.
- Corner frequency, which is where the magnitude leaves the flat region, compared against
  the W you predicted.

One extra thing the brief mentions: if your system is overdamped, the frequency response may
show a second bend that the step response could not reveal, which gives you the non-dominant
pole. Look for two separate corners in the magnitude plot. If you see them, say so, because
that is explicitly called out as going beyond the minimum.

## Before you leave

- Eleven CSV files on the stick: `step.csv`, `step2.csv`, `step3.csv`, and eight `sine_*.csv`.
- Open one of them in Excel on the lab machine and check it has three columns, a time, an
  input and an output displacement, with numbers under them. An empty or single column file
  means logging was not running and the test has to be redone.
- Write down the exact column headings as spelled in the file. The script matches them
  loosely, so spacing and brackets do not matter, but if it still cannot find a column it
  prints the headings it did find and they go straight into `_column()`.
- The handwritten numbers from Stage 1 and the frequency table from Stage 2.
- Photograph the graph on screen for at least one step run and one sine run. Not for the
  report, the script draws proper figures, but as a check that your CSV matches what you
  actually saw.
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

# Part 2: after the lab

## Run the analysis

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

The self-test builds a system with a known `k` and `b`, runs the whole analysis on it, and
checks the answers come back. If it passes, anything odd in your results is a property of
your suspension, not a bug in the code.

## Write the report

Open `report/`. Every value the lab produces is marked in red. Fill each one from what the
script printed. Anything still red is unfinished. Build with `pdflatex main.tex` twice from
inside `report/`, the second pass fills in the cross-references, then check the PDF has no
red left in it.

Several red gaps are sentences, not numbers, asking what your own data showed. Those are
where the marks are.

### The rules the brief actually grades on

The brief spends three pages on this. Short version:

**Every figure earns its place.** If you do not reference a figure in the text, delete it.

**Every important feature appears in three places.** The brief spells this out: on the plot
itself as an annotation, in the caption, and in the text where you use it as evidence. So
the overshoot gets marked on the step plot, mentioned in the caption, and used in the
paragraph where you calculate the damping ratio. This feels repetitive. Do it anyway, it is
directly worth marks.

**Captions describe, they do not label.** Not "Figure 2: Step Response". The brief's own good
example is "Figure 2: Step Response. The response is oscillatory, with an overshoot (OS) of
0.4 m. The final value (A) is 1 m and the period (T) is 1 s."

**No magic numbers.** Every value in a calculation says where it came from. Not "the gain is
0.26" but "the final value in Figure 2 is 0.26 m for a 1 V step, so the gain is 0.26 m/V".

**State the formula before you use it, and define every symbol**, even the obvious ones. Say
which assumptions you are making and why: zero initial conditions, mass normalised to 1,
second order model chosen because the response oscillates.

**Words and mathematics flow together.** Describe the step in a sentence, then show the
mathematics, then say what the result means. Do not stack equations with no text between.

### What separates a good mark from a full one

- Say what you expected before each test, then what you got. A prediction that turned out
  wrong, explained honestly, reads far better than numbers with no reasoning attached, and
  reasoning is exactly where marks are usually lost in this course.
- Compare the step answers against the frequency answers explicitly, number by number, and
  explain any disagreement. That comparison is the point of the whole lab.
- Quote the spread across your three step runs as the repeatability of the method, so the
  reader can tell measurement noise apart from model error.
- Say what the model cannot do. If you could not locate a second pole, say so and say why.
  The brief tells you outright this may happen and expects you to handle it honestly.
- Sketch what you think the physical system is like in plain words. A bouncy suspension with
  light damping, or a stiff sluggish one. It shows you understand what the numbers mean.

## Submit

`main.pdf` to Gradescope through the LTI link in the Lab 1 folder on Amathuba. Check the PDF
one last time for red text before you upload.
