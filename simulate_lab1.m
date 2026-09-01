function simulate_lab1
% Identify the suspension from the logs, simulate it, and draw the figures.
% Every number printed comes off a logged measurement. Paths are relative to
% this file, so it runs from anywhere.

here = fileparts(mfilename('fullpath'));
datadir = fullfile(here, 'data', 'lab-session-2026-08-31');
outdir = fullfile(here, 'report', 'figures');

steplog = fullfile(datadir, 'step1.csv');
[A, zeta, wn, Tp, overshoot, period, zeta2, dev] = identify_step(steplog);

k = wn^2;
b = 2*zeta*wn;
alpha = A*k;      % the input is scaled, so alpha is not 1

fprintf('Identified from the step test\n');
fprintf('  A = %.4f m/V   zeta = %.4f   wn = %.4f rad/s\n', A, zeta, wn);
fprintf('  overshoot %.4f at Tp = %.2f s, damped period %.2f s\n', ...
        overshoot, Tp, period);
if ~isnan(zeta2)
    fprintf('  zeta again %.4f, from the decay ratio of the %.0f and %.0f mm extrema\n', ...
            zeta2, 1000*dev(1), 1000*dev(2));
end
fprintf('  k = %.4f   b = %.4f   alpha = %.4f\n', k, b, alpha);
fprintf('  poles %.4f +/- %.4fj\n\n', -zeta*wn, wn*sqrt(1-zeta^2));

two_exponential_check(steplog, zeta, wn);
repeats_check(datadir);

G = tf(alpha, [1 b k]);

plot_input(steplog, outdir);
plot_step(steplog, outdir);
frequency_check(G, datadir, outdir);
step_check(G, steplog, outdir);
end


function [A, zeta, wn, Tp, overshoot, period, zeta2, dev] = identify_step(path)
% Gain, damping ratio and natural frequency from one step log.

[t, u, y] = read_log(path);
i0 = step_index(u);
tail = max(1, floor(0.1*numel(y)));

y0 = mean(y(1:i0-1));
yinf = mean(y(end-tail+1:end));
du = mean(u(end-tail+1:end)) - u(1);
dy = yinf - y0;

resp = y(i0:end) - y0;
tr = t(i0:end) - t(i0);

% One sample above the final value is quantisation noise, not a peak, and
% counting it wrecks the averaged period. Require several resolution steps.
levels = unique(y);
resolution = min(diff(levels));
[peaks, troughs] = turning_points(resp);
peaks = peaks(resp(peaks) > dy + 3*resolution);
if isempty(peaks)
    error('simulate_lab1:noOvershoot', ...
          'no peak clears the resolution, so the response is not underdamped');
end

% %OS = 100 exp(-pi zeta / sqrt(1 - zeta^2)), inverted for zeta
overshoot = (resp(peaks(1)) - dy) / dy;
ln = log(overshoot);
zeta = -ln / sqrt(pi^2 + ln^2);

% Average over the peaks when there are several. With only one, the time to
% it is half a damped period and that is all the timing there is.
if numel(peaks) >= 2
    period = mean(diff(tr(peaks)));
else
    period = 2 * tr(peaks(1));
end
wd = 2*pi / period;
wn = wd / sqrt(1 - zeta^2);

A = dy/du;
Tp = tr(peaks(1));

% Log decrement, a second reading of zeta. Peak to trough is one decay
% factor, peak to peak is two, and the second peak is under the resolution.
zeta2 = NaN;
dev = [NaN NaN];
after = troughs(troughs > peaks(1));
after = after(abs(resp(after) - dy) > 3*resolution);
if ~isempty(after)
    dev = abs(resp([peaks(1) after(1)]) - dy);
    lnr = log(dev(2) / dev(1));
    zeta2 = -lnr / sqrt(pi^2 + lnr^2);
end
end


function two_exponential_check(path, zeta, wn)
% Can two real poles fit this record? Fit c + P exp(-at) + Q exp(-bt) and see.
% P, Q and c are linear once the rates are fixed, so only the rates are searched.

[t, u, y] = read_log(path);
i0 = step_index(u);
tail = max(1, floor(0.1*numel(y)));
y0 = mean(y(1:i0-1));
dy = mean(y(end-tail+1:end)) - y0;
resp = y(i0:end) - y0;
tr = t(i0:end) - t(i0);

rates = logspace(log10(0.02), log10(20), 240);
best = Inf;
for i = 1:numel(rates)
    for j = i+1:numel(rates)
        basis = [ones(size(tr)), exp(-rates(i)*tr), exp(-rates(j)*tr)];
        c = basis \ resp;
        fit = basis * c;
        err = sqrt(mean((fit - resp).^2));
        if err < best
            best = err;
            bestfit = fit;
        end
    end
end

model = dy * step_response(tr, zeta, wn);
complexerr = sqrt(mean((model - resp).^2));

fprintf('Real or complex poles\n');
fprintf('  two real poles fit to %.1f mm RMS, complex pair to %.1f mm\n', ...
        1000*best, 1000*complexerr);
fprintf('  measured peak stands %.1f mm above the final value, the real fit reaches %.1f mm\n\n', ...
        1000*(max(resp) - dy), 1000*(max(bestfit) - bestfit(end)));
end


function repeats_check(datadir)
% The repeats say how much of the model error is measurement noise. The runs
% driven harder say whether the plant is linear, which everything here assumes.

runs = {'step1.csv', 'step2.csv', 'step3.csv', 'step_amp2.csv'};
fprintf('Repeats and linearity\n');
fprintf('  %14s %8s %8s %8s %8s\n', 'run', 'A', 'zeta', 'wn', 'Tp');
for n = 1:numel(runs)
    [A, zeta, wn, Tp] = identify_step(fullfile(datadir, runs{n}));
    fprintf('  %14s %8.4f %8.4f %8.4f %8.2f\n', runs{n}, A, zeta, wn, Tp);
end

% A linear plant scales its output with the input and leaves the ratio and
% the phase alone.
[m1, p1] = analyse_sine(fullfile(datadir, 'sine_1p0.csv'), 1.0);
[m3, p3] = analyse_sine(fullfile(datadir, 'sine_1p0_a3.csv'), 1.0);
fprintf('  1 rad/s at amplitude 1: ratio %.4f, phase %.2f deg\n', 10^(m1/20), p1);
fprintf('  1 rad/s at amplitude 3: ratio %.4f, phase %.2f deg\n\n', 10^(m3/20), p3);
end


function y = step_response(t, zeta, wn)
% Unit step response of wn^2 / (s^2 + 2 zeta wn s + wn^2), underdamped.

wd = wn * sqrt(1 - zeta^2);
phi = atan2(sqrt(1 - zeta^2), zeta);
env = exp(-zeta*wn*t) / sqrt(1 - zeta^2);
y = 1 - env .* sin(wd*t + phi);
end


function [peaks, troughs] = turning_points(y)
% Indices of local maxima and minima, ignoring flat runs.

peaks = []; troughs = [];
s = sign(diff(y(:)));
idx = find(s ~= 0);
if numel(idx) < 2
    return
end
ss = s(idx);
flips = find(ss(1:end-1) ~= ss(2:end));
peaks = idx(flips(ss(flips) > 0)) + 1;
troughs = idx(flips(ss(flips) < 0)) + 1;
end


function plot_input(path, outdir)
% The input as applied, so the step the analysis reads is on the record.

[t, u, ~] = read_log(path);

fig = figure('Visible', 'off', 'Position', [0 0 900 450]);
plot(t, u, 'LineWidth', 2, 'Color', [0.17 0.24 0.44]);
xlabel('Time [s]'); ylabel('Input [V]');
title('Applied input');
ylim([-0.15 1.35]); grid on;
save_figure(fig, outdir, 'input_step.png');
end


function plot_step(path, outdir)
% The measured step, with the baseline, the final value and the step marked.

[t, u, y] = read_log(path);
i0 = step_index(u);
tail = max(1, floor(0.1*numel(y)));
y0 = mean(y(1:i0-1));
yinf = mean(y(end-tail+1:end));
tstep = t(i0);

fig = figure('Visible', 'off', 'Position', [0 0 950 600]);
plot(t, y, 'LineWidth', 2, 'Color', [0.17 0.24 0.44]); hold on;
yline(y0, ':', 'Color', [0.7 0.1 0.1]);
yline(yinf, '--', 'Color', [0.7 0.1 0.1]);
xline(tstep, ':', 'Color', [0.35 0.35 0.35]);

xlabel('Time [s]'); ylabel('Displacement [m]');
title('Measured step response');
legend('measured displacement', 'Location', 'southeast');
grid on;
save_figure(fig, outdir, 'step_response.png');
end


function step_check(G, path, outdir)
% Simulate the step response, overlay it, and plot what is left over.

[t, u, y] = read_log(path);
i0 = step_index(u);
tr = t(i0:end) - t(i0);
du = u(end) - u(1);

model = y(1) + du * step(G, tr);
resid = 1000 * (y(i0:end) - model);
rms = sqrt(mean((y(i0:end) - model).^2));

% The model runs early by about one sample. Shifting it later is the time
% domain view of the delay the frequency test finds.
shifts = 0:0.001:0.3;
bestrms = Inf;
for n = 1:numel(shifts)
    shifted = interp1(tr, model, tr - shifts(n), 'linear', y(1));
    err = sqrt(mean((y(i0:end) - shifted).^2));
    if err < bestrms
        bestrms = err;
        bestshift = shifts(n);
    end
end
shiftresid = 1000 * (y(i0:end) - interp1(tr, model, tr - bestshift, 'linear', y(1)));

[worst, iworst] = max(abs(resid));
worst = resid(iworst);

levels = unique(y);
resolution = 1000 * min(diff(levels));

fig = figure('Visible', 'off', 'Position', [0 0 950 700]);
subplot(2,1,1);
plot(tr, y(i0:end), 'LineWidth', 2, 'Color', [0.17 0.24 0.44]); hold on;
plot(tr, model, '--', 'LineWidth', 1.8, 'Color', [0.7 0.1 0.1]);
ylabel('Displacement [m]');
title('Simulated against measured step response');
legend('measured', 'identified model', 'Location', 'southeast');
grid on;

subplot(2,1,2);
patch([0 max(tr) max(tr) 0], [-resolution -resolution resolution resolution], ...
      [0.85 0.85 0.85], 'EdgeColor', 'none'); hold on;
plot(tr, resid, 'LineWidth', 1.6, 'Color', [0.17 0.24 0.44]);
plot(tr, shiftresid, '--', 'LineWidth', 1.4, 'Color', [0.1 0.5 0.2]);
xlabel('Time since the step [s]'); ylabel('Residual [mm]');
legend('logger resolution', 'measured - model', ...
       sprintf('model shifted %.0f ms later', 1000*bestshift), 'Location', 'southeast');
grid on;
save_figure(fig, outdir, 'validation_step.png');

fprintf('Step response\n');
fprintf('  step of %.2f V at t = %.1f s\n', du, t(i0));
fprintf('  RMS error %.5f m over %d samples\n', rms, numel(tr));
fprintf('  worst residual %.1f mm at %.1f s\n', worst, tr(iworst));
fprintf('  best time shift %.0f ms, which takes the RMS to %.2f mm\n\n', ...
        1000*bestshift, 1000*bestrms);
end


function frequency_check(G, datadir, outdir)
% Fit each sine log at its drive frequency and compare against the model.

files = dir(fullfile(datadir, 'sine_*.csv'));
w = []; magdb = []; phase = [];
for n = 1:numel(files)
    drive = frequency_from_name(files(n).name);
    if isnan(drive)
        continue    % amplitude repeats are not sweep points
    end
    [m, p] = analyse_sine(fullfile(datadir, files(n).name), drive);
    w(end+1) = drive; %#ok<AGROW>
    magdb(end+1) = m; %#ok<AGROW>
    phase(end+1) = p; %#ok<AGROW>
end
[w, order] = sort(w);
magdb = magdb(order);
phase = phase(order);

[mm, mp] = bode_at(G, w);

% Second order cannot lag past 180 degrees, so the extra phase is either a
% third pole or a delay. A pole would bend the magnitude too; a delay fits
% -w*T and leaves it alone.
excess = (phase - mp) * pi/180;
T = -sum(excess .* w) / sum(w .* w);
Gd = tf(G.Numerator{1}, G.Denominator{1}, 'InputDelay', T);
[~, mpd] = bode_at(Gd, w);

fprintf('Frequency response\n');
fprintf('  %8s %10s %10s %10s %10s\n', 'w', '|G| meas', '|G| model', ...
        'ph meas', 'ph model');
for n = 1:numel(w)
    fprintf('  %8.2f %10.2f %10.2f %10.1f %10.1f\n', ...
            w(n), magdb(n), mm(n), phase(n), mpd(n));
end
fprintf('  magnitude %.2f dB RMS\n', sqrt(mean((magdb - mm).^2)));
fprintf('  fitted transport delay %.1f ms\n', 1000*T);
fprintf('  phase %.1f deg RMS without it, %.1f deg RMS with it\n', ...
        sqrt(mean((phase - mp).^2)), sqrt(mean((phase - mpd).^2)));

% The same four parameters from these points alone. The step answer only sets
% where the search starts, so this is a separate reading, not a restatement.
den = G.Denominator{1};
guess = [G.Numerator{1}(end)/den(end), den(2)/(2*sqrt(den(3))), sqrt(den(3)), 0.1];
wide = fit_frequency(w, magdb, phase, guess, 4);
low = fit_frequency(w(w <= 2), magdb(w <= 2), phase(w <= 2), guess, 4);

fprintf('  fitted to these points alone, using no step data:\n');
fprintf('    A %.4f m/V  zeta %.4f  wn %.4f rad/s  T %.1f ms  alpha %.4f\n', ...
        wide(1), wide(2), wide(3), 1000*wide(4), wide(1)*wide(3)^2);
fprintf('    below 2 rad/s, where the gain is the only thing measured:\n');
fprintf('    A %.4f m/V  zeta %.4f  wn %.4f rad/s  T %.1f ms  alpha %.4f\n', ...
        low(1), low(2), low(3), 1000*low(4), low(1)*low(3)^2);

% The same fit with the delay held fixed, to show what naming it is worth.
fprintf('  refitted with the delay held fixed:\n');
fprintf('    %8s %8s %8s %10s %10s\n', 'T [ms]', 'zeta', 'wn', 'mag RMS', 'ph RMS');
for held = [0 0.058 T]
    q = guess;
    q(4) = held;
    q = fit_frequency(w, magdb, phase, q, 3);
    r = freq_residual(q, w, magdb, phase);
    half = numel(w);
    fprintf('    %8.0f %8.3f %8.3f %8.2f dB %8.1f deg\n', 1000*held, q(2), q(3), ...
            sqrt(mean(r(1:half).^2)), sqrt(mean(r(half+1:end).^2)));
end
fprintf('\n');

grid_w = logspace(log10(min(w)/5), log10(max(w)*5), 400);
[gm, gp] = bode_at(G, grid_w);
[~, gpd] = bode_at(Gd, grid_w);

fig = figure('Visible', 'off', 'Position', [0 0 950 750]);
subplot(2,1,1);
semilogx(grid_w, gm, 'LineWidth', 1.6, 'Color', [0.3 0.3 0.3]); hold on;
semilogx(w, magdb, 'o', 'MarkerSize', 7, 'MarkerFaceColor', [0.17 0.24 0.44], ...
         'MarkerEdgeColor', 'none');
xline(sqrt(G.Denominator{1}(3)), '--', 'Color', [0.7 0.1 0.1]);
ylabel('Magnitude [dB]'); grid on;
title('Measured frequency response');
legend('identified model', 'measured', 'Location', 'southwest');

subplot(2,1,2);
semilogx(grid_w, gp, 'LineWidth', 1.6, 'Color', [0.3 0.3 0.3]); hold on;
semilogx(grid_w, gpd, '--', 'LineWidth', 1.6, 'Color', [0.1 0.5 0.2]);
semilogx(w, phase, 'o', 'MarkerSize', 7, 'MarkerFaceColor', [0.17 0.24 0.44], ...
         'MarkerEdgeColor', 'none');
yline(-90, ':', 'Color', [0.7 0.1 0.1]);
yline(-180, ':', 'Color', [0.5 0.5 0.5]);
xlabel('Frequency [rad/s]'); ylabel('Phase [deg]'); grid on;
legend('identified model', sprintf('with %.0f ms delay', 1000*T), 'measured', ...
       'Location', 'southwest');
save_figure(fig, outdir, 'bode_plot.png');

fig = figure('Visible', 'off', 'Position', [0 0 950 700]);
subplot(2,1,1);
semilogx(w, magdb - mm, 'o-', 'LineWidth', 1.4, 'MarkerSize', 6, ...
         'Color', [0.17 0.24 0.44], 'MarkerFaceColor', [0.17 0.24 0.44]); hold on;
yline(0, ':', 'Color', [0.5 0.5 0.5]);
ylabel('Magnitude error [dB]'); grid on;
title('Model against the measured frequency points');

subplot(2,1,2);
semilogx(w, phase - mp, 'o-', 'LineWidth', 1.4, 'MarkerSize', 6, ...
         'Color', [0.17 0.24 0.44], 'MarkerFaceColor', [0.17 0.24 0.44]); hold on;
semilogx(w, phase - mpd, 's--', 'LineWidth', 1.4, 'MarkerSize', 6, ...
         'Color', [0.1 0.5 0.2], 'MarkerFaceColor', [0.1 0.5 0.2]);
yline(0, ':', 'Color', [0.5 0.5 0.5]);
xlabel('Frequency [rad/s]'); ylabel('Phase error [deg]'); grid on;
legend('no delay', sprintf('with %.0f ms delay', 1000*T), 'Location', 'southwest');
save_figure(fig, outdir, 'validation_freq.png');
end


function p = fit_frequency(w, magdb, phase, guess, nfree)
% Gain, damping ratio, natural frequency and delay fitted to the sine points.
% Gauss-Newton with a numerical Jacobian, dB and degrees weighted equally.
% nfree = 3 holds the delay at its guess and fits the plant around it.

p = guess(:);
for iter = 1:200
    r = freq_residual(p, w, magdb, phase);
    J = zeros(numel(r), nfree);
    for i = 1:nfree
        h = 1e-6 * max(abs(p(i)), 1e-3);
        q = p;
        q(i) = q(i) + h;
        J(:,i) = (freq_residual(q, w, magdb, phase) - r) / h;
    end
    delta = J \ (-r);
    p(1:nfree) = p(1:nfree) + delta;
    if max(abs(delta ./ max(abs(p(1:nfree)), 1e-9))) < 1e-12
        break
    end
end
end


function r = freq_residual(p, w, magdb, phase)
% Magnitude and phase residuals of a delayed second order model, stacked.

s = 1j * w(:);
h = p(1)*p(3)^2 ./ (s.^2 + 2*p(2)*p(3)*s + p(3)^2);
mag = 20*log10(abs(h));
ph = unwrap(angle(h))*180/pi - p(4)*w(:)*180/pi;
r = [mag - magdb(:); ph - phase(:)];
end


function [magdb, phasedeg] = bode_at(sys, w)
% Magnitude in dB and phase in degrees at the given frequencies.

% .' not ', since ' conjugates and would flip the sign of every phase
h = squeeze(freqresp(sys, w));
h = h(:).';
magdb = 20*log10(abs(h));
phasedeg = unwrap(angle(h)) * 180/pi;
end


function [magdb, phasedeg] = analyse_sine(path, w)
% One frequency point, fitted by least squares at the known drive frequency.
% The transient is dropped by using only the second half of the record.

[t, u, y] = read_log(path);
start = floor(numel(t)/2) + 1;
t = t(start:end); u = u(start:end); y = y(start:end);

[au, pu] = fit_sinusoid(t, u, w);
[ay, py] = fit_sinusoid(t, y, w);

magdb = 20*log10(ay/au);
phasedeg = (py - pu) * 180/pi;
phasedeg = mod(phasedeg + 180, 360) - 180;
if phasedeg > 10
    % a lowpass only ever lags, so a large positive reading is a wrapped lag
    phasedeg = phasedeg - 360;
end
end


function [amp, ph] = fit_sinusoid(t, y, w)
% Least squares amplitude and phase of y at a known frequency.

basis = [cos(w*t), sin(w*t), ones(size(t))];
c = basis \ y;
amp = hypot(c(1), c(2));
ph = atan2(-c(2), c(1));
end


function [t, u, y] = read_log(path)
% Time, input and output displacement columns of a simulator log.

d = readmatrix(path);
t = d(:,1); u = d(:,2); y = d(:,3);
end


function i0 = step_index(u)
% First sample where the input leaves its resting value.

i0 = find(u ~= u(1), 1, 'first');
if isempty(i0)
    error('simulate_lab1:noStep', 'the input never changes in this log');
end
end


function w = frequency_from_name(name)
% Drive frequency from a name like sine_1p3.csv, where 1p3 means 1.3 rad/s.

[~, stem] = fileparts(name);
text = strrep(extractAfter(stem, 'sine_'), 'p', '.');
w = str2double(text);
end


function save_figure(fig, outdir, name)
% Write one figure into the report figures folder, on white for printing.
% MATLAB runs headless in a dark theme, which prints as a black page, so the
% colours are set back here rather than on every plot.

set(fig, 'Color', 'w');
axes_list = findall(fig, 'Type', 'axes');
for i = 1:numel(axes_list)
    set(axes_list(i), 'Color', 'w', 'XColor', 'k', 'YColor', 'k', ...
        'GridColor', [0.5 0.5 0.5]);
    set(axes_list(i).Title, 'Color', 'k');
end
legend_list = findall(fig, 'Type', 'legend');
for i = 1:numel(legend_list)
    set(legend_list(i), 'Color', 'w', 'TextColor', 'k', 'EdgeColor', [0.6 0.6 0.6]);
end

if ~exist(outdir, 'dir')
    mkdir(outdir);
end
exportgraphics(fig, fullfile(outdir, name), 'Resolution', 150, ...
               'BackgroundColor', 'white');
close(fig);
end
