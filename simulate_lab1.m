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

G = tf(alpha, [1 b k]);

repeats_check(G, datadir);

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

% Require several resolution steps; one sample above the final value is noise.
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

% Average over the peaks. A single peak still gives half a damped period.
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


function repeats_check(G, datadir)
% Repeats: measurement noise. Harder runs: linearity. The shift column checks
% whether the step onset lands between samples.

runs = {'step1.csv', 'step2.csv', 'step3.csv', 'step_amp2.csv'};
fprintf('Repeats and linearity\n');
fprintf('  %14s %8s %8s %8s %8s %10s %10s\n', ...
        'run', 'A', 'zeta', 'wn', 'Tp', 'step at', 'shift');
for n = 1:numel(runs)
    path = fullfile(datadir, runs{n});
    [A, zeta, wn, Tp] = identify_step(path);
    [t, u] = read_log(path);
    shift = best_shift(G, path);
    fprintf('  %14s %8.4f %8.4f %8.4f %8.2f %8.1f s %7.0f ms\n', ...
            runs{n}, A, zeta, wn, Tp, t(step_index(u)), 1000*shift);
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


function plot_step(path, outdir)
% The measured step, with the features the analysis reads labelled.

[t, u, y] = read_log(path);
i0 = step_index(u);
tail = max(1, floor(0.1*numel(y)));
y0 = mean(y(1:i0-1));
yinf = mean(y(end-tail+1:end));
tstep = t(i0);

resp = y(i0:end);
[ypeak, ipeak] = max(resp);
tpeak = t(i0 + ipeak - 1);

fig = figure('Visible', 'off', 'Position', [0 0 950 600]);
plot(t, y, 'LineWidth', 2, 'Color', [0.17 0.24 0.44]); hold on;
yline(y0, ':', 'Color', [0.7 0.1 0.1]);
yline(yinf, '--', 'Color', [0.7 0.1 0.1]);
xline(tstep, ':', 'Color', [0.35 0.35 0.35]);
plot(tpeak, ypeak, 'v', 'MarkerSize', 8, 'MarkerFaceColor', [0.7 0.1 0.1], ...
     'MarkerEdgeColor', 'none');

% The marking guide asks for the points of interest labelled on the graph, and
% the worked example in the brief puts the value next to each one.
span = ypeak + 0.06*(ypeak - y0);
plot([tstep tpeak], [span span], '-', 'LineWidth', 1.2, 'Color', [0.35 0.35 0.35]);
plot([tstep tstep], span + [-1 1]*0.015*(ypeak - y0), '-', 'Color', [0.35 0.35 0.35]);
plot([tpeak tpeak], span + [-1 1]*0.015*(ypeak - y0), '-', 'Color', [0.35 0.35 0.35]);
text(mean([tstep tpeak]), span, sprintf('T_p = %.2f s', tpeak - tstep), ...
     'Color', [0.35 0.35 0.35], 'FontSize', 11, ...
     'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');

% Overshoot drawn as the span it measures, from the settled value to the peak.
plot([tpeak tpeak], [yinf ypeak], '-', 'LineWidth', 1.2, 'Color', [0.7 0.1 0.1]);
text(tpeak, mean([yinf ypeak]), sprintf('  OS = %.0f mm', 1000*(ypeak - yinf)), ...
     'Color', [0.7 0.1 0.1], 'FontSize', 11);
% Gain drawn as the span between the levels, left of the step where it is flat.
tgain = tstep - 0.45*(tstep - t(1));
plot([tgain tgain], [y0 yinf], '-', 'LineWidth', 1.2, 'Color', [0.1 0.4 0.2]);
cap = 0.025*(tstep + 20 - t(1));
plot(tgain + [-1 1]*cap, [y0 y0], '-', 'LineWidth', 1.2, 'Color', [0.1 0.4 0.2]);
plot(tgain + [-1 1]*cap, [yinf yinf], '-', 'LineWidth', 1.2, 'Color', [0.1 0.4 0.2]);
text(tgain, mean([y0 yinf]), sprintf(' A = %.3f m/V', yinf - y0), ...
     'Color', [0.1 0.4 0.2], 'FontSize', 11);

tend = tstep + 20;
text(tend, y0, sprintf(' x_0 = %.3f m', y0), 'Color', [0.7 0.1 0.1], ...
     'FontSize', 11, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');
text(tend, yinf, sprintf(' x_\\infty = %.3f m', yinf), 'Color', [0.7 0.1 0.1], ...
     'FontSize', 11, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top');

% Stop at 20 s; past that the trace is flat to the resolution.
xlim([t(1) tend]);
ylim([y0 - 0.04*(ypeak - y0), span + 0.10*(ypeak - y0)]);
xlabel('Time [s]'); ylabel('Displacement [m]');
title('Measured step response');
grid on;

% The 6 mm undershoot vanishes on a 0.9 m axis; it gets an inset.
[trough, itrough] = min(resp(ipeak:end));
ttrough = t(i0 + ipeak + itrough - 2);
inset = axes('Position', [0.53 0.30 0.34 0.30]);
plot(inset, t, 1000*(y - yinf), 'LineWidth', 1.6, 'Color', [0.17 0.24 0.44]);
hold(inset, 'on');
yline(inset, 0, '--', 'Color', [0.7 0.1 0.1]);
plot(inset, [ttrough ttrough], [0 1000*(trough - yinf)], '-', ...
     'LineWidth', 1.2, 'Color', [0.7 0.1 0.1]);
text(ttrough, 500*(trough - yinf), sprintf('  %.0f mm', 1000*(yinf - trough)), ...
     'Parent', inset, 'Color', [0.7 0.1 0.1], 'FontSize', 10);
xlim(inset, [tpeak, tpeak + 8]);
ylim(inset, [1000*(trough - yinf) - 3, 4]);
xlabel(inset, 'Time [s]'); ylabel(inset, 'x - x_\infty [mm]');
title(inset, 'Undershoot, second reading of \zeta', 'FontSize', 10);
grid(inset, 'on');
box(inset, 'on');

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

[bestshift, bestrms] = best_shift(G, path);
shiftresid = 1000 * (y(i0:end) - interp1(tr, model, tr - bestshift, 'linear', y(1)));

[worst, iworst] = max(abs(resid));
worst = resid(iworst);

levels = unique(y);
resolution = 1000 * min(diff(levels));

% Same 20 s window as the step figure; past that both traces are flat.
tend = 20;

fig = figure('Visible', 'off', 'Position', [0 0 950 700]);
subplot(2,1,1);
plot(tr, y(i0:end), 'LineWidth', 2, 'Color', [0.17 0.24 0.44]); hold on;
plot(tr, model, '--', 'LineWidth', 1.8, 'Color', [0.7 0.1 0.1]);
ylabel('Displacement [m]');
title('Simulated against measured step response');
legend('measured', 'identified model', 'Location', 'southeast');
xlim([0 tend]);
grid on;

subplot(2,1,2);
patch([0 tend tend 0], [-resolution -resolution resolution resolution], ...
      [0.85 0.85 0.85], 'EdgeColor', 'none'); hold on;
plot(tr, resid, 'LineWidth', 1.6, 'Color', [0.17 0.24 0.44]);
plot(tr, shiftresid, '--', 'LineWidth', 1.4, 'Color', [0.1 0.5 0.2]);

% Mark the worst residual; the report quotes it.
plot(tr(iworst), worst, 'o', 'MarkerSize', 9, 'LineWidth', 1.5, 'Color', [0.7 0.1 0.1]);
text(tr(iworst), worst, sprintf('  %.1f mm at %.1f s', worst, tr(iworst)), ...
     'Color', [0.7 0.1 0.1], 'FontSize', 11, 'VerticalAlignment', 'middle');

xlabel('Time since the step [s]'); ylabel('Residual [mm]');
legend('logger resolution', 'measured - model', ...
       sprintf('model shifted %.0f ms later', 1000*bestshift), 'Location', 'southeast');
xlim([0 tend]);
grid on;
save_figure(fig, outdir, 'validation_step.png');

fprintf('Step response\n');
fprintf('  step of %.2f V at t = %.1f s\n', du, t(i0));
fprintf('  RMS error %.5f m over %d samples\n', rms, numel(tr));
fprintf('  worst residual %.1f mm at %.1f s\n', worst, tr(iworst));
fprintf('  best time shift %.0f ms, which takes the RMS to %.2f mm\n', ...
        1000*bestshift, 1000*bestrms);
fprintf('  worst sample after the shift %.1f mm, largest overshoot the other way %.1f mm\n\n', ...
        max(abs(shiftresid)), max(resid));
end


function frequency_features(G, Gd, w, magdb, phase, mp)
% The readings the report takes off the Bode plot by hand: the DC gain the
% low frequency points imply, the damping ratio from the resonant peak, the
% phase crossings and the high frequency slope.

den = G.Denominator{1};
A = G.Numerator{1}(end)/den(end);
zeta = den(2)/(2*sqrt(den(3)));
wn = sqrt(den(3));

fprintf('Read off the frequency response\n');

% Each low point sits above DC by a known lift; dividing it out recovers
% what the run would read at zero frequency.
low = w <= 0.5;
r = w(low)/wn;
lift = 1 ./ sqrt((1 - r.^2).^2 + (2*zeta*r).^2);
amp = 10.^(magdb/20);
fprintf('  %8s %10s %10s %10s\n', 'w', 'lift', '|G| meas', 'implied A');
for n = find(low)
    fprintf('  %8.2f %9.2f%% %10.4f %10.4f\n', ...
            w(n), 100*(lift(n) - 1), amp(n), amp(n)/lift(n));
end

% Resonant peak. Which gain it is measured against decides whether the
% reading is independent of the step test, so both are given.
[peak, ipeak] = max(amp);
adc = amp(1)/lift(1);
for ref = [amp(1), adc]
    mr = peak/ref;
    z = sqrt((1 - sqrt(1 - 1/mr^2))/2);
    fprintf('  peak %.4f at %.2f rad/s over %.4f: Mr %.4f, zeta %.4f, wn %.3f rad/s\n', ...
            peak, w(ipeak), ref, mr, z, w(ipeak)/sqrt(1 - 2*z^2));
end

fprintf('  phase crosses -90 deg at %.3f rad/s measured, %.3f model, %.3f with the delay\n', ...
        cross90(w, phase), crossing_of(G), crossing_of(Gd));

% Slope over the points well above the corner, against the model on the same
% points and against the model's own local slope on the way out to them.
top = w > 2*wn;
fprintf('  slope over the %d points above 2 wn: %.1f dB/decade measured, %.1f model\n', ...
        sum(top), slope(w(top), magdb(top)), slope(w(top), bode_at(G, w(top))));
for f = [1 1.5 2 9.33]
    fprintf('    model local slope at %.2f wn: %.1f dB/decade\n', f, ...
            slope(f*wn*[0.99 1.01], bode_at(G, f*wn*[0.99 1.01])));
end

% A pure delay is a constant number of degrees per rad/s. This one is not
% quite constant, which is the qualification the report makes.
excess = mp - phase;
for n = [1, find(abs(w - wn) == min(abs(w - wn)), 1), numel(w)]
    fprintf('    excess phase at %.2f rad/s: %.2f deg per rad/s\n', ...
            w(n), excess(n)/w(n));
end
fprintf('\n');
end


function s = slope(w, magdb)
% dB per decade through the points, by least squares on log frequency.
p = polyfit(log10(w(:)), magdb(:), 1);
s = p(1);
end


function wc = cross90(w, phase)
% Where the phase passes -90 degrees, interpolated between the two points
% either side of it.
wc = interp1(phase, w, -90);
end


function wc = crossing_of(sys)
% Same crossing for a model, found on a fine grid rather than between points.
grid = logspace(-1, 1.2, 4000);
[~, p] = bode_at(sys, grid);
wc = interp1(p, grid, -90);
end


function [bestshift, bestrms] = best_shift(G, path)
% The pure time shift that fits the model to a step log best.

[t, u, y] = read_log(path);
i0 = step_index(u);
tr = t(i0:end) - t(i0);
model = y(1) + (u(end) - u(1)) * step(G, tr);

shifts = 0:0.001:0.3;
bestrms = Inf;
bestshift = 0;
for n = 1:numel(shifts)
    shifted = interp1(tr, model, tr - shifts(n), 'linear', y(1));
    err = sqrt(mean((y(i0:end) - shifted).^2));
    if err < bestrms
        bestrms = err;
        bestshift = shifts(n);
    end
end
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

% The excess phase fits -w*T; a pole would bend the magnitude too.
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

frequency_features(G, Gd, w, magdb, phase, mp);

% Fit all four from these points alone; the step answer only seeds the search.
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

% Same fit with the delay pinned at each candidate value.
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

% Only a little past the measured points, or the data gets squashed.
grid_w = logspace(log10(min(w)/1.5), log10(max(w)*1.5), 400);
[gm, gp] = bode_at(G, grid_w);
[~, gpd] = bode_at(Gd, grid_w);

fig = figure('Visible', 'off', 'Position', [0 0 950 750]);
subplot(2,1,1);
semilogx(grid_w, gm, 'LineWidth', 1.6, 'Color', [0.3 0.3 0.3]); hold on;
semilogx(w, magdb, 'o', 'MarkerSize', 7, 'MarkerFaceColor', [0.17 0.24 0.44], ...
         'MarkerEdgeColor', 'none');
wn = sqrt(G.Denominator{1}(3));
xline(wn, '--', 'Color', [0.7 0.1 0.1]);
ylabel('Magnitude [dB]'); grid on;
title('Measured frequency response');
legend('identified model', 'measured', 'Location', 'southwest');

% Same as the step figure: the marking guide asks for these on the graph.
% The extra headroom is so the peak label has somewhere to sit.
[mpeak, ipeak] = max(magdb);
ylim([min(gm) - 2, mpeak + 8]);
text(wn, min(magdb), sprintf(' \\omega_n = %.2f rad/s', wn), ...
     'Color', [0.7 0.1 0.1], 'FontSize', 11, 'VerticalAlignment', 'bottom');

% Peak labelled as the ratio to the lowest measured point, as in the report.
text(w(ipeak), mpeak, sprintf(' M_r = %.3f', 10^((mpeak - magdb(1))/20)), ...
     'Color', [0.7 0.1 0.1], 'FontSize', 11, 'VerticalAlignment', 'bottom');

subplot(2,1,2);
semilogx(grid_w, gp, 'LineWidth', 1.6, 'Color', [0.3 0.3 0.3]); hold on;
semilogx(grid_w, gpd, '--', 'LineWidth', 1.6, 'Color', [0.1 0.5 0.2]);
semilogx(w, phase, 'o', 'MarkerSize', 7, 'MarkerFaceColor', [0.17 0.24 0.44], ...
         'MarkerEdgeColor', 'none');
yline(-90, ':', 'Color', [0.7 0.1 0.1]);
yline(-180, ':', 'Color', [0.5 0.5 0.5]);
% Bound by the measurements, not by where the delayed model runs off.
ylim([min(phase) - 25, 5]);

% Mark the -90 crossing with both frequencies; the gap is the delay argument.
w90 = cross90(w, phase);
semilogx(w90, -90, 'p', 'MarkerSize', 14, 'MarkerFaceColor', [0.7 0.1 0.1], ...
         'MarkerEdgeColor', 'none');
text(w90, -103, sprintf('  measured at %.2f rad/s, \\omega_n = %.2f', w90, wn), ...
     'Color', [0.7 0.1 0.1], 'FontSize', 11);

text(grid_w(1), -180, ' -180^\circ limit', 'Color', [0.4 0.4 0.4], ...
     'FontSize', 11, 'VerticalAlignment', 'bottom');
xlabel('Frequency [rad/s]'); ylabel('Phase [deg]'); grid on;
legend('identified model', sprintf('with %.0f ms delay', 1000*T), 'measured', ...
       'Location', 'southwest');
save_figure(fig, outdir, 'bode_plot.png');

fig = figure('Visible', 'off', 'Position', [0 0 950 700]);
% Draw the RMS magnitude error as a band the points can be read against.
magerr = magdb - mm;
magrms = sqrt(mean(magerr.^2));

subplot(2,1,1);
patch([min(w) max(w) max(w) min(w)], [-magrms -magrms magrms magrms], ...
      [0.87 0.87 0.87], 'EdgeColor', 'none'); hold on;
set(gca, 'XScale', 'log');
semilogx(w, magerr, 'o-', 'LineWidth', 1.4, 'MarkerSize', 6, ...
         'Color', [0.17 0.24 0.44], 'MarkerFaceColor', [0.17 0.24 0.44]);
yline(0, ':', 'Color', [0.5 0.5 0.5]);
text(min(w), magrms, sprintf('  %.2f dB RMS', magrms), 'FontSize', 11, ...
     'Color', [0.35 0.35 0.35], 'VerticalAlignment', 'bottom');
ylabel('Magnitude error [dB]'); grid on;
title('Model against the measured frequency points');
legend('RMS band', 'measured - model', 'Location', 'southwest');

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

basis = [cos(w*t), sin(w*t), ones(size(t))];
c = basis \ y;
amp = hypot(c(1), c(2));
ph = atan2(-c(2), c(1));
end


function [t, u, y] = read_log(path)

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
