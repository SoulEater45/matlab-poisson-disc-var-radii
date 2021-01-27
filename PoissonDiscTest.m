clear;
close all;
rng('default')
dimensions = [256, 256, 256] * 2;
radii = [10 20 30 40];
% dimensions = [1024, 786];
% radii = [5 10 20 30];
k = 50;
o = 10;
s = 0;
c = 1;
vis = 0;

sradii.name = 'normal';
sradii.mu = 15;
sradii.sigma = 5;
% sradii.name = 'unif';
% sradii.lower = 10;
% sradii.upper = 30;
% sradii.name = 'poisson';
% sradii.lambda = 15;
sradii.rmin = 10;
sradii.rmax = 30;

if logical(s)
    dimensions = dimensions / 2;
end

tic;
[p, r] = poissondisc(dimensions, sradii, ...
    'repetitions', k, ...
    'additionaloffset', o, ...
    'spherical', s, ...
    'centered', c, ...
    'visualise', vis);
t = toc;
fprintf('Time requiered: %f s\n', t);
fprintf('Resulting speed: %f pts/s\n', length(r)/t);
%% Plot results if not already done by the function
if all([~vis any(length(dimensions)==[2,3])])
    figure;
    cmap = @(sampleR) interp1(linspace(min(r), max(r), 256)', jet, sampleR);
    if length(dimensions) == 2
        axis([0, dimensions(1), 0, dimensions(2)].*(1+logical(s)) - ...
            dimensions([1 1 2 2]) * c * (1+logical(s)) / 2 + ...
            [-1, 1, -1, 1] * max(r) * 1.5);
        for j = 1:length(r)
            rectangle('Position', [p(j,:) 3*r(j) 3*r(j)]-r(j), ...
                'Curvature', [1 1], ...
                'FaceColor', cmap(r(j)));
        end
    elseif length(dimensions) == 3
        axis([0, dimensions(1), 0, dimensions(2), 0, dimensions(3)].*(1+logical(s)) - ...
            dimensions([1 1 2 2 3 3]) * c * (1+logical(s)) / 2 + ...
            [-1, 1, -1, 1, -1, 1] * max(r) * 1.5);
        [xS, yS, zS] = sphere;
        for j = 1:length(r)
            surf(...
                p(j,1) + xS .* r(j), ...  % Shift and scale x data
                p(j,2) + yS .* r(j), ...  % Shift and scale y data
                p(j,3) + zS .* r(j), ...  % Shift and scale z data
                'EdgeColor', 'none', ...
                'FaceColor', cmap(r(j)));
            hold on;
        end
        clear xS yS zS;
    end
    axis equal;
end
%% Plot Distribution of the radii
fld = ["name", "rmax", "rmin"];
params = string(fields(sradii));
params = params(~contains(params, fld));
eval(sprintf("pd = truncate(makedist('%s' %s), sradii.rmin, sradii.rmax);", sradii.name, sprintf(", '%s', sradii.%s", [params, params]')));
figure;

switch lower(sradii.name)
    case 'poisson'
        x = floor(sradii.rmin):ceil(sradii.rmax);
        plot(x, poisspdf(x, sradii.lambda));
    otherwise
        x = linspace(sradii.rmin, sradii.rmax, (sradii.rmax-sradii.rmin)*10);
        plot(x, pdf(pd, x));
end
hold on;
histogram(r,'Normalization','pdf');
clear x fld params;