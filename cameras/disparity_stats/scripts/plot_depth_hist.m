clear all; close all;

subj = 'bwsdrdp6';
flag = '_random';
task = 'walk_1';
save = 1;

dmat = load(['../data/' subj flag '/' subj '_task_' task '_depth.mat']);

dmat.upper = dmat.depth(1:81,:,:);
dmat.lower = dmat.depth(83:163,:,:);

dmat.depth = reshape(dmat.depth,1,163*163*size(dmat.depth,3));
dmat.upper = reshape(dmat.upper,1,81*163*size(dmat.upper,3));
dmat.lower = reshape(dmat.lower,1,81*163*size(dmat.lower,3));

%remove very large and very small values
dmat.depth(dmat.depth > quantile(dmat.depth,.99)) = NaN;
dmat.depth(dmat.depth < quantile(dmat.depth,.01)) = NaN;

dmat.upper(dmat.upper > quantile(dmat.depth,.99)) = NaN;
dmat.upper(dmat.upper < quantile(dmat.depth,.01)) = NaN;

dmat.lower(dmat.lower > quantile(dmat.depth,.99)) = NaN;
dmat.lower(dmat.lower < quantile(dmat.depth,.01)) = NaN;

dmat.depth_median = nanmedian(dmat.depth);
dmat.upper_median = nanmedian(dmat.upper);
dmat.lower_median = nanmedian(dmat.lower);

dmat.depth_mode = mode(round(dmat.depth));
dmat.upper_mode = mode(round(dmat.upper));
dmat.lower_mode = mode(round(dmat.lower));

f1 = figure(); hold on;
[f,xi] = ksdensity(dmat.depth);
plot(xi,f,'k');
plot([dmat.depth_median dmat.depth_median],[0 f(abs(xi-dmat.depth_median) == min(abs(xi-dmat.depth_median)))],'k:','LineWidth',1.5);
plot([xi(f == max(f)) xi(f == max(f))],[0 f(abs(xi-dmat.depth_mode) == min(abs(xi-dmat.depth_mode)))],'k--');

if save
    print(f1,'-depsc',['0bwsdata/depthhist_' subj flag '_' task '_1.eps']);
end

f2 = figure(); hold on;
[f,xi] = ksdensity(dmat.depth);
h(1) = plot(xi,f,'k');
plot([dmat.depth_median dmat.depth_median],[0 f(abs(xi-dmat.depth_median) == min(abs(xi-dmat.depth_median)))],'k:','LineWidth',1.5);
plot([xi(f == max(f)) xi(f == max(f))],[0 f(abs(xi-dmat.depth_mode) == min(abs(xi-dmat.depth_mode)))],'k--');

[f,xi] = ksdensity(dmat.upper);
h(2) = plot(xi,f,'r');
plot([dmat.upper_median dmat.upper_median],[0 f(abs(xi-dmat.upper_median) == min(abs(xi-dmat.upper_median)))],'r:','LineWidth',1.5);
plot([xi(f == max(f)) xi(f == max(f))],[0 f(abs(xi-dmat.upper_mode) == min(abs(xi-dmat.upper_mode)))],'r--');

[f,xi] = ksdensity(dmat.lower);
h(3) = plot(xi,f,'b');
plot([dmat.lower_median dmat.lower_median],[0 f(abs(xi-dmat.lower_median) == min(abs(xi-dmat.lower_median)))],'b:','LineWidth',1.5);
plot([xi(f == max(f)) xi(f == max(f))],[0 f(abs(xi-dmat.lower_mode) == min(abs(xi-dmat.lower_mode)))],'b--');

legend(h,'all','upper','lower');

if save
    print(f2,'-depsc',['0bwsdata/depthhist_' subj flag '_' task '_2.eps']);
end