%% Save the results

%% RES1 folder 
% pRMSE = [0.4833, 0.4829, 0.3997, 0.3409, 0.3773, 0.1304, 0.3478, 0.3211, 0.1915, 0.3922, 0.2957, 0.3592, 0.3503 ...
%          0.2895, 0.3632, 0.3613, 0.3817, 0.2613, 0.1221];
% nRMSE = [0.5625, 0.4409, NaN, 0.3821, 0.4427, 0.3897, 0.3869, 0.3067, 0.3804, 0.3067, 0.3037, 0.3510, 0.2800 ...
%         0.3233, 0.2800, 0.4224, 0.5001, 0.3763, 0.3912]; 
% pMAE = [0.8400, 1.1200, 0.8400, 0.8400, 0.8854, 0.5600, 0.7920, 0.8854, 0.7920, 0.8854, 0.8400, 1.1200, 0.8400 ...
%         0.8400, 1.1200, 0.8400, 0.8400, 0.8400, 0.5600]; 
% nMAE = [1.12, 0.6261, NaN, 0.5600, 0.5600, 0.5600, 0.5600, 0.3960, 0.5600, 0.3960, 0.5600, 0.5600, 0.2800, ...
%         0.5600, 0.2800, 0.5600, 1.1200, 0.7920, 0.7920];
% N_pos = 
% N_neg = 
%% RES1 figure
load('res1.mat');
pRMSE = []; 
nRMSE = []; 
pMAE = [];
nMAE = [];
N_pos = []; 
N_neg = []; 

for i = 1:length(CELL)
if isnan(CELL{i}.mSRE_positive_mm)
    pRMSE(i) = 0;
    pMAE(i) = 0; 
    N_pos(i) = 0;
else
    pRMSE(i) = CELL{i}.mSRE_positive_mm;
    pMAE(i) = CELL{i}.max_positive_mm;
    N_pos(i) = CELL{i}.N_positive;
end

if isnan(CELL{i}.mSRE_negative_mm) 
    nRMSE(i) = 0;
    nMAE(i) = 0;
    N_neg(i) = 0;
else
    nRMSE(i) = CELL{i}.mSRE_negative_mm;
    nMAE(i) = CELL{i}.max_negative_mm;
    N_neg(i) = CELL{i}.N_negative;
end

end

figure(1);clf;
plot(pRMSE, nRMSE, 'rx'); 
xlim([0,2]);
ylim([0,2]);
xlabel('pRMSE (mm)', 'Fontsize', 20);
ylabel('nRMSE (mm)', 'Fontsize', 20);
grid on;

figure(2);clf;
plot(pMAE, nMAE, 'rx'); 
xlim([0,2]);
ylim([0,2]);
xlabel('pMAE (mm)', 'Fontsize', 20);
ylabel('nMAE (mm)', 'Fontsize', 20);
grid on;

figure(3);clf
% nbins = 20;
% [histFreq, histXout] = hist(N_pos, 1:19);
bar(1:length(N_pos), N_pos);
% xlabel('Case index');
xlabel('Case index', 'Fontsize', 20);
ylabel('Number of positive (negative) pixels', 'Fontsize', 20);
hold on;
bar(1:length(N_pos), N_neg);
legend({'positive pixels', 'negative pixels'}, 'FontSize',20);
ylim([0,1000]);
% ylabel('Number of negative pixels', 'Fontsize', 20);
grid on;

% histogram(N_pos, 20); 
% h = histogram(N_pos, 1:19,'Normalization', 'probability');
% ytix = get(gca, 'YTick')
% set(gca, 'YTick',ytix, 'YTickLabel',ytix*100);

%% RES2 folde
% pRMSE = [0.1236, 0.4083, 0.1280, 0.4891, 0.2608, 0.2268, 0.3572, 0.2558, 0.4794, 0.2800];
% nRMSE = [0.5291, 0.3762, 0.2800, 0.6204, 0.2858, 0.3399, 0.4575, 0.2886, 0.3202, 0.2800]; 
% pMAE = [0.5600, 1.1545, 0.5600, 1.1200, 0.5600, 0.7920, 1.0096, 0.7920, 0.8400, 0.6261]; 
% nMAE = [1.2522, 0.6261, 0.2800, 1.1200, 0.3960, 0.6261, 1.0096, 0.3960, 0.5600, 0.2800];

%% RES2 figure

load('res2.mat');
pRMSE = []; 
nRMSE = []; 
pMAE = [];
nMAE = [];
N_pos = []; 
N_neg = []; 

for i = 1:length(CELL)
if isnan(CELL{i}.mSRE_positive_mm)
    pRMSE(i) = 0;
    pMAE(i) = 0; 
    N_pos(i) = 0;
else
    pRMSE(i) = CELL{i}.mSRE_positive_mm;
    pMAE(i) = CELL{i}.max_positive_mm;
    N_pos(i) = CELL{i}.N_positive;
end

if isnan(CELL{i}.mSRE_negative_mm) 
    nRMSE(i) = 0;
    nMAE(i) = 0;
    N_neg(i) = 0;
else
    nRMSE(i) = CELL{i}.mSRE_negative_mm;
    nMAE(i) = CELL{i}.max_negative_mm;
    N_neg(i) = CELL{i}.N_negative;
end

end

figure(1);hold on;
plot(pRMSE, nRMSE, 'bx'); 
xlim([0,2]);
ylim([0,2]);
xlabel('pRMSE (mm)', 'Fontsize', 20);
ylabel('nRMSE (mm)', 'Fontsize', 20);
axis equal;
grid on;

figure(2);hold on;
plot(pMAE, nMAE, 'bx'); 
xlim([0,2]);
ylim([0,2]);
xlabel('pMAE (mm)', 'Fontsize', 20);
ylabel('nMAE (mm)', 'Fontsize', 20);
grid on;
axis equal;

figure(4);clf
% nbins = 20;
% [histFreq, histXout] = hist(N_pos, 1:19);
bar(1:length(N_pos), N_pos);
% xlabel('Case index');
xlabel('Case index', 'Fontsize', 20);
ylabel('Number of positive (negative) pixels', 'Fontsize', 20);
hold on;
bar(1:length(N_pos), N_neg);
legend({'positive pixels', 'negative pixels'}, 'FontSize',20);
ylim([0,1000]);
% ylabel('Number of negative pixels', 'Fontsize', 20);
grid on;
