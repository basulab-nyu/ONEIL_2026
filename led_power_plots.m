
path = ['/Volumes/KO_Portable/MEC_LEC_CA3/analysis/mW']; cd(path); 
load('led_power.mat');

figure;
plot(blue(:,3),1, 'b|', 'LineWidth',2); hold on;
plot(red(:,3),1, 'r|', 'LineWidth',2);

figure;
plot(blue(1:9,3), 1,'b|', 'LineWidth',2); hold on;
plot(red(1:10,3), 1,'r|', 'LineWidth',2);

figure;
plot(blue(10:19,3),1, 'b|', 'LineWidth',2); hold on;
plot(red(11:19,3),1, 'r|', 'LineWidth',2);

figure; 
h = histogram(blue(1:9,3), 'BinWidth',0.05, 'BinLimits', [0.0326 0.4488]); hold on; 
h2 = histogram(red(1:10,3), 'BinWidth',0.05, 'BinLimits', [0.0326 0.4488], 'FaceAlpha', 0.5);

figure;
plot(blue(:,3),1, 'b|', 'LineWidth',2); hold on;
plot(red(:,3),1.1, 'r|', 'LineWidth',2);
plot(combo(:,3),1.2, 'm|', 'LineWidth',2);
ylim([0.5 1.5])
