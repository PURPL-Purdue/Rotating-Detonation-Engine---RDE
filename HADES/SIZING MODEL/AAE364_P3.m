clc; clear; close all;

% time vector
t = linspace(0,15,2000);

%% for part a
% transfer functions
G_a  = tf(10, [1 4.5 22 10]);
G1_a = tf(0.5, [1 0.5]);
G2_a = tf(20, [1 4 20]);

% step responses
[yG_a,  tG_a]  = step(G_a,  t);
[yG1_a, tG1_a] = step(G1_a, t);
[yG2_a, tG2_a] = step(G2_a, t);

% plot
figure();
plot(tG_a,  yG_a,  'LineWidth', 1.5); hold on;
plot(tG1_a, yG1_a, 'LineWidth', 1.5);
plot(tG2_a, yG2_a, 'LineWidth', 1.5);
grid on;
xlabel('Time (s)');
ylabel('Response');
title('Unit-Step Responses of G(s), G1(s), and G2(s)');
legend('G(s)', 'G1(s)', 'G2(s)', 'Location', 'Best');

% Step info
disp('Step Info for G(s):');
stepinfo(G_a)

%% for part c
% transfer functions
G_b  = tf(80, [1 10 26 80]);
G1_b = tf(8,  [1 8]);
G2_b = tf(10, [1 2 10]);

%step responses
[yG_b,  tG_b] = step(G_b,  t);
[yG1_b, tG1_b] = step(G1_b, t);
[yG2_b, tG2_b] = step(G2_b, t);

%plot
figure;
plot(tG_b,  yG_b,  'LineWidth', 1.5); hold on;
plot(tG1_b, yG1_b, 'LineWidth', 1.5);
plot(tG2_b, yG2_b, 'LineWidth', 1.5);
grid on;
xlabel('Time (s)'); ylabel('Response');
title('Unit-Step Responses of G(s), G1(s), and G2(s)');
legend('G(s)','G1(s)','G2(s)','Location','Best');

