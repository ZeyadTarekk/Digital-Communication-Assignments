clc;
clear all;
close all;

% Time axis
N = 20;
T = linspace(0, 1, N);

% Input signal s1(t)
s1 = ones(1, N);

% Input signal s2(t)
s2 = zeros(1, N);

% Setting from 0 to 0.75 as 1, else -1 
s2(1:15) = 1;
s2(16:end) = -1;

% Plot s1(t)
figure;
plot(T, s1);
title('Input Signal s1(t)');
xlabel('Time in seconds');
ylabel('Amplitude');

% Plot s2(t)
figure;
plot(T, s2);
title('Input Signal s2(t)');
xlabel('Time in seconds');
ylabel('Amplitude');

% Get phi1 & phi2 from the GM_Bases function
[phi1, phi2] = GM_Bases(s1, s2);

% Plot phi1
figure;
subplot(2,1,1);
plot(T, phi1);
title('Basis Function \phi_1');
xlabel('Time in seconds');
ylabel('Amplitude');

% Plot phi2
subplot(2,1,2);
plot(T, phi2);
title('Basis Function \phi_2');
xlabel('Time in seconds');
ylabel('Amplitude');

% Calculate the signal space representation using signal_space function
% For s1(t)
[s1_v1, s1_v2] = signal_space(s1, phi1, phi2);
% For s2(t)
[s2_v1, s2_v2] = signal_space(s2, phi1, phi2);

% Plot the signal space representation we just obtained
figure;
scatter(s1_v1, s1_v2, 'filled', 'MarkerFaceColor', 'red', 'DisplayName', 's1');
hold on;
scatter(s2_v1, s2_v2, 'filled', 'MarkerFaceColor', 'blue', 'DisplayName', 's2');
title('Signal Space Representation');
xlabel('Projection onto \phi_1(t)');
ylabel('Projection onto \phi_2(t)');
legend('Location', 'best');
axis([-0.1 1.1 -0.1 1]);
grid on;

% Line Plot
figure;
scatter(s1_v1, s1_v2, 'filled', 'MarkerFaceColor', 'red', 'DisplayName', 's1');
hold on;
scatter(s2_v1, s2_v2, 'filled', 'MarkerFaceColor', 'blue', 'DisplayName', 's2');
plot([0, s1_v1], [0, s1_v2], 'r', 'LineWidth', 1.5);
plot([0, s2_v1], [0, s2_v2], 'b', 'LineWidth', 1.5);
title('Signal Space Representation');
xlabel('Projection onto \phi_1(t)');
ylabel('Projection onto \phi_2(t)');
legend('S1', 'S2');
axis([-0.1 1.1 -0.1 1]);
grid on;

% Random Samples
% Set the random number generator seed for reproducibility
rng(0);

% Number of samples
numberOfSamples = 100;

% Variance values (dB) for different noise levels
values = [-5, 0, 10];

% Set the length of values
L = length(values);

% Get size of s1 & s2
sz_s1 = size(s1);
sz_s2 = size(s2);

% Get length of s1 & s2
len_s1 = length(s1);
len_s2 = length(s2);

for i = 1:L
    r1V1 = [];
    r1V2 = [];
    r2V1 = [];
    r2V2 = [];
    for sample = 1:numberOfSamples
        % Convert dB to linear scale
        val = 10^(values(i)/10);

        % Add noise to the original signal
        r1 = awgn(s1, val);
        r2 = awgn(s2, val);

        % Calculate the signal points of the input using signal_space function
        [r1_v1, r1_v2] = signal_space(r1, phi1, phi2);
        [r2_v1, r2_v2] = signal_space(r2, phi1, phi2);
        
        % Append to lists
        r1V1(end + 1) = r1_v1;
        r1V2(end + 1) = r1_v2;
        r2V1(end + 1) = r2_v1;
        r2V2(end + 1) = r2_v2;
    end
    disp(values(i));
    % Plot the signal points of generated samples of r1 & r2
    figure;
    % Plot r1 samples
    scatter(r1V1, r1V2, 'DisplayName', 'r1');
    hold on;
    % Plot r2 samples
    scatter(r2V1, r2V2, 'DisplayName', 'r2');
    hold on
    % Plot the input signal s1(t)
    scatter(s1_v1, s1_v2, 100, 'filled', 'MarkerFaceColor', 'red', 'DisplayName', 's1');
    hold on
    % Plot the input signal s2(t)
    scatter(s2_v1, s2_v2, 100, 'filled', 'MarkerFaceColor', 'blue', 'DisplayName', 's2');
    % Title/xLabel/yLabel
    title(sprintf('Signal Points (sigma^2 = %ddB)', values(i)));
    xlabel('Projection onto \phi_1');
    ylabel('Projection onto \phi_2');
    legend('r1', 'r2', 's1', 's2');
    axis([-0.1 1.1 -0.1 1.1]);
end

function [phi1, phi2] = GM_Bases(s1, s2)
    % Getting each signal length
    len_s1 = length(s1);
    len_s2 = length(s2);
    
    % Signal Energy
    sq_s1 = s1.^2;
    energy = sum(sq_s1);
    disp(energy);
    
    % Calculate the first basis function (phi1)
    phi1 = s1 / sqrt(energy);
    phi1 = phi1 * sqrt(len_s1);
    
    % Calculate s21
    s21 = sum(s2.*phi1) / len_s2;
    
    % g = s2 - (s21)(phi1)
    g2 = s2 - s21.*phi1;
    
    % Energy for s2 then we get phi2 using g2 like we did in the tutorial
    energy2 = sum(g2.^2);
    phi2 = g2 / sqrt(energy2);
    phi2 = phi2 * sqrt(len_s2);

    % Check if s2 is linearly independent from s1
    %if dot(s2, phi1) ~= 0
        % Calculate the second basis function (phi2)
       % phi2 = s2 - dot(s2, phi1) * phi1;
       % phi2 = phi2 / norm(phi2);
    %else
        % s2 is linearly dependent on s1, so phi2 is a zero vector
       % phi2 = zeros(size(s2));
    %end
    
    % Multiplying by the square root of signal length
    % phi1 = phi1 * sqrt(len_s1);
    % phi2 = phi2 * sqrt(len_s2);
end

function [v1, v2] = signal_space(s, phi1, phi2)
    % Get the length of each basis function
    phi1_len = length(phi1);
    phi2_len = length(phi2);
    
    % Set new phi1 after dividing by the square root of their lengths
    phi1 = phi1 / sqrt(phi1_len);
    phi2 = phi2 / sqrt(phi2_len);
    
    % Calculating the dot product between the signal and phi just like the
    % lecture
    dotProduct_1 = dot(s, phi1);
    dotProduct_2 = dot(s, phi2);
    
    % Calculate the projections (correlations) of s over phi1 and phi2
    % Then divide by the square root of each phi
    v1 = dotProduct_1 / sqrt(phi1_len);
    v2 = dotProduct_2 / sqrt(phi2_len);
end