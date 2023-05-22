clear all;
close all;

% Time axis
N = 100;
T = linspace(0, 1, N);

% Input signal s1(t)
s1 = ones(1, N);

% Input signal s2(t)
s2 = zeros(1, N);

% Setting from 0 to 0.75 as 1, else -1 
s2(1:75) = 1;
s2(76:end) = -1;

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
legend('Location', 'best');
axis([-0.1 1.1 -0.1 1]);
grid on;

% Random Samples
% Set the random number generator seed for reproducibility
rng('default');

% Number of samples
numberOfSamples = 100;

% Variance values (dB) for different noise levels
values = [-5, 0, 10];

% Set the length of values
L = length(values);

% Get size of s1 & s2
len_s1 = size(s1);
len_s2 = size(s2);

for i = 1:L
    r1V1 = [];
    r1V2 = [];
    r2V1 = [];
    r2V2 = [];
    for sample = 1:numberOfSamples
        % Convert dB to linear scale
        val = 10^(values(i)/10);

        % Calculate standard deviation
        sigma = 1 / sqrt(val);

        % Noise Signals
        w1 = sqrt(sigma) * randn(len_s1);
        w2 = sigma * randn(len_s2);

        % Generate samples of r1(t) and r2(t) using s1(t) and s2(t)
        r1 = s1 + w1;
        r2 = s2 + w2;

        % Calculate the signal points of the input using signal_space function
        [r1_v1, r1_v2] = signal_space(r1, phi1, phi2);
        [r2_v1, r2_v2] = signal_space(r2, phi1, phi2);
        
        % Append to lists
        r1V1(end + 1) = r1_v1;
        r1V2(end + 1) = r1_v2;
        r2V1(end + 1) = r2_v1;
        r2V2(end + 1) = r2_v2;
    end
    
    % Plot the signal points of generated samples of r1 & r2
    figure;
    scatter(r1V1, r1V2, 'DisplayName', 'r1');
    hold on;
    scatter(r2V1, r2V2, 'DisplayName', 'r2');
    hold on
    scatter(s1_v1, s1_v2, 150, 'filled', 'MarkerFaceColor', 'red', 'DisplayName', 's1');
    hold on
    scatter(s2_v1, s2_v2, 150, 'filled', 'MarkerFaceColor', 'blue', 'DisplayName', 's2');
    title(sprintf('Signal Points (sigma^2 = %ddB)', values(i)));
    xlabel('Projection onto \phi_1');
    ylabel('Projection onto \phi_2');
end

function [phi1, phi2] = GM_Bases(s1, s2)
    % Calculate the first basis function (phi1)
    phi1 = s1 / norm(s1);

    % Check if s2 is linearly independent from s1
    if dot(s2, phi1) ~= 0
        % Calculate the second basis function (phi2)
        phi2 = s2 - dot(s2, phi1) * phi1;
        phi2 = phi2 / norm(phi2);
    else
        % s2 is linearly dependent on s1, so phi2 is a zero vector
        phi2 = zeros(size(s2));
    end
end

function [v1, v2] = signal_space(s, phi1, phi2)
    % Get the length of each basis function
    phi1_len = length(phi1);
    phi2_len = length(phi2);
    
    % Set new phi1 after dividing by the square root of their lengths
    phi1 = phi1 / sqrt(phi1_len);
    phi2 = phi2 / sqrt(phi2_len);
    
    % Calculate the projections (correlations) of s over phi1 and phi2
    v1 = dot(s, phi1);
    v2 = dot(s, phi2);
end