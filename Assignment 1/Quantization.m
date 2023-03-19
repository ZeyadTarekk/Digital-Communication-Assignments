close all;
clear all;

% ---------------- Req 3 ----------------
ramp = -6:0.01:6;

% m = 0
q = UniformQuantizer(ramp, 3, 6, 0);
dq = UniformDequantizer(q, 3, 6, 0);
figure();
hold on
plot(ramp, q);
plot(ramp, ramp);
plot(ramp, dq);
title('m = 0 - Midrise (Ramp)');
legend('Quantizer', 'Ramp Signal', 'Dequantizer');
hold off

% m = 1
q = UniformQuantizer(ramp, 3, 6, 1);
dq = UniformDequantizer(q, 3, 6, 1);
figure();
hold on
plot(ramp, q);
plot(ramp, ramp);
plot(ramp, dq);
title('m = 1 - Midtread (Ramp)');
legend('Quantizer', 'Ramp Signal', 'Dequantizer');
hold off

% ---------------- Req 4 ----------------
% Continuous Uniform Random Variables between -5 & 5
N = 10000;
lowerBound = -5;
upperBound = 5;
input = unifrnd(lowerBound, upperBound, 1, N);

% Set xmax = 5 & m to midrise as given in the requirement
xmax = 5;
m = 0;

% Define the number of bits array & get its length
n_bits = [2, 3, 4, 5, 6, 7, 8];
len = length(n_bits);

% Initialize Theoretical and Simulated SNR
theoreticalSNR = zeros(1, len);
simulatedSNR = zeros(1, len);

% Define the levels array (2 ^ n) & calculate P
Levels = 2 .^ n_bits;
P = mean(input .^ 2);

for i = 1:len
    % Get the dequantized signal with corresponding inputs
    q = UniformQuantizer(input, n_bits(i), xmax, m);
    dq = UniformDequantizer(q, n_bits(i), xmax, m);
    
    % Calculate the Theoretical SNR for each number of bits
    theoreticalSNR(i) = 10 * log10(P / (((xmax) .^ 2) / (3 * ((Levels(i) .^ 2)))));
    
    % Get the quantization error by subtracting input values & the
    % dequantized signal
    quantizedError = input - dq;
    
    % Calculate the simulated SNR using the quantized error
    simulatedSNR(i) = 10 * log10(P / mean(quantizedError .^ 2));
end

% Sketching the simulation & Theoretical SNR
figure();
hold on
plot(n_bits, theoreticalSNR, 'B-');
plot(n_bits, simulatedSNR, 'Ro');
xlabel('Number of Bits');
ylabel('SNR (in dB)');
title('Uniform Random Variables');
legend('Theoretical SNR', 'Simulation');
hold off

% ---------------- Req 5 ----------------
% Define polarity & give it a value +/- with probability 0.5
polarity = 2 * randi([0 1], 1, N) - 1;

% Sample magnitude following an exponential distribution
magnitude = exprnd(1, 1, N);

% Apply the random polarity to magnitude
input = magnitude .* polarity;

% Re-initialize Theoretical & Simulated SNR
theoreticalSNR = zeros(1, len);
simulatedSNR = zeros(1, len);

% Calculate signal power as the variance of the input
signalPower = var(input);

% Re-define P & xmax accordingly (m is still 0)
P = mean(input .^ 2);
xmax = max(abs(input));

for i = 1:len
    % Get the dequantized signal with corresponding inputs
    q = UniformQuantizer(input, n_bits(i), xmax, m);
    dq = UniformDequantizer(q, n_bits(i), xmax, m);
    
    % Calculate the error power
    errorPower = mean((dq - input) .^ 2);
    
    % Get theoratical SNR (similar to Req 4)
    theoreticalSNR(i) = 10 * log10(P / (((xmax) .^ 2) / (3 * ((Levels(i) .^ 2)))));
    
    % Define simulated SNR as signal power divided by error power
    simulatedSNR(i) = 10 * log10(signalPower / errorPower);
end

% Plot
figure();
hold on
plot(n_bits, theoreticalSNR, 'B-');
plot(n_bits, simulatedSNR, 'R--');
xlabel('Number of Bits');
ylabel('SNR (in dB)');
title('Non-uniform Random Input');
legend('Theoretical SNR', 'Simulation');
hold off

% ---------------- Req 6 ----------------
% Non-uniform mu quantization
mu = [0 ,5, 100, 200];
colors = ['r' , 'b' , 'g' , 'k'];

% Normalize the Input
input_N = input / xmax;

figure();
hold on
for j = 1:length(mu)
    % Initialize Theoretical and Simulated SNR
    simulatedSNR = n_bits;
    theoreticalSNR = zeros(1, len);
   
    for i = 1:len
        % Compress the Signal if mu is greater than zero
        if (mu(j) > 0)
            y = polarity .* (log(1+mu(j) * abs(input_N)) / log(1+mu(j)));
        else
            y = input_N;
        end
        
        ymax = max(abs(y));
        
        q = UniformQuantizer(y, n_bits(i), ymax, m);
        dq = UniformDequantizer(q, n_bits(i), ymax, m);

        if (mu(j) > 0)
            z = polarity .*(((1+mu(j)) .^ abs(dq)-1) / mu(j));
        else
            z =  dq;
        end
        
        ex = z * xmax;
        
        error = abs(input - ex);
        simulatedSNR(i) = 10*log10(mean(input .^ 2) / mean(error .^ 2));
        
        if (mu(j) > 0)
            theoreticalSNR(i) = 10*log10((3*(Levels(i) .^ 2))/((log(1 + mu(j))) .^ 2));
        else
            theoreticalSNR(i) = 10*log10(P / (((xmax) .^ 2) / (3 * ((Levels(i) .^ 2)))));
        end
    end
    
    % Plot
    plot(n_bits, theoreticalSNR, sprintf('%s-' , colors(j))  , 'LineWidth', 1);
    plot(n_bits, simulatedSNR, sprintf('%s--' , colors(j)) , 'LineWidth', 1);
end    
hold off

legend('mu=0 (T)', 'mu=0 (S)', 'mu=5 (T)', 'mu=5 (S)', 'mu=100 (T)', 'mu=100 (S)', ...
       'mu=200 (T)', 'mu=200 (S)');

% ---------------- Req 1 ----------------
function q_ind = UniformQuantizer(in_val, n_bits, xmax, m)
    % Calculate the number of quantization levels
    numberOfLevels = 2^n_bits;
    
    % Get the width of each interval
    delta = 2*xmax / numberOfLevels;
    
    % Define the range of the quantizer reconstruction levels
    d = ((1 - m)*delta) / 2;
    levels = (d - xmax):delta:(d + xmax);
    
    % Get size of the input
    size = length(in_val);
    
    % Initialize quantized indices
    q_ind = zeros(1, size);
    
    % Find the closest reconstruction level to the input sample
    for i = 1:size
        M = abs(in_val(i) - levels);
        q_ind(i) = find(M == min(M), 1);
    end
end

% ---------------- Req 2 ----------------
function deq_val = UniformDequantizer(q_ind, n_bits, xmax, m)
    % Calculate the number of levels
    numberOfLevels = 2 ^ n_bits;
    
    % Calculate Delta
    delta = 2 * xmax / numberOfLevels;
    
    % Get the Levels
    levels = (((1 - m)*delta / 2) - xmax):delta:(((1 - m)*delta / 2) + xmax);
    
    % Restore each level to original amplitude
    deq_val = levels(q_ind);
end