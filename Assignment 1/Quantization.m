close all;

% ---------------- Req 3 ----------------
ramp = -6:0.01:6;

% m = 0
q0 = UniformQuantizer(ramp, 3, 6, 0);
dq0 = UniformDequantizer(q0, 3, 6, 0);
figure(1);
hold on
plot(ramp, q0);
plot(ramp, ramp);
plot(ramp, dq0);
title('m = 0 (Midrise)');
legend('Quantizer', 'Ramp Signal', 'Dequantizer');
hold off

% m = 1
q1 = UniformQuantizer(ramp, 3, 6, 1);
dq1 = UniformDequantizer(q1, 3, 6, 1);
figure(2);
hold on
plot(ramp, q1);
plot(ramp, ramp);
plot(ramp, dq1);
title('m = 1 (Midtread)');
legend('Quantizer', 'Ramp Signal', 'Dequantizer');
hold off

% ---------------- Req 4 ----------------
randomSequence = unifrnd(-5,5,1,10000);

xmax = 5;
m = 0;
n_bits = [2, 3, 4, 5, 6, 7, 8];
theoreticalSNR = zeros(1, length(n_bits));
simulatedSNR = zeros(1, length(n_bits));
Levels = 2 .^ n_bits;
P = mean(randomSequence .^ 2);

for i = 1:length(n_bits)
    q2 = UniformQuantizer(randomSequence, n_bits(i), xmax, m);
    dq2 = UniformDequantizer(q2, n_bits(i), xmax, m);
    theoreticalSNR(i) = 10 * log10(P / (((xmax) .^ 2) / (3 * ((Levels(i) .^ 2)))));
    quantizedError = randomSequence - dq2;
    simulatedSNR(i) = 10 * log10(P / mean(quantizedError .^ 2));
end

figure(3);
hold on
plot(n_bits, theoreticalSNR, 'B-');
plot(n_bits, simulatedSNR, 'Ro');
xlabel('Number of Bits');
ylabel('Theoretical SNR (in dB)');
title('Uniform Random Variables');
legend('Theoretical SNR', 'Simulation');
hold off

% ---------------- Req 5 ----------------
polarity = 2 * randi([0 1], 1, 10000) - 1; % +/- with probability 0.5
randomSequence_E = exprnd(0, 1, 10000) .* polarity;

theoreticalSNR_E = zeros(1, length(n_bits));
simulatedSNR_E = zeros(1, length(n_bits));
signalPower = var(randomSequence_E);
P = mean(randomSequence_E .^ 2);

for i = 1:length(n_bits)
    q = UniformQuantizer(randomSequence_E, n_bits(i), xmax, m);
    dq = UniformDequantizer(q, n_bits(i), xmax, m);
    errorPower = mean(dq - randomSequence_E) ^ 2;
    theoreticalSNR_E(i) = 10 * log10(P / (((xmax) .^ 2) / (3 * ((Levels(i) .^ 2)))));
    simulatedSNR_E(i) = signalPower / errorPower;
end

figure(4);
hold on
plot(n_bits, theoreticalSNR_E, 'B-');
plot(n_bits, simulatedSNR_E, 'Ro');
xlabel('Number of Bits');
ylabel('Theoretical SNR (in dB)');
title('Non-uniform Random Input');
legend('Theoretical SNR', 'Simulation');
hold off

% ---------------- Req 6 ----------------

% ---------------- Req 1 ----------------
function q_ind = UniformQuantizer(in_val, n_bits, xmax, m)
    % Calculate the number of quantization levels
    numberOfLevels = 2^n_bits;
    
    % Get the width of each interval
    delta = 2*xmax / numberOfLevels;
    
    % Define the range of the quantizer reconstruction levels
    d = (m*delta) / 2;
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
    levels = ((m*delta / 2) - xmax):delta:((m*delta / 2) + xmax);
    
    % Restore each level to original amplitude
    deq_val = levels(q_ind);
end