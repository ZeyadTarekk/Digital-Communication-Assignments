% Req 3
ramp = -6:0.01:6;
% m = 0
q0 = UniformQuantizer(ramp, 3, 6, 0);
dq0 = UniformDequantizer(q0, 3, 6, 0);
figure(1);
hold on
plot(ramp, q0);
plot(ramp,ramp);
plot(ramp,dq0);
title('m = 0');
hold off
% m = 1
q1 = UniformQuantizer(ramp, 3, 6, 1);
dq1 = UniformDequantizer(q1, 3, 6, 1);
figure(2);
hold on
plot(ramp, q1);
plot(ramp,dq1);
title('m = 1');
hold off

% Testing the Quantizer
in_val = [0.5, 1.2, -0.3, 2.1];
n_bits = 3;
xmax = 2;
m = 0;
q_ind = UniformQuantizer(in_val, n_bits, xmax, m);
disp(q_ind);

% Req 1
function q_ind = UniformQuantizer(in_val, n_bits, xmax, m)
    % Calculate the number of quantization levels
    numberOfLevels = 2^n_bits;
    
    % Get the width of each interval
    delta = 2*xmax / numberOfLevels;
    
    % Define the range of the quantizer reconstruction levels
    d = (m*delta) / 2;
    levels = linspace(-xmax + delta/2, xmax - delta/2, numberOfLevels) + d;
    disp('Levels');
    disp(levels);
    
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



function deq_val = UniformDequantizer(q_ind, n_bits, xmax, m)

    % Calculate the number of levels
    numberOfLevels = 2 ^ n_bits;
    % Get the width of each interval
    delta = 2 * xmax / numberOfLevels;
    deq_val = ((q_ind) * delta) + ((m + 1) * (delta / 2) - xmax);
end