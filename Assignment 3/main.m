% Parameters
size = 100;
T = linspace(0, 1, size);

% Signal s1
s1 = ones(1, size);
figure(1);
plot(T, s1);

% Signal s2
s2 = zeros(1, size);
s2(1:75) = 1;
s2(76:end) = -1;
figure(2);
plot(T, s2);



% Calculate the basis functions using GM_Bases function
[phi1, phi2] = GM_Bases(s1, s2);

% Plot the obtained basis functions
figure;
subplot(2,1,1);
stem(phi1);
title('Basis Function \phi_1');
subplot(2,1,2);
stem(phi2);
title('Basis Function \phi_2');

% Calculate the signal space representation using signal_space function
[v1, v2] = signal_space(s1, phi1, phi2);

% Plot the signal space representation
figure;
scatter(v1, v2);
title('Signal Space Representation');
xlabel('Projection onto \phi_1');
ylabel('Projection onto \phi_2');

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
    % Calculate the projections (correlations) of s over phi1 and phi2
    v1 = dot(s, phi1);
    v2 = dot(s, phi2);
end

