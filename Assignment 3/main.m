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