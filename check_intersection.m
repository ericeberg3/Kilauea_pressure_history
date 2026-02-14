function is_intersecting = check_intersection(m, n_sources)
% CHECK_INTERSECTION Checks if bounding spheres of sources intersect
%
% Inputs:
%   m:                  Current model vector [src1_params; src2_params; ...]
%   n_sources:          Number of spheroid sources (e.g., 2)
%   n_params_per_source: Number of parameters per source (e.g., 8)
%
% Output:
%   is_intersecting:    true if any sources overlap, false otherwise

    % Define indices for geometry (Modify these if your m structure differs)
    % Assuming typical structure: [x, y, z, a, b/a, strike, dip, P]
    idx_xyz = 5:7; % Indices for coordinates (relative to source start)
             % Index for semi-major axis (bounding radius)

    is_intersecting = false;

    % Loop through pairs of sources to check overlap
    for i = 1:n_sources
        % Extract parameters for source i
        start_i = (i-1) * 8;
        [~, idx_a] = max(m(start_i+1:start_i + 2));
        center_i = m(start_i + idx_xyz);
        radius_i = m(start_i + idx_a);

        for j = i+1:n_sources
            % Extract parameters for source j
            start_j = (j-1) * 8;
            center_j = m(start_j + idx_xyz);
            radius_j = m(start_j + idx_a);

            % Calculate Euclidean distance between centers
            dist = norm(center_i - center_j);

            % Check if distance is less than sum of radii (Overlap condition)
            if dist < (radius_i + radius_j)
                is_intersecting = true;
                return; % Exit early if intersection found
            end
        end
    end
end