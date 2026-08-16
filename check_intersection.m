function is_intersecting = check_intersection(m, n_sources)
% CHECK_INTERSECTION Checks if bounding spheres of sources intersect
%
% Inputs:
%   m:          Model vector [src1_params; src2_params; ...]
%   n_sources:  Number of spheroid sources (e.g., 2)
%
% Output:
%   is_intersecting: true if any two bounding spheres overlap
%
% Per-source layout (8 params): [a, b, ..., x, y, z, P]
%   semi-axes at local indices 1:2, centroid at local indices 5:7.
% Bounding radius = larger of the two semi-axes (the semi-major axis).

    idx_xyz = 5:7;
    is_intersecting = false;

    % Precompute centroid and bounding radius for every source
    centers = zeros(n_sources, 3);
    radii   = zeros(n_sources, 1);
    for s = 1:n_sources
        start_s      = (s-1) * 8;
        centers(s,:) = m(start_s + idx_xyz);
        radii(s)     = max(m(start_s + 1 : start_s + 2));  % own argmax
    end

    % Check all unique pairs
    for i = 1:n_sources
        for j = i+1:n_sources
            dist = norm(centers(i,:) - centers(j,:));
            if dist < (radii(i) + radii(j))
                is_intersecting = true;
                return;
            end
        end
    end
end