function data = crg_view(file, opencrgDir)
% CRG_VIEW Open an OpenCRG file and show it in 3D with the OpenCRG MATLAB tools.
%   DATA = CRG_VIEW(FILE) reads FILE with crg_read and shows the road surface
%   in world (x, y, z) coordinates.
%   DATA = CRG_VIEW(FILE, OPENCRGDIR) first adds the toolbox of a local clone
%   of https://github.com/asam-ev/OpenCRG to the path (OPENCRGDIR is the
%   repository root, the one containing the 'matlab' folder).
%
%   Example:
%   data = crg_view('belgian_block.crg', 'C:\src\OpenCRG');

if nargin > 1
    run(fullfile(opencrgDir, 'matlab', 'crg_init.m'));
end
if ~exist('crg_read', 'file')
    error('OpenCRG MATLAB toolbox not on path: call crg_view(file, opencrgDir) or run crg_init');
end

data = crg_read(file);          % read + consistency check

figure('Name', file, 'NumberTitle', 'off');
data = crg_plot_road_xyz_map(data);   % 3D surface over the curved grid
axis tight; view(3); rotate3d on;
title(strrep(file, '_', '\_'));

% for the full set of standard plots (u/v map, x/y/z map, cuts) use instead:
%   crg_show_road_surface(data);
end
