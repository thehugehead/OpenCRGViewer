function data = crg_view(file, varargin)
% CRG_VIEW Open an OpenCRG file, show it in 3D and optionally export it.
%   DATA = CRG_VIEW(FILE) reads FILE with the OpenCRG MATLAB toolbox and
%   shows the road surface in world (x, y, z) coordinates.
%   DATA = CRG_VIEW(FILE, OPENCRGDIR, ...) first adds the toolbox of a local
%   clone of https://github.com/asam-ev/OpenCRG to the path (OPENCRGDIR is
%   the repository root, the one containing the 'matlab' folder).
%
%   Name/value options:
%   'ZScale'     vertical exaggeration of the plot (default 1)
%   'MaxPoints'  max. plotted points per grid direction (default 300)
%   'Export'     write the surface to a file instead of plotting:
%                *.ply  binary point cloud (x, y, z double; u, v scalar fields)
%                *.xyz / *.txt / *.csv  ASCII columns x y z u v
%   'Step'       with 'Export': every STEP-th point in u and v (default 1)
%   'Show'       with 'Export': also show the plot (default false)
%
%   Examples:
%   data = crg_view('belgian_block.crg', 'C:\src\OpenCRG');
%   data = crg_view('belgian_block.crg', 'ZScale', 10);
%   data = crg_view('country_road.crg', 'Export', 'country_road.ply');
%
%   The plot uses a local origin (first point of the reference line) and
%   'axis vis3d' so that rotating and zooming keep the road inside the axes.

if mod(numel(varargin), 2) == 1
    run(fullfile(varargin{1}, 'matlab', 'crg_init.m'));
    varargin(1) = [];
end
opt = struct('zscale', 1, 'maxpoints', 300, 'export', '', 'step', 1, 'show', false);
for k = 1:2:numel(varargin)
    name = lower(varargin{k});
    if ~isfield(opt, name)
        error('crg_view: unknown option ''%s''', varargin{k});
    end
    opt.(name) = varargin{k+1};
end
if ~exist('crg_read', 'file')
    error('OpenCRG MATLAB toolbox not on path: call crg_view(file, opencrgDir) or run crg_init');
end

data = crg_read(file);
if ~isfield(data, 'ok')
    data = crg_check(data);
    if ~isfield(data, 'ok')
        error('crg_view: check of CRG data was not completely successful');
    end
end
[nu, nv] = size(data.z);

% grid axes, as in crg_plot_road_xyz_map
u = data.head.ubeg + (0:nu-1) * data.head.uinc;
if isfield(data.head, 'vinc')
    v = data.head.vmin + (0:nv-1) * data.head.vinc;
else
    v = double(data.v(:)');
end

if ~isempty(opt.export)
    [data, n] = export_surface(data, u, v, opt.export, max(1, round(opt.step)));
    fprintf('wrote %d points to %s\n', n, opt.export);
    if ~opt.show
        return
    end
end

% decimate for plotting only
iu = unique([1:max(1, ceil(nu / opt.maxpoints)):nu, nu]);
iv = unique([1:max(1, ceil(nv / opt.maxpoints)):nv, nv]);
[x, y, z, data] = eval_grid(data, u(iu), v(iv));
[xr, yr, zr, data] = eval_grid(data, u(iu), 0);   % reference line

% local origin avoids precision problems with large world coordinates
x0 = xr(1);
y0 = yr(1);

fig = figure('Name', file, 'NumberTitle', 'off');
ax = axes('Parent', fig);
surf(ax, x - x0, y - y0, z, z, 'EdgeColor', 'none');
hold(ax, 'on');
plot3(ax, xr - x0, yr - y0, zr, 'r-', 'LineWidth', 1);
hold(ax, 'off');

daspect(ax, [1 1 1 / opt.zscale]);
axis(ax, 'tight');
axis(ax, 'vis3d');                     % fixed aspect: no re-fitting while rotating
set(ax, 'ClippingStyle', 'rectangle');  % clip to the axes area when zoomed in
view(ax, 3);
grid(ax, 'on');
colormap(ax, parula);
cb = colorbar(ax);
ylabel(cb, 'z [m]');

xlabel(ax, sprintf('x - %.3f [m]', x0));
ylabel(ax, sprintf('y - %.3f [m]', y0));
if opt.zscale == 1
    zlabel(ax, 'z [m]');
else
    zlabel(ax, sprintf('z [m] (x%g)', opt.zscale));
end
title(ax, sprintf('%s  (%d x %d grid)', strrep(file, '_', '\_'), nu, nv));
rotate3d(fig, 'on');   % double-click restores the initial view
end


function [x, y, z, data] = eval_grid(data, u, v)
% evaluate x, y, z on the (u, v) grid; returns (numel(u) x numel(v)) arrays
[U, V] = ndgrid(u, v);
puv = [U(:) V(:)];
[pxy, data] = crg_eval_uv2xy(data, puv);
[pz, data] = crg_eval_uv2z(data, puv);
x = reshape(pxy(:, 1), size(U));
y = reshape(pxy(:, 2), size(U));
z = reshape(pz, size(U));
end


function [data, n] = export_surface(data, u, v, fname, step)
% write the surface point by point, in blocks of rows to limit memory
u = u(1:step:end);
v = v(1:step:end);
[~, ~, ext] = fileparts(fname);
ext = lower(ext);
if ~any(strcmp(ext, {'.ply', '.xyz', '.txt', '.csv'}))
    error('crg_view: unsupported export format ''%s'' (use .ply, .xyz, .txt or .csv)', ext);
end

% count valid points first: the PLY header needs the vertex count
block = max(1, floor(2e6 / numel(v)));
n = 0;
for i = 1:block:numel(u)
    [x, y, z, data] = eval_grid(data, u(i:min(i+block-1, end)), v);
    n = n + nnz(isfinite(x) & isfinite(y) & isfinite(z));
end

fid = fopen(fname, 'w');
if fid < 0
    error('crg_view: cannot open ''%s'' for writing', fname);
end
cleanup = onCleanup(@() fclose(fid));
if strcmp(ext, '.ply')
    fprintf(fid, ['ply\nformat binary_little_endian 1.0\n' ...
        'comment exported from OpenCRG by crg_view.m\n' ...
        'element vertex %d\nproperty double x\nproperty double y\nproperty double z\n' ...
        'property float u\nproperty float v\nend_header\n'], n);
elseif strcmp(ext, '.csv')
    fprintf(fid, 'x,y,z,u,v\n');
else
    fprintf(fid, '//x y z u v\n');
end

for i = 1:block:numel(u)
    ub = u(i:min(i+block-1, end));
    [x, y, z, data] = eval_grid(data, ub, v);
    [U, V] = ndgrid(ub, v);
    ok = isfinite(x) & isfinite(y) & isfinite(z);
    if strcmp(ext, '.ply')
        % interleave 3 doubles + 2 singles per vertex as raw bytes
        xyz = typecast(reshape([x(ok) y(ok) z(ok)]', [], 1), 'uint8');
        uv  = typecast(single(reshape([U(ok) V(ok)]', [], 1)), 'uint8');
        rec = [reshape(xyz, 24, []); reshape(uv, 8, [])];
        fwrite(fid, rec, 'uint8', 0, 'l');
    else
        sep = ' ';
        if strcmp(ext, '.csv'), sep = ','; end
        fmt = [strjoin(repmat({'%.6f'}, 1, 5), sep) '\n'];
        fprintf(fid, fmt, [x(ok) y(ok) z(ok) U(ok) V(ok)]');
    end
end
end
