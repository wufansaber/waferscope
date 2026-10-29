%% Independent planar near-field to far-field transform (MATLAB R2016b+)
% 独立程序：不调用 SCATMECH、不连接任何 UI。
% 适用条件：固定 z 的平面位于散射体上方；其上方为均匀、无源、
% 无损介质，待变换场只向 +z 出射。输出为该上半空间的远场。
% 输入：8 列 x y Ex_re Ex_im Ey_re Ey_im Ez_re Ez_im，或
%       9 列 x y z Ex_re Ex_im Ey_re Ey_im Ez_re Ez_im。
% 场单位须为 V/m；坐标单位在下面显式指定。每个文件只含单一频率。
% 这是近远场后处理，不会修正原始电磁求解、有限窗口或周期边界误差。

clear; clc; close all;

%% ---------------- USER INPUT / 用户设置 ------------------------------
file_name = '100nm_P_NP_2um.txt';
reference_file = '';            % 可选：匹配的无颗粒参考文件；留空则不相减
coordinate_unit = 'm';          % 'm'、'mm'、'um'、'nm'；必须与 TXT 坐标一致
lam = 266e-9;                   % 真空波长 [m]
n_out = 1.0;                    % 上方均匀、无损、非磁介质的实折射率
input_exp_plus_iwt = true;      % COMSOL exp(+i*w*t): true；exp(-i*w*t): false
z_plane_m = 0;                  % 仅 8 列文件使用：采样平面的真实 z 坐标 [m]
                               % 9 列文件直接读取 z 并检查是否为同一平面
dx = 'auto'; dy = 'auto';       % 'auto' 或手动指定步长 [m]
coordinate_tolerance_m = 1e-15; % 合并同一坐标的导出舍入误差 [m]
grid_relative_tolerance = 1e-6; % 等间距网格容差，相对于 dx/dy

pad_factor = 4;                 % FFT 零填充；只增加频谱采样，不增加物理信息
max_fft_memory_gib = 1.5;       % FFT 工作数组的保守内存规划上限
transform_method = 'fft';       % 'fft': FFT+线性插值；'direct': 逐方向直接积分
                               % direct 较慢，但不需要频谱插值，也不使用 pad_factor
enable_window = false;
alpha_tukey = 0.2;              % 开窗后不使用所谓“精确能量修正”

theta_max_deg = 85;             % 严格小于 90；有限平面不宜直接外推至精确掠射
theta_step_deg = 0.5;
phi_samples = 360;              % 周期方位角网格，不重复 -180/+180 端点
incident_intensity_w_m2 = NaN;  % 入射强度；已知时填正数，才计算 DSC [m^2/sr]
observation_distance_m = NaN;   % 可选：远场距离；填正数才输出该距离处 E [V/m]
output_root = '';              % 留空：在近场文件旁创建 farfield_output
save_results = true;
make_plots = true;

%% ---------------- 参数和输入文件 ------------------------------------
assert(isscalar(lam) && isreal(lam) && isfinite(lam) && lam > 0, 'lam 必须为有限正标量。');
assert(isscalar(n_out) && isreal(n_out) && isfinite(n_out) && n_out > 0, 'n_out 必须为正实数；本程序不处理吸收介质远场。');
assert(isscalar(theta_max_deg) && isfinite(theta_max_deg) && theta_max_deg > 0 && theta_max_deg < 90, 'theta_max_deg 必须在 (0,90) 内。');
assert(isscalar(theta_step_deg) && isfinite(theta_step_deg) && theta_step_deg > 0 && theta_step_deg <= theta_max_deg, 'theta_step_deg 无效。');
assert(isscalar(phi_samples) && isfinite(phi_samples) && phi_samples >= 8 && phi_samples == fix(phi_samples), 'phi_samples 必须为不小于 8 的整数。');
assert(isscalar(pad_factor) && isfinite(pad_factor) && pad_factor >= 1, 'pad_factor 必须 >= 1。');
assert(isscalar(max_fft_memory_gib) && isfinite(max_fft_memory_gib) && max_fft_memory_gib > 0, 'FFT 内存上限必须为正数。');
assert(isscalar(coordinate_tolerance_m) && isfinite(coordinate_tolerance_m) && coordinate_tolerance_m >= 0, '坐标容差必须为有限非负数。');
assert(isscalar(grid_relative_tolerance) && isfinite(grid_relative_tolerance) && grid_relative_tolerance > 0, '网格容差必须为正数。');
assert(isscalar(alpha_tukey) && isfinite(alpha_tukey) && alpha_tukey >= 0 && alpha_tukey <= 1, 'alpha_tukey 必须在 [0,1] 内。');
assert(isscalar(incident_intensity_w_m2) && (isnan(incident_intensity_w_m2) || (isfinite(incident_intensity_w_m2) && incident_intensity_w_m2 > 0)), '入射强度必须为正数或 NaN。');
assert(isscalar(observation_distance_m) && (isnan(observation_distance_m) || (isfinite(observation_distance_m) && observation_distance_m > 0)), '观察距离必须为正数或 NaN。');
assert(any(strcmpi(transform_method, {'fft','direct'})), 'transform_method 只能是 fft 或 direct。');

script_folder = fileparts(mfilename('fullpath'));
file_name = ff_resolve_file(file_name, script_folder);
if isempty(file_name)
    [selected_name, selected_path] = uigetfile({'*.txt;*.csv;*.dat','Near-field data'}, '选择近场文件');
    if isequal(selected_name, 0), error('未选择近场文件。'); end
    file_name = fullfile(selected_path, selected_name);
end
near = ff_read_plane(file_name, coordinate_unit, z_plane_m, dx, dy, ...
    coordinate_tolerance_m, grid_relative_tolerance);

if ~isempty(reference_file)
    reference_file = ff_resolve_file(reference_file, fileparts(file_name));
    assert(~isempty(reference_file), '找不到指定的参考文件。');
    ref = ff_read_plane(reference_file, coordinate_unit, z_plane_m, dx, dy, ...
        coordinate_tolerance_m, grid_relative_tolerance);
    tol = max(coordinate_tolerance_m, grid_relative_tolerance * min(near.dx, near.dy));
    assert(isequal(size(near.Ex), size(ref.Ex)), '样品与参考网格尺寸不匹配。');
    assert(max(abs(near.x-ref.x)) <= tol && max(abs(near.y-ref.y)) <= tol && abs(near.z-ref.z) <= tol, '样品与参考必须使用相同的采样平面和坐标。');
    near.Ex = near.Ex-ref.Ex;
    near.Ey = near.Ey-ref.Ey;
    near.Ez = near.Ez-ref.Ez;
    fprintf('处理复差场：样品 - 参考（没有拟合幅度或整体相位）。\n');
else
    fprintf('处理导出的原始复场；未自动减去入射场或晶圆背景场。\n');
end
if input_exp_plus_iwt
    near.Ex = conj(near.Ex); near.Ey = conj(near.Ey); near.Ez = conj(near.Ez);
end

[Nx, Ny] = size(near.Ex);
lambda_medium = lam/n_out;
k = 2*pi/lambda_medium;
eta = 376.730313668/n_out; % 非磁、无损介质，峰值相量约定
fprintf('网格 %d x %d，dx=%.6g m，dy=%.6g m，z=%.6g m。\n', Nx, Ny, near.dx, near.dy, near.z);
fprintf('窗口离散积分宽度：%.6g x %.6g um。\n', Nx*near.dx*1e6, Ny*near.dy*1e6);
fprintf('假定采样面上方均匀、无源，且场向 +z 出射；本程序不从文件推断几何条件。\n');

%% ---------------- 有限窗口及可选加窗 --------------------------------
W = ones(Nx, Ny);
if enable_window
    W = ff_tukey(Nx, alpha_tukey) * ff_tukey(Ny, alpha_tukey).';
    fprintf('启用 Tukey 窗 alpha=%.4g。输出属于窗口化场，不进行能量或峰值补偿。\n', alpha_tukey);
end
Ex = near.Ex .* W; Ey = near.Ey .* W; Ez = near.Ez .* W;
near_norm = abs(near.Ex).^2+abs(near.Ey).^2+abs(near.Ez).^2;
assert(sum(near_norm(:)) > 0, '输入场或样品-参考差场全为零。');
edge = false(Nx, Ny);
bx = max(1, ceil(.05*Nx)); by = max(1, ceil(.05*Ny));
edge(1:bx,:) = true; edge(end-bx+1:end,:) = true;
edge(:,1:by) = true; edge(:,end-by+1:end) = true;
edge_norm_fraction = sum(near_norm(edge))/sum(near_norm(:));
window_norm = abs(Ex).^2+abs(Ey).^2+abs(Ez).^2;
window_norm_ratio = sum(window_norm(:))/sum(near_norm(:));
fprintf('边缘 5%% 区域占平面 |E|^2 积分 %.4g%%（不是截断误差界）。\n', 100*edge_norm_fraction);

%% ---------------- 观察方向 ------------------------------------------
theta_deg = (0:theta_step_deg:theta_max_deg).';
if theta_deg(end) < theta_max_deg, theta_deg(end+1,1) = theta_max_deg; end
phi_deg = -180+(0:phi_samples-1)*360/phi_samples;
theta = theta_deg*pi/180; phi = phi_deg*pi/180;
[Th, Ph] = ndgrid(theta, phi);
KXq = k*sin(Th).*cos(Ph); KYq = k*sin(Th).*sin(Ph);
assert(k*sin(theta(end)) < min(pi/near.dx, pi/near.dy), '采样间距不能支持指定角度范围，请减小 dx/dy 或 theta_max_deg。');

%% ---------------- 平面傅里叶积分 ------------------------------------
% 围绕窗口中心插值，避免数组起点带来的快速线性相位；
% 插值之后在查询方向上精确恢复全局坐标及 z 平面的相位。
xc = (near.x(1)+near.x(end))/2;
yc = (near.y(1)+near.y(end))/2;
x_local = near.x-xc; y_local = near.y-yc;
transverse_residual = NaN;
Nfx = NaN; Nfy = NaN;
if strcmpi(transform_method, 'fft')
    Nfx = 2^nextpow2(Nx*pad_factor);
    Nfy = 2^nextpow2(Ny*pad_factor);
    planning_gib = Nfx*Nfy*16*12/2^30;
    assert(planning_gib <= max_fft_memory_gib, 'FFT 工作数组规划 %.3g GiB 超过上限；请减少 pad_factor 或使用 direct。', planning_gib);
    kx = 2*pi*(-floor(Nfx/2):ceil(Nfx/2)-1).'/(Nfx*near.dx);
    ky = 2*pi*(-floor(Nfy/2):ceil(Nfy/2)-1).'/(Nfy*near.dy);
    phase = exp(-1i*kx*x_local(1)) * exp(-1i*ky.'*y_local(1));
    Ax = fftshift(fft2(Ex, Nfx, Nfy))*near.dx*near.dy .* phase;
    Ay = fftshift(fft2(Ey, Nfx, Nfy))*near.dx*near.dy .* phase;
    Az = fftshift(fft2(Ez, Nfx, Nfy))*near.dx*near.dy .* phase;
    clear phase;
    fprintf('FFT 零填充：%d x %d -> %d x %d。\n', Nx, Ny, Nfx, Nfy);
    supported_k = min([max(kx), -min(kx), max(ky), -min(ky)]);
    assert(k*sin(theta(end)) <= supported_k, '查询角度超出实际 FFT 网格，请增加填充或减小角度范围。');
    Axq = interpn(kx, ky, Ax, KXq, KYq, 'linear', NaN);
    Ayq = interpn(kx, ky, Ay, KXq, KYq, 'linear', NaN);

    % 使用测量的 Ez 做加权横波残差诊断；不开启“通过即证明精度”的判定。
    [KX, KY] = ndgrid(kx, ky);
    prop = KX.^2+KY.^2 < k^2*(1-1e-10);
    kz_prop = sqrt(k^2-KX(prop).^2-KY(prop).^2);
    residual = KX(prop).*Ax(prop)+KY(prop).*Ay(prop)+kz_prop.*Az(prop);
    denominator = k^2*sum(abs(Ax(prop)).^2+abs(Ay(prop)).^2+abs(Az(prop)).^2);
    if denominator > 0
        transverse_residual = sqrt(sum(abs(residual).^2)/denominator);
    end
    fprintf('传播区加权横波残差：%.6g；有限窗口及加窗也会影响它。\n', transverse_residual);
else
    fprintf('直接计算每个方向的离散平面积分，无 FFT 插值；方向较多时耗时较长。\n');
    Axq = ff_direct(Ex, x_local, y_local, KXq, KYq, near.dx, near.dy);
    Ayq = ff_direct(Ey, x_local, y_local, KXq, KYq, near.dx, near.dy);
end
assert(all(isfinite(Axq(:))) && all(isfinite(Ayq(:))), '角谱存在非有限值或越界查询，停止输出。');

%% ---------------- 远场振幅与有单位的功率 ------------------------------
phase_global = exp(-1i*(KXq*xc+KYq*yc+k*cos(Th)*near.z));
Ftheta = -1i*k/(2*pi) .* (Axq.*cos(Ph)+Ayq.*sin(Ph)) .* phase_global;
Fphi = -1i*k/(2*pi) .* cos(Th) .* (-Axq.*sin(Ph)+Ayq.*cos(Ph)) .* phase_global;
% E_far(r) = exp(i*k*r)/r * (Ftheta*e_theta + Fphi*e_phi).
% 若近场 E 的单位为 V/m，则 F 的单位为 V。
F_squared = abs(Ftheta).^2+abs(Fphi).^2; % V^2
peak_F_squared = max(F_squared(:));
assert(isfinite(peak_F_squared) && peak_F_squared > 0, '远场全为零或出现非有限值。');
U = F_squared/peak_F_squared; % 只用于相对方向图
dP_dOmega = F_squared/(2*eta); % W/sr，不是采样面能量密度
P_cap = trapz(theta, sum(dP_dOmega,2)*(2*pi/phi_samples).*sin(theta));
DSC = NaN(size(U));
if ~isnan(incident_intensity_w_m2)
    DSC = dP_dOmega/incident_intensity_w_m2; % m^2/sr
end

[~, peak_index] = max(U(:));
[peak_t, peak_p] = ind2sub(size(U), peak_index);
peak_theta_deg = theta_deg(peak_t); peak_phi_deg = phi_deg(peak_p);
if peak_theta_deg == 0, peak_phi_deg = NaN; end % 法线方向没有唯一方位角
i0 = ff_phi_index(phi, 0); i90 = ff_phi_index(phi, pi/2);
i180 = ff_phi_index(phi, pi); i270 = ff_phi_index(phi, -pi/2);
beam_xz = ff_beamwidth(theta_deg, U(:,i0), U(:,i180));
beam_yz = ff_beamwidth(theta_deg, U(:,i90), U(:,i270));
fprintf('\n峰值方向 theta=%.5g deg，phi=%.5g deg（法线时 phi 无定义）。\n', peak_theta_deg, peak_phi_deg);
fprintf('远场振幅平方峰值：%.8g V^2。\n', peak_F_squared);
fprintf('角功率密度峰值：%.8g W/sr。\n', max(dP_dOmega(:)));
fprintf('0 <= theta <= %.5g deg 的角积分功率：%.8g W。\n', theta_max_deg, P_cap);
fprintf('该积分只覆盖上述角帽；有限窗口/加窗结果不能自动视为真实总散射功率。\n');
ff_print_beamwidth('XZ', beam_xz); ff_print_beamwidth('YZ', beam_yz);
if isnan(incident_intensity_w_m2)
    fprintf('未指定入射强度，DSC 未计算（导出为 NaN）。\n');
end

%% ---------------- 保留完整复场结果 ----------------------------------
settings = struct('wavelength_vacuum_m',lam, 'n_out',n_out, ...
    'input_exp_plus_iwt',input_exp_plus_iwt, 'coordinate_unit',coordinate_unit, ...
    'transform_method',transform_method, 'pad_factor',pad_factor, ...
    'enable_window',enable_window, 'alpha_tukey',alpha_tukey, ...
    'theta_max_deg',theta_max_deg, 'theta_step_deg',theta_step_deg, ...
    'phi_samples',phi_samples, 'incident_intensity_w_m2',incident_intensity_w_m2, ...
    'observation_distance_m',observation_distance_m);
result = struct('settings',settings, 'source_file',file_name, ...
    'reference_file',reference_file, 'theta_deg',theta_deg, 'phi_deg',phi_deg, ...
    'Ftheta_V',Ftheta, 'Fphi_V',Fphi, 'F_squared_V2',F_squared, ...
    'normalized_intensity',U, 'dP_dOmega_W_sr',dP_dOmega, 'DSC_m2_sr',DSC, ...
    'integrated_cap_power_W',P_cap, 'peak_theta_deg',peak_theta_deg, ...
    'peak_phi_deg',peak_phi_deg, 'beam_xz',beam_xz, 'beam_yz',beam_yz, ...
    'x_m',near.x, 'y_m',near.y, 'z_plane_m',near.z, 'dx_m',near.dx, 'dy_m',near.dy, ...
    'FFT_size',[Nfx,Nfy], 'edge_field_norm_fraction',edge_norm_fraction, ...
    'window_field_norm_ratio',window_norm_ratio, 'transverse_weighted_residual',transverse_residual);
result.conventions = 'exp(-i*omega*t); global phase origin (0,0,0); theta from +Z; phi from +X toward +Y; peak-amplitude phasors';
result.limitations = {'One-way field in the homogeneous upper half-space required.', ...
    'Finite aperture and spatial sampling need independent convergence checks.', ...
    'Windowing changes the modeled field; no energy compensation is applied.', ...
    'No automatic incident-field subtraction or periodic-copy removal.', ...
    'Output covers an angular cap, not the full sphere.'};
if ~isnan(observation_distance_m)
    result.Etheta_V_m = Ftheta*exp(1i*k*observation_distance_m)/observation_distance_m;
    result.Ephi_V_m = Fphi*exp(1i*k*observation_distance_m)/observation_distance_m;
end

%% ---------------- 绘图 ----------------------------------------------
figures = gobjects(0);
if make_plots
    if strcmpi(transform_method, 'fft')
        figures(end+1) = figure('Name','Plane spectrum including evanescent components','Color','w');
        sx = abs(kx/k) <= 1.5; sy = abs(ky/k) <= 1.5;
        spectrum_power = abs(Ax(sx,sy)).^2+abs(Ay(sx,sy)).^2;
        spectrum_max = max(spectrum_power(:));
        if spectrum_max > 0, spectrum_power = spectrum_power/spectrum_max; end
        imagesc(kx(sx)/k, ky(sy)/k, spectrum_power.');
        axis xy image; colorbar; hold on;
        circle = linspace(0,2*pi,361); plot(cos(circle),sin(circle),'w--','LineWidth',1.2);
        xlabel('k_x/k'); ylabel('k_y/k');
        title('Normalized tangential spectrum |A_x|^2+|A_y|^2; circle K=1');
    end
    figures(end+1) = figure('Name','Continuous far-field angles','Color','w');
    closed_phi = [Ph, Ph(:,1)+2*pi]; closed_th = [Th,Th(:,1)]; closed_u = [U,U(:,1)];
    pcolor(sin(closed_th).*cos(closed_phi), sin(closed_th).*sin(closed_phi), closed_u); shading flat;
    axis equal tight; colorbar; caxis([0,1]);
    xlabel('sin(theta) cos(phi)'); ylabel('sin(theta) sin(phi)');
    title(sprintf('Normalized far-field intensity; theta <= %.3g deg',theta_max_deg));

    figures(end+1) = figure('Name','Azimuth cuts','Color','w');
    plot(theta_deg,U(:,i0),'LineWidth',1.5); hold on;
    plot(theta_deg,U(:,i90),'LineWidth',1.5);
    plot(theta_deg,U(:,i180),'--','LineWidth',1.5);
    xlabel('theta from +Z (deg)'); ylabel('Intensity / global peak');
    legend(sprintf('phi=%.3g deg',phi_deg(i0)),sprintf('phi=%.3g deg',phi_deg(i90)), ...
        sprintf('phi=%.3g deg',phi_deg(i180)),'Location','best');
    grid on; xlim([0,theta_max_deg]); ylim([0,1.05]);

    figures(end+1) = figure('Name','Signed principal-plane cuts','Color','w');
    signed_theta = [-flipud(theta_deg(2:end)); theta_deg];
    ux = [flipud(U(2:end,i180));U(:,i0)];
    uy = [flipud(U(2:end,i270));U(:,i90)];
    if max(ux) > 0, ux = ux/max(ux); end
    if max(uy) > 0, uy = uy/max(uy); end
    plot(signed_theta,ux,'LineWidth',1.5); hold on;
    plot(signed_theta,uy,'LineWidth',1.5);
    plot([-theta_max_deg,theta_max_deg],[.5,.5],'k:');
    xlabel('Signed angle in each plane (deg)'); ylabel('Intensity / each cut peak');
    legend('XZ','YZ','Half-power level','Location','best'); grid on;

    figures(end+1) = figure('Name','3D upper-half-space intensity','Color','w');
    % 只为绘图补齐周期接缝；积分和 CSV 不重复端点。
    closed_phi = [Ph, Ph(:,1)+2*pi]; closed_th = [Th,Th(:,1)]; closed_u = [U,U(:,1)];
    [xx,yy,zz] = sph2cart(closed_phi,pi/2-closed_th,closed_u);
    surf(xx,yy,zz,closed_u,'EdgeColor','none'); axis equal; colorbar; caxis([0,1]);
    xlabel('x'); ylabel('y'); zlabel('z'); view(45,25);
    title('Normalized upper-half-space intensity; radius = intensity');
end

%% ---------------- MAT、CSV 及图片导出 --------------------------------
if save_results
    if isempty(output_root), output_root = fullfile(fileparts(file_name),'farfield_output'); end
    output_dir = fullfile(output_root,datestr(now,'yyyymmdd_HHMMSS'));
    suffix = 0; base_dir = output_dir;
    while exist(output_dir,'dir')
        suffix = suffix+1; output_dir = [base_dir,'_',num2str(suffix)];
    end
    mkdir(output_dir);
    save(fullfile(output_dir,'farfield_result.mat'),'result','-v7');
    Theta_deg = Th(:)*180/pi; Phi_deg = Ph(:)*180/pi;
    Fx_re = real(Ftheta(:)); Fx_im = imag(Ftheta(:));
    Fy_re = real(Fphi(:)); Fy_im = imag(Fphi(:));
    data_table = table(Theta_deg,Phi_deg,Fx_re,Fx_im,Fy_re,Fy_im, ...
        F_squared(:),U(:),dP_dOmega(:),DSC(:), 'VariableNames', ...
        {'theta_deg','phi_deg','Ftheta_real_V','Ftheta_imag_V','Fphi_real_V', ...
        'Fphi_imag_V','F_squared_V2','normalized_intensity','dP_dOmega_W_sr','DSC_m2_sr'});
    writetable(data_table,fullfile(output_dir,'farfield_angles.csv'));
    for fi = 1:numel(figures)
        saveas(figures(fi),fullfile(output_dir,sprintf('figure_%02d.png',fi)));
    end
    fprintf('\n结果已保存到：%s\n',output_dir);
end

%% ---------------- LOCAL FUNCTIONS -----------------------------------
function path = ff_resolve_file(name, base)
    path = '';
    if exist(name,'file') == 2
        [ok, attr] = fileattrib(name);
        if ok, path = attr.Name; else, path = name; end
    elseif exist(fullfile(base,name),'file') == 2
        path = fullfile(base,name);
    end
end

function near = ff_read_plane(path, unit, z_fallback, dx_request, dy_request, abs_tol, rel_tol)
    switch lower(unit)
        case 'm', scale = 1;
        case 'mm', scale = 1e-3;
        case 'um', scale = 1e-6;
        case 'nm', scale = 1e-9;
        otherwise, error('坐标单位只支持 m、mm、um、nm。');
    end
    fid = fopen(path,'r');
    assert(fid > 0,'无法打开文件：%s',path);
    closer = onCleanup(@() fclose(fid)); %#ok<NASGU>
    rows = {}; column_count = 0; line_no = 0;
    while true
        text_line = fgetl(fid);
        if ~ischar(text_line), break; end
        line_no = line_no+1; text_line = strtrim(text_line);
        if isempty(text_line) || any(text_line(1) == '%#'), continue; end
        if isempty(regexp(text_line,'^[+\-.0-9]','once')), continue; end
        values = sscanf(strrep(text_line,',',' '),'%f').';
        assert(numel(values) == 8 || numel(values) == 9, '第 %d 行必须恰为 8 或 9 列实数；请分别导出复场实部、虚部。',line_no);
        assert(all(isfinite(values)),'第 %d 行含 NaN/Inf；不自动置零。',line_no);
        if column_count == 0, column_count = numel(values); end
        assert(numel(values) == column_count,'第 %d 行与前面的列数不一致。',line_no);
        rows{end+1,1} = values; %#ok<AGROW>
    end
    assert(~isempty(rows),'没有读取到有效的近场数据。');
    matrix = vertcat(rows{:});
    xr = matrix(:,1)*scale; yr = matrix(:,2)*scale;
    if column_count == 9
        zr = matrix(:,3)*scale; z = mean(zr); first_field = 4;
        z_tol = max(abs_tol,64*eps(max(abs(zr))));
        assert(max(abs(zr-z)) <= z_tol, 'z 坐标不恒定；单平面角谱法不能用于曲面或多个 z 平面的混合数据。');
    else
        assert(isscalar(z_fallback) && isfinite(z_fallback) && isreal(z_fallback),'8 列文件必须指定有限的 z_plane_m。');
        z = z_fallback; first_field = 3;
    end
    [x,ix,dx_value] = ff_axis(xr,dx_request,abs_tol,rel_tol,'x');
    [y,iy,dy_value] = ff_axis(yr,dy_request,abs_tol,rel_tol,'y');
    nx = numel(x); ny = numel(y);
    linear = sub2ind([nx,ny],ix,iy);
    assert(numel(unique(linear)) == size(matrix,1), '存在重复的 (x,y) 点；即使总行数正确也不能继续。');
    assert(size(matrix,1) == nx*ny,'采样点不能组成完整矩形网格，存在缺失点。');
    fields = complex(matrix(:,first_field:2:first_field+4),matrix(:,first_field+1:2:first_field+5));
    near = struct('x',x,'y',y,'z',z,'dx',dx_value,'dy',dy_value);
    names = {'Ex','Ey','Ez'};
    for j = 1:3
        field = zeros(nx,ny); field(linear) = fields(:,j); near.(names{j}) = field;
    end
end

function [axis_values,index,step] = ff_axis(raw, requested, abs_tol, rel_tol, name)
    tolerance = max(abs_tol,64*eps(max(abs(raw))));
    [sorted, permutation] = sort(raw);
    groups = cumsum([true;diff(sorted)>tolerance]);
    axis_values = accumarray(groups,sorted,[],@mean);
    assert(numel(axis_values) >= 2,'至少需要两个不同的 %s 坐标。',name);
    index = zeros(size(raw)); index(permutation) = groups;
    if ischar(requested) || isstring(requested)
        assert(strcmpi(requested,'auto'),'步长必须为 auto 或正数。');
        step = median(diff(axis_values));
    else
        step = requested;
    end
    assert(isscalar(step) && isreal(step) && isfinite(step) && step > 0,'步长必须为有限正数。');
    expected = axis_values(1)+(0:numel(axis_values)-1).'*step;
    allowed = max(5*tolerance,rel_tol*step);
    assert(max(abs(axis_values-expected)) <= allowed, '%s 坐标不是完整等间距网格，或手动步长错误；不压缩缺失列/行。',name);
    assert(max(abs(raw-axis_values(index))) <= allowed,'%s 坐标误差超过容差。',name);
    axis_values = expected;
end

function w = ff_tukey(n, alpha)
    % 自包含实现，不依赖 Signal Processing Toolbox。
    w = ones(n,1);
    if alpha == 0 || n == 1, return; end
    t = (0:n-1).'/(n-1);
    left = t < alpha/2; right = t > 1-alpha/2;
    w(left) = .5*(1+cos(pi*(2*t(left)/alpha-1)));
    w(right) = .5*(1+cos(pi*(2*t(right)/alpha-2/alpha+1)));
end

function spectrum = ff_direct(field, x, y, qx, qy, dx, dy)
    output_size = size(qx); qx = qx(:); qy = qy(:);
    spectrum = complex(zeros(size(qx)));
    for first = 1:64:numel(qx)
        take = first:min(first+63,numel(qx));
        px = exp(-1i*x(:)*qx(take).');
        py = exp(-1i*y(:)*qy(take).');
        spectrum(take) = sum((px.'*field).*py.',2)*dx*dy;
    end
    spectrum = reshape(spectrum,output_size);
end

function index = ff_phi_index(phi, target)
    [~,index] = min(abs(atan2(sin(phi-target),cos(phi-target))));
end

function beam = ff_beamwidth(theta, positive, negative)
    angle = [-flipud(theta(2:end));theta];
    power = [flipud(negative(2:end));positive];
    [peak,ip] = max(power);
    beam = struct('peak_signed_deg',angle(ip),'cut_peak',peak, ...
        'half_level',peak/2,'left_deg',NaN,'right_deg',NaN,'HPBW_deg',NaN);
    if peak <= 0, return; end
    il = find(power(1:ip) <= peak/2,1,'last');
    ir = find(power(ip:end) <= peak/2,1,'first');
    if ~isempty(il) && il < ip
        beam.left_deg = ff_crossing(angle(il),angle(il+1),power(il),power(il+1),peak/2);
    end
    if ~isempty(ir)
        ir = ip+ir-1;
        if ir > ip
            beam.right_deg = ff_crossing(angle(ir-1),angle(ir),power(ir-1),power(ir),peak/2);
        end
    end
    if isfinite(beam.left_deg) && isfinite(beam.right_deg)
        beam.HPBW_deg = beam.right_deg-beam.left_deg;
    end
end

function angle = ff_crossing(x1,x2,y1,y2,level)
    if y1 == y2, angle = (x1+x2)/2;
    else, angle = x1+(level-y1)*(x2-x1)/(y2-y1);
    end
end

function ff_print_beamwidth(name, beam)
    if isfinite(beam.HPBW_deg)
        fprintf('%s 切面最强瓣 HPBW=%.6g deg，左右交点 %.6g / %.6g deg。\n', name,beam.HPBW_deg,beam.left_deg,beam.right_deg);
    else
        fprintf('%s 切面未找到完整的两侧半功率交点，HPBW=NaN。\n',name);
    end
end
