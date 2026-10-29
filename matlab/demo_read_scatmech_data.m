% Select a SCATMECH NPZ export; no Python or extra toolbox is required.
[result, source, info] = read_scatmech_data();
fprintf('Archive mode: %s\n', info.mode);
if isequal(info.manifest_verified, true)
    fprintf('SHA256 verified.\n');
else
    fprintf('SHA256 not verified (manifest missing or verification disabled).\n');
end

native = result;
if isfield(result, 'plot_kind') && strcmp(result.plot_kind, 'objective_collection')
    channels = result.channel_pupils;
    sharedMaximum = 0;
    for ch = 1:numel(channels)
        p = channels{ch};
        values = p.sampled_dsc_um2_sr(p.aperture_mask);
        values = values(isfinite(values));
        if ~isempty(values)
            sharedMaximum = max(sharedMaximum, max(values));
        end
    end
    if sharedMaximum <= 0
        sharedMaximum = 1;
    end
    figure('Name', 'SCATMECH objective channels');
    for ch = 1:numel(channels)
        p = channels{ch};
        dsc = p.sampled_dsc_um2_sr;
        dsc(~p.aperture_mask) = NaN;
        subplot(1, numel(channels), ch);
        surf(p.x, p.y, dsc, 'EdgeColor', 'none');
        view(2); axis image; colorbar; caxis([0 sharedMaximum]);
        xlabel('Pupil x (mm)'); ylabel('Pupil y (mm)');
        title(sprintf('CH%d DSC (um^2/sr)', ch));
    end
    native = source;
elseif isfield(result, 'plot_kind') && strcmp(result.plot_kind, 'particle_size_sweep')
    figure('Name', 'SCATMECH size sweep');
    plot(result.diameters_nm, result.channel_collected_dsc_um2, '-o');
    xlabel('Particle diameter (nm)'); ylabel('Collected DSC (um^2)');
    legend('CH1', 'CH2', 'CH3'); grid on;
    native = [];
elseif isfield(result, 'plot_kind') && strcmp(result.plot_kind, 'coherent_imaging')
    channels = result.channels;
    figure('Name', 'SCATMECH image plane');
    for ch = 1:numel(channels)
        p = channels{ch};
        subplot(1, numel(channels), ch);
        imagesc(p.image_axis_um, p.image_axis_um, p.irradiance_w_um2);
        axis image; axis xy; colorbar;
        xlabel('Image x (um)'); ylabel('Image y (um)');
        title(sprintf('CH%d (W/um^2)', ch));
    end
    native = [];
end

if isstruct(native) && isfield(native, 'theta') && isfield(native, 'intensity')
    figure('Name', 'SCATMECH native far field');
    surf(native.x, native.y, native.intensity, 'EdgeColor', 'none');
    view(2); axis image; colorbar;
    xlabel('k_x/k'); ylabel('k_y/k'); title('Absolute DSC (um^2/sr)');
    thetaAxis = native.theta(:, 1);
    phiAxis = native.phi(1, :);
    integrand = bsxfun(@times, native.intensity, sin(thetaAxis));
    sampledDSC = trapz(thetaAxis, trapz(phiAxis, integrand, 2), 1);
    fprintf('DSC integrated over saved angular range: %.9g um^2\n', sampledDSC);
end
