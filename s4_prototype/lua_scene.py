"""Generate numerical-only Lua input for the pinned S4 command-line runtime."""
import math


OUTPUT = r'''
local function row(tag, ...)
    io.write('WS4\t', tag)
    for _, value in ipairs({...}) do io.write('\t', string.format('%.17g', value)) end
    io.write('\n')
end
local orders = S:GetGList()
local f, b = S:GetAmplitudes('top', 0)
local power = S:GetPowerFluxByOrder('top', 0)
row('COUNT', #orders)
row('TOP', S:GetPowerFlux('top', 0))
row('BOTTOM', S:GetPowerFlux('bottom', 0))
for i, g in ipairs(orders) do
    local n = #orders
    row('MODE', i, g[1], g[2],
        f[i][1], f[i][2], f[n+i][1], f[n+i][2],
        b[i][1], b[i][2], b[n+i][1], b[n+i][2],
        power[i][1], power[i][2], power[i][3], power[i][4])
end
'''


def number(value):
    value = float(value)
    if not math.isfinite(value):
        raise ValueError('Lua inputs must be finite')
    return format(value, '.17g')


def material(name, epsilon):
    return "S:AddMaterial('{}', {{{},{}}})".format(name, number(epsilon.real), number(epsilon.imag))


def script(scene, particle=False, polarization='p', one_dimensional=False,
           factorization='default', resolution=8, near_points=(), epsilon_method='analytic',
           stripe_rule=False, finite_patch=None):
    scene.validate()
    if finite_patch is not None:
        finite_patch.validate(scene)
        if stripe_rule or one_dimensional:
            raise ValueError('Finite line ends require full 2D material convolution, not the invariant stripe rule')
    if polarization not in ('s', 'p') or factorization not in ('normal', 'default', 'vector', 'jones'):
        raise ValueError('Unsupported S4 illumination/factorization')
    if epsilon_method not in ('analytic', 'fft'):
        raise ValueError('Unsupported S4 material Fourier method')
    if scene.harmonics_x is not None and one_dimensional:
        raise ValueError('Rectangular truncation requires a 2D supercell')
    if stripe_rule and epsilon_method != 'analytic':
        raise ValueError('The stripe adapter uses analytic material coefficients')
    if particle and (one_dimensional or min(scene.lx_um, scene.ly_um) <= scene.diameter_um):
        raise ValueError('A sphere requires a 2D supercell without touching copies')
    if isinstance(resolution, bool) or not isinstance(resolution, int) or resolution < 2:
        raise ValueError('Normal-vector resolution must be an integer >= 2')
    commands = ['local S = S4.NewSimulation()']
    if one_dimensional:
        commands.append('S:SetLattice('+number(scene.pitch_um)+')')
    else:
        commands.append('S:SetLattice({'+number(scene.lx_um)+',0},{0,'+number(scene.ly_um)+'})')
    commands += ['S:SetNumG({})'.format(scene.modes), 'S:SetResolution({})'.format(resolution)]
    if scene.harmonics_x is not None:
        commands.append('S:SetRectangularBasis({},{})'.format(scene.harmonics_x, scene.harmonics_y))
    if stripe_rule:
        commands.append('S:UseAnalyticStripeRule()')
    if epsilon_method == 'fft':
        commands.append('S:UseDiscretizedEpsilon()')
    if factorization != 'default':
        commands.append('S:UsePolarizationDecomposition()')
        if factorization == 'normal':
            commands.append('S:UseNormalVectorBasis()')
        elif factorization == 'jones':
            commands.append('S:UseJonesVectorBasis()')
    commands += [material('Air', 1+0j), material('Wafer', scene.eps_substrate),
                 material('Particle', scene.eps_particle), "S:AddLayer('top',0,'Air')"]
    volume = 0.
    if particle:
        dz, radius = scene.diameter_um/scene.slices, scene.diameter_um/2
        for i in range(scene.slices):
            r = math.sqrt(radius**2-((i+.5)*dz-radius)**2)
            name = 'particle_'+str(i)
            commands.append("S:AddLayer('{}',{},'Air')".format(name, number(dz)))
            if scene.eps_particle != 1.:
                commands.append("S:SetLayerPatternCircle('{}','Particle',{{0,0}},{})".format(name, number(r)))
            volume += math.pi*r*r*dz
        if scene.gap_um:
            commands.append("S:AddLayer('gap',{},'Air')".format(number(scene.gap_um)))
    else:
        commands.append("S:AddLayer('cap',{},'Air')".format(number(scene.cap_top_um)))
    base = 'Air' if scene.kind == 'grating' and finite_patch is None else 'Wafer'
    commands.append("S:AddLayer('grating',{},'{}')".format(number(scene.depth_um), base))
    if finite_patch is not None:
        width = scene.pitch_um-scene.ridge_width_um
        for center in finite_patch.centers(scene):
            commands.append("S:SetLayerPatternRectangle('grating','Air',{{0,{}}},0,{{{},{}}})".format(
                number(center), number(finite_patch.length_um/2), number(width/2)))
    elif scene.kind == 'grating':
        for j in range(1 if one_dimensional else scene.periods_y):
            y = (j*scene.pitch_um+scene.ly_um/2) % scene.ly_um-scene.ly_um/2
            if one_dimensional:
                commands.append("S:SetLayerPatternRectangle('grating','Wafer',{{0,0}},0,{{{},{}}})".format(
                    number(scene.ridge_width_um/2), number(scene.lx_um/2)))
            else:
                commands.append("S:SetLayerPatternRectangle('grating','Wafer',{{0,{}}},0,{{{},{}}})".format(
                    number(y), number(scene.lx_um/2), number(scene.ridge_width_um/2)))
    azimuth = 0 if one_dimensional else 90
    commands += ["S:AddLayer('bottom',0,'Wafer')", 'S:SetFrequency({})'.format(number(1/scene.wavelength_um)),
        'S:SetExcitationPlanewave({'+number(scene.incidence_deg)+','+str(azimuth)+'},{'+
        ('1,0},{0,0})' if polarization == 's' else '0,0},{1,0})'), OUTPUT]
    for point in near_points:
        x, y, z = map(float, point)
        sx, sy = (y, -x) if one_dimensional else (x, y)
        commands.append("row('FIELD', {}, {}, {}, S:GetFields({{{},{},{}}}))".format(
            number(x), number(y), number(z), number(sx), number(sy), number(scene.cap_top_um-z)))
    commands.append("row('DONE',1)")
    return '\n'.join(commands)+'\n', volume
