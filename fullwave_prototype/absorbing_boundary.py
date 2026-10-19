"""Radial complex stretching that preserves the planar substrate interface."""


def stretch_factors(radius, p):
    """Tangential/radial Jacobian eigenvalues, for radius >= physical radius.

    r_tilde = r + i * strength * T * ((r-R)/T)**power.
    Equal outer attenuation distance for every power; powers >= 2 start
    with the identity Jacobian at the physical/PML interface.
    Works with both NumPy arrays and NGSolve coefficient functions.
    """
    depth = (radius-p.domain_radius_um)/p.pml_thickness_um
    shift = 1j*p.pml_strength*p.pml_thickness_um*depth**p.pml_power
    tangential = 1+shift/radius
    radial = 1+1j*p.pml_strength*p.pml_power*depth**(p.pml_power-1)
    return tangential, radial


def set_boundary(mesh, p):
    import ngsolve as ng
    if p.pml_power == 1:
        transformation = ng.pml.Radial(rad=p.domain_radius_um,
            alpha=1j*p.pml_strength, origin=(0, 0, 0))
    else:
        position = (ng.x, ng.y, ng.z)
        radius = ng.sqrt(sum(x*x for x in position))
        tangential, radial = stretch_factors(radius, p)
        transformed = ng.CF(tuple(tangential*x for x in position))
        jacobian = ng.CF(tuple(
            tangential*(i == j)+(radial-tangential)*position[i]*position[j]/radius**2
            for i in range(3) for j in range(3)), dims=(3, 3))
        transformation = ng.pml.Custom(transformed, jacobian)
    mesh.SetPML(transformation, 'pml_.*')
