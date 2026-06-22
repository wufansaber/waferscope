"""Reference optical constants used by the SCATMECH user interface."""

# Reference values near 355 nm. Every n and k field remains user-editable.
WAFER_MATERIALS = {
    "Silicon (Si)": 5.57 + 0.387j,
    "Fused silica (SiO2)": 1.476 + 0.0j,
    "Sapphire (Al2O3)": 1.79 + 0.0j,
    "Silicon carbide (SiC)": 2.75 + 0.02j,
    "Gallium arsenide (GaAs)": 3.65 + 1.80j,
    "Custom": 1.50 + 0.0j,
}

PARTICLE_MATERIALS = {
    "PSL": 1.59 + 0.0j,
    "Silica (SiO2)": 1.476 + 0.0j,
    "Silicon (Si)": 5.57 + 0.387j,
    "Alumina (Al2O3)": 1.79 + 0.0j,
    "Aluminum (Al)": 0.49 + 4.84j,
    "Gold (Au)": 1.48 + 1.86j,
    "Custom": 1.50 + 0.0j,
}
