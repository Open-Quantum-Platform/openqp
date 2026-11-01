"""Physical constants"""
ANGSTROM_TO_BOHR = 0.52917721090299996

"""
1 Hartree = 4.3597482e-18 Nm
1 Bohr = 5.2917721067121e-11 m
1 Hartree/Bohr ** 2 = 1556.8943461281847 N/m
1 N/m = 1 kg/s ** 2
1 Hartree/Bohr ** 2 = 1556.8943461281847 kg/s ** 2
1 g/mol = 1.6605390671738466e-27 kg
1 Hartree/Bohr ** 2/(g/mol) = 9.37583690083209e+29 1/s ** 2
1 [Hartree/Bohr ** 2/(g/mol)]/2π = 154107989120439.12 1/s
1 c = 29979245800 cm / s
1 [Hartree/Bohr ** 2/(g/mol)]/2πc = 5140.489195376594 cm -1
"""
SPEED_OF_LIGHT = 2.99792458e10  # cm / s
ATMOS = 101.325  # kPa
BOHR = 0.52917721090299996e-10  # m
FREQ_TO_INV_CM = 5140.489195376594  # cm -1
AMU_to_KG = 1.66053886E-27  # kg
J_TO_AU = 1 / (4.184 * 627.509541 * 1000.0)  # Hartree
J_TO_INV_CM = 5.03412e22  # cm-1
GAS_CONSTANT = 8.3144621  # J / K / mol
PLANCK_CONSTANT = 6.62606957e-34  # J * s
BOLTZMANN_CONSTANT = 1.3806488e-23  # J / K
AVOGADRO_CONSTANT = 6.0221415e23  # 1 / mol

HA_TO_WAVENUM = 219474.6313708  # cm-1, source/physical_constants.F90
FINE_STRUCTURE = 7.2973525693e-3  # source/physical_constants.F90

#: Multiplier for the raw spin-orbits of ``OQP::soc_hsoc_re/_im``, which the
#: native kernel stores in atomic units WITHOUT the alpha^2/2 prefactor
#: (source/modules/soc_mrsf.F90).  Matching the value the SOC log table prints.
SOC_AU_TO_WAVENUM = 0.5 * FINE_STRUCTURE ** 2 * HA_TO_WAVENUM  # ~5.8438 cm-1/a.u.
