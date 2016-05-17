from numpy import sqrt

# Evolutionary Factor
def Evolution_factors(Cosmology,zred):
    Ez   = sqrt(Cosmology.Omega_M*(1.0+zred)**3 + Cosmology.Omega_DE)
    EzRo = sqrt(Cosmology.Omega_M*(1.0+0.23)**3 + Cosmology.Omega_DE)
    return (Ez,EzRo)

# Evolutionary Factor Lx
def Evolution_factors_Lx(Cosmology,zred):
    Ez   = sqrt(Cosmology.Omega_M*(1.0+zred)**3 + Cosmology.Omega_DE)
    EzRo = sqrt(Cosmology.Omega_M*(1.0+0.23)**3 + Cosmology.Omega_DE)
    return (Ez,EzRo)

# Evolutionary Factor Lx
def Evolution_factors_Tx(Cosmology,zred):
    Ez   = sqrt(Cosmology.Omega_M*(1.0+zred)**3 + Cosmology.Omega_DE)
    EzRo = sqrt(Cosmology.Omega_M*(1.0+0.23)**3 + Cosmology.Omega_DE)
    return (Ez,EzRo)
