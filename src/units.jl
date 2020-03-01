module Units

import Unitful
using Unitful: @unit, @derived_dimension,𝐋

# Re-export units
import Unitful: nm, m, cm, s, kg, g, eV, keV, MeV, c, u, b, C, J, A, H, Hz, rad, T, V, Ω, W, K, Pa, angstrom, Å, cal, atm, Torr, mol, lm
import Unitful: k, h, c0

# Define New Units
@unit Gauss "Gauss"    Gauss       (1//10_000)T   true
@unit kB     "kB"    Boltzmann          k         false
@unit hc     "hc"   PlanckLight        h*c0       false

@derived_dimension NumberDensity 𝐋^-3

export nm, m, cm, s, kg, g, eV, keV, MeV, c, u, b, C, J, A, Hz, H, rad, T, V, Ω, W, K, Pa, angstrom, Å, cal, atm, Torr, mol, lm
export Gauss, kB, hc, NumberDensity

#Register Units
const localunits = Unitful.basefactors
function __init__()
    merge!(Unitful.basefactors, localunits)
    Unitful.register(Units)
end

end
