# Source CODATA2018 https://physics.nist.gov/cuu/Constants/
module Constants

using ..Units: m, s, c, u, kg, eV, keV, MeV, lm, W, K, m, atm, Torr, Hz, Ω

# import constants defined in Unitful
using Unitful
import Unitful: c0, k, μ0, q, ϵ0, ε0, Z0, me, mp, mn, G, gn, h, ħ, μB, Na, R, σ, R∞, μB

# Masses and Charges
const AtomicMassUnit = uconvert(kg, 1u)
const mu = AtomicMassUnit
export AtomicMassUnit

const ElectronMass = me
const ElectronAtomicMass = 5.48579909065e-4u
const Ae = ElectronAtomicMass
const ElectronRestMass = 0.51099895000MeV/c^2
const me0 = ElectronRestMass
const ElectronChargeNumber = -1
const Ze = ElectronChargeNumber
export ElectronMass, ElectronAtomicMass, ElectronRestMass, ElectronChargeNumber

const ProtonMass = mp
const ProtonAtomicMass = 1.007276466621u
const Ap = ProtonAtomicMass
const ProtonRestMass = 938.27208816MeV/c^2
const mp0 = ProtonRestMass
const ProtonChargeNumber = 1
const Zp = ProtonChargeNumber
export ProtonMass, ProtonAtomicMass, ProtonRestMass, ProtonChargeNumber

const NeutronMass = mn
const NeutronAtomicMass = 1.00866491595u
const An = NeutronAtomicMass
const NeutronRestMass = 939.56542052MeV/c^2
const mn0 = NeutronRestMass
const NeutronChargeNumber = 0
const Zn = NeutronChargeNumber
export NeutronMass, NeutronAtomicMass, NeutronRestMass, NeutronChargeNumber

const DeuteronMass = 3.3435837724e-27kg
const md = DeuteronMass
const DeuteronAtomicMass = 2.013553212745u
const Ad = DeuteronAtomicMass
const DeuteronRestMass = 1875.61294257MeV/c^2
const md0 = DeuteronRestMass
const DeuteronChargeNumber = 1
const Zd = DeuteronChargeNumber
export DeuteronMass, DeuteronAtomicMass, DeuteronRestMass, DeuteronChargeNumber

const TritonMass = 5.0073567446e-27kg
const mt = TritonMass
const TritonAtomicMass = 3.01550071621u
const At = TritonAtomicMass
const TritonRestMass = 2808.92113298MeV/c^2
const mt0 = TritonRestMass
const TritonChargeNumber = 1
const Zt = TritonChargeNumber
export TritonMass, TritonAtomicMass, TritonRestMass, TritonChargeNumber

const HelionMass = 5.0064127796e-27kg
const mh = HelionMass
const HelionAtomicMass = 3.014932247175u
const Ah = HelionAtomicMass
const HelionRestMass = 2808.39160743MeV/c^2
const mh0 = HelionRestMass
const HelionChargeNumber = 2
const Zh = HelionChargeNumber
export HelionMass, HelionAtomicMass, HelionRestMass, HelionChargeNumber

const AlphaMass = 6.6446573357e-27kg
const mα = AlphaMass
const ma = mα
const AlphaAtomicMass = 4.001506179127u
const Aα = AlphaAtomicMass
const Aa = Aα
const AlphaRestMass = 3727.3794066MeV/c^2
const mα0 = AlphaRestMass
const ma0 = mα0
const AlphaChargeNumber = 2
const Zα = AlphaChargeNumber
const Za = Zα
export AlphaMass, AlphaAtomicMass, AlphaRestMass, AlphaChargeNumber

# Impurity Atomic Masses
const LithiumAtomicMass = 6.941u
const Lithium6AtomicMass = 6.938u
const Lithium7AtomicMass = 6.997u
const LithiumMass = uconvert(kg, LithiumAtomicMass)
const Lithium6Mass = uconvert(kg, Lithium6AtomicMass)
const Lithium7Mass = uconvert(kg, Lithium7AtomicMass)
const LithiumChargeNumber = 1
export LithiumAtomicMass, Lithium6AtomicMass, Lithium7AtomicMass
export LithiumMass, Lithium6Mass, Lithium7Mass
export LithiumChargeNumber

const BerylliumAtomicMass = 9.0121831u
const BerylliumMass = uconvert(kg, BerylliumAtomicMass)
const BerylliumChargeNumber = 2
export BerylliumAtomicMass, BerylliumMass, BerylliumChargeNumber

const BoronAtomicMass = 10.811u
const Boron10AtomicMass = 10.806u
const Boron11AtomicMass = 10.821u
const BoronMass = uconvert(kg, BoronAtomicMass)
const Boron10Mass = uconvert(kg, Boron10AtomicMass)
const Boron11Mass = uconvert(kg, Boron11AtomicMass)
const BoronChargeNumber = 5
export BoronAtomicMass, Boron10AtomicMass, Boron11AtomicMass
export BoronMass, Boron10Mass, Boron11Mass
export BoronChargeNumber

const CarbonAtomicMass = 12.011u
const Carbon12AtomicMass = 12.0096u
const Carbon13AtomicMass = 12.0116u
const CarbonMass = uconvert(kg, CarbonAtomicMass)
const Carbon12Mass = uconvert(kg, Carbon12AtomicMass)
const Carbon13Mass = uconvert(kg, Carbon13AtomicMass)
const CarbonChargeNumber = 6
export CarbonAtomicMass, Carbon12AtomicMass, Carbon13AtomicMass
export CarbonMass, Carbon12Mass, Carbon13Mass
export CarbonChargeNumber

const NitrogenAtomicMass = 14.007u
const Nitrogen14AtomicMass = 14.00643u
const Nitrogen15AtomicMass = 14.00728u
const NitrogenMass = uconvert(kg, NitrogenAtomicMass)
const Nitrogen14Mass = uconvert(kg, Nitrogen14AtomicMass)
const Nitrogen15Mass = uconvert(kg, Nitrogen15AtomicMass)
const NitrogenChargeNumber = 7
export NitrogenAtomicMass, Nitrogen14AtomicMass, Nitrogen15AtomicMass
export NitrogenMass, Nitrogen14Mass, Nitrogen15Mass
export NitrogenChargeNumber

const OxygenAtomicMass = 15.999u
const OxygenMass = uconvert(kg, OxygenAtomicMass)
const OxygenChargeNumber = 8
export OxygenAtomicMass, OxygenMass, OxygenChargeNumber


# Physical Constants
const ProtonElectronMassRatio = 1836.15267343
const ElectronProtonMassRatio = 5.44617021487e-4
export ProtonElectronMassRatio, ElectronProtonMassRatio

const GravitationalConstant = G
const GravitationalAcceleration = gn
export GravitationalConstant, GravitationalAcceleration

const PlanckConstant = h
const ReducedPlanckConstant = ħ
export PlanckConstant, ReducedPlanckConstant

const ElementaryCharge = q
const ElectronChargeMassRatio = q/me
export ElementaryCharge, ElectronChargeMassRatio

const BoltzmannConstant = k
const StephanBoltzmannConstant = σ
export BoltzmannConstant, StephanBoltzmannConstant

const SpeedOfLightInVacuum = c0
const PermittivityOfFreeSpace = ε0
const PermeabilityOfFreeSpace = μ0
const ImpedanceOfVacuum = Z0
export SpeedOfLightInVacuum, PermittivityOfFreeSpace, PermeabilityOfFreeSpace, ImpedanceOfVacuum

const RydbergConstant = R∞
const Ry = 13.605693122994eV
const RydbergEnergy = Ry
const RydbergFrequency = c0*R∞
const cR∞ = RydbergFrequency
export RydbergConstant, RydbergEnergy, RydbergFrequency

const BohrRadius = ϵ0*h^2/(pi*me*q^2)
const a0 = BohrRadius
const AtomicCrossSection = pi*a0^2
const r_e = (q^2)/(4pi*ϵ0*me*c0^2)
const ClassicalElectronRadius = r_e
const ThomsonCrossSection = (8pi/3)*r_e^2
const σTh = ThomsonCrossSection
const BohrMagneton = μB
export BohrRadius, AtomicCrossSection, ClassicalElectronRadius, ThomsonCrossSection, BohrMagneton

const ComptonWavelength = h/(me*c0)
const λC  = ComptonWavelength
const ReducedComptonWavelength = ħ/(me*c0)
const ƛC = ReducedComptonWavelength
export ComptonWavelength, ReducedComptonWavelength

const FineStructureConstant = 7.2973525693e-3
const α = FineStructureConstant
const InverseFineStructureConstant = 137.035999084
const inv_α = InverseFineStructureConstant
export FineStructureConstant, InverseFineStructureConstant

const FirstRadiationConstant = 2π*h*c0
const c1 = FirstRadiationConstant
const SecondRadiationConstant = h*c0/k
const c2 = SecondRadiationConstant
export FirstRadiationConstant, SecondRadiationConstant

const λ0 = h*c0/q
const ElectronVoltWavelength = λ0
const ν0 = q/h
const ElectronVoltFrequency = ν0
const k0 = q/(h*c0)
const ElectronVoltWaveNumber = k0
const ElectronVoltTemperature = k/q
const TemperatureElectronVolt = q/k
export ElectronVoltWavelength, ElectronVoltFrequency, ElectronVoltWaveNumber, ElectronVoltTemperature, TemperatureElectronVolt

const AvogadroNumber = Na
const FaradayConstant = Na*q
const F = FaradayConstant
const MolarGasConstant = R
const LoschmidtsNumber = 2.6868e25m^-3
const n0 = LoschmidtsNumber
const StandardTemperature = 273.15K
const T0 = StandardTemperature
const StandardPressure = 1atm
const p0 = 1atm
const mmHg = 1Torr
const StandardVolume = R*T0/p0
const V0 = StandardVolume
const M_air = 2.8971e-2kg
const MolarWeightOfAir = M_air
export AvogadroNumber, FaradayConstant, MolarGasConstant, LoschmidtsNumber
export StandardTemperature, StandardPressure, StandardVolume, MolarWeightOfAir

# Miscellaneous Constants
const Cs133HyperfineTransitionFrequency = 9192631770Hz
const ΔνCs133 = Cs133HyperfineTransitionFrequency
const LuminousEfficacy = 683lm/W
const Kcd = LuminousEfficacy
const JosephsonConstant = 2q/h
const KJ = JosephsonConstant
const VonKlitzingConstant = uconvert(Ω, h/q^2)
const RK = VonKlitzingConstant
export Cs133HyperfineTransitionFrequency, LuminousEfficacy, JosephsonConstant, VonKlitzingConstant

end
