# Re-export units
@reexport using Unitful: m, cm, s, kg, g, eV, keV, c, u, b, C, J, A, Hz, rad, T, V, Ω, W, K, Pa, angstrom, Å, cal, atm, Torr, mol, lm
const MeV = u"MeV"
export MeV

# Define Gauss Units
@unit Gs "Gs" Gauss (1//10_000)T true
export Gs

import Unitful: uconvert, TemperatureUnits, EnergyUnits, MassUnits, Temperature, Energy, Mass, k, c0

# Temperature <--> Energy Conversion 1 J = k_b * 1 K
uconvert(a::T, x::S) where {T<:EnergyUnits, S<:Temperature} = uconvert(a, uconvert(K,x)*k)
uconvert(a::T, x::S) where {T<:TemperatureUnits, S<:Energy} = uconvert(a, uconvert(J,x)/k)

# import constants defined in Unitful
import Unitful: c0, μ0, k, q, ϵ0, ε0, Z0, me, mp, mn, G, gn, h, ħ, μB, Na, R, σ, R∞, μB

# MFE Formulary http://library.psfc.mit.edu/catalog/online_pubs/MFE_formulary_2014.pdf
# Values from CODATA

# Masses
const AtomicMassUnit = uconvert(kg, 1u)
const mu = AtomicMassUnit
export AtomicMassUnit

const ElectronMass = me
const ElectronAtomicMass = 5.48579909065e-4u
const Ae = ElectronAtomicMass
const ElectronRestMass = 0.51099895000MeV/c^2
const me0 = ElectronRestMass
export ElectronMass, ElectronAtomicMass, ElectronRestMass

const ProtonMass = mp
const ProtonAtomicMass = 1.007276466621u
const Ap = ProtonAtomicMass
const ProtonRestMass = 938.27208816MeV/c^2
const mp0 = ProtonRestMass
export ProtonMass, ProtonAtomicMass, ProtonRestMass

const NeutronMass = mn
const NeutronAtomicMass = 1.00866491595u
const An = NeutronAtomicMass
const NeutronRestMass = 939.56542052MeV/c^2
const mn0 = NeutronRestMass
export NeutronMass, NeutronAtomicMass, NeutronRestMass

const DeuteronMass = 3.3435837724e-27kg
const md = DeuteronMass
const DeuteronAtomicMass = 2.013553212745u
const Ad = DeuteronAtomicMass
const DeuteronRestMass = 1875.61294257MeV/c^2
const md0 = DeuteronRestMass
export DeuteronMass, DeuteronAtomicMass, DeuteronRestMass

const TritonMass = 5.0073567446e-27kg
const mt = TritonMass
const TritonAtomicMass = 3.01550071621u
const At = TritonAtomicMass
const TritonRestMass = 2808.92113298MeV/c^2
const mt0 = TritonRestMass
export TritonMass, TritonAtomicMass, TritonRestMass

const HelionMass = 5.0064127796e-27kg
const mh = HelionMass
const HelionAtomicMass = 3.014932247175u
const Ah = HelionAtomicMass
const HelionRestMass = 2808.39160743MeV/c^2
const mh0 = HelionRestMass
export HelionMass, HelionAtomicMass, HelionRestMass

const AlphaMass = 6.6446573357e-27kg
const mα = AlphaMass
const ma = mα
const AlphaAtomicMass = 4.001506179127u
const Aα = AlphaAtomicMass
const Aa = Aα
const AlphaRestMass = 3727.3794066MeV/c^2
const mα0 = AlphaRestMass
const ma0 = mα0
export AlphaMass, AlphaAtomicMass, AlphaRestMass

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
const kB = k
const StephanBoltzmannConstant = σ
export BolzmannConstant, StephanBoltzmannConstant

const SpeedOfLightInVacuum = c0
const PermittivityOfFreeSpace = ε0
const PermeabilityOfFreeSpace = μ0
export SpeedofLightInVacuum, PermittivityOfFreeSpace, PermeabilityOfFreeSpace

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

# CODATA Wallet https://physics.nist.gov/cuu/pdf/wallet_2018.pdf
const Cs133HyperfineTransitionFrequency = 9192631770Hz
const ΔνCs = Cs133HyperfineTransitionFrequency
const LuminousEfficacy = 683lm/W
const Kcd = LuminousEfficacy
const JosephsonConstant = 2q/h
const KJ = JosephsonConstant
const VonKlitzingConstant = uconvert(Ω, h/q^2)
const RK = VonKlitzingConstant
export Cs133HyperfineTransitionFrequency, LuminousEfficacy, JosephsonConstant, VonKlitzingConstant
