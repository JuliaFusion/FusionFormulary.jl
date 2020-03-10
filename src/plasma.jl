module PlasmaParameters

using Unitful
using ..Units
using ..Constants: q, me, ϵ0

function electron_cyclotron_frequency(B::Unitful.BField)
    return uconvert(rad/s, q*B/me)
end

function electron_cyclotron_frequency(B::U) where U
    B_units = isa(B, Unitful.Quantity) ? B : B*T
    return ustrip(electron_cyclotron_frequency(B_units))
end

function ion_cyclotron_frequency(B::Unitful.BField, Zi::Int, mi::Unitful.Mass)
    return uconvert(rad/s, Zi*q*B/mi)
end

function ion_cyclotron_frequency(B, Zi, mi)
    B_units = isa(B, Unitful.Quantity) ? B : B*T
    mi_units = isa(mi, Unitful.Quantity) ? mi : mi*u
    return ustrip(ion_cyclotron_frequency(B_units, Zi, mi_units))
end

export electron_cyclotron_frequency, ion_cyclotron_frequency

function electron_plasma_frequency(ne::NumberDensity)
    return uconvert(rad/s,sqrt((ne*q^2)/(me*ϵ0)))
end

function electron_plasma_frequency(ne)
    ne_units = isa(ne, Unitful.Quantity) ? ne : ne*cm^-3
    return ustrip(electron_plasma_frequency(ne_units))
end

function ion_plasma_frequency(ni::NumberDensity, Zi::Int, mi::Unitful.Mass)
    return uconvert(rad/s,sqrt((ni*(Zi*q)^2)/(mi*ϵ0)))
end

function ion_plasma_frequency(ni, Zi, mi)
    ni_units = isa(ni, Unitful.Quantity) ? ni : ni*cm^-3
    mi_units = isa(mi, Unitful.Quantity) ? mi : mi*u
    return ustrip(ion_plasma_frequency(ni_units, Zi, mi_units))
end

export electron_plasma_frequency, ion_plasma_frequency

function electron_thermal_speed(Te::Unitful.Energy)
    return uconvert(m/s,sqrt(Te/me))
end

function electron_thermal_speed(Te)
    Te_units = isa(Te, Unitful.Quantity) ? Te : Te*keV
    return ustrip(electron_thermal_speed(Te_units))
end

function ion_thermal_speed(Ti::Unitful.Energy, mi::Unitful.Mass)
    return uconvert(m/s,sqrt(Ti/mi))
end

function ion_thermal_speed(Ti, mi)
    Ti_units = isa(Ti, Unitful.Quantity) ? Ti : Ti*keV
    mi_units = isa(mi, Unitful.Quantity) ? mi : mi*u
    return ustrip(ion_thermal_speed(Ti_units, mi_units))
end

export electron_thermal_speed, ion_thermal_speed

function electron_gyro_radius(Te::Unitful.Energy, B::Unitful.BField)
    return uconvert(m, sqrt(2*me*Te)/(q*B))
end

function electron_gyro_radius(Te, B)
    Te_units = isa(Te, Unitful.Quantity) ? Te : Te*keV
    B_units = isa(B, Unitful.Quantity) ? B : B*T
    return ustrip(electron_gyro_radius(Te_units,B_units))
end

function ion_gyro_radius(Ti::Unitful.Energy, B::Unitful.BField, Zi::Int, mi::Unitful.Mass)
    return uconvert(m, sqrt(2*mi*Ti)/(Zi*q*B))
end

function ion_gyro_radius(Ti, B, Zi, mi)
    Ti_units = isa(Ti, Unitful.Quantity) ? Ti : Ti*keV
    B_units = isa(B, Unitful.Quantity) ? B : B*T
    mi_units = isa(mi, Unitful.Quantity) ? mi : mi*u
    return ustrip(ion_gyro_radius(Ti_units,B_units, Zi, mi_units))
end

export electron_gyro_radius, ion_gyro_radius

end #module
