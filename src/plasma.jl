module PlasmaParameters

using Unitful
using ..Units
using ..Constants: q, me, ϵ0

function electron_cyclotron_frequency(B::Unitful.BField)
    return uconvert(rad/s, q*B/me)
end

function ion_cyclotron_frequency(B::Unitful.BField, Zi::Int, mi::Unitful.Mass)
    return uconvert(rad/s, Zi*q*B/mi)
end

export electron_cyclotron_frequency, ion_cyclotron_frequency

function electron_plasma_frequency(ne::NumberDensity)
    return uconvert(rad/s,sqrt((ne*q^2)/(me*ϵ0)))
end

function ion_plasma_frequency(ni::NumberDensity, Zi::Int, mi::Unitful.Mass)
    return uconvert(rad/s,sqrt((ni*(Zi*q)^2)/(mi*ϵ0)))
end

export electron_plasma_frequency, ion_plasma_frequency

function electron_thermal_speed(Te::Unitful.Energy)
    return uconvert(m/s,sqrt(Te/me))
end

function ion_thermal_speed(Ti::Unitful.Energy, mi::Unitful.Mass)
    return uconvert(m/s,sqrt(Ti/mi))
end

export electron_thermal_speed, ion_thermal_speed

function electron_gyro_radius(Te::Unitful.Energy, B::Unitful.BField)
    return uconvert(m, sqrt(2*me*Te)/(q*B))
end

function ion_gyro_radius(Ti::Unitful.Energy, B::Unitful.BField, Zi::Int, mi::Unitful.Mass)
    return uconvert(m, sqrt(2*mi*Ti)/(Zi*q*B))
end

export electron_gyro_radius, ion_gyro_radius
end #module
