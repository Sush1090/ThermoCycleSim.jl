
begin 
pt_entropy(model::EoSModel,p,T,z) = Clapeyron.entropy(model::EoSModel,p,T,z,phase = "unknown")
@register_symbolic pt_entropy(model::EoSModel,p,T,z::Array)

pt_enthalpy(model::EoSModel,p,T,z) = Clapeyron.enthalpy(model::EoSModel,p,T,z,phase = "unknown")
@register_symbolic pt_enthalpy(model::EoSModel,p,T,z::Array)

ph_temperature(model::EoSModel,p,h,z) = PH.temperature(model::EoSModel,p,h,z,phase = "unknown")
@register_symbolic ph_temperature(model::EoSModel,p,h,z::Array)

ph_entropy(model::EoSModel,p,h,z) = PH.entropy(model::EoSModel,p,h,z,phase = "unkown")
@register_symbolic ph_entropy(model::EoSModel,p,h,z::Array)

ph_mass_density(model::EoSModel,p,h,z) = PH.mass_density(model::EoSModel,p,h,z,phase = "unkown")
@register_symbolic ph_mass_density(model::EoSModel,p,h,z::Array)

ph_volume(model::EoSModel,p,h,z) = PH.volume(model::EoSModel,p,h,z,phase = "unkown")
@register_symbolic ph_volume(model::EoSModel,p,h,z::Array)

ps_temperature(model::EoSModel,p,s,z) = PS.temperature(model::EoSModel,p,s,z,phase = "unknown")
@register_symbolic ps_temperature(model::EoSModel,p,s,z::Array)

ps_enthalpy(model::EoSModel,p,s,z) = PS.enthalpy(model::EoSModel,p,s,z,phase = "unknown")
@register_symbolic ps_enthalpy(model::EoSModel,p,s,z::Array)

Tproperty_S(model::EoSModel,p,s,z) = Clapeyron.Tproperty(model::EoSModel,p,s,z,entropy,phase = "unkown")
@register_symbolic Tproperty_S(model::EoSModel,p,s,z::Array)

Tproperty_H(model::EoSModel,p,s,z) = Clapeyron.Tproperty(model::EoSModel,p,s,z,enthalpy,phase = "unkown",verbose = false)
@register_symbolic Tproperty_S(model::EoSModel,p,s,z::Array)

Bubble_temperature(model::EoSModel,p,z) = Clapeyron.bubble_temperature(model,p,z,FugBubbleTemperature())[1]
@register_symbolic Bubble_temperature(model::EoSModel,p,z::Array)

Dew_temperature(model::EoSModel,p,z) = Clapeyron.dew_temperature(model,p,z,FugDewTemperature())[1]
@register_symbolic Dew_temperature(model::EoSModel,p,z::Array)

Bubble_pressure(model::EoSModel,T,z) = Clapeyron.bubble_pressure(model,T,z,FugBubblePressure())[1]
@register_symbolic Bubble_pressure(model::EoSModel,T,z::Array)

Dew_pressure(model::EoSModel,T,z) = Clapeyron.dew_pressure(model,T,z,FugDewPressure())[1]
@register_symbolic Dew_pressure(model::EoSModel,T,z::Array)




function CritcalTemperature(model::EoSModel,z)
    if size(z,1) == 1
        return Clapeyron.crit_pure(model)[1]
    end
    if size(z,1) != 1
        return Clapeyron.crit_mix(model,z)[1]
    end
end
@register_symbolic CritcalTemperature(model::EoSModel,z::Array)

function CriticalPressure(model::EoSModel,z)
    if size(z,1) == 1
        return Clapeyron.crit_pure(model)[2]
    end
    if size(z,1) != 1
        return Clapeyron.crit_mix(model,z)[2]
    end
end
@register_symbolic CriticalPressure(model::EoSModel,z::Array)

end