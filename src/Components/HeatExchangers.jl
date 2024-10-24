


function HEX_type(;T_in = nothing,T_out = nothing,Qdot = nothing, mdot = nothing,p = nothing)

end

"""
`IsobaricHeatSource(;name,Q_dot,fluid)`
   A heat source independent of temperature and no pressure drop
*    Arguments: 
    1. `Q_dot`     : Total heat supplied
    
*    Local Variables:
    1. `P`      : Power  
    2. `s_in`   : Inlet Entropy
    3. `p_in`   : Inlet Pressure
    4. `T_in`   : Inlet Temperature
    5. `h_in`   : Inlet Enthalpy
    6. `ρ_in`   : Inlet Density
    7. `s_out`  : Outlet Entropy
    8. `p_out`  : Outlet Pressure
    9. `T_out`  : Outlet Temperature
    10. `h_out` : Outlet Enthalpy
    11. `ρ_out` : Outlet Density

*    Port Variables:
    1. `inport`         : `p` and `h`
    2. `outport`        : `p` and `h`
"""
function IsobaricHeatSource(;name,Q_dot,fluid::AbstractString = set_fluid)
    if isnothing(fluid)
        throw(error("Fluid not selected"))
    end
    @assert Q_dot >= 0
    @named inport = CoolantPort()
    @named outport = CoolantPort()
    para = @parameters begin
        
    end
    vars = @variables begin
        P(t)
        s_in(t)
        p_in(t)
        T_in(t)
        h_in(t)
        ρ_in(t)

        s_out(t)
        p_out(t)
        T_out(t)
        h_out(t)
        ρ_out(t)
     end
   eqs = [  outport.mdot ~ abs(inport.mdot) 
            outport.p ~ inport.p
            outport.h ~ inport.h + Q_dot/outport.mdot
            P ~ abs(inport.mdot)*(outport.h - inport.h)
            s_in ~ PropsSI("S","H",inport.h,"P",inport.p,fluid)
            p_in ~ inport.p
            T_in ~ PropsSI("T","H",inport.h,"P",inport.p,fluid)
            h_in ~ inport.h
            ρ_in ~ PropsSI("D","H",inport.h,"P",inport.p,fluid)
            s_out ~ PropsSI("S","H",outport.h,"P",outport.p,fluid)
            p_out ~ outport.p
            T_out ~ PropsSI("T","H",outport.h,"P",outport.p,fluid)
            h_out ~ outport.h
            ρ_out ~ PropsSI("D","H",outport.h,"P",outport.p,fluid)
   ]
   compose(ODESystem(eqs, t, vars, para;name), inport, outport)
end

export IsobaricHeatSource

"""
`IsobaricHeatSink(;name,Q_dot,fluid)`
   A heat sink independent of temperature and no pressure drop
*    Arguments: 
    1. `Q_dot`     : Total heat supplied
    
*    Local Variables:
    1. `P`      : Power  
    2. `s_in`   : Inlet Entropy
    3. `p_in`   : Inlet Pressure
    4. `T_in`   : Inlet Temperature
    5. `h_in`   : Inlet Enthalpy
    6. `ρ_in`   : Inlet Density
    7. `s_out`  : Outlet Entropy
    8. `p_out`  : Outlet Pressure
    9. `T_out`  : Outlet Temperature
    10. `h_out` : Outlet Enthalpy
    11. `ρ_out` : Outlet Density

*    Port Variables:
    1. `inport`         : `p` and `h`
    2. `outport`        : `p` and `h`
"""
function IsobaricHeatSink(;name,Q_dot,fluid::AbstractString = set_fluid)
    if isnothing(fluid)
        throw(error("Fluid not selected"))
    end
    @assert Q_dot <= 0 
    @named inport = CoolantPort()
    @named outport = CoolantPort()
    para = @parameters begin
        
    end
    vars = @variables begin
        P(t)
        s_in(t)
        p_in(t)
        T_in(t)
        h_in(t)
        ρ_in(t)

        s_out(t)
        p_out(t)
        T_out(t)
        h_out(t)
        ρ_out(t)
     end
   eqs = [  outport.mdot ~ abs(inport.mdot) 
            outport.p ~ inport.p
            outport.h ~ inport.h + Q_dot/outport.mdot
            P ~ abs(inport.mdot)*(outport.h - inport.h)
            s_in ~ PropsSI("S","H",inport.h,"P",inport.p,fluid)
            p_in ~ inport.p
            T_in ~ PropsSI("T","H",inport.h,"P",inport.p,fluid)
            h_in ~ inport.h
            ρ_in ~ PropsSI("D","H",inport.h,"P",inport.p,fluid)
            s_out ~ PropsSI("S","H",outport.h,"P",outport.p,fluid)
            p_out ~ outport.p
            T_out ~ PropsSI("T","H",outport.h,"P",outport.p,fluid)
            h_out ~ outport.h
            ρ_out ~ PropsSI("D","H",outport.h,"P",outport.p,fluid)
   ]
   compose(ODESystem(eqs, t, vars, para;name), inport, outport)
end

export IsobaricHeatSink


"""
`Preheater(;name,fluid,Δp)`
    Preheater for the Evaporator.     
*    Arguments: 
    1. `Δp`     : Pressure Drop across Evaporator
    
*    Local Variables:
    1. `P`      : Power  
    2. `s_in`   : Inlet Entropy
    3. `p_in`   : Inlet Pressure
    4. `T_in`   : Inlet Temperature
    5. `h_in`   : Inlet Enthalpy
    6. `ρ_in`   : Inlet Density
    7. `s_out`  : Outlet Entropy
    8. `p_out`  : Outlet Pressure
    9. `T_out`  : Outlet Temperature
    10. `h_out` : Outlet Enthalpy
    11. `ρ_out` : Outlet Density
    12. `T_sat` : Saturation Temperature

*    Port Variables:
    1. `inport`         : `p` and `h`
    2. `outport`        : `p` and `h`
"""
function Preheater(;name,Δp,fluid::AbstractString = set_fluid)
    if isnothing(fluid)
        throw(error("Fluid not selected"))
    end
    @named inport = CoolantPort()
    @named outport = CoolantPort()
    para = @parameters begin

    end
    vars =  @variables begin
        P(t)
        s_in(t)
        p_in(t)
        T_in(t)
        h_in(t)
        ρ_in(t)

        s_out(t)
        p_out(t)
        T_out(t)
        h_out(t)
        ρ_out(t)

        T_sat(t)
    end
    eqs = [
        outport.mdot ~ abs(inport.mdot)
        outport.p ~ inport.p - Δp
        outport.h ~ PropsSI("H","Q",0,"P",outport.p,fluid)

        s_in ~ PropsSI("S","H",inport.h,"P",inport.p,fluid)
        p_in ~ inport.p
        T_in ~ PropsSI("T","H",inport.h,"P",inport.p,fluid)
        h_in ~ inport.h
        ρ_in ~ PropsSI("D","H",inport.h,"P",inport.p,fluid)
            
        s_out ~ PropsSI("S","H",outport.h,"P",outport.p,fluid)
        p_out ~ outport.p
        T_out ~ PropsSI("T","H",outport.h,"P",outport.p,fluid)
        h_out ~ outport.h
        ρ_out ~ PropsSI("D","H",outport.h,"P",outport.p,fluid)
        
        P ~ outport.mdot*(outport.h - inport.h)
        
        T_sat ~ PropsSI("T","Q",0,"P",outport.p,fluid)
    ]
    compose(ODESystem(eqs, t, vars, para;name), inport, outport)
end
"""
`TwoPhaseEvap(;name,fluid,Δp)`
    Twophase part for the Evaporator.     
*    Arguments: 
    1. `Δp`     : Pressure Drop across Evaporator
    
*    Local Variables:
    1. `P`      : Power  
    2. `s_in`   : Inlet Entropy
    3. `p_in`   : Inlet Pressure
    4. `T_in`   : Inlet Temperature
    5. `h_in`   : Inlet Enthalpy
    6. `ρ_in`   : Inlet Density
    7. `s_out`  : Outlet Entropy
    8. `p_out`  : Outlet Pressure
    9. `T_out`  : Outlet Temperature
    10. `h_out` : Outlet Enthalpy
    11. `ρ_out` : Outlet Density
    12. `T_sat` : Saturation Temperature

*    Port Variables:
    1. `inport`         : `p` and `h`
    2. `outport`        : `p` and `h`
"""
function TwoPhaseEvap(;name,Δp,fluid::AbstractString = set_fluid)
    if isnothing(fluid)
        throw(error("Fluid not selected"))
    end
    @named inport = CoolantPort()
    @named outport = CoolantPort()
    para = @parameters begin

    end
    vars =  @variables begin
        P(t)
        s_in(t)
        p_in(t)
        T_in(t)
        h_in(t)
        ρ_in(t)

        s_out(t)
        p_out(t)
        T_out(t)
        h_out(t)
        ρ_out(t)

        T_sat(t)
    end
    eqs = [
        outport.mdot ~ abs(inport.mdot)
        outport.p ~ inport.p - Δp
        outport.h ~ PropsSI("H","Q",1,"P",outport.p,fluid)
        
        s_in ~ PropsSI("S","H",inport.h,"P",inport.p,fluid)
        p_in ~ inport.p
        T_in ~ PropsSI("T","H",inport.h,"P",inport.p,fluid)
        h_in ~ inport.h
        ρ_in ~ PropsSI("D","H",inport.h,"P",inport.p,fluid)
            
        s_out ~ PropsSI("S","H",outport.h,"P",outport.p,fluid)
        p_out ~ outport.p
        T_out ~ PropsSI("T","H",outport.h,"P",outport.p,fluid)
        h_out ~ outport.h
        ρ_out ~ PropsSI("D","H",outport.h,"P",outport.p,fluid)
        
        P ~ outport.mdot*(outport.h - inport.h)
        
        T_sat ~ PropsSI("T","Q",0,"P",outport.p,fluid)
    ]
    compose(ODESystem(eqs, t, vars, para;name), inport, outport)
end

"""
`SuperHeat(;name,fluid,ΔT_sh,Δp)`
    Superheated part of the evaporator.
*    Arguments: 
    1. `Δp`     : Pressure Drop across Evaporator
    2. `ΔT_sh`  : Super heated Temperature
    
*    Local Variables:
    1. `P`      : Power  
    2. `s_in`   : Inlet Entropy
    3. `p_in`   : Inlet Pressure
    4. `T_in`   : Inlet Temperature
    5. `h_in`   : Inlet Enthalpy
    6. `ρ_in`   : Inlet Density
    7. `s_out`  : Outlet Entropy
    8. `p_out`  : Outlet Pressure
    9. `T_out`  : Outlet Temperature
    10. `h_out` : Outlet Enthalpy
    11. `ρ_out` : Outlet Density
    12. `T_sat` : Saturation Temperature

*    Port Variables:
    1. `inport`         : `p` and `h`
    2. `outport`        : `p` and `h`
"""
function SuperHeat(;name,ΔT_sh,Δp,fluid::AbstractString = set_fluid)
    if isnothing(fluid)
        throw(error("Fluid not selected"))
    end
    @named inport = CoolantPort()
    @named outport = CoolantPort()
    para = @parameters begin

    end
    vars =  @variables begin
        P(t)
        s_in(t)
        p_in(t)
        T_in(t)
        h_in(t)
        ρ_in(t)

        s_out(t)
        p_out(t)
        T_out(t)
        h_out(t)
        ρ_out(t)

        T_sat(t)
    end
    eqs = [
        outport.mdot ~ abs(inport.mdot)
        outport.p ~ inport.p - Δp
        T_in ~ PropsSI("T","H",inport.h,"P",outport.p,fluid)
        T_out ~ T_in + ΔT_sh
        outport.h ~ PropsSI("H","T",T_out,"P",outport.p,fluid)

        s_in ~ PropsSI("S","H",inport.h,"P",inport.p,fluid)
        p_in ~ inport.p
        h_in ~ inport.h
        ρ_in ~ PropsSI("D","H",inport.h,"P",inport.p,fluid)
            
        s_out ~ PropsSI("S","H",outport.h,"P",outport.p,fluid)
        p_out ~ outport.p
        h_out ~ outport.h
        ρ_out ~ PropsSI("D","H",outport.h,"P",outport.p,fluid)
        
        P ~ outport.mdot*(outport.h - inport.h)
        
        T_sat ~ PropsSI("T","Q",1,"P",outport.p,fluid)
    ]
    compose(ODESystem(eqs, t, vars, para;name), inport, outport)
end


"""
`SimpleEvaporator(;name,fluid,Δp,ΔT_sh)`
    Composed of multiple `ODESystem` - `Preheater`, `TwoPhaseEvap`, and `SuperHeat`
*    Arguments: 
    1. `Δp`     : Pressure Drop across Evaporator
    2. `ΔT_sh`  : Super heated Temperature
    
*    Local Variables:
    1. `P`      : Power  
    2. `s_in`   : Inlet Entropy
    3. `p_in`   : Inlet Pressure
    4. `T_in`   : Inlet Temperature
    5. `h_in`   : Inlet Enthalpy
    6. `ρ_in`   : Inlet Density
    7. `s_out`  : Outlet Entropy
    8. `p_out`  : Outlet Pressure
    9. `T_out`  : Outlet Temperature
    10. `h_out` : Outlet Enthalpy
    11. `ρ_out` : Outlet Density
    12. `T_sat` : Saturation Temperature

*    Port Variables:
    1. `inport`         : `p` and `h`
    2. `outport`        : `p` and `h`
    3. `Preheater`      : Component with internal `CoolantPort` --> `inport` and `outport`
    4. `TwoPhaseEvap`   : Component with internal `CoolantPort` --> `inport` and `outport`
    5. `SuperHeat`      : Component with internal `CoolantPort` --> `inport` and `outport`
"""
function SimpleEvaporator(;name,Δp::AbstractVector = [0,0,0],ΔT_sh,fluid::AbstractString = set_fluid)
    if isnothing(fluid)
        throw(error("Fluid not selected"))
    end
    @assert size(Δp,1) ==3 "pressure drop vector has to be of size 3 for preheater,twophase and superheatear"
    @assert ΔT_sh > 1e-3 "Keep subcooling temperature away from Saturation curve to avoid CoolProp assertion errors"
    @named inport = CoolantPort()
    @named preheater = Preheater(fluid = fluid,Δp= Δp[1]) 
    @named twophase = TwoPhaseEvap(fluid = fluid,Δp= Δp[2])
    @named superheater = SuperHeat(fluid = fluid,Δp= Δp[3],ΔT_sh = ΔT_sh)
    @named outport = CoolantPort()

    para = @parameters begin

    end
    vars =  @variables begin
        P(t)
        s_in(t)
        p_in(t)
        T_in(t)
        h_in(t)
        ρ_in(t)

        s_out(t)
        p_out(t)
        T_out(t)
        h_out(t)
        ρ_out(t)

        T_sat(t)
    end
    eqs = [
    connect(inport,preheater.inport)
    connect(preheater.outport,twophase.inport)
    connect(twophase.outport,superheater.inport)
    connect(superheater.outport,outport)
    T_sat ~ twophase.T_sat
    T_in ~ preheater.T_in
    T_out ~ superheater.T_out
    P ~ preheater.P + twophase.P + superheater.P

    s_in ~ preheater.s_in
    p_in ~ preheater.p_in
    h_in ~ preheater.h_in
    ρ_in ~ preheater.ρ_in

    s_out ~ superheater.s_out
    p_out ~ superheater.p_out
    h_out ~ superheater.h_out
    ρ_out ~ superheater.ρ_out
]

compose(ODESystem(eqs, t, vars, para;name), inport,preheater,twophase,superheater, outport)

end


export SimpleEvaporator, SuperHeat,Preheater,TwoPhaseEvap


"""
`Precooler(;name,fluid,Δp)`
    Precooler for the Condensor.     
*    Arguments: 
    1. `Δp`     : Pressure Drop across Evaporator
    
*    Local Variables:
    1. `P`      : Power  
    2. `s_in`   : Inlet Entropy
    3. `p_in`   : Inlet Pressure
    4. `T_in`   : Inlet Temperature
    5. `h_in`   : Inlet Enthalpy
    6. `ρ_in`   : Inlet Density
    7. `s_out`  : Outlet Entropy
    8. `p_out`  : Outlet Pressure
    9. `T_out`  : Outlet Temperature
    10. `h_out` : Outlet Enthalpy
    11. `ρ_out` : Outlet Density
    12. `T_sat` : Saturation Temperature

*    Port Variables:
    1. `inport`         : `p` and `h`
    2. `outport`        : `p` and `h`
"""
function Precooler(;name,Δp,fluid::AbstractString = set_fluid)
    if isnothing(fluid)
        throw(error("Fluid not selected"))
    end
    @named inport = CoolantPort()
    @named outport = CoolantPort()
    para = @parameters begin

    end
    vars =  @variables begin
        P(t)
        s_in(t)
        p_in(t)
        T_in(t)
        h_in(t)
        ρ_in(t)

        s_out(t)
        p_out(t)
        T_out(t)
        h_out(t)
        ρ_out(t)

        T_sat(t)
    end
    eqs = [
        outport.mdot ~ abs(inport.mdot)
        outport.p ~ inport.p - Δp
        outport.h ~ PropsSI("H","Q",1,"P",outport.p,fluid)

        s_in ~ PropsSI("S","H",inport.h,"P",inport.p,fluid)
        p_in ~ inport.p
        T_in ~ PropsSI("T","H",inport.h,"P",inport.p,fluid)
        h_in ~ inport.h
        ρ_in ~ PropsSI("D","H",inport.h,"P",inport.p,fluid)
            
        s_out ~ PropsSI("S","H",outport.h,"P",outport.p,fluid)
        p_out ~ outport.p
        T_out ~ PropsSI("T","H",outport.h,"P",outport.p,fluid)
        h_out ~ outport.h
        ρ_out ~ PropsSI("D","H",outport.h,"P",outport.p,fluid)
        
        P ~ outport.mdot*(outport.h - inport.h)
        
        T_sat ~ PropsSI("T","Q",1,"P",outport.p,fluid) 
    ]
    compose(ODESystem(eqs, t, vars, para;name), inport, outport)
end


"""
`TwoPhaseCond(;name,fluid,Δp)`
    Twophase part for the Condensor.     
*    Arguments: 
    1. `Δp`     : Pressure Drop across Evaporator
    
*    Local Variables:
    1. `P`      : Power  
    2. `s_in`   : Inlet Entropy
    3. `p_in`   : Inlet Pressure
    4. `T_in`   : Inlet Temperature
    5. `h_in`   : Inlet Enthalpy
    6. `ρ_in`   : Inlet Density
    7. `s_out`  : Outlet Entropy
    8. `p_out`  : Outlet Pressure
    9. `T_out`  : Outlet Temperature
    10. `h_out` : Outlet Enthalpy
    11. `ρ_out` : Outlet Density
    12. `T_sat` : Saturation Temperature

*    Port Variables:
    1. `inport`         : `p` and `h`
    2. `outport`        : `p` and `h`  
"""
function TwoPhaseCond(;name,Δp,fluid::AbstractString = set_fluid)
    if isnothing(fluid)
        throw(error("Fluid not selected"))
    end
    @named inport = CoolantPort()
    @named outport = CoolantPort()
    para = @parameters begin

    end
    vars =  @variables begin
        P(t)
        s_in(t)
        p_in(t)
        T_in(t)
        h_in(t)
        ρ_in(t)

        s_out(t)
        p_out(t)
        T_out(t)
        h_out(t)
        ρ_out(t) 

        T_sat(t)
    end
    eqs = [
        outport.mdot ~ abs(inport.mdot)
        outport.p ~ inport.p - Δp
        outport.h ~ PropsSI("H","Q",0,"P",outport.p,fluid)
        
        s_in ~ PropsSI("S","H",inport.h,"P",inport.p,fluid)
        p_in ~ inport.p
        T_in ~ PropsSI("T","H",inport.h,"P",inport.p,fluid)
        h_in ~ inport.h
        ρ_in ~ PropsSI("D","H",inport.h,"P",inport.p,fluid)
            
        s_out ~ PropsSI("S","H",outport.h,"P",outport.p,fluid)
        p_out ~ outport.p
        T_out ~ PropsSI("T","H",outport.h,"P",outport.p,fluid)
        h_out ~ outport.h
        ρ_out ~ PropsSI("D","H",outport.h,"P",outport.p,fluid)
        
        P ~ outport.mdot*(outport.h - inport.h)
        
        T_sat ~ PropsSI("T","Q",0,"P",outport.p,fluid)
    ]
    compose(ODESystem(eqs, t, vars, para;name), inport, outport)
end



"""
`SuperCooler(;name,fluid,ΔT_sh,Δp)`
    SuperCooler part of the Condensor.
*    Arguments: 
    1. `Δp`     : Pressure Drop across Evaporator
    2. `ΔT_sh`  : Super heated Temperature
    
*    Local Variables:
    1. `P`      : Power  
    2. `s_in`   : Inlet Entropy
    3. `p_in`   : Inlet Pressure
    4. `T_in`   : Inlet Temperature
    5. `h_in`   : Inlet Enthalpy
    6. `ρ_in`   : Inlet Density
    7. `s_out`  : Outlet Entropy
    8. `p_out`  : Outlet Pressure
    9. `T_out`  : Outlet Temperature
    10. `h_out` : Outlet Enthalpy
    11. `ρ_out` : Outlet Density
    12. `T_sat` : Saturation Temperature

*    Port Variables:
    1. `inport`         : `p` and `h`
    2. `outport`        : `p` and `h`
"""
function SuperCooler(;name,ΔT_sc,Δp,fluid::AbstractString = set_fluid)
    if isnothing(fluid)
        throw(error("Fluid not selected"))
    end
    @named inport = CoolantPort()
    @named outport = CoolantPort()
    para = @parameters begin

    end
    vars =  @variables begin
        P(t)
        s_in(t)
        p_in(t)
        T_in(t)
        h_in(t)
        ρ_in(t)

        s_out(t)
        p_out(t)
        T_out(t)
        h_out(t)
        ρ_out(t)

        T_sat(t)
    end
    eqs = [
        outport.mdot ~ abs(inport.mdot)
        outport.p ~ inport.p - Δp
        T_in ~ PropsSI("T","H",inport.h,"P",outport.p,fluid)
        T_out ~ T_in - ΔT_sc
        outport.h ~ PropsSI("H","T",T_out,"P",outport.p,fluid)

        s_in ~ PropsSI("S","H",inport.h,"P",inport.p,fluid)
        p_in ~ inport.p
        h_in ~ inport.h
        ρ_in ~ PropsSI("D","H",inport.h,"P",inport.p,fluid)
            
        s_out ~ PropsSI("S","H",outport.h,"P",outport.p,fluid)
        p_out ~ outport.p
        h_out ~ outport.h
        ρ_out ~ PropsSI("D","H",outport.h,"P",outport.p,fluid)
        
        P ~ outport.mdot*(outport.h - inport.h)
        
        T_sat ~ PropsSI("T","Q",1,"P",outport.p,fluid)
    ]
    compose(ODESystem(eqs, t, vars, para;name), inport, outport)
end



"""
`SimpleCondensor(;name,fluid,Δp,ΔT_sh)`
    Composed of multiple `ODESystem` - `Precooler`, `TwoPhaseCond`, and `SuperCooler`
*    Arguments: 
    1. `Δp`     : Pressure Drop across Evaporator
    2. `ΔT_sh`  : Super heated Temperature
    
*    Local Variables:
    1. `P`      : Power  
    2. `s_in`   : Inlet Entropy
    3. `p_in`   : Inlet Pressure
    4. `T_in`   : Inlet Temperature
    5. `h_in`   : Inlet Enthalpy
    6. `ρ_in`   : Inlet Density
    7. `s_out`  : Outlet Entropy
    8. `p_out`  : Outlet Pressure
    9. `T_out`  : Outlet Temperature
    10. `h_out` : Outlet Enthalpy
    11. `ρ_out` : Outlet Density
    12. `T_sat` : Saturation Temperature

*    Port Variables:
    1. `inport`         : `p` and `h`
    2. `outport`        : `p` and `h`
    3. `Precooler`      : Component with internal `CoolantPort` --> `inport` and `outport`
    4. `TwoPhaseCond`   : Component with internal `CoolantPort` --> `inport` and `outport`
    5. `SuperCooler`   : Component with internal `CoolantPort` --> `inport` and `outport`
"""
function SimpleCondensor(;name,Δp::AbstractVector = [0,0,0],ΔT_sc,fluid::AbstractString = set_fluid)
    if isnothing(fluid)
        throw(error("Fluid not selected"))
    end
    @assert size(Δp,1) ==3 "pressure drop vector has to be of size 3 for Precooler,TwoPhaseCond and SuperCooler"
    @assert ΔT_sc > 1e-3 "Keep subcooling temperature away from Saturation curve to avoid CoolProp assertion errors"
    @named inport = CoolantPort()
    @named precooler = Precooler(fluid = fluid,Δp= Δp[1]) 
    @named twophasecond = TwoPhaseCond(fluid = fluid,Δp= Δp[2])
    @named supercooler = SuperCooler(fluid = fluid,Δp= Δp[3],ΔT_sc = ΔT_sc)
    @named outport = CoolantPort()

    para = @parameters begin

    end
    vars =  @variables begin
        P(t)
        s_in(t)
        p_in(t)
        T_in(t)
        h_in(t)
        ρ_in(t)

        s_out(t)
        p_out(t)
        T_out(t)
        h_out(t)
        ρ_out(t)

        T_sat(t)
    end
    eqs = [
    connect(inport,precooler.inport)
    connect(precooler.outport,twophasecond.inport)
    connect(twophasecond.outport,supercooler.inport)
    connect(supercooler.outport,outport)
    T_sat ~ twophasecond.T_sat
    T_in ~ precooler.T_in
    T_out ~ supercooler.T_out
    P ~ precooler.P + twophasecond.P + supercooler.P

    s_in ~ precooler.s_in
    p_in ~ precooler.p_in
    h_in ~ precooler.h_in
    ρ_in ~ precooler.ρ_in

    s_out ~ supercooler.s_out
    p_out ~ supercooler.p_out
    h_out ~ supercooler.h_out
    ρ_out ~ supercooler.ρ_out
]

compose(ODESystem(eqs, t, vars, para;name), inport,precooler,twophasecond,supercooler, outport)

end

export SimpleCondensor, Precooler, TwoPhaseCond, SuperCooler