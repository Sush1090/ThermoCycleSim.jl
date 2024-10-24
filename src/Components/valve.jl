
struct IsenthalpicExpansionValve
    πc
end

"""
`Valve(;name,πc,fluid): Isenthalpic valve`
"""
function Valve(;name,πc,fluid::AbstractString = set_fluid) 
    if isnothing(fluid)
        throw(error("Fluid not selected"))
    end
    @named inport = CoolantPort()
    @named outport = CoolantPort()
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
    para = @parameters begin
        
    end
    eqs = [
        outport.mdot ~ abs(inport.mdot) 
        outport.p ~  inport.p/πc
        outport.h ~ inport.h
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
"""
`Valve(type::IsenthalpicExpansionValve;name,fluid): Isenthalpic valve`

*    Arguments: 
    1. `type` : `IsenthalpicExpansionValve` contains --> pressure ratio 
    
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
    1. `inport`  : `p` and `h`
    2. `outport` : `p` and `h`
"""
function Valve(type::IsenthalpicExpansionValve;name,fluid = set_fluid) 
    if isnothing(fluid)
        throw(error("Fluid not selected"))
    end
    @unpack πc = type
    @named inport = CoolantPort()
    @named outport = CoolantPort()
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
    para = @parameters begin
        
    end
    eqs = [
        outport.mdot ~ abs(inport.mdot) 
        outport.p ~  inport.p/πc
        outport.h ~ inport.h
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

export IsenthalpicExpansionValve, Valve
