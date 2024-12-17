
"""
`Isentropic_η(;η = 0.9,πc = 5)`
A struct that passes to `Compressor` and `Expander`
"""
struct Isentropic_η
    η
    πc
    function Isentropic_η(;η = 0.9,πc = 5)
        @assert η <=1
        @assert πc >= 1
        new(η,πc)
    end
end
export Isentropic_η

"""
`Compressor(type::Isentropic_η=Isentropic_η();name,fluid = set_fluid)`

*    Arguments: 
    1. `type` : `Isentropic_η` contains --> isentropic effeiciency and pressure ratio parameters
    
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
function Compressor(type::Isentropic_η=Isentropic_η();name,fluid = set_fluid) 
    if isnothing(fluid)
        throw(error("Fluid not selected"))
    end
    # @unpack η,πc = type
    @named inport = CoolantPort()
    @named outport = CoolantPort()
    para = @parameters begin
        η=0.8, [description = "Isentropic Effeciency",bounds = (0.0,1.0)]
        πc, [description = "Pressure ratio"]
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
            outport.p ~ πc * inport.p
            outport.h ~ IsentropicCompression(πc, inport.h, inport.p,fluid,η)
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



struct Isothermal_comp
    πc
    function Isothermal_comp(;πc = 5)
        @assert πc >=1
        new(πc)
    end
end
export Isothermal_comp

"""
`Compressor(type::Isothermal_comp;name,fluid = set_fluid)`

*    Arguments: 
    1. `type` : `Isothermal_comp` contains --> pressure ratio 
    
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
function Compressor(type::Isothermal_comp;name,fluid = set_fluid) 
    if isnothing(fluid)
        throw(error("Fluid not selected"))
    end
    # @unpack πc = type
    @named inport = CoolantPort()
    @named outport = CoolantPort()
    para = @parameters begin
        η, [description = "Isentropic Effeciency"]
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
            outport.p ~ πc * inport.p
            outport.h ~ IsothermalCompression(πc, inport.h, inport.p,fluid)
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


struct Isochoric_comp
πc
    function Isochoric_comp(;πc=5)
        @assert πc >=1
        new(πc)
    end
end
export Isochoric_comp
"""
`Compressor(type::Isochoric_comp;name,fluid = set_fluid)`

*    Arguments: 
    1. `type` : `Isothermal_comp` contains --> pressure ratio 
    
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
function Compressor(type::Isochoric_comp;name,fluid = set_fluid) 
    if isnothing(fluid)
        throw(error("Fluid not selected"))
    end
    # @unpack πc = type
    @named inport = CoolantPort()
    @named outport = CoolantPort()
    para = @parameters begin
        πc, [description = "Pressure ratio"]
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
            outport.p ~ πc * inport.p
            outport.h ~ IsochoricCompression(πc, inport.h, inport.p,fluid)
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


export Compressor


