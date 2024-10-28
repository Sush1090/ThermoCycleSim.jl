
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
        πc ~  inport.p/outport.p
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





struct ThreeFacedValveSplit
ratio::AbstractVector
function ThreeFacedValveSplit(;ratio::AbstractVector)
    @assert isapprox(sum(ratio),1) "Sum of the ratio of mass flow rates should be 1"
    @assert size(ratio,1) == 2 "Only splits in two streams hence provide only two ratios"
    new(ratio)
end
end


"""
`Valve(type::ThreePhaseValveSplit;name,fluid=set_fluid)`

"""
function Valve(type::ThreeFacedValveSplit;name,fluid=set_fluid)
    if isnothing(fluid)
        throw(error("Fluid not selected"))
    end
    @unpack ratio= type
    @named inport = CoolantPort()
    @named outport1 = CoolantPort()
    @named outport2 = CoolantPort()
    vars = @variables begin
        P1(t)
        P2(t)
        P(t)

        s_in(t)
        p_in(t)
        T_in(t)
        h_in(t)
        ρ_in(t)

        s_out1(t)
        p_out1(t)
        T_out1(t)
        h_out1(t)
        ρ_out1(t)

        s_out2(t)
        p_out2(t)
        T_out2(t)
        h_out2(t)
        ρ_out2(t)
    end
    para = @parameters begin
        
    end
    eqs = [
        outport1.mdot ~ abs(inport.mdot)*ratio[1]
        outport2.mdot ~ abs(inport.mdot)*ratio[2]
        outport1.p ~  inport.p
        outport2.p ~  inport.p
        outport1.h ~ inport.h
        outport2.h ~ inport.h
        P1 ~ abs(outport1.mdot)*(outport1.h - inport.h)
        P2 ~ abs(outport2.mdot)*(outport2.h - inport.h)
        P ~ P1+ P2
        s_in ~ PropsSI("S","H",inport.h,"P",inport.p,fluid)
        p_in ~ inport.p
        T_in ~ PropsSI("T","H",inport.h,"P",inport.p,fluid)
        h_in ~ inport.h
        ρ_in ~ PropsSI("D","H",inport.h,"P",inport.p,fluid)

        s_out1 ~ PropsSI("S","H",outport1.h,"P",outport1.p,fluid)
        p_out1 ~ outport1.p
        T_out1 ~ PropsSI("T","H",outport1.h,"P",outport1.p,fluid)
        h_out1 ~ outport1.h
        ρ_out1 ~ PropsSI("D","H",outport1.h,"P",outport1.p,fluid)

        s_out2 ~ PropsSI("S","H",outport2.h,"P",outport2.p,fluid)
        p_out2 ~ outport2.p
        T_out2 ~ PropsSI("T","H",outport2.h,"P",outport2.p,fluid)
        h_out2 ~ outport2.h
        ρ_out2 ~ PropsSI("D","H",outport2.h,"P",outport2.p,fluid)
    ]
    compose(ODESystem(eqs, t, vars, para;name), inport, outport1,outport2)
end



struct ThreeFacedValveCombine

end


"""
`Valve(type::ThreePhaseValveSplit;name,fluid=set_fluid)`

"""
function Valve(type::ThreeFacedValveCombine;name,fluid=set_fluid)
    if isnothing(fluid)
        throw(error("Fluid not selected"))
    end
    @named inport1 = CoolantPort()
    @named inport2 = CoolantPort()
    @named outport = CoolantPort()
    vars = @variables begin
        r(t)

        P1(t)
        P2(t)
        P(t)

        s_in2(t)
        p_in2(t)
        T_in2(t)
        h_in2(t)
        ρ_in2(t)

        s_in1(t)
        p_in1(t)
        T_in1(t)
        h_in1(t)
        ρ_in1(t)

        s_out(t)
        p_out(t)
        T_out(t)
        h_out(t)
        ρ_out(t)
    end
    para = @parameters begin
        
    end
    eqs = [
        outport.mdot ~ abs(inport1.mdot) + abs(inport2.mdot) 
        outport.p ~  inport1.p*(abs(inport1.mdot)/outport.mdot) + inport2.p*(abs(inport2.mdot)/outport.mdot) 
        outport.h ~ inport1.h*(abs(inport1.mdot)/outport.mdot) + inport2.h*(abs(inport2.mdot)/outport.mdot)
        
        P1 ~ (abs(inport1.mdot)/outport.mdot)*(outport.h - inport1.h)
        P2 ~ (abs(inport2.mdot)/outport.mdot)*(outport.h - inport2.h)
        P ~ P1+ P2
        
        s_in1 ~ PropsSI("S","H",inport1.h,"P",inport1.p,fluid)
        p_in1 ~ inport1.p
        T_in1 ~ PropsSI("T","H",inport1.h,"P",inport1.p,fluid)
        h_in1 ~ inport1.h
        ρ_in1 ~ PropsSI("D","H",inport1.h,"P",inport1.p,fluid)

        s_in2 ~ PropsSI("S","H",inport2.h,"P",inport2.p,fluid)
        p_in2 ~ inport2.p
        T_in2 ~ PropsSI("T","H",inport2.h,"P",inport2.p,fluid)
        h_in2 ~ inport2.h
        ρ_in2 ~ PropsSI("D","H",inport2.h,"P",inport2.p,fluid)


        s_out ~ PropsSI("S","H",outport.h,"P",outport.p,fluid)
        p_out ~ outport.p
        T_out ~ PropsSI("T","H",outport.h,"P",outport.p,fluid)
        h_out ~ outport.h
        ρ_out ~ PropsSI("D","H",outport.h,"P",outport.p,fluid)
    ]
    compose(ODESystem(eqs, t, vars, para;name), inport1, inport2,outport)
end



export IsenthalpicExpansionValve, Valve, ThreeFacedValveSplit,ThreeFacedValveCombine