using ModelingToolkit

ModelingToolkit.@independent_variables t
D = Differential(t)
const AtmosphericPressure = 101305 #Pa
const AmbientTemperature = 300 #K

PropsSI(out::AbstractString, name1::AbstractString, value1::Real, name2::AbstractString, value2::Real, fluid::AbstractString) = CoolProp.PropsSI(out, name1, value1, name2, value2, fluid)
@register_symbolic PropsSI(out::AbstractString, name1::AbstractString, value1::Real, name2::AbstractString, value2::Real, fluid::AbstractString)

PhaseSI(name1::AbstractString, value1::Real, name2::AbstractString, value2::Real, fluid::AbstractString) = CoolProp.PhaseSI(name1, value1, name2, value2, fluid)
@register_symbolic PhaseSI(name1::AbstractString, value1::Real, name2::AbstractString, value2::Real, fluid::AbstractString)


"""
Makes single node at ports. This node is Pressure,Enthalpy and Massflowrate
"""
@connector  function CoolantPort(;name) 
    vars = @variables begin 
        p(t),  [input = true]
        h(t)
        mdot(t), [connect = Flow] # alternates sign of variable after every component.
    end
    ODESystem(Equation[], t, vars, [];name=name)
end


"""
This is double node at ports. Inport and outport. 
"""
function CoolComponent(;name) 
    @named inport = CoolantPort()
    @named outport = CoolantPort()
    vars = @variables begin
        Δp(t) 
        Δmdot(t) 
        Δh(t) 
    end
    eqs = [
        Δp ~ outport.p - inport.p
        Δh ~ outport.h - inport.h
        Δmdot ~ outport.mdot - inport.mdot
    ]
    compose(ODESystem(eqs, t, vars, [];name=name), inport, outport)
end


"""
Mass source -  Use when the cycle needs a start point. Requires initial temperature,pressure and Massflowrate
"""
function MassSource(;name,source_pressure = 101305,source_enthalpy=1e6,source_mdot=5,fluid) 
    @named port = CoolantPort()
    para = @parameters begin

    end
    vars = @variables begin
        mdot(t)
        s(t)
        p(t)
        T(t)
        h(t)
        ρ(t)
     end

    eqs = [
        port.mdot ~ source_mdot # Outflow is negative
        port.p ~ source_pressure
        port.h ~ source_enthalpy
        mdot ~ port.mdot
        s ~ PropsSI("S","H",port.h,"P",port.p,fluid)
        p ~ port.p
        T ~ PropsSI("T","H",port.h,"P",port.p,fluid)
        h ~ port.h
        ρ ~ PropsSI("D","H",port.h,"P",port.p,fluid)
    ]
    compose(ODESystem(eqs, t, vars, para;name),port)
end

"""
Mass sink -  Use when the cycle needs a end point. Sets the final port input values to the variables
"""
function MassSink(;name,fluid) 
    @named    port = CoolantPort()
    para = @parameters begin
        
    end
    vars = @variables begin
        mdot(t)
        s(t)
        p(t)
        T(t)
        h(t)
        ρ(t)
     end

   eqs = [
    port.p ~ p
    port.h ~ h
    mdot ~ port.mdot
    s ~ PropsSI("S","H",port.h,"P",port.p,fluid)
    p ~ port.p
    T ~ PropsSI("T","H",port.h,"P",port.p,fluid)
    h ~ port.h
    ρ ~ PropsSI("D","H",port.h,"P",port.p,fluid)
   ]
   compose(ODESystem(eqs, t, vars, para;name),port)
end

# """
# CycleSource: An attempt to close the cycle
# """
# function CycleSource(;name,source_pressure = 101305,source_enthalpy=1e6,source_mdot=5)
#     @named inport = CoolantPort()
#     @named outport = CoolantPort()
#     para = @parameters begin
        
#     end
#     vars = @variables begin
#         mdot(t)
#         s(t)
#         p(t)
#         T(t)
#         h(t)
#         ρ(t)
#      end
#      eqs = [
#         inport.mdot ~ abs(source_mdot) # Outflow is negative
#         inport.p ~ source_pressure
#         inport.h ~ source_enthalpy
#     s ~ PropsSI("S","H",port.h,"P",port.p,fluid)
#     p ~ port.p
#     T ~ PropsSI("T","H",port.h,"P",port.p,fluid)
#     h ~ port.h
#     ρ ~ PropsSI("D","H",port.h,"P",port.p,fluid)
#             outport.h ~ source_enthalpy
#             outport.p ~ source_pressure
#             outport.mdot ~ source_mdot 
#    ]
#    compose(ODESystem(eqs, t, vars, para;name),inport,outport)

# end



"""
ComputeSpecificLatentHeat: Computes the specific latent heat of the give fluid at a particular varliable value. var1 should not be enthalpy or vapour quality
*Arguments:
-'var1'     : Variable String. Uses CoolProp variable strings
-'value1'   : Value of the variale chosen
-'fluid'    : Fluid name string
"""
function ComputeSpecificLatentHeat(var1::AbstractString,value1,fluid::AbstractString)
    @assert var1 != "Q"
    H_L = PropsSI("H",var1,value1,"Q",0,fluid)
    H_V = PropsSI("H",var1,value1,"Q",1,fluid)
    return H_V - H_L
end
@register_symbolic ComputeSpecificLatentHeat(var1::AbstractString,value1,fluid::AbstractString)

export CoolantPort,CoolComponent,MassSource,MassSink,AmbientTemperature,AtmosphericPressure,ComputeSpecificLatentHeat, CycleSource