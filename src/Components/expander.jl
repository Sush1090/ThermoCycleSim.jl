
"""
`Expander(type::Isentropic_η;name)`

*    Arguments: 
    1. `type` : `Isentropic_η` contains --> isentropic effeiciency and pressure ratio parameters
    
*    Local Variables:
    1. `P` : Power required 

*    Port Variables:
    1. `inport`  : `p` and `h`
    2. `outport` : `p` and `h`
"""
function Expander(type::Isentropic_η;name)
    @unpack η,πc = type
    @named inport = CoolantPort()
    @named outport = CoolantPort()
    para = @parameters begin
        
    end
    vars = @variables begin
        P(t)
     end
   eqs = [  outport.mdot ~ abs(inport.mdot) 
            outport.p ~ πc * inport.p
            outport.h ~ IsentropicExpansion(πc, inport.h, inport.p,fluid,η)
            P ~ abs(inport.mdot)*(outport.h - inport.h)
   ]
   compose(ODESystem(eqs, t, vars, para;name), inport, outport)
end




export Expander