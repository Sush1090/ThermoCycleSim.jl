# CoolPropCycles

[![Build Status](https://github.com/Sush1090/CoolPropCycles.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/Sush1090/CoolPropCycles.jl/actions/workflows/CI.yml?query=branch%3Amain)


This package combines ModelingToolkit.jl with CoolProp.jl to make basic thermodynamic cycles:
    1. Organic Rakine Cycle
    2. Vapour Compression Cycle

It is for steady state calculations, hence the evaporators and compressors are  not only modeled using energy conservation. 

It is extendible as users can add their own component in the following manner 

```
function MyComp(type::abc;name,...)
    @named inport = CoolantPort()
    @named outport = CoolantPort()
    para = @parameters begin
        MyParas ...
    end
    vars = @variables begin
        MyVars ...
     end
   eqs = [  outport.mdot ~ abs(inport.mdot) 
            outport.p ~ eq1 ...
            outport.h ~ eq2 ...
            ..
   ]
   compose(ODESystem(eqs, t, vars, para;name), inport, outport)
end
```

The `CoolantPort()` is fixed for the pressure `p`, enthalpy `h`, and mass flow rate `mdot`. The mass flow rate of the system is conserved throught the system. 
