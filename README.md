# CoolPropCycles.jl

[![Build Status](https://github.com/Sush1090/CoolPropCycles.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/Sush1090/CoolPropCycles.jl/actions/workflows/CI.yml?query=branch%3Amain)


This package combines ModelingToolkit.jl with CoolProp.jl to make basic thermodynamic cycles:

    1. Organic Rakine Cycle
    2. Vapour Compression Cycle
    3. Brayton Cycle

It is for steady state calculations, hence the evaporators and compressors are  not only modeled using energy conservation. The principle is state-point manipulations.  The system works better for subcritical cycle parameters

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


It provides the following components:

    1. Compressor using Isentropic Efficiency 
    2. Expander using Isentropic Efficiency
    3. Isenthalpic Expander 
    4. Evaporator
    5. Condensor
    6. Heat Exchanger -  No phase change - Recuperator
    7. Electric Heating

It also provides basic functions that find the pressure to match the pitch points. 


## Basic Usage

## Examples 
### Organic Rankine Cycle
Example of Organic Rankine Cycle

### Vapour Compression cycle
Example of Vapour Compression

## Phase Diagram Plotting

A phase ploting routine is provided for the given cycle. 

## Basic Thermodynamic Cycle Optimization
The cycles created can be wrapped with functions and sent to optimization routines. Most of the optimal solutions of purel theromodynamic systems lie at the boundary of constrains or saturation curve. Hence the initial box of constrain chosen has to be robust enough to have decent volume of feasible solutions.

The most trusted algorithms for thermodynamic optimizations are Genetic Algorithms. It is well integrated with
[Optimization.jl](https://docs.sciml.ai/Optimization/stable/) and [Metaheuristics.jl](https://github.com/jmejia8/Metaheuristics.jl)

Usage Framework 

```
ORC(x,p) = ...
using Optimization, OptimizationMetaheuristics

x0 = [...]
p = [...]
f = OptimizationFunction(ORC)
prob = Optimization.OptimizationProblem(f, x0, p, lb = [...], ub = [...])
sol = solve(prob, PSO(), maxiters = 100000, maxtime = 100.0)
```
The thing to note here is that the `lb` and `ub` give to the optimizer has to be in limits such that the interanal functions of CoolProp are computable. 

A working example is in [Examples\OptimizationExample.jl](https://github.com/Sush1090/CoolPropCycles.jl/blob/main/Examples/OptimizationExample.jl)
