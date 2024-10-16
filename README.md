# CoolPropCycles.jl

[![Build Status](https://github.com/Sush1090/CoolPropCycles.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/Sush1090/CoolPropCycles.jl/actions/workflows/CI.yml?query=branch%3Amain)


This package combines [ModelingToolkit.jl](https://github.com/SciML/ModelingToolkit.jl) "Acausal Modeling" with [CoolProp.jl](https://github.com/CoolProp/CoolProp.jl) to make thermodynamic cycles:

    1. Organic Rakine Cycle
    2. Vapour Compression Cycle
    3. Brayton Cycle

It is for steady state calculations, hence the evaporators and compressors are modeled using energy conservation. The principle is state-point manipulations.  The system works better for subcritical cycle parameters

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
    6. Heat Source - Temperature independent

 <!-- It also provides basic functions that find the pressure to match the pitch points.  -->
## Installation

## Basic Usage
Every default component has the following variables : `P(t)`,`s_in(t)`,`p_in(t)`,`T_in(t)`,`h_in(t)`,`ρ_in(t)`,`s_out(t)`,`p_out(t)`,`T_out(t)`,`h_out(t)`,`ρ_out(t)` and two ports `inport` and `outport`. 
Follow the [ModelingToolkit.jl](https://github.com/SciML/ModelingToolkit.jl) documentation for further understanding variable handling. 

Note: sign convetion - Power supplied to the system is +ve while from thee system is -ve

## Examples 
### Organic Rankine Cycle
Example of Organic Rankine Cycle using R134A

```
using CoolPropCycles, ModelingToolkit, DifferentialEquations, CoolProp, Plots


@independent_variables t
fluid = "R134A"
_system = Isentropic_η(η = 0.8,πc = 5) # fix the isentropic Efficiency of compressor and pressre ratio

start_T =    290; # Temperature at source 
start_p = PropsSI("P","Q",0,"T",start_T,fluid) + 1e3 # pressure at source.
# As it is ORC the inlet state is liquid and bit away from saturation curv. Hence 1e3Pa of pressure is added
ΔT_subcool = PropsSI("T","P",start_p,"Q",0,fluid) - start_T; # ensure the subcoolin temperature to reach bck to starting state.
@assert ΔT_subcool > 1e-3 # stay away from saturaton curve to aviod coolprop assertion
start_h = PropsSI("H","T",start_T,"P",start_p,fluid); start_mdot = 0.2 #kg/s


@named source = MassSource(source_enthalpy = start_h,source_pressure = start_p,source_mdot = start_mdot,fluid = fluid)
@named comp = Compressor(_system, fluid =fluid)
@named evap = SimpleEvaporator(Δp = [0,0,0],ΔT_sh = 19.4058339236015846,fluid = fluid)
@named exp = Expander(_system,fluid= fluid)
@named cond = SimpleCondensor(ΔT_sc = ΔT_subcool,Δp = [0,0,0],fluid = fluid)
@named sink = MassSink(fluid = fluid)

# Define equations
eqs = [
    connect(source.port,comp.inport)
    connect(comp.outport,evap.inport)
    connect(evap.outport,exp.inport)
    connect(exp.outport,cond.inport)
    connect(cond.outport,sink.port)
]
systems=[source,comp,evap,exp,cond,sink] # Define system

@named dis_test = ODESystem(eqs, t, systems=systems)

u0 = []
tspan = (0.0, 10.0)
sys = structural_simplify(dis_test)
prob = ODEProblem(sys,u0,tspan,guesses = [])
sol = solve(prob)


#compute Efficiency of the cycle.
#Note: sign convetion: Power supplied to the system is +ve while from thee system is -ve
@show η = (sol[exp.P][1] + sol[comp.P][1])/sol[evap.P][1]

#Check if the final state is close to the inital state. 
Compute_cycle_error(sol,systems)
```


## Phase Diagram Plotting

As of now the direct plotting of T-S phase digram of the cycle is provided. It requires the `system` to have the first as `MassSource(...)` and last variable as `MassSink(...)`.   
```
CoolPropCycles.PhasePlot(PhasePlotType_TS(),sol,systems,fluid)
```

Insert Diagram

## Thermodynamic Cycle Optimization
The cycles created can be wrapped with functions and sent to optimization routines. Most of the optimal solutions of purely theromodynamic systems lie at the boundary of constrains or saturation curve. Hence the initial box of constrain chosen has to be robust enough to have decent volume of feasible solutions.

The most trusted algorithms for thermodynamic optimizations are Genetic Algorithms. It is well integrated with
[Optimization.jl](https://docs.sciml.ai/Optimization/stable/) and [Metaheuristics.jl](https://github.com/jmejia8/Metaheuristics.jl)

Note: The heuristic approach will not always converge. But will return the best answer in the possible iterations.   

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
The thing to note here is that the `lb` and `ub` given to the optimizer has to be in limits such that the internal functions of CoolProp are computable. 

A working example is in [Examples\OptimizationExample.jl](https://github.com/Sush1090/CoolPropCycles.jl/blob/main/Examples/OptimizationExample.jl)
