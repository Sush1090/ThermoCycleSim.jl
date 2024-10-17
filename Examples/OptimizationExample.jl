using CoolPropCycles, ModelingToolkit, DifferentialEquations, CoolProp

global fluid = "R245fa"
function ORC(x,p)

    @independent_variables t
  
    _system = Isentropic_η(η = 0.9,πc = x[2]) 
    start_T = x[1]; # Temperature at source 
    start_p = PropsSI("P","Q",0,"T",start_T,fluid) + 1e3 
    ΔT_subcool = PropsSI("T","P",start_p,"Q",0,fluid) - start_T;
    @assert ΔT_subcool > 1e-3   
    start_h = PropsSI("H","T",start_T,"P",start_p,fluid); start_mdot = 0.2 #kg/s

    @named source = MassSource(source_enthalpy = start_h,source_pressure = start_p,source_mdot = start_mdot,fluid = fluid)
    @named comp = Compressor(_system, fluid =fluid)
    @named evap = SimpleEvaporator(Δp = [0,0,0],ΔT_sh = x[3],fluid = fluid)
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

    #exp_phase = PhaseSI("H",sol[exp.h_out][1],"P",sol[exp.p_out][1],fluid)
    

    #@show
     η = (sol[exp.P][1] + sol[comp.P][1])/sol[evap.P][1] 

    penalty1 = 0     # for outlet to be gas phase
    # if exp_phase != "gas"
    #     penalty1 = abs(η)
    # end

    if sol[exp.T_out][1] <= sol[cond.T_sat][1] + 3
        penalty1 = abs(η)
    end
    @assert p[1]>p[2]
    penalty2 = 0 # evaporator outlte temperature is desired to be less than 380
    if (sol[evap.T_out][1] > p[1])
        penalty2 = abs(η)
    end

     penalty3 = 0
    if (sol[evap.T_sat][1] > p[2])
        penalty3 = abs(η)
    end
    cost = η + penalty1 + penalty2  +penalty3

    #@show cost
    Compute_cycle_error(sol,systems)

    return cost
end



"""
NOTE:  Generally the themodynamically optimial solutions lie at the boundary of the of phase changes. 
Hence choose the parameters box of variables carefully. 
The major constraint that are faced:
1) The limits on CoolProp computations.
2) As optimial solutions many times lie close to the boundary of phase changes, CoolProp assertion might start to fail.

hence Choose the box of variables carefully. 
It is recommended to use Genetic Algorithms instead of Line search Algorithms.
"""

using Optimization, OptimizationMetaheuristics

#[start_T,πc,T_sh]
x0 = [300,5,100]
p = [365,355]

p_max_start = PropsSI("P","Q",1,"T",300,fluid) + 1e3
p_crit = PropsSI("PCRIT",fluid)
πc_min = p_crit/p_max_start

f = OptimizationFunction(ORC)
# prob = Optimization.OptimizationProblem(f, x0, p, lb = [290, 1.1,2], ub = [300, 10,100])
# sol = solve(prob, PSO(), maxiters = 100000, maxtime = 100.0)

prob = Optimization.OptimizationProblem(f, x0, p, lb = [280, 1.1,2], ub = [300, πc_min,100])
sol = solve(prob, DE(), maxiters = 50, maxtime = 1000.0)

if (f(sol.u,p) >=0)
    @warn "No Feasible point found"
end