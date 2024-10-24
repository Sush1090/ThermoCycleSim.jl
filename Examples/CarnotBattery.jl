

using CoolPropCycles, ModelingToolkit, DifferentialEquations, CoolProp




function HP(x,p)
    @independent_variables t
    fluid = "R134A"
    _system = Isentropic_η(η = 0.59,πc = x[1]) # fix the isentropic Efficiency of compressor and pressre ratio
    valve_system  = IsenthalpicExpansionValve(x[1])
    start_T = 273; # Temperature at source 
    start_p = PropsSI("P","Q",1,"T",start_T,fluid) - 1e2 # pressure at source. For HP we need gas at source
    @assert PhaseSI("T",start_T,"P",start_p,fluid) == "gas"
    ΔT_superheat = 10#start_T - PropsSI("T","P",start_p,"Q",0,fluid) ; # ensure the superheat temperature to reach bck to starting state.
    start_h = PropsSI("H","T",start_T,"P",start_p,fluid); start_mdot = 0.2 #kg/s


    @named source = MassSource(source_enthalpy = start_h,source_pressure = start_p,source_mdot = start_mdot,fluid = fluid)
    @named comp = Compressor(_system, fluid =fluid)
    @named cond = SimpleCondensor(ΔT_sc = 2,Δp = [0,0,0],fluid = fluid)
    @named exp = Valve(valve_system,fluid= fluid)
    @named evap = SimpleEvaporator(ΔT_sh = ΔT_superheat,Δp = [0,0,0],fluid = fluid)
    @named sink = MassSink(fluid = fluid)

    eqs = [
        connect(source.port,comp.inport)
        connect(comp.outport,cond.inport)
        connect(cond.outport,exp.inport)
        connect(exp.outport,evap.inport)
        connect(evap.outport,sink.port)
    ]
    systems=[source,comp,cond,exp,evap,sink] # Define system

    @named HP = ODESystem(eqs, t, systems=systems)

    u0 = []
    tspan = (0.0, 10.0)
    sys = structural_simplify(HP)
    prob = ODEProblem(sys,u0,tspan,guesses = [])
    sol = solve(prob)

    COP = sol[cond.P][1]/sol[comp.P][1]

    if COP > 0
        COP = 0
    end
    if sol[cond.T_sat][1] <  p[1] # minimum temperature needed. Else might have ∞*0 case for CB
        COP = 0
    end

    return COP,sol[cond.T_sat][1]

end

function ORC(x,p)
    @independent_variables t
    fluid = "R601A"
    _system = Isentropic_η(η =0.75,πc =x[1]) # fix the isentropic Efficiency of compressor and pressre ratio

    start_T =     273; # Temperature at source 
    start_p = PropsSI("P","Q",0,"T",start_T,fluid) + 1e3 # pressure at source.
    # As it is ORC the inlet state is liquid and bit away from saturation curv. Hence 1e3Pa of pressure is added
    ΔT_subcool = PropsSI("T","P",start_p,"Q",0,fluid) - start_T; # ensure the subcoolin temperature to reach bck to starting state.
    @assert ΔT_subcool > 1e-3 # stay away from saturaton curve to aviod coolprop assertion
    start_h = PropsSI("H","T",start_T,"P",start_p,fluid); start_mdot = 0.2 #kg/s


    @named source = MassSource(source_enthalpy = start_h,source_pressure = start_p,source_mdot = start_mdot,fluid = fluid)
    @named comp = Compressor(_system, fluid =fluid)
    @named evap = SimpleEvaporator(Δp = [0,0,0],ΔT_sh = x[2],fluid = fluid)
    @named exp = Expander(_system,fluid= fluid)
    @named cond = SimpleCondensor(ΔT_sc = ΔT_subcool,Δp = [0,0,0],fluid = fluid)
    @named sink = MassSink(fluid = fluid)

    eqs = [
        connect(source.port,comp.inport)
        connect(comp.outport,evap.inport)
        connect(evap.outport,exp.inport)
        connect(exp.outport,cond.inport)
        connect(cond.outport,sink.port)
    ]
    systems=[source,comp,evap,exp,cond,sink] # Define system

    @named orc = ODESystem(eqs, t, systems=systems)

    u0 = []
    tspan = (0.0, 100.0)
    sys = structural_simplify(orc)
    prob = ODEProblem(sys,u0,tspan,guesses = [])
    sol = solve(prob)
    
    η = (sol[exp.P][1] .+ sol[comp.P][1])./sol[evap.P][1]
    
    if sol[evap.T_out][1] > p[1]
        η = 0;
    end
    if η > 0
        η = 0
    end
    return η
end


function CB(x,p_hp)
    x1 = [x[1]];x2 = x[2:end]
    COP,T_sat_hp = HP(x1,p_hp)
   
    p = [T_sat_hp]
    η = ORC(x2,p)
    @show -COP*η
    return -COP*η
    
end
# hp_pr,orc_pr,orc_suph
# x0 = [ 2,2.2,2.0]
# pp = [330]



# using Optimization, OptimizationMetaheuristics
# f = OptimizationFunction(CB)
# prob = Optimization.OptimizationProblem(f, x0, pp, lb = [1.1, 2,2], ub = [9.5, 15,100])
# sol = solve(prob, DE(), maxiters = 100, maxtime = 6000.0)