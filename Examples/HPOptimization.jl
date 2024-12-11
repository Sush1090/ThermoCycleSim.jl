using CarnotCycles, ModelingToolkit, DifferentialEquations, CoolProp



function HP(x,p)
    @independent_variables t
    fluid = "R1233zdE"
    @load_fluid "R1233zdE"
    _system = Isentropic_η(η = 0.55,πc = x[1]) # fix the isentropic Efficiency of compressor and pressre ratio
    valve_system  = IsenthalpicExpansionValve(x[1])    
    start_T = 273.15+20; # Temperature at source 
    start_p = PropsSI("P","Q",1,"T",start_T-10,fluid) - 1e2 #2*101325 #PropsSI("P","Q",1,"T",start_T-10,fluid) - 1e2 # pressure at source. For HP we need gas at source
    p_crit = PropsSI("PCRIT",fluid)
    # if (πc*start_p >= p_crit)
    #     throw(error("Pressure is in super-critical zone. Put parameters for sub-critical cycle"))
    # end
    @assert PhaseSI("T",start_T,"P",start_p,fluid) == "gas"
    ΔT_superheat = start_T - PropsSI("T","P",start_p,"Q",0,fluid) ; # ensure the superheat temperature to reach bck to starting state.
    start_mdot = 0.034 #kg/s

    state_src = initialize_state(T_start = start_T,p_start=start_p,mdot = start_mdot)

    @named source = MassSource(state_src)
    @named comp = Compressor(_system)
    @named cond = SimpleCondensor(ΔT_sc = x[2],Δp = [0,0,0])
    @named expander = Valve(valve_system)
    @named evap = SimpleEvaporator(ΔT_sh = ΔT_superheat,Δp = [0,0,0])
    @named sink = MassSink()

    eqs = [
        connect(source.port,comp.inport)
        connect(comp.outport,cond.inport)
        connect(cond.outport,expander.inport)
        connect(expander.outport,evap.inport)
        connect(evap.outport,sink.port)
    ]
    systems=[source,comp,cond,expander,evap,sink] # Define system

    @named hp = ODESystem(eqs, t, systems=systems)

    u0 = []
    tspan = (0.0, 10.0)
    sys = structural_simplify(hp)
    prob = ODEProblem(sys,u0,tspan)
    sol = solve(prob)

   @show COP = sol[cond.P][1]/sol[comp.P][1]
    
    Q = abs(sol[cond.P][1])
    @show Q
    @show sol[comp.P][1]
    #Tsat = sol[cond.T_sat][1]
    T_out = sol[cond.T_out][1]
    pen1 = 0;pen2 = 0
    if Q < p[1]
        pen1 = abs(COP)
    end
    if T_out < p[2]
        pen2 = abs(COP)
    end
    @show cost = COP + pen1 + pen2
#Check if the final state is close to the inital state. 
    Compute_cycle_error(sol,systems)       
    return cost 
end

# x0 = [14,2]
# p = [35000,340]


#HP(x0,p)
# using Optimization, OptimizationMetaheuristics
# f = OptimizationFunction(HP)
# prob = Optimization.OptimizationProblem(f, x0, p, lb = [2.0,1], ub = [14.0,30])
# sol = solve(prob, DE(), maxiters = 500, maxtime = 1000.0,x_tol = 1e-2)