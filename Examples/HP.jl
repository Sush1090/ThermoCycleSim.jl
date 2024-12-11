

using CarnotCycles, ModelingToolkit, DifferentialEquations, CoolProp


@independent_variables t
#fluid = "R1234YF"
@load_fluid "Propane"
πc = 1.8
_system = Isentropic_η(η = 0.7,πc = πc) # fix the isentropic Efficiency of compressor and pressre ratio
valve_system  = IsenthalpicExpansionValve(πc)

    start_T = 273+60; # Temperature at source 
    # start_p = PropsSI("P","Q",1,"T",start_T-10,fluid) - 1e2
    # p_crit = PropsSI("PCRIT",fluid)
    # if (πc*start_p >= p_crit)
    #     throw(error("Pressure is in super-critical zone. Put parameters for sub-critical cycle"))
    # end
    # @assert PhaseSI("T",start_T,"P",start_p,fluid) == "gas"
    ΔT_superheat = 6.713304480203021 #start_T - PropsSI("T","P",start_p,"Q",0,fluid) ; # ensure the superheat temperature to reach bck to starting state.
    start_mdot = 0.0407 #kg/s

    #state_src = initialize_state(T_start = start_T,p_start=start_p,mdot = start_mdot)

    @named source = MassSource(:gas,Δp_sub = 2e5,source_mdot = start_mdot,source_temperature = start_T)
    @named comp = Compressor(_system)
    @named cond = SimpleCondensor(ΔT_sc = 2,Δp = [0,0,0])
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
    tspan = (0,10)
    u0 = []
    sys = structural_simplify(hp)
    prob = ODEProblem(sys,u0,tspan)
    sol = solve(prob)

    @show COP = sol[cond.P][1]/sol[comp.P][1]
#Check if the final state is close to the inital state. 
Compute_cycle_error(sol,systems,reltol = 1e-3)
#CarnotCycles.PhasePlot(PhasePlotType_TS(),sol,systems)
fig1 = CarnotCycles.PhasePlot(PhasePlotType_TS(),sol,systems)
fig2 = CarnotCycles.PhasePlot(PhasePlotType_PH(),sol,systems)
@show sol[evap.T_sat]
@show sol[cond.T_sat]