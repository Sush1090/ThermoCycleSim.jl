

using ThermoCycleSim, ModelingToolkit, DifferentialEquations, CoolProp


@independent_variables t
@load_fluid "R134A"
_system = Isentropic_η(η = 1,πc = 3.35) # fix the isentropic Efficiency of compressor and pressre ratio
valve_system  = IsenthalpicExpansionValve(3.35)
start_T = 273.15; # Temperature at source 
    start_p = PropsSI("P","Q",1,"T",start_T-10,fluid) - 1e2 #2*101325 #PropsSI("P","Q",1,"T",start_T-10,fluid) - 1e2 # pressure at source. For HP we need gas at source
    @assert PhaseSI("T",start_T,"P",start_p,"R134A") == "gas"
    ΔT_superheat = start_T - PropsSI("T","P",start_p,"Q",0,"R134A") ; # ensure the superheat temperature to reach bck to starting state.
    start_h = PropsSI("H","T",start_T,"P",start_p,"R134A"); start_mdot = 0.2 #kg/s


    @named source = MassSource(source_enthalpy = start_h,source_pressure = start_p,source_mdot = start_mdot)
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

    u0 = []
    tspan = (0.0, 10.0)
    sys = structural_simplify(hp)
    prob = ODEProblem(sys,u0,tspan)
    sol = solve(prob)

    @show COP = sol[cond.P][1]/sol[comp.P][1]
#Check if the final state is close to the inital state. 
Compute_cycle_error(sol,systems)
#ThermoCycleSim.PhasePlot(PhasePlotType_TS(),sol,systems)
ThermoCycleSim.PhasePlot(PhasePlotType_TS(),sol,systems)