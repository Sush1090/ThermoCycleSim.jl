

using ThermodynamicCycleSim, ModelingToolkit, DifferentialEquations, CoolProp


@independent_variables t
@load_fluid "R134A"
_system = Isentropic_η(η = 0.8,πc = 7.5) # fix the isentropic Efficiency of compressor and pressre ratio
valve_system  = IsenthalpicExpansionValve(7.5)
start_T = 260; # Temperature at source 
start_p = PropsSI("P","Q",1,"T",start_T-10,"R134A") - 1e2 # pressure at source. For HP we need gas at source
@assert PhaseSI("T",start_T,"P",start_p,"R134A") == "gas"
ΔT_superheat = start_T - PropsSI("T","P",start_p,"Q",0,"R134A") ; # ensure the superheat temperature to reach bck to starting state.
start_mdot = 0.2 #kg/s

state_init = initialize_state(T_start = start_T,p_start=start_p,mdot = start_mdot) 

@named source = MassSource(state_init)
@named comp = Compressor(_system)
@named cond = SimpleCondensor(ΔT_sc = 1e-2,Δp = [0,0,0])
@named expander = Valve(valve_system)
@named evap = SimpleEvaporator(ΔT_sh = ΔT_superheat,Δp = [0,0,0])
@named sink = MassSink()

# Define equations
eqs = [
    connect(source.port,comp.inport)
    connect(comp.outport,cond.inport)
    connect(cond.outport,expander.inport)
    connect(expander.outport,evap.inport)
    connect(evap.outport,sink.port)
]
systems=[source,comp,cond,expander,evap,sink] # Define system

@named HP = ODESystem(eqs, t, systems=systems)

u0 = []
tspan = (0.0, 10.0)
sys = structural_simplify(HP)
prob = ODEProblem(sys,u0,tspan)
sol = solve(prob)

#compute COP
#Note: sign convetion: Power supplied to the system is +ve while from thee system is -ve
@show COP = sol[cond.P][1]/sol[comp.P][1]
#Check if the final state is close to the inital state. 
Compute_cycle_error(sol,systems)
ThermodynamicCycleSim.PhasePlot(PhasePlotType_TS(),sol,systems)