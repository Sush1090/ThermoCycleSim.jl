using CarnotCycles, ModelingToolkit, DifferentialEquations, CoolProp


@independent_variables t
fluid = "Argon"
@load_fluid "Argon"
isentropic_system = Isentropic_η(η =1.0,πc =10)
isothermal_system = Isothermal_comp(πc = 5)
start_T = 450; start_p = 2*101325;mdot = 0.2
init_state = initialize_state(T_start = start_T,p_start = start_p,mdot = mdot)

@named source = MassSource(init_state)
@named isothermal_exp = Expander(isothermal_system)
@named isentropic_exp = Expander(isentropic_system)
@named isothermal_comp = Compressor(isothermal_system)
@named isentropic_comp = Compressor(isentropic_system)
@named sink = MassSink()

eqs = [
    connect(source.port,isothermal_exp.inport)
    connect(isothermal_exp.outport,isentropic_exp.inport)
    connect(isentropic_exp.outport,isothermal_comp.inport)
    connect(isothermal_comp.outport,isentropic_comp.inport)
    connect(isentropic_comp.outport,sink.port)
]

systems = [source,isothermal_exp,isentropic_exp,isothermal_comp,isentropic_comp,sink]
@named CarnotCycle = ODESystem(eqs,t,systems=systems)
u0 = []
tspan = (0.0, 10.0)
sys = structural_simplify(CarnotCycle)
prob = ODEProblem(sys,u0,tspan)
sol = solve(prob)

Compute_cycle_error(sol,systems)
CarnotCycles.CyclePlot(PhasePlotType_TS(),sol,systems)