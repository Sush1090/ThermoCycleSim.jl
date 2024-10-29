using CarnotCycles, ModelingToolkit, DifferentialEquations, CoolProp



fluid = "Argon"
@load_fluid "Argon"
_system = Isentropic_η(η = 1,πc = 5)
_isochoric = Isochoric_comp(πc = 5)
_isothermal = Isothermal_comp(πc =5)
@independent_variables t
start_T = 300;
start_p = 101325
start_h = PropsSI("H","T",start_T,"P",start_p,fluid);start_mdot = 0.2

@named source = MassSource(source_enthalpy = start_h,source_pressure = start_p,source_mdot = start_mdot)
@named comp_isen = Compressor(_system)
@named exp_isen = Expander(_system)
@named comp_ther = Compressor(_isothermal)
@named exp_ther = Expander(_isothermal)
@named comp_chor = Compressor(_isochoric)
@named exp_chor = Expander(_isochoric)
@named sink = MassSink()

eqs = [
    connect(source.port,comp_isen.inport)
    connect(comp_isen.outport,exp_isen.inport)
    connect(exp_isen.outport,comp_ther.inport)
    connect(comp_ther.outport,exp_ther.inport)
    connect(exp_ther.outport,comp_chor.inport)
    connect(comp_chor.outport,exp_chor.inport)
    connect(exp_chor.outport,sink.port)
]
systems = [source,comp_isen,exp_isen,comp_ther,exp_ther,comp_chor,exp_chor,sink]
@named test_processes = ODESystem(eqs, t, systems=systems)
u0 = []
tspan = (0.0, 1.0)
sys = structural_simplify(test_processes)
prob = ODEProblem(sys,u0,tspan)
sol = solve(prob)

fig1 = CarnotCycles.CyclePlot(PhasePlotType_TS(),sol,systems)
fig2 = CarnotCycles.CyclePlot(PhasePlotType_PH(),sol,systems)