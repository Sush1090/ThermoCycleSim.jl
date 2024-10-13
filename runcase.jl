using CoolPropCycles, ModelingToolkit, DifferentialEquations, CoolProp


@independent_variables t
fluid = "Argon"
comp_system = Isentropic_η(η = 1.0,πc = 5)

start_p = 101325; start_T = 300;
start_h = PropsSI("H","T",start_T,"P",start_p,fluid); start_mdot = 2 #kg/s


@named source = MassSource(source_enthalpy = start_h,source_pressure = start_p,source_mdot = start_mdot,fluid = fluid)
@named comp = Compressor(comp_system, fluid =fluid)
@named exp = Expander(comp_system,fluid= fluid)
@named sink = MassSink(fluid = fluid)

eqs = [
    connect(source.port,comp.inport)
    connect(comp.outport,exp.inport)
    connect(exp.outport,sink.port)
]

@named dis_test = ODESystem(eqs, t, systems=[source,comp,exp,sink])

u0 = []
tspan = (0.0, 10.0)
sys = structural_simplify(dis_test)
prob = ODEProblem(sys,u0,tspan,guesses = [])
sol = solve(prob)