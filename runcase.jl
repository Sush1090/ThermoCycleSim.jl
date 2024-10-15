using CoolPropCycles, ModelingToolkit, DifferentialEquations, CoolProp


@independent_variables t
fluid = "R134A"
_system = Isentropic_η(η = 1,πc = 5)

start_T = 300;
start_p = PropsSI("P","Q",0,"T",start_T,fluid) + 1e3
ΔT_subcool = PropsSI("T","P",start_p,"Q",0,fluid) - start_T;
@assert ΔT_subcool > 1e-3
start_h = PropsSI("H","T",start_T,"P",start_p,fluid); start_mdot = 0.2 #kg/s


@named source = MassSource(source_enthalpy = start_h,source_pressure = start_p,source_mdot = start_mdot,fluid = fluid)
@named comp = Compressor(_system, fluid =fluid)
@named evap = SimpleEvaporator(Δp = 0,ΔT_sh = 5,fluid = fluid)
@named exp = Expander(_system,fluid= fluid)
@named cond = SimpleCondensor(ΔT_sc = ΔT_subcool,Δp = 0,fluid = fluid)
@named sink = MassSink(fluid = fluid)

eqs = [
    connect(source.port,comp.inport)
    connect(comp.outport,evap.inport)
    connect(evap.outport,exp.inport)
    connect(exp.outport,cond.inport)
    connect(cond.outport,sink.port)
]
systems=[source,comp,evap,exp,cond,sink]
@named dis_test = ODESystem(eqs, t, systems=systems)

u0 = []
tspan = (0.0, 10.0)
sys = structural_simplify(dis_test)
prob = ODEProblem(sys,u0,tspan,guesses = [])
sol = solve(prob)

@show η = (sol[exp.P][1] + sol[comp.P][1])/sol[evap.P][1]

Compute_cycle_error(sol,systems)