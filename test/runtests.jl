using CoolPropCycles, ModelingToolkit, DifferentialEquations, CoolProp
using Test

@testset "Isentropic Process" begin
    # Write your tests here.

    fluid = "R134A"
    _system = Isentropic_η(η = 1,πc = 5)
@independent_variables t
start_T = 300;
start_p = PropsSI("P","Q",0,"T",start_T,fluid) + 1e3
ΔT_subcool = PropsSI("T","P",start_p,"Q",0,fluid) - start_T;
@assert ΔT_subcool > 1e-3
start_h = PropsSI("H","T",start_T,"P",start_p,fluid); start_mdot = 0.2 #kg/s


@named source = MassSource(source_enthalpy = start_h,source_pressure = start_p,source_mdot = start_mdot,fluid = fluid)
@named comp = Compressor(_system, fluid =fluid)
@named exp = Expander(_system,fluid= fluid)
@named sink = MassSink(fluid = fluid)

eqs = [
    connect(source.port,comp.inport)
    connect(comp.outport,exp.inport)
    connect(exp.outport,sink.port)
]
systems = [source,comp,exp,sink]
@named test_isentropic = ODESystem(eqs, t, systems=systems)
u0 = []
tspan = (0.0, 1.0)
sys = structural_simplify(test_isentropic)
prob = ODEProblem(sys,u0,tspan,guesses = [])
sol = solve(prob)
bool1 = isapprox(sol[comp.s_in][1],sol[comp.s_out][1])
bool2 = isapprox(sol[exp.s_in][1],sol[exp.s_out][1])
@test bool1 == true
@test bool2 == true
@test isapprox(sol[source.p][1],sol[sink.p][1])
@test isapprox(sol[source.h][1],sol[sink.h][1])
end


@testset "Evaporator phase" begin
    out_phase = "gas"

fluid = "R134A"
@independent_variables t
start_T = 300;
start_p = PropsSI("P","Q",0,"T",start_T,fluid) + 1e3
start_h = PropsSI("H","T",start_T,"P",start_p,fluid); start_mdot = 0.2 #kg/s

@named source = MassSource(source_enthalpy = start_h,source_pressure = start_p,source_mdot = start_mdot,fluid = fluid)
@named evporator = SimpleEvaporator(ΔT_sh = 5,fluid = fluid)
@named sink = MassSink(fluid = fluid)

eqs = [
    connect(source.port,evporator.inport)
    connect(evporator.outport,sink.port)
]
systems = [source,evporator,sink]
@named test_evap = ODESystem(eqs, t, systems=systems)
u0 = []
tspan = (0.0, 1.0)
sys = structural_simplify(test_evap)
prob = ODEProblem(sys,u0,tspan,guesses = [])
sol = solve(prob)

phase_test = PhaseSI("T",sol[evporator.T_out][1],"P",sol[evporator.p_out][1],fluid)

@test phase_test == out_phase
end
#@testset 