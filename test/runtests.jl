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
if _system.η == 1
    @test isapprox(sol[comp.s_in][1],sol[comp.s_out][1])
    @test isapprox(sol[exp.s_in][1],sol[exp.s_out][1])
end
@test isapprox(sol[comp.p_out][1]/sol[comp.p_in][1],_system.πc)
@test isapprox(sol[exp.p_in][1]/sol[exp.p_out][1],_system.πc)
end


@testset "Evaporator" begin
    out_phase = "gas"
    in_phase_liquid = "liquid"
    #in_phase_supercrit
    ΔT_sh = 5
fluid = "R134A"
@independent_variables t
start_T = 300;
start_p = PropsSI("P","Q",0,"T",start_T,fluid) + 1e3
start_h = PropsSI("H","T",start_T,"P",start_p,fluid); start_mdot = 0.2 #kg/s

@named source = MassSource(source_enthalpy = start_h,source_pressure = start_p,source_mdot = start_mdot,fluid = fluid)
@named evporator = SimpleEvaporator(ΔT_sh = ΔT_sh,fluid = fluid)
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

out_phase_test = PhaseSI("T",sol[evporator.T_out][1],"P",sol[evporator.p_out][1],fluid)
in_phase_test = PhaseSI("T",sol[evporator.T_in][1],"P",sol[evporator.p_in][1],fluid)

@test out_phase_test == out_phase
@test isapprox(sol[evporator.T_out][1]-ΔT_sh,sol[evporator.T_sat][1])
@test in_phase_test == in_phase_liquid
end
#@testset 

@testset "Condensor" begin
    out_phase = "liquid"
    in_phase = "gas"
    #in_phase_supercrit
    ΔT_sc = 3
fluid = "R134A"
@independent_variables t
start_T = 300;
start_p = PropsSI("P","Q",0,"T",start_T,fluid) - 101325
start_h = PropsSI("H","T",start_T,"P",start_p,fluid); start_mdot = 0.2 #kg/s

@named source = MassSource(source_enthalpy = start_h,source_pressure = start_p,source_mdot = start_mdot,fluid = fluid)
@named condensor = SimpleCondensor(ΔT_sc = ΔT_sc,fluid = fluid)
@named sink = MassSink(fluid = fluid)

eqs = [
    connect(source.port,condensor.inport)
    connect(condensor.outport,sink.port)
]
systems = [source,condensor,sink]
@named test_condensor = ODESystem(eqs, t, systems=systems)
u0 = []
tspan = (0.0, 1.0)
sys = structural_simplify(test_condensor)
prob = ODEProblem(sys,u0,tspan,guesses = [])
sol = solve(prob)

out_phase_test = PhaseSI("T",sol[condensor.T_out][1],"P",sol[condensor.p_out][1],fluid)
in_phase_test = PhaseSI("T",sol[condensor.T_in][1],"P",sol[condensor.p_in][1],fluid)

@test out_phase_test == out_phase
@test isapprox(sol[condensor.T_out][1]+ΔT_sc,sol[condensor.T_sat][1])
@test in_phase_test == in_phase

end


@testset "Valve" begin

fluid = "R134A"
@independent_variables t
start_T = 300;
start_p = PropsSI("P","Q",0,"T",start_T,fluid) 
start_h = PropsSI("H","Q",0,"P",start_p,fluid); start_mdot = 0.2 #kg/s

valve_system  = IsenthalpicExpansionValve(4.5)

@named source = MassSource(source_enthalpy = start_h,source_pressure = start_p,source_mdot = start_mdot,fluid = fluid)
@named exp = Valve(valve_system,fluid= fluid)
@named sink = MassSink(fluid = fluid)

eqs = [
    connect(source.port,exp.inport)
    connect(exp.outport,sink.port)
]
systems = [source,exp,sink]
@named test_condensor = ODESystem(eqs, t, systems=systems)
u0 = []
tspan = (0.0, 1.0)
sys = structural_simplify(test_condensor)
prob = ODEProblem(sys,u0,tspan,guesses = [])
sol = solve(prob)
@test isapprox(sol[source.h][1],sol[sink.h][1]) 
end