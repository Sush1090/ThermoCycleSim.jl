using CarnotCycles, ModelingToolkit, DifferentialEquations

@independent_variables t
@load_fluid "Water"

@named src = MassSource(:liquid,Δp_super = 1e5, source_temperature = 300,source_mdot = 1)
@named comp = Compressor()
@named sink = MassSink()
eqs = [
    connect(src.port,comp.inport)
    connect(comp.outport,sink.port)
]

system=[src,comp,sink]
@named model = ODESystem(eqs,t,systems= system)
sys = structural_simplify(model)
para = [comp.πc => 3,comp.η => 1]
u0 = []
tspan = (0,1)

prob = ODEProblem(sys,u0,tspan,para)
sol = solve(prob)