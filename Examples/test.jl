using CarnotCycles, ModelingToolkit, DifferentialEquations, Clapeyron

@independent_variables t
model = cPR(["ethane","methane"],idealmodel = ReidIdeal);
load_fluid(model)

@named src = MassSource()
@named sink = MassSink()
eqs = [
    connect(src.port,sink.port)
]

system=[src,sink]
@named cycle = ODESystem(eqs,t,systems= system)
sys = structural_simplify(cycle)
para = [sys.src.z=>[1.0,1.0]]
u0 = []

prob = SteadyStateProblem(sys,u0,para)
sol = solve(prob)