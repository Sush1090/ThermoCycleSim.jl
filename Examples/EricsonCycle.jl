using CoolPropCycles, ModelingToolkit, DifferentialEquations, CoolProp


@independent_variables t
fluid = "Argon"
_system = Isothermal_Δp(20) # fix the isentropic Efficiency of compressor and pressre ratio

start_T = 300; # Temperature at source 
start_p = 101325# pressure at source.
start_h = PropsSI("H","T",start_T,"P",start_p,fluid); start_mdot = 0.2 #kg/s

@named source = MassSource(source_enthalpy = start_h,source_pressure = start_p,source_mdot = start_mdot,fluid = fluid)
@named comp = Compressor2(_system, fluid =fluid)
@named heatsrc = IsobaricHeatSource(Q_dot=1e5,fluid=fluid)
@named exp = Expander2(_system,fluid= fluid)
@named heatsink = IsobaricHeatSink(Q_dot=-1e5+867.478,fluid=fluid)
@named sink = MassSink(fluid = fluid)


# Define equations
eqs = [
    connect(source.port,comp.inport)
    connect(comp.outport,heatsrc.inport)
    connect(heatsrc.outport,exp.inport)
    connect(exp.outport,heatsink.inport)
    connect(heatsink.outport,sink.port)
]
systems=[source,comp,heatsrc,exp,heatsink,sink] # Define system

@named dis_test = ODESystem(eqs, t, systems=systems)

u0 = []
tspan = (0.0, 100.0)
sys = structural_simplify(dis_test)
prob = ODEProblem(sys,u0,tspan,guesses = [])
sol = solve(prob)

@show η = (sol[exp.P] .+ sol[comp.P])./sol[heatsrc.P]

#Check if the final state is close to the inital state. 
Compute_cycle_error(sol,systems)

propx =:T_in
propy =:s_in
propp =:p_in
proph =:h_in
TT = [sol[getproperty(i, propx)][1] for i in plot_sys]
ss = [sol[getproperty(i, propy)][1] for i in plot_sys]
pp = [sol[getproperty(i, propp)][1] for i in plot_sys]
hh = [sol[getproperty(i, proph)][1] for i in plot_sys]
scatter(ss,TT)
