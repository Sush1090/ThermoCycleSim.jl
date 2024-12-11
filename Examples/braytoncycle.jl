using CarnotCycles, ModelingToolkit, DifferentialEquations, CoolProp


@independent_variables t
@load_fluid "Argon"
_system = Isentropic_η(η =0.7,πc =2) # fix the isentropic Efficiency of compressor and pressre ratio
_systemcomp = Isentropic_η(η =0.7,πc =2.0)
start_T =  280; # Temperature at source 
start_mdot = 0.2 #kg/s

source_h = PropsSI("H","P",101325,"T",280,"Argon")

@named source = MassSource(source_pressure = 101325,source_mdot = start_mdot,source_enthalpy = source_h)
@named comp = Compressor(_systemcomp)
@named superheater = IsobaricHeatSource(Q_dot=2e5)
@named expander = Expander(_system)
@named sink = MassSink()

# Define equations
eqs = [
    connect(source.port,comp.inport)
    connect(comp.outport,superheater.inport)
    connect(superheater.outport,expander.inport)
    connect(expander.outport,sink.port)
]
system=[source,comp,superheater,expander,sink] # Define system

@named dis_test = ODESystem(eqs, t, systems=system)

u0 = []
tspan = (0.0, 100.0)
sys = structural_simplify(dis_test)
prob = ODEProblem(sys,u0,tspan)
sol = solve(prob)


#compute Efficiency of the cycle.
#Note: sign convetion: Power supplied to the system is +ve while from thee system is -ve
@show η = (sol[expander.P].+ sol[comp.P])./sol[superheater.P]

#Check if the final state is close to the inital state. 
Compute_cycle_error(sol,system,reltol = 1e-2)
CarnotCycles.CyclePlot(PhasePlotType_PH(),sol,system)