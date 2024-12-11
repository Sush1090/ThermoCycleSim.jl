using CarnotCycles, ModelingToolkit, DifferentialEquations, CoolProp


@independent_variables t
@load_fluid "R134A"
_system = Isentropic_η(η =0.7,πc =10) # fix the isentropic Efficiency of compressor and pressre ratio
_systempump = Isentropic_η(η =0.7,πc =10)
_systemcomp = Isentropic_η(η =0.7,πc =1.0)
start_T =  280; # Temperature at source 
start_mdot = 0.2 #kg/s


@named source = MassSource(:liquid,Δp_super = 1e1,source_mdot = start_mdot,source_temperature = start_T)
@named pump = Compressor(_systempump)
@named evap = SimpleEvaporator(Δp = [0,0,0],ΔT_sh = 10)
@named comp = Compressor(_systemcomp)
@named superheater = IsobaricHeatSource(Q_dot=1e4)
@named expander = Expander(_system)
@named cond = SimpleCondensor(ΔT_sc = 1e-2,Δp = [0,0,0])
@named sink = MassSink()

# Define equations
eqs = [
    connect(source.port,pump.inport)
    connect(pump.outport,evap.inport)
    connect(evap.outport,comp.inport)
    connect(comp.outport,superheater.inport)
    connect(superheater.outport,expander.inport)
    connect(expander.outport,cond.inport)
    connect(cond.outport,sink.port)
]
system=[source,pump,evap,comp,superheater,expander,cond,sink] # Define system

@named dis_test = ODESystem(eqs, t, systems=system)

u0 = []
tspan = (0.0, 100.0)
sys = structural_simplify(dis_test)
prob = ODEProblem(sys,u0,tspan)
sol = solve(prob)
@show sol[evap.T_sat][1]


#compute Efficiency of the cycle.
#Note: sign convetion: Power supplied to the system is +ve while from thee system is -ve
@show η = (sol[expander.P] .+ sol[pump.P] .+ sol[comp.P])./sol[evap.P]

#Check if the final state is close to the inital state. 
Compute_cycle_error(sol,system,reltol = 1e-2)
CarnotCycles.PhasePlot(PhasePlotType_TS(),sol,system)