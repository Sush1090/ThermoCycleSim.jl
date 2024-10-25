using ThermoCycleSim, ModelingToolkit, DifferentialEquations, CoolProp


@independent_variables t
fluid = "R134A"
@load_fluid "R134A"
_system = Isentropic_η(η =0.5,πc =2.8) # fix the isentropic Efficiency of compressor and pressre ratio

start_T =     290; # Temperature at source 
start_p = PropsSI("P","Q",0,"T",start_T,fluid) + 1e5 # pressure at source.
# As it is ORC the inlet state is liquid and bit away from saturation curv. Hence 1e3Pa of pressure is added
ΔT_subcool = PropsSI("T","P",start_p,"Q",0,fluid) - start_T; # ensure the subcoolin temperature to reach bck to starting state.
@assert ΔT_subcool > 1e-3 # stay away from saturaton curve to aviod coolprop assertion
start_h = PropsSI("H","T",start_T,"P",start_p,fluid); start_mdot = 0.2 #kg/s


@named source = MassSource(source_enthalpy = start_h,source_pressure = start_p,source_mdot = start_mdot,)
@named comp = Compressor(_system, fluid =fluid)
@named evap = SimpleEvaporator(Δp = [0,0,0],ΔT_sh = 2,)
@named exp = Expander(_system,)
@named cond = SimpleCondensor(ΔT_sc = ΔT_subcool,Δp = [0,0,0])
@named sink = MassSink()

# Define equations
eqs = [
    connect(source.port,comp.inport)
    connect(comp.outport,evap.inport)
    connect(evap.outport,exp.inport)
    connect(exp.outport,cond.inport)
    connect(cond.outport,sink.port)
]
systems=[source,comp,evap,exp,cond,sink] # Define system

@named dis_test = ODESystem(eqs, t, systems=systems)

u0 = []
tspan = (0.0, 100.0)
sys = structural_simplify(dis_test)
prob = ODEProblem(sys,u0,tspan)
sol = solve(prob)



#compute Efficiency of the cycle.
#Note: sign convetion: Power supplied to the system is +ve while from thee system is -ve
@show η = (sol[exp.P] .+ sol[comp.P])./sol[evap.P]

#Check if the final state is close to the inital state. 
Compute_cycle_error(sol,systems)
ThermoCycleSim.PhasePlot(PhasePlotType_PH(),sol,systems,fluid)