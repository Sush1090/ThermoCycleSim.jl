module CoolPropCycles

# Write your package code here.

using ModelingToolkit, CoolProp#, DifferentialEquations

include("Utils.jl")
include("Processes.jl")
include("ThemodynamicClaculations.jl")

include("Components/compressor.jl")
include("Components/expander.jl")
include("Components/HeatExchangers.jl")

end
