module CoolPropCycles

# Write your package code here.

using ModelingToolkit, CoolProp#, DifferentialEquations

include("Utils.jl")
include("Processes.jl")
include("ThemodynamicClaculations.jl")

include("Components/Compressor.jl")
include("Components/HeatExchangers.jl")

end
