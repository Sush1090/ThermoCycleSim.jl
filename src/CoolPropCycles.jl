module CoolPropCycles

# Write your package code here.

using ModelingToolkit, CoolProp, DifferentialEquations, Plots

include("Utils.jl")
include("Processes.jl")
include("ThemodynamicClaculations.jl")

include("Components/compressor.jl")
include("Components/expander.jl")
include("Components/HeatExchangers.jl")
include("Components/valve.jl")

include("Plotting/PhasePlot.jl")
end
