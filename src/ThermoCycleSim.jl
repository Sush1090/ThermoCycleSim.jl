module ThermoCycleSim

# Write your package code here.

using ModelingToolkit, CoolProp, DifferentialEquations, Plots
using QML
include("Utils.jl")
include("Processes.jl")
include("ThemodynamicClaculations.jl")

include("Components/compressor.jl")
include("Components/expander.jl")
include("Components/HeatExchangers.jl")
include("Components/valve.jl")
include("Components/recuperator.jl")

include("Plotting/PhasePlot.jl")
end
