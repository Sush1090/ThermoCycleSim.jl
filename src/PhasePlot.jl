

abstract type PhasePlotType end

struct PhasePlotType_PH <: PhasePlotType end
struct PhasePlotType_TS <: PhasePlotType end


function PhasePlot(type::PhasePlotType = PhasePlotType_TS(),sol::ODESolution)
    
end

