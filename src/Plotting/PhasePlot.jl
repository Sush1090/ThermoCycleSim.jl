

abstract type PhasePlotType end

struct PhasePlotType_PH <: PhasePlotType end
struct PhasePlotType_TS <: PhasePlotType end

export PhasePlotType_PH,PhasePlotType_TS

T_line(p,h,fluid) = PropsSI("T","P",p,"H",h,fluid)
S_line(p,h,fluid) = PropsSI("S","P",p,"H",h,fluid)
ρ_line(p,h,fluid) = PropsSI("D","P",p,"H",h,fluid)

function LineData(f::Function,x1,y1,x2,y2,fluid::AbstractString;points = 100)
    x = collect(range(x1,x2,points))
    y = collect(range(y1,y2,points))
    data = f.(x,y,fluid)
    return x,y,data
end


function CollectPhaseData(sol::ODESolution,system::Vector{ODESystem},fluid::AbstractString;points = 100,p_width::Tuple = (1e3,1e3), h_width::Tuple = (1e3,1e3))
    plot_sys = system[2:end-1];
    propx = :p
    propy = :h
    pp = [sol[getproperty(i.inport, propx)][1] for i in plot_sys]
    hh = [sol[getproperty(i.inport, propy)][1] for i in plot_sys]
    pmin = minimum(pp); pmax = maximum(pp)
    hmin = minimum(hh); hmax = maximum(hh)
    p_range = collect(range(pmin-p_width[1],pmax + p_width[2],points));
    h_range = collect(range(hmin-h_width[1],hmax + h_width[2],points));

    T_sat = PropsSI.("T","Q",1,"P",p_range,fluid)
    s_sat_liquid = PropsSI.("S","Q",0,"P",p_range,fluid)
    s_sat_gas = PropsSI.("S","Q",1,"P",p_range,fluid)
    h_sat_liquid = PropsSI.("H","Q",0,"P",p_range,fluid)
    h_sat_gas = PropsSI.("H","Q",1,"P",p_range,fluid)
    return p_range,h_range,T_sat,s_sat_liquid,s_sat_gas,h_sat_liquid,h_sat_gas
end

function PhasePlot(type::PhasePlotType_TS,sol::ODESolution,system::Vector{ODESystem},fluid)
    plot_sys = system[2:end-1];
    propx = :p
    propy = :h
    pp = [sol[getproperty(i.inport, propx)][1] for i in plot_sys]
    hh = [sol[getproperty(i.inport, propy)][1] for i in plot_sys]

    p_range,h_range,T_sat,s_sat_liquid,s_sat_gas,h_sat_liquid,h_sat_gas = CollectPhaseData(sol,system,fluid)
    T_cycle_points = PropsSI.("T","P",pp,"H",hh,fluid) 
    s_cycle_points = PropsSI.("S","P",pp,"H",hh,fluid)
    T_data = []; s_data = []
    cycle_data_pp = [pp;pp[1]]; cycle_data_hh = [hh;hh[1]]
    for i = 1:size(pp,1)
        x1 = cycle_data_pp[i]; y1 = cycle_data_hh[i];
        x2 = cycle_data_pp[i+1]; y2 = cycle_data_hh[i+1];
        _,_,TT = LineData(T_line,x1,y1,x2,y2,fluid)
        T_data = append!(T_data,TT)
        _,_,ss = LineData(S_line,x1,y1,x2,y2,fluid)
        s_data = append!(s_data,ss)
    end
    fig1 = plot(s_data,T_data,label = "Cycle")
     scatter!(s_cycle_points,T_cycle_points,label="Cycle Points")
     plot!(s_sat_liquid,T_sat,label ="Saturation Curve liquid")
     plot!(s_sat_gas,T_sat,label ="Saturation Curve gas")
    xlabel!("Entropy")
    ylabel!("Temperature")
    return fig1
end

export PhasePlot