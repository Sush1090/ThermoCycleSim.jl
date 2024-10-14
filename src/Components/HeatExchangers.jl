


function HEX_type(;T_in = nothing,T_out = nothing,Qdot = nothing, mdot = nothing,p = nothing)

end


function IsobaricHeatSource(;name,Q_dot,fluid)
    @named inport = CoolantPort()
    @named outport = CoolantPort()
    para = @parameters begin
        
    end
    vars = @variables begin
        P(t)
        s_in(t)
        p_in(t)
        T_in(t)
        h_in(t)
        ρ_in(t)

        s_out(t)
        p_out(t)
        T_out(t)
        h_out(t)
        ρ_out(t)
     end
   eqs = [  outport.mdot ~ abs(inport.mdot) 
            outport.p ~ inport.p
            outport.h ~ inport.h + Q_dot/outport.mdot
            P ~ abs(inport.mdot)*(outport.h - inport.h)
            s_in ~ PropsSI("S","H",inport.h,"P",inport.p,fluid)
            p_in ~ inport.p
            T_in ~ PropsSI("T","H",inport.h,"P",inport.p,fluid)
            h_in ~ inport.h
            ρ_in ~ PropsSI("D","H",inport.h,"P",inport.p,fluid)
            s_out ~ PropsSI("S","H",outport.h,"P",outport.p,fluid)
            p_out ~ outport.p
            T_out ~ PropsSI("T","H",outport.h,"P",outport.p,fluid)
            h_out ~ outport.h
            ρ_out ~ PropsSI("D","H",outport.h,"P",outport.p,fluid)
   ]
   compose(ODESystem(eqs, t, vars, para;name), inport, outport)
end

export IsobaricHeatSource


function IsobaricHeatSink(;name,Q_dot,fluid)
    @named inport = CoolantPort()
    @named outport = CoolantPort()
    para = @parameters begin
        
    end
    vars = @variables begin
        P(t)
        s_in(t)
        p_in(t)
        T_in(t)
        h_in(t)
        ρ_in(t)

        s_out(t)
        p_out(t)
        T_out(t)
        h_out(t)
        ρ_out(t)
     end
   eqs = [  outport.mdot ~ abs(inport.mdot) 
            outport.p ~ inport.p
            outport.h ~ inport.h - Q_dot/outport.mdot
            P ~ abs(inport.mdot)*(outport.h - inport.h)
            s_in ~ PropsSI("S","H",inport.h,"P",inport.p,fluid)
            p_in ~ inport.p
            T_in ~ PropsSI("T","H",inport.h,"P",inport.p,fluid)
            h_in ~ inport.h
            ρ_in ~ PropsSI("D","H",inport.h,"P",inport.p,fluid)
            s_out ~ PropsSI("S","H",outport.h,"P",outport.p,fluid)
            p_out ~ outport.p
            T_out ~ PropsSI("T","H",outport.h,"P",outport.p,fluid)
            h_out ~ outport.h
            ρ_out ~ PropsSI("D","H",outport.h,"P",outport.p,fluid)
   ]
   compose(ODESystem(eqs, t, vars, para;name), inport, outport)
end

export IsobaricHeatSink



function Preheater(;name,fluid,Δp)
    @named inport = CoolantPort()
    @named outport = CoolantPort()
    para = @parameters begin

    end
    vars =  @variables begin
        P(t)
        s_in(t)
        p_in(t)
        T_in(t)
        h_in(t)
        ρ_in(t)

        s_out(t)
        p_out(t)
        T_out(t)
        h_out(t)
        ρ_out(t)

        T_sat(t)
    end
    eqs = [
        outport.mdot ~ abs(inport.mdot)
        outport.p ~ inport.p - Δp
        outport.h ~ PropsSI("H","Q",0,"P",outport.p,fluid)

        s_in ~ PropsSI("S","H",inport.h,"P",inport.p,fluid)
        p_in ~ inport.p
        T_in ~ PropsSI("T","H",inport.h,"P",inport.p,fluid)
        h_in ~ inport.h
        ρ_in ~ PropsSI("D","H",inport.h,"P",inport.p,fluid)
            
        s_out ~ PropsSI("S","H",outport.h,"P",outport.p,fluid)
        p_out ~ outport.p
        T_out ~ PropsSI("T","H",outport.h,"P",outport.p,fluid)
        h_out ~ outport.h
        ρ_out ~ PropsSI("D","H",outport.h,"P",outport.p,fluid)
        
        P ~ outport.mdot*(outport.h - inport.h)
        
        T_sat ~ PropsSI("T","Q",0,"P",outport.p,fluid)
    ]
    compose(ODESystem(eqs, t, vars, para;name), inport, outport)
end

function TwoPhaseEvap(;name,fluid,Δp)
    @named inport = CoolantPort()
    @named outport = CoolantPort()
    para = @parameters begin

    end
    vars =  @variables begin
        P(t)
        s_in(t)
        p_in(t)
        T_in(t)
        h_in(t)
        ρ_in(t)

        s_out(t)
        p_out(t)
        T_out(t)
        h_out(t)
        ρ_out(t)

        T_sat(t)
    end
    eqs = [
        outport.mdot ~ abs(inport.mdot)
        outport.p ~ inport.p - Δp
        outport.h ~ PropsSI("H","Q",1,"P",outport.p,fluid)
        
        s_in ~ PropsSI("S","H",inport.h,"P",inport.p,fluid)
        p_in ~ inport.p
        T_in ~ PropsSI("T","H",inport.h,"P",inport.p,fluid)
        h_in ~ inport.h
        ρ_in ~ PropsSI("D","H",inport.h,"P",inport.p,fluid)
            
        s_out ~ PropsSI("S","H",outport.h,"P",outport.p,fluid)
        p_out ~ outport.p
        T_out ~ PropsSI("T","H",outport.h,"P",outport.p,fluid)
        h_out ~ outport.h
        ρ_out ~ PropsSI("D","H",outport.h,"P",outport.p,fluid)
        
        P ~ outport.mdot*(outport.h - inport.h)
        
        T_sat ~ PropsSI("T","Q",0,"P",outport.p,fluid)
    ]
    compose(ODESystem(eqs, t, vars, para;name), inport, outport)
end

function SuperHeat(;name,fluid,ΔT_sh,Δp)
    @named inport = CoolantPort()
    @named outport = CoolantPort()
    para = @parameters begin

    end
    vars =  @variables begin
        P(t)
        s_in(t)
        p_in(t)
        T_in(t)
        h_in(t)
        ρ_in(t)

        s_out(t)
        p_out(t)
        T_out(t)
        h_out(t)
        ρ_out(t)

        T_sat(t)
    end
    eqs = [
        outport.mdot ~ abs(inport.mdot)
        outport.p ~ inport.p - Δp
        T_in ~ PropsSI("T","H",inport.h,"P",outport.p,fluid)
        T_out ~ T_in + ΔT_sh
        outport.h ~ PropsSI("H","T",T_out,"P",outport.p,fluid)

        s_in ~ PropsSI("S","H",inport.h,"P",inport.p,fluid)
        p_in ~ inport.p
        h_in ~ inport.h
        ρ_in ~ PropsSI("D","H",inport.h,"P",inport.p,fluid)
            
        s_out ~ PropsSI("S","H",outport.h,"P",outport.p,fluid)
        p_out ~ outport.p
        h_out ~ outport.h
        ρ_out ~ PropsSI("D","H",outport.h,"P",outport.p,fluid)
        
        P ~ outport.mdot*(outport.h - inport.h)
        
        T_sat ~ PropsSI("T","Q",1,"P",outport.p,fluid)
    ]
    compose(ODESystem(eqs, t, vars, para;name), inport, outport)
end


"""
`SimpleEvaporator(;name,fluid,Δp,ΔT_sh)`
    Composed of multiple `ODESystem` - `Preheater`, `TwoPhaseEvap`, and `SuperHeat`
*    Arguments: 
    1. `Δp`     : Pressure Drop across Evaporator
    2. `ΔT_sh`  : Super heated Temperature
    
*    Local Variables:
    1. `P`      : Power  
    2. `s_in`   : Inlet Entropy
    3. `p_in`   : Inlet Pressure
    4. `T_in`   : Inlet Temperature
    5. `h_in`   : Inlet Enthalpy
    6. `ρ_in`   : Inlet Density
    7. `s_out`  : Outlet Entropy
    8. `p_out`  : Outlet Pressure
    9. `T_out`  : Outlet Temperature
    10. `h_out` : Outlet Enthalpy
    11. `ρ_out` : Outlet Density
    12. `T_sat` : Saturation Temperature

*    Port Variables:
    1. `inport`         : `p` and `h`
    2. `outport`        : `p` and `h`
    3. `preheater`      : Component with internal `CoolantPort` --> `inport` and `outport`
    4. `TwoPhaseEvap`   : Component with internal `CoolantPort` --> `inport` and `outport`
    5. `superheater`   : Component with internal `CoolantPort` --> `inport` and `outport`
"""
function SimpleEvaporator(;name,fluid,Δp,ΔT_sh)
    @named inport = CoolantPort()
    @named preheater = Preheater(fluid = fluid,Δp= Δp/3) 
    @named TwoPhaseEvap = TwoPhaseEvap(fluid = fluid,Δp= Δp/3)
    @named superheater = SuperHeat(fluid = fluid,Δp= Δp/3,ΔT_sh = ΔT_sh)
    @named outport = CoolantPort()

    para = @parameters begin

    end
    vars =  @variables begin
        P(t)
        s_in(t)
        p_in(t)
        T_in(t)
        h_in(t)
        ρ_in(t)

        s_out(t)
        p_out(t)
        T_out(t)
        h_out(t)
        ρ_out(t)

        T_sat(t)
    end
    eqs = [
    connect(inport,preheater.inport)
    connect(preheater.outport,TwoPhaseEvap.inport)
    connect(TwoPhaseEvap.outport,superheater.inport)
    connect(superheater.outport,outport)
    T_sat ~ TwoPhaseEvap.T_sat
    T_in ~ preheater.T_in
    T_out ~ superheater.T_out
    P ~ preheater.P + TwoPhaseEvap.P + superheater.P

    s_in ~ preheater.inport.s
    p_in ~ preheater.inport.p
    h_in ~ preheater.inport.h
    ρ_in ~ preheater.inport.ρ

    s_out ~ superheater.outport.s
    p_out ~ superheater.outport.p
    h_out ~ superheater.outport.h
    ρ_out ~ superheater.outport.ρ
]

compose(ODESystem(eqs, t, vars, para;name), inport,preheater,TwoPhaseEvap,superheater, outport)

end


export SimpleEvaporator, SuperHeat,Preheater,TwoPhaseEvap



function Precooler(;name,fluid,Δp)
    @named inport = CoolantPort()
    @named outport = CoolantPort()
    para = @parameters begin

    end
    vars =  @variables begin
        P(t)
        s_in(t)
        p_in(t)
        T_in(t)
        h_in(t)
        ρ_in(t)

        s_out(t)
        p_out(t)
        T_out(t)
        h_out(t)
        ρ_out(t)

        T_sat(t)
    end
    eqs = [
        outport.mdot ~ abs(inport.mdot)
        outport.p ~ inport.p - Δp
        outport.h ~ PropsSI("H","Q",1,"P",outport.p,fluid)

        s_in ~ PropsSI("S","H",inport.h,"P",inport.p,fluid)
        p_in ~ inport.p
        T_in ~ PropsSI("T","H",inport.h,"P",inport.p,fluid)
        h_in ~ inport.h
        ρ_in ~ PropsSI("D","H",inport.h,"P",inport.p,fluid)
            
        s_out ~ PropsSI("S","H",outport.h,"P",outport.p,fluid)
        p_out ~ outport.p
        T_out ~ PropsSI("T","H",outport.h,"P",outport.p,fluid)
        h_out ~ outport.h
        ρ_out ~ PropsSI("D","H",outport.h,"P",outport.p,fluid)
        
        P ~ outport.mdot*(outport.h - inport.h)
        
        T_sat ~ PropsSI("T","Q",1,"P",outport.p,fluid)
    ]
    compose(ODESystem(eqs, t, vars, para;name), inport, outport)
end


function TwoPhaseCond(;name,fluid,Δp)
    @named inport = CoolantPort()
    @named outport = CoolantPort()
    para = @parameters begin

    end
    vars =  @variables begin
        P(t)
        s_in(t)
        p_in(t)
        T_in(t)
        h_in(t)
        ρ_in(t)

        s_out(t)
        p_out(t)
        T_out(t)
        h_out(t)
        ρ_out(t)

        T_sat(t)
    end
    eqs = [
        outport.mdot ~ abs(inport.mdot)
        outport.p ~ inport.p - Δp
        outport.h ~ PropsSI("H","Q",0,"P",outport.p,fluid)
        
        s_in ~ PropsSI("S","H",inport.h,"P",inport.p,fluid)
        p_in ~ inport.p
        T_in ~ PropsSI("T","H",inport.h,"P",inport.p,fluid)
        h_in ~ inport.h
        ρ_in ~ PropsSI("D","H",inport.h,"P",inport.p,fluid)
            
        s_out ~ PropsSI("S","H",outport.h,"P",outport.p,fluid)
        p_out ~ outport.p
        T_out ~ PropsSI("T","H",outport.h,"P",outport.p,fluid)
        h_out ~ outport.h
        ρ_out ~ PropsSI("D","H",outport.h,"P",outport.p,fluid)
        
        P ~ outport.mdot*(outport.h - inport.h)
        
        T_sat ~ PropsSI("T","Q",0,"P",outport.p,fluid)
    ]
    compose(ODESystem(eqs, t, vars, para;name), inport, outport)
end



function SuperCooler(;name,fluid,ΔT_sc,Δp)
    @named inport = CoolantPort()
    @named outport = CoolantPort()
    para = @parameters begin

    end
    vars =  @variables begin
        P(t)
        s_in(t)
        p_in(t)
        T_in(t)
        h_in(t)
        ρ_in(t)

        s_out(t)
        p_out(t)
        T_out(t)
        h_out(t)
        ρ_out(t)

        T_sat(t)
    end
    eqs = [
        outport.mdot ~ abs(inport.mdot)
        outport.p ~ inport.p - Δp
        T_in ~ PropsSI("T","H",inport.h,"P",outport.p,fluid)
        T_out ~ T_in - ΔT_sc
        outport.h ~ PropsSI("H","T",T_out,"P",outport.p,fluid)

        s_in ~ PropsSI("S","H",inport.h,"P",inport.p,fluid)
        p_in ~ inport.p
        h_in ~ inport.h
        ρ_in ~ PropsSI("D","H",inport.h,"P",inport.p,fluid)
            
        s_out ~ PropsSI("S","H",outport.h,"P",outport.p,fluid)
        p_out ~ outport.p
        h_out ~ outport.h
        ρ_out ~ PropsSI("D","H",outport.h,"P",outport.p,fluid)
        
        P ~ outport.mdot*(outport.h - inport.h)
        
        T_sat ~ PropsSI("T","Q",1,"P",outport.p,fluid)
    ]
    compose(ODESystem(eqs, t, vars, para;name), inport, outport)
end




function SimpleCondensor(;name,fluid,Δp,ΔT_sc)
    @named inport = CoolantPort()
    @named precooler = Precooler(fluid = fluid,Δp= Δp/3) 
    @named twophasecond = TwoPhaseCond(fluid = fluid,Δp= Δp/3)
    @named supercooler = SuperCooler(fluid = fluid,Δp= Δp/3,ΔT_sc = ΔT_sc)
    @named outport = CoolantPort()

    para = @parameters begin

    end
    vars =  @variables begin
        P(t)
        s_in(t)
        p_in(t)
        T_in(t)
        h_in(t)
        ρ_in(t)

        s_out(t)
        p_out(t)
        T_out(t)
        h_out(t)
        ρ_out(t)

        T_sat(t)
    end
    eqs = [
    connect(inport,precooler.inport)
    connect(precooler.outport,twophasecond.inport)
    connect(twophasecond.outport,supercooler.inport)
    connect(supercooler.outport,outport)
    T_sat ~ twophasecond.T_sat
    T_in ~ precooler.T_in
    T_out ~ supercooler.T_out
    P ~ precooler.P + twophasecond.P + supercooler.P

    s_in ~ precooler.inport.s
    p_in ~ precooler.inport.p
    h_in ~ precooler.inport.h
    ρ_in ~ precooler.inport.ρ

    s_out ~ supercooler.outport.s
    p_out ~ supercooler.outport.p
    h_out ~ supercooler.outport.h
    ρ_out ~ supercooler.outport.ρ
]

compose(ODESystem(eqs, t, vars, para;name), inport,precooler,twophasecond,supercooler, outport)

end

export SimpleCondensor, Precooler, TwoPhaseCond, SuperCooler