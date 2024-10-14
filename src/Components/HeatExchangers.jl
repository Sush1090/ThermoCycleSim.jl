


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

function TwoPhase(;name,fluid,Δp)
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

function SimpleEvaporator(;name,fluid,Δp,ΔT_sh)
    @named inport = CoolantPort()
    @named pre = Preheater(fluid = fluid,Δp= Δp/3) 
    @named twophase = TwoPhase(fluid = fluid,Δp= Δp/3)
    @named super = SuperHeat(fluid = fluid,Δp= Δp/3,ΔT_sh = ΔT_sh)
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
    connect(inport,pre.inport)
    connect(pre.outport,twophase.inport)
    connect(twophase.outport,super.inport)
    connect(super.outport,outport)
    T_sat ~ twophase.T_sat
    T_in ~ pre.T_in
    T_out ~ super.T_out
    P ~ pre.P + twophase.P + super.P

    s_in ~ pre.inport.s
    p_in ~ pre.inport.p
    h_in ~ pre.inport.h
    ρ_in ~ pre.inport.ρ

    s_out ~ super.outport.s
    p_out ~ super.outport.p
    h_out ~ super.outport.h
    ρ_out ~ super.outport.ρ
]

compose(ODESystem(eqs, t, vars, para;name), inport,pre,twophase,super, outport)

end


export SimpleEvaporator, SuperHeat,Preheater,TwoPhase