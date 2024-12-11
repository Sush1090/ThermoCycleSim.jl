

struct DimensionlessAnalytic end
@component function storage(type::DimensionlessAnalytic;name,N=10)
    @named inport = CoolantPort()
    @named outport =  CoolantPort()
    para = @parameters begin

    end

    vars = @variables begin
        T(t)[1:N]
        x(t)[1:N]
        Δx(t)
    end

    eqs = [ 
        Δx ~ 25/N
        [x[i] ~ i*Δx for i = 1:N]
        [T[i] ~ 0.5*erfc((x[i] - t)/(2sqrt(t+eps()))) + sqrt(t/π)exp((-(x[i] - t)^2)/(4t)) - 0.5*(1 + x[i] + t)exp(x[i])erfc((x[i] + t)/(2sqrt(t+eps())))  for i = 1:N] 
    ]
    ODESystem(eqs,t,vars,para;name =name)

end

# I still need to figure out the relevant connections to the storage as a condensor and evaporator. This is dimensionless hence Temperature is between 0 and 1 but this temperature is not true for the wf. 