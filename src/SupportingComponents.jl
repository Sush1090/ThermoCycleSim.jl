


@component function AmbientAirSource(f::Function;name)
    @named   air_port = AmbientNodePort()
    para = @parameters begin
        
    end
    vars = @variables begin
        T(t)
     end

   eqs = [
     D(T) ~ f(t)
     air_port.T ~ T
   ]
   compose(ODESystem(eqs, t, vars, para;name),air_port)
end