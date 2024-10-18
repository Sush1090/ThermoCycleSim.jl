
using CoolProp, ModelingToolkit

"""
`IsentropicCompression(πc, h_in, p_in,fluid,η)`

* Arguments:
    1. `πc`   : Pressure Ratio
    2. `h_in` : Inlet Enthalpy
    3. `p_in` : Inlet Pressure
    4. `fluid`: Fluid
    5. `η`    : Isentropic Efficiency

* Output -> Outlet enthalpy after isentropic compression
"""
function IsentropicCompression(πc, h_in, p_in,fluid,η)
    @assert η <= 1 "Efficiency more than 1"
    s_in = PropsSI("S", "H", h_in, "P", p_in, fluid)
    h_is = PropsSI("H", "S", s_in, "P",πc*p_in, fluid)
    h_out = h_in  + (h_is -h_in)/η
    return h_out
end
@register_symbolic IsentropicCompression(πc, h_in, p_in,fluid,η)
export IsentropicCompression


"""
`IsentropicExpansion(πc, h_in, p_in,fluid,η)`

* Arguments:
    1. `πc`   : Pressure Ratio
    2. `h_in` : Inlet Enthalpy
    3. `p_in` : Inlet Pressure
    4. `fluid`: Fluid
    5. `η`    : Isentropic Efficiency

* Output -> Outlet enthalpy after isentropic expansion
"""
function IsentropicExpansion(πc, h_in, p_in,fluid,η)
    @assert η <= 1 "Efficiency more than 1"
    s_in = PropsSI("S", "H", h_in, "P", p_in, fluid)
    h_is = PropsSI("H", "S", s_in, "P",p_in/πc, fluid)
    h_out = h_in  - η*(h_in - h_is)
    return h_out
end
@register_symbolic IsentropicExpansion(πc, h_in, p_in,fluid,η)
export IsentropicExpansion

"""
`IsothermalCompression(πc, h_in, p_in,fluid)`

* Arguments:
    1. `πc`   : Pressure Ratio
    2. `h_in` : Inlet Enthalpy
    3. `p_in` : Inlet Pressure
    4. `fluid`: Fluid

* Output -> Outlet enthalpy after Isothermal Compression
"""
function IsothermalCompression(πc, h_in, p_in,fluid)
    T_in = PropsSI("T", "H", h_in, "P", p_in, fluid)
    h_out = PropsSI("H", "T", T_in, "P",πc*p_in, fluid)
    return h_out
end
@register_symbolic IsothermalCompression(πc, h_in, p_in,fluid)
export IsothermalCompression



"""
`IsothermalExpansion(πc, h_in, p_in,fluid)`

* Arguments:
    1. `πc`   : Pressure Ratio
    2. `h_in` : Inlet Enthalpy
    3. `p_in` : Inlet Pressure
    4. `fluid`: Fluid

* Output -> Outlet enthalpy after Isothermal Expansion
"""
function IsothermalExpansion(πc, h_in, p_in,fluid)
    T_in = PropsSI("T", "H", h_in, "P", p_in, fluid)
    h_out = PropsSI("H", "T", T_in, "P", p_in/πc, fluid)
    return h_out
end
@register_symbolic IsothermalExpansion(πc, h_in, p_in,fluid)
export IsothermalExpansion
