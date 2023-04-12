# helper functions for dealing with units in initial conditions and parameter dicts
# see: https://docs.sciml.ai/ModelingToolkit/dev/basics/Validation/#Parameter-and-Initial-Condition-Values
function remove_units(p::Dict)
    Dict(k => Unitful.ustrip(ModelingToolkit.get_unit(k), v) for (k, v) in p)
end

add_units(p::Dict) = Dict(k => v * ModelingToolkit.get_unit(k) for (k, v) in p)