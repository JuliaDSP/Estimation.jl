module Estimation

using Docile
@docstrings

export
    apes,
    quin_fernandes

include("frequency.jl")
include("APES.jl")

end # module
