module RRCFilters

using Random: AbstractRNG, default_rng, randn
using Statistics: mean

export agc, awgn, biterr, carriersync, convenc, crcconfig, crcdetect, crcgenerate, equalize, evm, freqoffset, poly2trellis, qamdemod, qammod, pskdemod, pskmod, rcosdesign, symbolsync, timingoffset, upfirdn, vitdec

include("common.jl")
include("metrics.jl")
include("qam.jl")
include("psk.jl")
include("awgn.jl")
include("rcosdesign.jl")
include("upfirdn.jl")
include("carriersync.jl")
include("freqoffset.jl")
include("symbolsync.jl")
include("timingoffset.jl")
include("agc.jl")
include("fec.jl")
include("crc.jl")
include("equalizer.jl")

end # module RRCFilters
