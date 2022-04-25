module STFT

using FFTW


_fft(x::AbstractMatrix{<:Real}, d) = rfft(x, d)
_fft(x::AbstractMatrix{<:Complex}, d) = fft(x, d)



function analysis(
    x::AbstractVector{T},
    w::AbstractVector{T},
    L::Integer = 0,
    N::Integer = length(w);
)::AbstractMatrix{<:Complex} where {T<:Number}
    X = length(x)       # Length of the signal in samples
    W = length(w)       # Length of the window in samples
    H = W - L           # Hop
    S = (X-L) ÷ H       # Number of segments
    N = N < W ? W : N   # DFT size
    sc = zeros(T, N, S) # Allocate container for signal segments

    for s ∈ 1:S, n ∈ 1:W # Slice the signal
        sc[n, s] = w[n] * x[(s-1)*H+n]
    end
    _fft(sc, 1) # Convert segments to frequency-domain
end

function analysis(
    xs::AbstractArray{<:AbstractVector{T}},
    w::AbstractVector{T},
    L::Integer = 0,
    N::Integer = length(w);
)::AbstractArray{<:AbstractMatrix{<:Complex}} where {T<:Number}
    [analysis(x, w, L, N) for x ∈ xs]
end


end # module
