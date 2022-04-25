module STFT

using FFTW


_fft(x::AbstractMatrix{<:Real}, d) = rfft(x, d)
_fft(x::AbstractMatrix{<:Complex}, d) = fft(x, d)



"""
    analysis(x::Vector, w::Vector, L=0, N=length(w)) -> Matrix
    analysis(x::Array{Vector}, w::Vector, L=0, N=length(w)) -> Array{Matrix}


Analyse discrete time-domain signal ``\\mathrm{x}[n]``
using Short-Time Fourier Transform given by

```math
\\mathrm{X}[sH, \\omega] =
    \\sum_{n = -\\infty}^{+\\infty}
        \\mathrm{w}[n - sH] \\ \\mathrm{x}[n] e^{-j\\omega n},
```

where ``s`` and ``\\omega`` denotes segment index and angular frequency
respectively, ``\\mathrm{w}[n]`` is a discrete time-domain signal of analysis
window, ``H`` is nonnegative integer value that determine number of
samples between two consecutive signal segments (also known as `hop`).


# Parameters

- `x` - An array containing samples of a discrete time-domain signal.
- `w` - An array containing samples of a discrete time-domain window.
- `L` - An overlap in samples between two consecutive segments.
        Default value is `0`.
- `N` - A number of discrete frequency bins used in DFT computation.
        Default value is `length(w)`.
        If `N < length(w)`, then `N=length(w)` is enforced to avoid loss of
        information.


# Returns

- `X` - A complex matrix containing STFT-domain signal.


# Note

1. Relation between ``H`` and ``L`` is given as ``H = W - L`` where ``W`` is
   a length of a window. 
2. For real-valued (`x isa Real`) input signals function returns matrix is of
   size `(N÷2+1, S)` where `S` is a number of segments; i.e., one-sided
   spectrum.
3. For complex-valued (`x isa Complex`) input signals function returns matrix
   is of size `(N, S)` where `S` is a number of segments; i.e., two-sided
   spectrum.

# Examples
```julia
import STFT

x = rand(10000) # Generate mock signal
W = 64          # Window length
w = ones(W)     # Rectangular analysis window
H = 10          # Hop
L = W - H       # Overlap

X = STFT.analysis(x, w, L)  # Analysis
```
"""
function analysis()
end

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




function synthesis(
    X::AbstractMatrix{<:Complex},
    w::AbstractVector{<:Real},
    L::Integer = 0,
    N::Integer = length(w);
)::AbstractVector{<:Real}
    S = size(X, 2) # Number of segments
    W = length(w)  # Length of the window in samples
    H = W - L      # Hop
    K = H*(S-1)+W  # Expected length of synthesised signal
    w² = w.^2      # Squred window
    xn = zeros(K)  # Allocate memory for time-domain signal; numerator
    xd = zeros(K)  # Allocate memory for time-domain signal; denominator

    xs = irfft(X, N, 1) # Convert segments to time-domain

    for s ∈ 1:S
        ss = (s-1)*H # Segment start
        for k = 1:W
            xn[ss+k] += xs[k, s] * w[k]
            xd[ss+k] += w²[k]
        end
    end
    xn ./ xd # Normalize
end



end # module
