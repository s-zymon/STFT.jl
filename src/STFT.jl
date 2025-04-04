module STFT

using LoopVectorization
using FFTW

export stft, istft

_fft(x::AbstractArray{<:Real}, d) = rfft(x, d)
_fft(x::AbstractArray{<:Complex}, d) = fft(x, d)

function malloc_stft(
    x::M,
    w::V,
    L::I = zero(I),
    N::I = length(w)
) where {T<:Number,I<:Integer,V<:AbstractVector{T},M<:AbstractMatrix{T}}
    X, K = size(x)      # Length of the signal in samples
    W = length(w)       # Length of the window in samples
    S = (X - L) ÷ (W - L) # Number of segments
    N = N < W ? W : N   # DFT size
    return zeros(T, (N, S, K)) # Allocate container for signal segments
end


doc_analysis = """
    analysis(x::Vector, w::Vector, L=0, N=length(w)) -> Matrix
    analysis(x::Array{Vector}, w::Vector, L=0, N=length(w)) -> Array{Matrix}

    stft(x::Vector, w::Vector, L=0, N=length(w)) -> Matrix
    stft(x::Array{Vector}, w::Vector, L=0, N=length(w)) -> Array{Matrix}


Analyse discrete time-domain signal ``\\mathrm{x}[n]``
using Short-Time Fourier Transform given by

```math
\\mathrm{X}[sH, \\omega] =
    \\sum_{n = -\\infty}^{+\\infty}
        \\mathrm{w}[n - sH] \\ \\mathrm{x}[n] e^{-j\\omega n},
```

where ``s`` and ``\\omega`` denotes segment index and angular frequency
respectively, ``\\mathrm{w}[n]`` is a discrete time-domain signal of analysis
window, and ``H`` is nonnegative integer value that determine number of
samples between two consecutive signal segments (also known as `hop`).


# Parameters

- `x` - An array containing samples of a discrete time-domain signal.
- `w` - An array containing samples of a discrete time-domain analysis window.
- `L` - An overlap in samples between two consecutive segments.
        Default value is `0`.
- `N` - A number of discrete frequency bins to be used in the DFT computation.
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

```julia
using STFT

x = rand(100)   # Generate mock signal
W = 64          # Window length
w = ones(W)     # Rectangular analysis window
H = 4           # Hop
L = W - H       # Overlap

X  = stft(x, w, L)  # Analysis
```

"""

"$doc_analysis"
function analysis() end

"$doc_analysis"
stft(x, w, L=0, N=length(w)) = analysis(x, w, L, N)

function analysis(
    x::M,
    w::V,
    L::I=zero(I),
    N::I=length(w)
) where {T<:Number,I<:Integer,V<:AbstractVector{T},M<:AbstractMatrix{T}}
    sc = malloc_stft(x, w, L, N)
    N, S, K = sc |> size
    W = w |> length

    @turbo for k in 1:K, s in 1:S, n in 1:W
        sc[n, s, k] = w[n] * x[(s-1)*(W-L)+n, k]
    end
    return _fft(sc, 1) # Convert segments to frequency-domain
end

function analysis(
    x::V,
    w::V,
    L::I=zero(I),
    N::I=length(w)
)::Matrix{T |> complex} where {T<:Number,I<:Integer,V<:AbstractVector{T}}
    xx = @view x[:, :]
    return @view analysis(xx, w, L, N)[:, :, 1]
end

function analysis(
    xs::AbstractArray{<:AbstractVector{T}},
    w::AbstractVector{T},
    L::Integer=0,
    N::Integer=length(w)
)::AbstractArray{<:AbstractMatrix{<:Complex}} where {T<:Number}
    return [analysis(x, w, L, N) for x in xs]
end



doc_synthesis = """
    synthesis(X::Matrix, w::Vector, L=0, N=length(w)) -> Vector

    istft(X::Matrix, w::Vector, L=0, N=length(w)) -> Vector

Syntesise discrete time-domain signal ``y[n]`` from STFT-domain signal
``Y_w[sH, n]``.
An arbitrary STFT-domain signal ``Y_w[sH, n]``, in general, is not a valid STFT
in the sense that there is no discrete time-domian signal whoes STFT is given
by ``Y_w[sH, n]`` [1].
As result, time-domain signal must be estimated using following formula [1]:
```math
y[n] = \\frac{
    \\sum\\limits_{s=-\\infty}^{+\\infty} w[n - sH] \\ y_w[sH, n]
}{
    \\sum\\limits_{s=-\\infty}^{+\\infty} w^2[n - sH]
},
```
where ``w[n]`` is time-domain signal of analysis window,
``y_w[sH, n]`` is time-domain representation of ``Y_w[sH, n]``,
and ``H`` is nonnegative integer value that determine number of
samples between two consecutive signal segments (also known as `hop`).

# Parameters

- `X` - An matrix containing samples of a discrete STFT-domain signal.
- `w` - An array containing samples of a discrete time-domain analysis window.
- `L` - An overlap in samples between two consecutive segments.
        Default value is `0`.
- `N` - A number of discrete frequency bins to be used in the inverse DFT
        computation. Default value is `length(w)`.
        If `N < length(w)`, then `N=length(w)` is enforced to avoid loss of
        information.

# Returns

- `x` - A real-valued vector containing estimated time-domain signal.

# Note

1. Relation between ``H`` and ``L`` is given as ``H = W - L`` where ``W`` is
   a length of a window.
2. This function supports only synthesis of real-valued signals from
   one-sided STFT-domain signal.
3. Synthesised time-domain signal might be shorter than analysed one,
   since only whole segments are analysed.


# Example
```julia
import STFT

x = rand(100)   # Generate mock signal
W = 64          # Window length
w = ones(W)     # Rectangular analysis window
H = 4           # Hop
L = W - H       # Overlap

X  = STFT.analysis(x, w, L)  # Analysis
xr = STFT.synthesis(X, w, L) # Synthesis
```

```julia
using STFT

x = rand(100)   # Generate mock signal
W = 64          # Window length
w = ones(W)     # Rectangular analysis window
H = 4           # Hop
L = W - H       # Overlap

X  = stft(x, w, L)  # Analysis
xr = istft(X, w, L) # Synthesis
```


# References
1. D. Griffin and J. Lim, “Signal estimation from modified short-time
   Fourier transform,” IEEE Transactions on Acoustics, Speech, and
   Signal Processing, vol. 32, no. 2, pp. 236–243, Apr. 1984,
   doi: 10.1109/TASSP.1984.1164317.
   \\[[IEEE Xplore](https://ieeexplore.ieee.org/abstract/document/1164317)\\]
"""

"$doc_synthesis"
function synthesis() end

"$doc_synthesis"
istft(X, w, L=0, N=length(w)) = synthesis(X, w, L, N)

function synthesis(
    X::AbstractMatrix{<:Complex},
    w::AbstractVector{<:Real},
    L::Integer=0,
    N::Integer=length(w)
)::AbstractVector{<:Real}
    S = size(X, 2) # Number of segments
    W = length(w)  # Length of the window in samples
    H = W - L      # Hop
    K = H * (S - 1) + W  # Expected length of synthesised signal
    w² = w .^ 2      # Squred window
    xn = zeros(K)  # Allocate memory for time-domain signal; numerator
    xd = zeros(K)  # Allocate memory for time-domain signal; denominator

    xs = irfft(X, N, 1) # Convert segments to time-domain

    @turbo for s in 0:(S-1)
        ss = s * H # Segment shift
        for n in 1:W
            xn[ss+n] += xs[n, s+1] * w[n]
            xd[ss+n] += w²[n]
        end
    end
    return xn ./ xd # Normalize
end

"""
Real-value signal STFT with constant input size.
"""
function rSTFTm(A, w, L, N=length(w))
    mem = STFT.malloc_stft(A, w, L, N)
    N, S, K = mem |> size
    W = w |> length
    P = plan_rfft(mem, 1)
    function f(x)
        @turbo for k in 1:K, s in 1:S, n in 1:W
            mem[n, s, k] = w[n] * x[(s-1)*(W-L)+n, k]
        end
        return P * mem
    end
    return f
end


end # module
