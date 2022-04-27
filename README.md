# STFT.jl


[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://docs.zymon.org/STFT.jl/)

`STFT.jl` is a julia package implementing Short-Time Fourier Transform routines.
It provides signal analysis (time-domain signal to STFT-domain signal; stft)
and signal synthesis (STFT-domain siganl to time-domain signal; istft).


# Examples

```julia
import STFT

x = rand(10000) # Generate mock signal
W = 64          # Window length
w = ones(W)     # Rectangular analysis window
H = 10          # Hop
L = W - H       # Overlap

X = STFT.analysis(x, w, L)  # Analysis

# Compute spectogram of the signal
spectogram = abs2.(X)

# X = f(X) # Modify STFT-domain signal

# Reconsturction
xr = STFT.synthesis(X, w, L) # Synthesis
```
