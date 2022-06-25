# `STFT.jl` - Short-Time Fourier Transform


[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://docs.zymon.org/STFT.jl/)

`STFT.jl` is a small Julia package implementing just Short-Time Fourier Transform (STFT) routines.
It provides the following core functionality:
- **_signal analysis_**; transform time-domain signal to STFT-domain signal.
- **_signal synthesis_**; transform STFT-domain signal to time-domain signal.

Check the [documentation](https://docs.zymon.org/STFT.jl/) for more insights.

## Installation

The package is currently available in General, the default Julia package registry.
To install this package from General registry, use the following command in Julia REPL:
```julia
] add STFT
```
Alternatively, directly via repository:
```julia
pkg> add https://github.com/s-zymon/STFT.jl
```


## Examples

Below you can find a few standalone examples with basic usage of the package.

### Show spectrogram

```julia
using STFT
using Plots

x = randn(10000)  # Generate mock signal
W = 64            # Window length
w = ones(W)       # Rectangular analysis window
H = 10            # Hop
L = W - H         # Overlap

X = stft(x, w, L)    # Analysis
s = abs2.(X)         # Compute spectrogram
heatmap(10log10.(s)) # Display spectrogram
```

### Analyse signal, modify, and synthesise
```julia
using STFT

x = randn(10000)   # Generate mock signal
W = 64             # Window length
w = ones(W)        # Rectangular analysis window
H = 10             # Hop
L = W - H          # Overlap

X = stft(x, w, L)  # Analysis
X = f(X)           # Modify STFT-domain signal
y = istft(X, w, L) # Synthesis
```

Alternatively, instead of `using STFT`, you can `import STFT`,
and use an alternative API, i.e., `analysis` and `synthesis`.

```julia
import STFT

x = randm(10000) # Generate mock signal
W = 64           # Window length
w = ones(W)      # Rectangular analysis window
H = 10           # Hop
L = W - H        # Overlap

X = STFT.analysis(x, w, L)  # Analysis
X = f(X)                    # Modify STFT-domain signal
y = STFT.synthesis(X, w, L) # Synthesis
```

