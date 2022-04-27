using Documenter
using STFT

makedocs(
    sitename = "STFT",
    format = Documenter.HTML(),
    modules = [STFT]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
