# # Energy Distributions at Generator level

# In this example, we will generate 1 million events for ``2\nu`` and ``0\nu`` decays of ``^{76}\textrm{Ge}`` and plot summed electron energy spectra as sanity checks.

using DoubleBetaDecayGenerators
using FHist
using CairoMakie

dat_0nu = ZeroNuDBDData(:Ge76)
dat_2nu = TwoNuDBDData(:Ge76);

# Note that the data files needed for the generators are handled automatically via Artifacts API.

energies_0nu = [
    let (E1, E2, cosθ12) = rand(dat_0nu)
            E1 + E2
    end for i in 1:1_000_000
]
energies_2nu = [
    let (E1, E2, cosθ12) = rand(dat_2nu)
            E1 + E2
    end for i in 1:1_000_000
]

stephist(
    Hist1D(energies_0nu, binedges = 0:2040),
    axis = (;
        yscale = log10, limits = (nothing, (0.5, 10^7)),
        xlabel = L"$E_1$ + $E_2$ [keV]", ylabel = "Counts/(1 keV)",
        title = L"0\nu\beta\beta\;(N=10^6)",
    )
)

# As we expected, the ``0\nu\beta\beta`` spectrum is a peak at the ``Q``-value (2039 keV for ``^{76}\rm{Ge}``), while the ``2\nu\beta\beta`` spectrum is a broad continuum from 0 to the ``Q``-value.

stephist(
    Hist1D(energies_2nu, binedges = 0:2040),
    axis = (;
        yscale = log10, limits = (nothing, (0.5, 10^4)),
        xlabel = L"$E_1$ + $E_2$ [keV]", ylabel = "Counts/(1 keV)",
        title = L"2\nu\beta\beta\;(N=10^6)"
    )
)
