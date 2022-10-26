### A Pluto.jl notebook ###
# v0.19.9

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 97130b40-e41e-46c1-ac48-7fc734e05d95
begin
using Pkg
Pkg.activate("/Users/hernando/work/investigacion/NEXT/software/julias/Clouds")
import Clouds as cl
end

# ╔═╡ 82c494e0-b75e-474e-b08c-cde3791a121e
begin
Pkg.activate("/Users/hernando/work/investigacion/NEXT/software/julias/Nana")
import Nana as na
end

# ╔═╡ 8aaabf9e-08f7-11ed-1880-e798d62d8295
# NEXT Julia with Pluto

begin
using HDF5
using DataFrames
using Plots
using Statistics
using PlutoUI
import GraphPlot  as GP
#using PlotlyJS
end

# ╔═╡ 03bfbaeb-f6f2-4a8b-8815-66d293cbbb66
begin
#ng =ingredients("../src/julne.jl")
PlutoUI.TableOfContents(title="NEXT-100 Event Gallery", indent=true)
end

# ╔═╡ cdc50171-b288-40b6-9d0d-9511901218e0
md"""

## Description


This NB displays NEXT-100 MC events using Beersheba-reconstructed and MC hits of $^{214}$Bi and $\beta\beta0\nu$ tracks in the $Q_{\beta\beta}$ energy range


J.A. Hernado, JJ. Gomez-Cadenas

Donostia, July 2022

---
"""

# ╔═╡ f8dbfa77-c88d-42f6-bc2a-600b49f8f98d
# NEXT-100 214Bi
begin
datadir   = "/Users/hernando/work/investigacion/NEXT/data/NEXT100/"
filenames = Dict(:Bi214 => "Bi/label_beersheba_554mm_214Bi_ICS.h5",
                 :bb0nu => "bb0nu/v1/label_beersheba_554mm_0nubb.h5");
#filename = "/Users/hernando/work/investigacion/NEXT/data/NEXT100/Bi/label_beersheba_554mm_214Bi_ICS.h5"
end;

# ╔═╡ 8fe43246-a006-4ff9-bb64-3601d4e03938
begin
biso = @bind isotope Select([:Bi214, :bb0nu])

md"""

## Load the data

Select an isotope: $(biso)
"""
end

# ╔═╡ 6e7e284e-88a1-47f9-9da8-4a79797d8ebd
begin
filename      = string(datadir, filenames[isotope])
df, mc, steps = na.load_data(filename)
end;

# ╔═╡ 03d81ddf-105e-4ea6-a882-e1b40b7ecbfc
begin
events = sort(collect(Set(df.dataset_id)))

md"""

Selected isotope  : **$(isotope)**

Datafile          : $(filename)


Number of events : **$(length(events))**

Total number of entries in Beersheba's hit table: $(size(df)[1])

Voxel size:

| x (mm) | y (mm) | z (mm)|
| :-- | :-- | :-- |
| $(steps[1])| $(steps[2]) | $(steps[3]) |


"""
end

# ╔═╡ 98ab81f8-f630-4d7e-9624-829ec179dd10
begin
selevts = Dict()
bievts  = Dict()
bievts["number of tracks > 1"] = [34, 43, 112, 142, 157, 166, 191, 264, 296, 317, 329, 335, 348, 354, 364, 409, 421, 474, 494, 503]
bievts["graph no in roi"] = [2, 10, 23, 25, 32, 34, 40, 43, 86, 93, 98, 100, 101, 107, 112, 121, 134, 142, 157, 166]
bievts["graph has no blobs"] =  [12, 95, 138, 188, 237, 244, 261, 288, 302, 393]

bbevts = Dict()
bbevts["number of tracks > 1"]  = [34, 43, 112, 142, 157, 166, 191, 264, 296, 317, 329, 335, 348, 354, 364, 409, 421, 474, 494, 503]
bbevts["graph no in roi"] = [2, 10, 23, 25, 32, 34, 40, 43, 86, 93, 98, 100, 101, 107, 112, 121, 134, 142, 157, 166]
bbevts["graph has no blobs"] = [12, 95, 138, 188, 237, 244, 261, 288, 302, 393]


if (isotope == :Bi214)
	bievts["events"] = events[1:20]
else (isotope == :bb0nu)
	bbevts["events"] = events[1:20]
end

isoevts = isotope == :Bi214 ? bievts : bbevts

end;

# ╔═╡ 7913b885-4f25-40ae-ab84-3b5a4d50967c
begin


btypeevt = @bind typeevt Select(collect(keys(isoevts)))

md"""
---

Select event number : $(btypeevt)
"""
end

# ╔═╡ 44b35007-912d-49e3-90bb-09aa1360cbe9
begin

bevent = @bind event Select(isoevts[typeevt])

md"""
Select event number : $(bevent)
"""
end

# ╔═╡ 14160a98-07a0-4efc-9e71-9aab36ed01b6
begin
idf = na.get_event(df, event)
imc = na.get_event(mc, event)

md"""

Event $(event) has $(size(idf)[1]) Beersheba and $(size(imc)[1]) MC hits.

----
"""
end

# ╔═╡ 99a11ca4-b453-4e0b-a737-d9733b0e8c59
md"## Event Gallery"

# ╔═╡ 718fbf98-1725-429b-9a1c-609cb38430ac
begin
xwidth, ywidth, zwidth = steps
zrange = range(minimum(idf.z), maximum(idf.z) + zwidth; step = zwidth)
xrange = range(minimum(idf.x), maximum(idf.x) + xwidth; step = xwidth)
yrange = range(minimum(idf.y), maximum(idf.y) + ywidth; step = ywidth)
dx = maximum(idf.x) - minimum(idf.x)
dy = maximum(idf.y) - minimum(idf.y)
dz = maximum(idf.z) - minimum(idf.z)

md"""

Event window:

| coordinate | min (mm) | max (mm)| width (mm)
| :-- | :-- | :-- | :-- |
| x | $(minimum(idf.x)) | $(maximum(idf.x))| $(dx) |
| y | $(minimum(idf.y)) | $(maximum(idf.y))| $(dy) |
| z | $(minimum(idf.z)) | $(maximum(idf.z))| $(dz) |

"""
end

# ╔═╡ 65597330-86b4-4d17-abae-5df638603cfc
plotly();
#gr()

# ╔═╡ ca22eb02-282a-412b-9297-15a1e379a7b2
md"""

### MC: S2 signal, hit proyections
"""

# ╔═╡ 5f913a6e-1798-42da-8fbe-3dad64fc7bb7
begin

hezmc_ = histogram(imc.z, bins = zrange, weights = imc.energy,
	xlabel = "z (mm)", ylabel = "energy (MeV)")
hxzmc_ = histogram2d(imc.x, imc.z, bins = (xrange, zrange),  xlabel = "x (mm)", zlabel = "x (mm)")
hzymc_ = histogram2d(imc.z, imc.y, bins = (zrange, yrange),  xlabel = "z (mm)", ylabel = "y (mm)")
hxymc_ = histogram2d(imc.x, imc.y, bins = (xrange, yrange),  xlabel = "x (mm)", ylabel = "y (MeV)")
plot(hezmc_, hxzmc_, hzymc_, hxymc_; layout = grid(2, 2), size = (700, 600))
end

# ╔═╡ bd05aad2-f55a-424d-a467-f17d32286017
md"""

### RECO: S2 signal hit proyections
"""

# ╔═╡ aa55ee52-f3d0-4425-ba88-cafc4743ecf4
begin

hez_ = histogram(idf.z, bins = zrange, weights = idf.energy,
	xlabel = "z (mm)", ylabel = "energy (MeV)")
hxz_ = histogram2d(idf.x, idf.z, bins = (xrange, zrange),  xlabel = "x (mm)", zlabel = "x (mm)")
hzy_ = histogram2d(idf.z, idf.y, bins = (zrange, yrange),  xlabel = "z (mm)", ylabel = "y (mm)")
hxy_ = histogram2d(idf.x, idf.y, bins = (xrange, yrange),  xlabel = "x (mm)", ylabel = "y (MeV)")
plot(hez_, hxz_, hzy_, hxy_; layout = grid(2, 2), size = (700, 600))
end

# ╔═╡ d13e0a1a-8eb4-49eb-a8b7-b732cf9ffe31
md"""
### Scatter MC and RECO
"""

# ╔═╡ f547f778-8a74-4923-9373-cb96b2254c9b
begin
scxz = scatter(idf.x, idf.z,  markersize = 1e2*idf.energy,
	xlabel = "x (mm)", ylabel = "z (mm)")
scatter!(imc.x, imc.z, markercolor = "white", alpha = 0.1)
scyz = scatter(idf.z, idf.y,  markersize = 1e2*idf.energy,
	xlabel = "z (mm)", ylabel = "y (mm)")
scatter!(imc.z, imc.y, markercolor = "white", alpha = 0.1)
scxy = scatter(idf.x, idf.y,  markersize = 1e2*idf.energy,
		xlabel = "x (mm)", ylabel = "y (mm)")
scatter!(imc.x, imc.y, markercolor = "white", alpha = 0.1)
plot(scxy, scyz, scxz, grid = (2, 2), size = (700, 600), legend = false)
end

# ╔═╡ 420900c2-74f9-4972-bddd-982e766fc22f
md"""
### 3D Reco and MC hits
"""

# ╔═╡ 5e7771ee-372b-4fa9-8f8a-888c6e46796c
md"""

### MC 3D Hits  segmentation and energy
"""

# ╔═╡ 4f6f8003-b9f0-4aba-a42c-c79184bec77d
md"""
## Code
----

"""

# ╔═╡ 69920bad-6d3e-4e7a-a539-f08bc63460c4
kdf = na.get_event(df, 1)

# ╔═╡ e1c6ef24-fe3a-41b5-95a7-f0ddb7b3bf6c
begin
coors    = Tuple((kdf[!, :x], kdf[!, :y], kdf[!, :z]))
contents = kdf[!, :energy]
xsteps   = Tuple(steps)
end

# ╔═╡ 59f1b12f-3224-4ea1-8f55-98f228cdcfb6
xsteps

# ╔═╡ 7397c3e4-d057-42d2-947f-800dab52ee6f
xcl, xnd, graph, edges = cl.clouds(coors, contents, xsteps)

# ╔═╡ 59d569ec-c2a4-485b-88ac-86495ba1d6ab
xnd

# ╔═╡ 41e79a3a-fea2-449b-b42c-50ba020e5482
GP.gplot(graph)

# ╔═╡ aa5c46fc-ccf6-49e1-b20b-7714bf2c7630
"""
linear scale a vector of values between a minimum and a maximum

Parameters:
	var : Vector{Real}
	emin: Real
	emax: Real

Return:
	vvar: Vector{Real}

"""
function vscale(var, emin = 1, emax = 4)
	vvar = (var .- minimum(var)) ./(maximum(var) - minimum(var)) .*(emax-emin) .+ emin
	return vvar
end

# ╔═╡ 30970210-d5f6-4a1b-b36a-601789b0fd8b
begin
theme(:dark)
sc2 = scatter3d(idf.x, idf.y, idf.z, marker_z = idf.energy,
	markersize = vscale(idf.energy),
	markertype = "circle", label = false, alpha = 0.4, c = :inferno,
	xlabel = "x (mm)", ylabel = "y (mm)", zlabel = "z (mm)")
scatter3d!(imc.x, imc.y, imc.z, color = "white", alpha = 0.4, markersize = 0.5)
plot(sc2, size = (700, 600))
end

# ╔═╡ 5e61e1c2-a1d7-49a1-bd3f-649d4c6041cf
begin
theme(:dark)
sc3 = scatter3d(imc.x, imc.y, imc.z, marker_z = imc.segclass,
		markersize = vscale(imc.energy, 2, 6),
		markertype = "circle", label = false, alpha = 0.4, c = :inferno,
		xlabel = "x (mm)", ylabel = "y (mm)", zlabel = "z (mm)")
plot(sc3, size = (700, 600))
end

# ╔═╡ ed0c4b7e-cf34-4b20-9100-43ba5428e096
function load_data(filename)

	#dfs           = ng.julne.get_dfs(filename)
	dfs           = ng.julne.get_dfs(filename)
	x0, x1, steps = _get_dimensions(dfs)

	df            = dfs["BeershebaVoxels"]
	df[!, "x"] = x0[1] .+ (df[!, "xbin"] .+ 0.5) * steps[1]
	df[!, "y"] = x0[2] .+ (df[!, "ybin"] .+ 0.5) * steps[2]
	df[!, "z"] = x0[3] .+ (df[!, "zbin"] .+ 0.5) * steps[3]

	mc            = dfs["MCHits"]
	return df, mc, steps
end;

# ╔═╡ 0b5acafe-fc8f-42f6-823a-bc3e783ceab4
md"""

### System utilities
"""

# ╔═╡ 3d81d403-fbcc-4bc9-85f7-c1f257dc7695
function ingredients(path::String)
# this is from the Julia source code (evalfile in base/loading.jl)
# but with the modification that it returns the module instead of the last object
	name = Symbol(basename(path))
m = Module(name)
Core.eval(m,
        Expr(:toplevel,
             :(eval(x) = $(Expr(:core, :eval))($name, x)),
             :(include(x) = $(Expr(:top, :include))($name, x)),
             :(include(mapexpr::Function, x) = $(Expr(:top, :include))(mapexpr, $name, x)),
             :(include($path))))
m
end

# ╔═╡ Cell order:
# ╠═97130b40-e41e-46c1-ac48-7fc734e05d95
# ╠═82c494e0-b75e-474e-b08c-cde3791a121e
# ╠═8aaabf9e-08f7-11ed-1880-e798d62d8295
# ╟─03bfbaeb-f6f2-4a8b-8815-66d293cbbb66
# ╟─cdc50171-b288-40b6-9d0d-9511901218e0
# ╟─f8dbfa77-c88d-42f6-bc2a-600b49f8f98d
# ╟─8fe43246-a006-4ff9-bb64-3601d4e03938
# ╟─6e7e284e-88a1-47f9-9da8-4a79797d8ebd
# ╟─03d81ddf-105e-4ea6-a882-e1b40b7ecbfc
# ╟─98ab81f8-f630-4d7e-9624-829ec179dd10
# ╟─7913b885-4f25-40ae-ab84-3b5a4d50967c
# ╟─44b35007-912d-49e3-90bb-09aa1360cbe9
# ╟─14160a98-07a0-4efc-9e71-9aab36ed01b6
# ╟─99a11ca4-b453-4e0b-a737-d9733b0e8c59
# ╟─718fbf98-1725-429b-9a1c-609cb38430ac
# ╟─65597330-86b4-4d17-abae-5df638603cfc
# ╟─ca22eb02-282a-412b-9297-15a1e379a7b2
# ╟─5f913a6e-1798-42da-8fbe-3dad64fc7bb7
# ╟─bd05aad2-f55a-424d-a467-f17d32286017
# ╟─aa55ee52-f3d0-4425-ba88-cafc4743ecf4
# ╟─d13e0a1a-8eb4-49eb-a8b7-b732cf9ffe31
# ╟─f547f778-8a74-4923-9373-cb96b2254c9b
# ╟─420900c2-74f9-4972-bddd-982e766fc22f
# ╟─30970210-d5f6-4a1b-b36a-601789b0fd8b
# ╟─5e7771ee-372b-4fa9-8f8a-888c6e46796c
# ╟─5e61e1c2-a1d7-49a1-bd3f-649d4c6041cf
# ╠═4f6f8003-b9f0-4aba-a42c-c79184bec77d
# ╠═69920bad-6d3e-4e7a-a539-f08bc63460c4
# ╠═e1c6ef24-fe3a-41b5-95a7-f0ddb7b3bf6c
# ╠═59f1b12f-3224-4ea1-8f55-98f228cdcfb6
# ╠═7397c3e4-d057-42d2-947f-800dab52ee6f
# ╠═59d569ec-c2a4-485b-88ac-86495ba1d6ab
# ╠═41e79a3a-fea2-449b-b42c-50ba020e5482
# ╟─aa5c46fc-ccf6-49e1-b20b-7714bf2c7630
# ╟─ed0c4b7e-cf34-4b20-9100-43ba5428e096
# ╟─0b5acafe-fc8f-42f6-823a-bc3e783ceab4
# ╠═3d81d403-fbcc-4bc9-85f7-c1f257dc7695
