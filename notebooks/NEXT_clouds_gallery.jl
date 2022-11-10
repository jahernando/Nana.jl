### A Pluto.jl notebook ###
# v0.19.14

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

# ╔═╡ 82c494e0-b75e-474e-b08c-cde3791a121e
begin
using Pkg
Pkg.activate("/Users/hernando/work/investigacion/NEXT/software/julias/Nana")
import Nana as na
end

# ╔═╡ 8aaabf9e-08f7-11ed-1880-e798d62d8295
begin
using HDF5
using LinearAlgebra
using DataFrames
using Plots
using Statistics
import StatsBase as SB
import Clouds as cl
import Graphs as GG
import PlutoUI as PUI
import Clouds as cl
import GraphPlot  as GP
import Colors as Colors
#using PlotlyJS
end

# ╔═╡ cdc50171-b288-40b6-9d0d-9511901218e0
md"""

## Description


This NB displays: 

	1) NEXT-100 MC events using Beersheba-reconstructed and MC hits.

    2) Clouds Cells, nodes and the Spline

	3) A summary of the umber of nodes and blobs


J.A. Hernando,

Santiago de Compostela, October 2022

---
"""

# ╔═╡ 03bfbaeb-f6f2-4a8b-8815-66d293cbbb66
PUI.TableOfContents(title="NEXT-100 Event Gallery", indent=true)

# ╔═╡ f8dbfa77-c88d-42f6-bc2a-600b49f8f98d
# NEXT-100 data
begin
datadir   = "/Users/hernando/work/investigacion/NEXT/data/NEXT100/"
filename  = "bb0nu/v2/beersheba_fixed_label_1_0nubb.h5"
end;

# ╔═╡ 6e7e284e-88a1-47f9-9da8-4a79797d8ebd
df, mc, steps = na.load_data(string(datadir, filename));

# ╔═╡ 03d81ddf-105e-4ea6-a882-e1b40b7ecbfc
begin
events = sort(collect(Set(df.dataset_id)))

md"""


Datafile          : **$(filename)**


Number of events : **$(length(events))**

Total number of entries in Beersheba's hit table: $(size(df)[1])

Voxel size:

| x (mm) | y (mm) | z (mm)|
| :-- | :-- | :-- |
| $(steps[1])| $(steps[2]) | $(steps[3]) |


"""
end

# ╔═╡ 6da8778c-1b25-414d-b269-2d85978fb72d
md"""

-------

### Configuration
"""

# ╔═╡ beb55efd-4b0e-49a3-9dc8-3bc0de573fbf
begin
	
breco = @bind reco PUI.Select([:reco, :mc])

md"""
Select reco or mc voxels : $(breco)
"""
end

# ╔═╡ 4611084e-712c-4b50-8ab9-bc3be16bc231
begin
minv, maxv = 0.0, 0.02
	
btrheshold = @bind ethreshold PUI.Slider(minv:0.0001:maxv)

md"""

Beersheba threshold : $(btrheshold)
"""
end

# ╔═╡ 6fb542df-43a5-46f6-bf89-e9e88bd563aa
md"""

Energy Beersheba's hits min energy threshold : $(ethreshold)

"""

# ╔═╡ 76dc4710-0771-41ff-a69d-1a23be120502
begin

btype = @bind type PUI.Select([:clouds, :paulina])

md"""

Nodes via Clouds or Paulina: $(btype)
"""
end

# ╔═╡ 897c8506-b03e-44fc-9e37-884e8d28684c
begin

bnsteps = @bind nsteps PUI.Select([1, 2, 3, 4])

md"""
Select factor of steps: $(bnsteps)
"""
end

# ╔═╡ 37b9f270-57b5-4a43-9957-bc0ff4f0bd98
cellnode = type === :paulina;

# ╔═╡ 44b35007-912d-49e3-90bb-09aa1360cbe9
begin

bevent = @bind event PUI.Select(events)

md"""

-------

Select event number : $(bevent)
"""
end

# ╔═╡ 14160a98-07a0-4efc-9e71-9aab36ed01b6
begin
idf = na.get_event(df, event)
imc = na.get_event(mc, event)

md"""

----------

##### Event $(event) has $(size(idf)[1]) Beersheba and $(size(imc)[1]) MC hits.

----
"""
end

# ╔═╡ a1791e77-fa3f-4b2a-b582-71e1a6737f06
md"""

## Run Event $(event), $(reco), $(type)

#### Configuration : *$(reco), cellnode = $(cellnode), nsteps = $(nsteps)*

"""

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

# ╔═╡ 34f80b89-82ae-442f-a148-69885ef7b304
#draw_mchits(imc)

# ╔═╡ 9e1ee657-a53e-4c1e-933e-d103562db116
md"""

### Diplay MC and Beersheba Hits

"""

# ╔═╡ 4853fff3-dd52-4b6a-b7d1-44fef36b1cef
begin
xdf = reco == :reco ? idf : imc
xcoors    = Tuple((xdf[!, :x], xdf[!, :y], xdf[!, :z]))
xcontents = copy(xdf[!, :energy])
xsteps    = Tuple(nsteps * steps)
end;

# ╔═╡ 561c837c-0c9f-4291-9ae9-5d174a8e0b48
histogram(xcontents, nbins = 100)

# ╔═╡ bf401c22-da4e-4b72-8fd2-f68c54a11978
md"""

Total energy $(sum(xcontents)) MeV

Number of empty fixed voxels $(sum(xcontents .== 0.0))

"""

# ╔═╡ 2ce7a34b-ca5c-4e12-a700-3cb49d563175
begin
xcontents[xcontents .== 0.0] .= 1e-3
md"""
Total energy $(sum(xcontents)) MeV (after addition of 1 keV per empty cell)
"""
end

# ╔═╡ e46aef6a-3618-4b2a-8a83-8501e75c3e52
begin
sel       = xcontents .>= ethreshold
ycontents = xcontents[sel]
ycoors    = Tuple((x[sel] for x in xcoors))
end;

# ╔═╡ 0b3712ad-e3a9-405d-b173-a54910f2d15a
xcl, xnd, xspine, xedges = cl.clouds(ycoors, ycontents, xsteps; cellnode = cellnode);

# ╔═╡ edf5589e-f940-4e4c-bba0-c7dfe592f54d
md"""

### Display the Cells, Nodes and Spine

"""

# ╔═╡ 8cd9d7d2-0e6c-4f1e-a574-9ddc59137412
begin
clabel = na.label_cell(xedges, xcl.cells, xcoors, xdf.segclass)
nlabel = na.label_node(xcl.node, clabel)
end;

# ╔═╡ cd9e43f2-a81f-41e0-b726-9f58d8f3592b
md"""

#### Spine

Nodes at the extremes: $(xspine.extremes_maxcontents)
"""

# ╔═╡ eef8dfc6-b4a6-4b9a-ba3d-a65dcf093cdb
begin
ncells = length(xcl.cells)
nnodes = length(xnd.size)
nblobs = sum(nlabel .== 3)
nother = sum(nlabel .== 1)
nextr  = sum(xnd.extremes .== 1)
nbextr = sum(nlabel[xnd.extremes .== 1] .== 3)
sum(nlabel.==3)
md"""

### Event Summary

| cells | nodes | nodes-blob | nodes-other | extremes | extremes-blob| 
| :-- | :-- | :-- | :-- | :-- | :-- |
| $(ncells) | $(nnodes) | $(nblobs) | $(nother) |$(nextr) | $(nbextr) |
"""
end

# ╔═╡ 3d7b4885-8f3d-4449-aa10-5b481a1f6fba
md"""

voxel sizes: $(xsteps) mm

"""

# ╔═╡ 5ad5a6a2-2967-4448-8d3e-4bb1c7015a42
md"""

## DEV
"""

# ╔═╡ 5e2b6643-b5ba-4b7b-8aa7-7c9ec24c5df9
tree = GG.bfs_tree(xspine.spine, xspine.extremes_maxcontents[1][1]);

# ╔═╡ ce50597b-7959-40b2-acdf-773670ca96ff
md"""
### Display Non-Cycle Spine

Nodes at the extremes: $(xspine.extremes_maxcontents)

"""

# ╔═╡ a024311d-e35f-472d-94c7-2c894f1d1dca
md"""

## Code

"""

# ╔═╡ d543b6b1-dcee-4dd5-8613-2105ecb74888
"""
linear scale a vector of values between a minimum and a maximum

Parameters:
	var : Vector{Real}
	emin: Real
	emax: Real

Return:
	vvar: Vector{Real}

"""
function vscale(var, emin = 1.5, emax = 4.)
	vvar = (var .- minimum(var)) ./(maximum(var) - minimum(var)) .*(emax-emin) .+ emin
	return vvar
end

# ╔═╡ e04b9d5c-64d2-4da3-90fb-e79bbf169cc6
function draw_hits(df, mc; title = "")
	theme(:dark)
	sc3 = scatter3d(df.x, df.y, df.z, marker_z = df.energy,
	markersize = 2., #vscale(df.energy),
	markertype = :circle, label = false, alpha = 0.4, c = :inferno,
	xlabel = "x (mm)", ylabel = "y (mm)", zlabel = "z (mm)", title = title)
	scatter3d!(mc.x, mc.y, mc.z, alpha = 0.4, markersize = 0.5)
	plot(sc3, size = (500, 400))
end

# ╔═╡ 96dc3227-bdc4-42b9-a0e7-c8f25995058b
draw_hits(idf, imc; title =  "hits")

# ╔═╡ 709da913-34aa-4e87-9fa4-d85eb902cec0
function draw(rc, mc; title = "")
	theme(:dark)
	sc3 = scatter3d(rc.coors..., marker_z = rc.contents, 
		markersize = vscale(rc.contents), markertype = :circle, label = false, alpha = 0.4, c = :inferno,
		xlabel = "x (mm)", ylabel = "y (mm)", zlabel = "z (mm)", title = title)
	scatter3d!(mc.x, mc.y, mc.z, alpha = 0.4, markersize = 0.5)
	plot(sc3, size = (500, 400))
end

# ╔═╡ 29a5823d-46a7-42c6-b573-43793f53a365
draw(xcl, imc; title = "Cells")

# ╔═╡ 2b269e4f-a829-4f82-a8c2-b32b38139e0f
draw(xnd, imc; title = "Nodes")

# ╔═╡ 8afff426-08e7-49e6-9345-ec727d9110b8
function draw_mchits(mc, title = "MC hits")
	theme(:dark)
	mcseg = cgrad(:matter, 3, categorical = true)
	sc3 = scatter3d(mc.x, mc.y, mc.z, alpha = 0.4, marker_z = mc.segclass, c = mcseg, markersize = 0.5)
	plot(sc3, size = (500, 400))
end


# ╔═╡ d8659fa1-c140-43c9-b2d2-d4ee4c73534e
begin
function draw_graph(graph, nd, labels)
	nvertices = GG.nv(graph)
	reds      = 0.3 .* labels
	nodefillc = [Colors.RGBA(red, red , 1.0-red, 0.6) for red in reds]
	GP.gplot(graph, nodelabel=1:nvertices, nodesize = nd.contents,
		nodefillc = nodefillc)
end

function draw_graph(graph, nd)
	nvertices = GG.nv(graph)
	GP.gplot(graph, nodelabel=1:nvertices, nodesize = nd.contents)
end
end

# ╔═╡ 9d29ba39-ffb9-474d-b34f-0c493e28b09b
draw_graph(xspine.spine, xnd, nlabel)

# ╔═╡ 3e6c7109-ada7-4e2b-b78e-b77f352ba407
draw_graph(tree, xnd, nlabel)

# ╔═╡ aa5c46fc-ccf6-49e1-b20b-7714bf2c7630


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

# ╔═╡ Cell order:
# ╟─82c494e0-b75e-474e-b08c-cde3791a121e
# ╟─cdc50171-b288-40b6-9d0d-9511901218e0
# ╟─8aaabf9e-08f7-11ed-1880-e798d62d8295
# ╟─03bfbaeb-f6f2-4a8b-8815-66d293cbbb66
# ╟─f8dbfa77-c88d-42f6-bc2a-600b49f8f98d
# ╟─6e7e284e-88a1-47f9-9da8-4a79797d8ebd
# ╟─03d81ddf-105e-4ea6-a882-e1b40b7ecbfc
# ╟─6da8778c-1b25-414d-b269-2d85978fb72d
# ╟─beb55efd-4b0e-49a3-9dc8-3bc0de573fbf
# ╟─4611084e-712c-4b50-8ab9-bc3be16bc231
# ╟─6fb542df-43a5-46f6-bf89-e9e88bd563aa
# ╟─76dc4710-0771-41ff-a69d-1a23be120502
# ╟─897c8506-b03e-44fc-9e37-884e8d28684c
# ╟─37b9f270-57b5-4a43-9957-bc0ff4f0bd98
# ╟─44b35007-912d-49e3-90bb-09aa1360cbe9
# ╟─14160a98-07a0-4efc-9e71-9aab36ed01b6
# ╟─a1791e77-fa3f-4b2a-b582-71e1a6737f06
# ╟─718fbf98-1725-429b-9a1c-609cb38430ac
# ╟─65597330-86b4-4d17-abae-5df638603cfc
# ╟─34f80b89-82ae-442f-a148-69885ef7b304
# ╟─9e1ee657-a53e-4c1e-933e-d103562db116
# ╟─96dc3227-bdc4-42b9-a0e7-c8f25995058b
# ╟─4853fff3-dd52-4b6a-b7d1-44fef36b1cef
# ╠═561c837c-0c9f-4291-9ae9-5d174a8e0b48
# ╟─bf401c22-da4e-4b72-8fd2-f68c54a11978
# ╟─2ce7a34b-ca5c-4e12-a700-3cb49d563175
# ╠═e46aef6a-3618-4b2a-8a83-8501e75c3e52
# ╠═0b3712ad-e3a9-405d-b173-a54910f2d15a
# ╟─edf5589e-f940-4e4c-bba0-c7dfe592f54d
# ╠═29a5823d-46a7-42c6-b573-43793f53a365
# ╠═2b269e4f-a829-4f82-a8c2-b32b38139e0f
# ╟─8cd9d7d2-0e6c-4f1e-a574-9ddc59137412
# ╟─cd9e43f2-a81f-41e0-b726-9f58d8f3592b
# ╟─9d29ba39-ffb9-474d-b34f-0c493e28b09b
# ╟─eef8dfc6-b4a6-4b9a-ba3d-a65dcf093cdb
# ╟─3d7b4885-8f3d-4449-aa10-5b481a1f6fba
# ╟─5ad5a6a2-2967-4448-8d3e-4bb1c7015a42
# ╟─5e2b6643-b5ba-4b7b-8aa7-7c9ec24c5df9
# ╟─ce50597b-7959-40b2-acdf-773670ca96ff
# ╟─3e6c7109-ada7-4e2b-b78e-b77f352ba407
# ╟─a024311d-e35f-472d-94c7-2c894f1d1dca
# ╠═d543b6b1-dcee-4dd5-8613-2105ecb74888
# ╠═e04b9d5c-64d2-4da3-90fb-e79bbf169cc6
# ╠═709da913-34aa-4e87-9fa4-d85eb902cec0
# ╠═8afff426-08e7-49e6-9345-ec727d9110b8
# ╠═d8659fa1-c140-43c9-b2d2-d4ee4c73534e
# ╟─aa5c46fc-ccf6-49e1-b20b-7714bf2c7630
# ╟─ed0c4b7e-cf34-4b20-9100-43ba5428e096
