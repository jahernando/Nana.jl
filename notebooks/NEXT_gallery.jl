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

# ╔═╡ 03bfbaeb-f6f2-4a8b-8815-66d293cbbb66
PUI.TableOfContents(title="NEXT-100 Event Gallery", indent=true)

# ╔═╡ cdc50171-b288-40b6-9d0d-9511901218e0
md"""

## Description


This NB displays NEXT-100 MC events using Beersheba-reconstructed and MC hits.

J.A. Hernando,

Santiago de Compostela, October 2022

---
"""

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

# ╔═╡ 44b35007-912d-49e3-90bb-09aa1360cbe9
begin

bevent = @bind event PUI.Select(events)

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
md"""

## Event Display

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

# ╔═╡ d85c31dc-4768-407d-bc5b-4b3992002b09
md"""

## Reco Clouds

"""

# ╔═╡ 4853fff3-dd52-4b6a-b7d1-44fef36b1cef
begin
xcoors    = Tuple((idf[!, :x], idf[!, :y], idf[!, :z]))
xcontents = copy(idf[!, :energy])
xsteps    = Tuple(steps)
end;

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

# ╔═╡ c44c6d1c-ae49-45ba-836a-e5f67994749d
xcl, xnd, graph, edges = cl.clouds(xcoors, xcontents, xsteps);

# ╔═╡ 8cd9d7d2-0e6c-4f1e-a574-9ddc59137412
begin
clabel = na.label_cell(edges, xcl.cells, xcoors, idf.segclass)
nlabel = na.label_node(xcl.node, clabel)
end

# ╔═╡ a1e1e2a5-f524-4581-b4d3-c59160ca5608
md"""

## MC Clouds

"""

# ╔═╡ c662109d-8c93-44ce-837e-22c73c64d800
begin
mccoors    = Tuple((imc[!, :x], imc[!, :y], imc[!, :z]))
mccontents = copy(imc[!, :energy])
mcsteps    = (2.0, 2.0, 2.0);
end;

# ╔═╡ c793bb48-6fed-46c6-a68a-0d26c5c57d50
mccl, mcnd, mcgraph, mcedges = cl.clouds(mccoors, mccontents, xsteps);

# ╔═╡ 3d7b4885-8f3d-4449-aa10-5b481a1f6fba
md"""

MC steps $(xsteps) mm in each direction

"""

# ╔═╡ c075866b-0b1e-44ef-ac7a-718d1c1d1a4b
begin
mc_clabel = na.label_cell(mcedges, mccl.cells, mccoors, imc.segclass)
mc_nlabel = na.label_node(mccl.node, mc_clabel)
end

# ╔═╡ 348bc84f-ad16-4cec-a3df-e1ad7063754e
md"""
MC graph
"""

# ╔═╡ 5ad5a6a2-2967-4448-8d3e-4bb1c7015a42
md"""

## DEV
"""

# ╔═╡ 476e7039-1199-4aa8-a7f3-7e8a14215ade
xnd

# ╔═╡ 5e2b6643-b5ba-4b7b-8aa7-7c9ec24c5df9
tree = GG.bfs_tree(graph, extr[1][1])

# ╔═╡ a9e159a0-9aa1-41d0-93e8-f565bbe4565f
begin
mcextrs, mcdists = na.graph_extremes(mcgraph)
mcextr, _, _     = na.graph_extremes_maxcontents(mcgraph, mcnd)
end

# ╔═╡ 1f8fe651-2ee4-4b69-a21c-e7773b0f0283
xx = GG.bfs_tree(mcgraph, mcextr[1][1])

# ╔═╡ ad8a7d63-401a-40c9-9520-7306ae976f92
begin
mcnextr = zeros(length(mcnd.contents))
for (i, j) in mcextr
	mcnextr[i] = 1
	mcnextr[j] = 1
end
end

# ╔═╡ b9ea99e0-4d9f-46a6-9fb1-a342934487fb
mcnextr

# ╔═╡ b40e10f4-01e8-4c98-a07d-689df32a5706
mcnd

# ╔═╡ fc5ee335-bb67-4de2-9448-a286edf8c126
function dfnodes(xnd, nlabel)

	dd = Dict()
	dd[:contents] = xnd.contents
	dd[:size]     = xnd.size
	for i in 1:length(xnd.coors)
		dd[Symbol(string("x", i, "_mean"))] = xnd.coors[i]
		dd[Symbol(string("x", i, "_std")) ] = xnd.coors_std[i]
	end
	
	dd[:maxgrad] = xnd.maxgrad
	dd[:maxlap]  = xnd.maxlap
	dd[:minlap]  = xnd.minlap
	dd[:maxcur]  = xnd.curmax
	dd[:mincur]  = xnd.curmin

	dd[:nedges] = xnd.nedges
	dd[:cloud]  = xnd.cloud
	dd[:ecc]    = xnd.ecc
	dd[:extreme] = xnd.extremes

	dd[:label]   = nlabel
	
	return DataFrame(dd)
	
end

# ╔═╡ 63bb90d4-178a-490a-ae8c-8d21e6c10d1d
dfnd = dfnodes(mcnd, mc_nlabel)

# ╔═╡ dc39374d-52a4-460d-83d8-4701e74b5454
dfnd[!, :event] .= event

# ╔═╡ 7b881784-ecaa-4a2a-99be-ee85a2dd9cec
dfnd

# ╔═╡ ea2c540a-3c15-4689-95fb-a4d896ffb20a


# ╔═╡ 5c7ba791-9302-41e9-b7f7-b5f0237b581a
string("x", 1, "_mean")

# ╔═╡ ffbb28d9-bb4c-4708-8fd8-3633124173a5
mcnd

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
	plot(sc3, size = (400, 300))
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
	plot(sc3, size = (400, 300))
end

# ╔═╡ 29a5823d-46a7-42c6-b573-43793f53a365
draw(xcl, imc; title = "reco cells")

# ╔═╡ 2b269e4f-a829-4f82-a8c2-b32b38139e0f
draw(xnd, imc; title = "reco nodes")

# ╔═╡ 83853cc3-5ec9-4f37-a27e-c7a8563237b8
draw(mccl, imc; title =  "mc-cells")

# ╔═╡ 1338cf59-de08-4d7e-a192-b8b2c42efb24
draw(mcnd, imc; title = "mc-nodes")

# ╔═╡ 8afff426-08e7-49e6-9345-ec727d9110b8
function draw_mchits(mc, title = "MC hits")
	theme(:dark)
	mcseg = cgrad(:matter, 3, categorical = true)
	sc3 = scatter3d(mc.x, mc.y, mc.z, alpha = 0.4, marker_z = mc.segclass, c = mcseg, markersize = 0.5)
	plot(sc3, size = (400, 300))
end


# ╔═╡ 34f80b89-82ae-442f-a148-69885ef7b304
draw_mchits(imc)

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
draw_graph(graph.spine, xnd, nlabel)

# ╔═╡ 05848241-42ea-4dab-aff2-a1cdecffaf8d
draw_graph(mcgraph.spine, mcnd, mc_nlabel)

# ╔═╡ 3e6c7109-ada7-4e2b-b78e-b77f352ba407
draw_graph(tree, xnd, nlabel)

# ╔═╡ c1ce0881-420b-4a55-a091-26a5d62f051c
draw_graph(xx, mcnd, mc_nlabel)

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
# ╠═82c494e0-b75e-474e-b08c-cde3791a121e
# ╠═8aaabf9e-08f7-11ed-1880-e798d62d8295
# ╟─03bfbaeb-f6f2-4a8b-8815-66d293cbbb66
# ╟─cdc50171-b288-40b6-9d0d-9511901218e0
# ╠═f8dbfa77-c88d-42f6-bc2a-600b49f8f98d
# ╠═6e7e284e-88a1-47f9-9da8-4a79797d8ebd
# ╠═03d81ddf-105e-4ea6-a882-e1b40b7ecbfc
# ╟─44b35007-912d-49e3-90bb-09aa1360cbe9
# ╠═14160a98-07a0-4efc-9e71-9aab36ed01b6
# ╟─99a11ca4-b453-4e0b-a737-d9733b0e8c59
# ╟─718fbf98-1725-429b-9a1c-609cb38430ac
# ╟─65597330-86b4-4d17-abae-5df638603cfc
# ╟─d85c31dc-4768-407d-bc5b-4b3992002b09
# ╠═34f80b89-82ae-442f-a148-69885ef7b304
# ╟─96dc3227-bdc4-42b9-a0e7-c8f25995058b
# ╠═4853fff3-dd52-4b6a-b7d1-44fef36b1cef
# ╟─bf401c22-da4e-4b72-8fd2-f68c54a11978
# ╟─2ce7a34b-ca5c-4e12-a700-3cb49d563175
# ╠═c44c6d1c-ae49-45ba-836a-e5f67994749d
# ╠═29a5823d-46a7-42c6-b573-43793f53a365
# ╠═2b269e4f-a829-4f82-a8c2-b32b38139e0f
# ╠═8cd9d7d2-0e6c-4f1e-a574-9ddc59137412
# ╠═9d29ba39-ffb9-474d-b34f-0c493e28b09b
# ╟─a1e1e2a5-f524-4581-b4d3-c59160ca5608
# ╠═c662109d-8c93-44ce-837e-22c73c64d800
# ╟─c793bb48-6fed-46c6-a68a-0d26c5c57d50
# ╟─3d7b4885-8f3d-4449-aa10-5b481a1f6fba
# ╟─83853cc3-5ec9-4f37-a27e-c7a8563237b8
# ╠═1338cf59-de08-4d7e-a192-b8b2c42efb24
# ╠═c075866b-0b1e-44ef-ac7a-718d1c1d1a4b
# ╠═348bc84f-ad16-4cec-a3df-e1ad7063754e
# ╠═05848241-42ea-4dab-aff2-a1cdecffaf8d
# ╠═5ad5a6a2-2967-4448-8d3e-4bb1c7015a42
# ╠═476e7039-1199-4aa8-a7f3-7e8a14215ade
# ╠═5e2b6643-b5ba-4b7b-8aa7-7c9ec24c5df9
# ╠═3e6c7109-ada7-4e2b-b78e-b77f352ba407
# ╠═a9e159a0-9aa1-41d0-93e8-f565bbe4565f
# ╠═1f8fe651-2ee4-4b69-a21c-e7773b0f0283
# ╠═c1ce0881-420b-4a55-a091-26a5d62f051c
# ╠═ad8a7d63-401a-40c9-9520-7306ae976f92
# ╠═b9ea99e0-4d9f-46a6-9fb1-a342934487fb
# ╠═b40e10f4-01e8-4c98-a07d-689df32a5706
# ╠═fc5ee335-bb67-4de2-9448-a286edf8c126
# ╠═63bb90d4-178a-490a-ae8c-8d21e6c10d1d
# ╠═dc39374d-52a4-460d-83d8-4701e74b5454
# ╠═7b881784-ecaa-4a2a-99be-ee85a2dd9cec
# ╠═ea2c540a-3c15-4689-95fb-a4d896ffb20a
# ╠═5c7ba791-9302-41e9-b7f7-b5f0237b581a
# ╠═ffbb28d9-bb4c-4708-8fd8-3633124173a5
# ╟─a024311d-e35f-472d-94c7-2c894f1d1dca
# ╠═d543b6b1-dcee-4dd5-8613-2105ecb74888
# ╠═e04b9d5c-64d2-4da3-90fb-e79bbf169cc6
# ╠═709da913-34aa-4e87-9fa4-d85eb902cec0
# ╠═8afff426-08e7-49e6-9345-ec727d9110b8
# ╠═d8659fa1-c140-43c9-b2d2-d4ee4c73534e
# ╟─aa5c46fc-ccf6-49e1-b20b-7714bf2c7630
# ╟─ed0c4b7e-cf34-4b20-9100-43ba5428e096
