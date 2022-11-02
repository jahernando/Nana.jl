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

# ╔═╡ 348bc84f-ad16-4cec-a3df-e1ad7063754e
md"""
MC graph
"""

# ╔═╡ 5ad5a6a2-2967-4448-8d3e-4bb1c7015a42
md"""

## DEV
"""

# ╔═╡ fb26d9d2-6496-4388-a8fa-8f8157ecf6a3
xcl.cells

# ╔═╡ f359135b-5a9c-4640-9026-92f65b1fd525
begin
function bins_ids(coors, edges)
    h = SB.fit(SB.Histogram, coors, edges)
    ucoors = reduce(hcat, coors)
    bin_ids = map(xi -> SB.binindex(h, Tuple(xi)), eachrow(ucoors))
    return bin_ids
end

function _label(labs)
	for ll in (3, 1, 2, 6, 4, 5, 7)
		if ll in labs
			return ll
		end
	end
end

function cell_label(edges, cells, coors, labels)
	bins = bins_ids(coors, edges)
	xlabel = zeros(Int64, length(cells))
	for (i, icell) in enumerate(cells)
		sel = [b == icell for b in bins]
		labs = Set(labels[sel])
		xlabel[i] = _label(labs)
	end
	return xlabel
end

function node_label(cell_node, cell_label)
	nnodes = maximum(cell_node)
	xlabel = zeros(Int64, nnodes)
	for i in 1:nnodes
		sel = cell_node .== i
		labs = Set(cell_label[sel])
		xlabel[i] = _label(labs)
	end
	return xlabel
end
end

# ╔═╡ 8cd9d7d2-0e6c-4f1e-a574-9ddc59137412
begin
clabel = cell_label(edges, xcl.cells, xcoors, idf.segclass)
nlabel = node_label(xcl.node, clabel)
end

# ╔═╡ c075866b-0b1e-44ef-ac7a-718d1c1d1a4b
begin
mc_clabel = cell_label(mcedges, mccl.cells, mccoors, imc.segclass)
mc_nlabel = node_label(mccl.node, mc_clabel)
end

# ╔═╡ 8e77df59-8a20-4972-a657-f73f25a45d27
idbis = bins_ids(xcoors, edges)

# ╔═╡ d059f2b3-41a4-4824-90b7-d062e2df92a4
length(idf.segclass)

# ╔═╡ 4d506bb2-52d8-4eeb-b0d5-04ade6f42d06
length(xcoors[1])

# ╔═╡ 55d1569e-94f3-43ff-9db4-3531d3964625
xnd

# ╔═╡ 4fb921cd-9a6b-4514-aeed-c0d0efa72ab0
labs = Set(idf.segclass)

# ╔═╡ 70d073c5-7186-43a8-81ee-a816eab2c75e
xcl.cells

# ╔═╡ 2931bf98-f582-44ad-a815-dbb2616d8f17
bins = bins_ids(xcoors, edges)

# ╔═╡ 580ed7e6-8633-4241-9a54-24d614238d85
xlabels = get_label(xcl.cells, bins, idf.segclass)

# ╔═╡ 5e2aa018-24c0-45f1-8722-48edc915caac


# ╔═╡ 4334854e-5e8e-4dab-83e9-5338650b5be1
Set(xlabels)

# ╔═╡ ee14b333-1c4f-4dc2-9cfc-89ed7064b519
begin
sel = [b == xcl.cells[1] for b in bins]
labels = Set(idf.segclass[sel])
end

# ╔═╡ 476e7039-1199-4aa8-a7f3-7e8a14215ade


# ╔═╡ ff69699e-0929-44eb-8fb3-f4ea252a2d91
xcl.cells[1] .== bins

# ╔═╡ 4a469534-0380-41e1-bf6f-805918d1f556
xcl.cells[1]

# ╔═╡ 944df99e-f941-4454-9ca7-e4c2938d5148
get_label(xcl.cells, bins, idf.segclass)

# ╔═╡ 333173ae-6bae-4901-988a-289ca89987f9
length(xcl.cells)

# ╔═╡ 2224f228-be41-427f-8c90-995acf5a4041
length(idf.segclass)

# ╔═╡ 72ca32c1-d67b-472b-a92d-55a741b0913b
_label(labs)

# ╔═╡ ec6d446f-dd80-4832-a833-21e2cfe7ef8e


# ╔═╡ 8622dc4d-9f9f-42fe-94a2-9845991ad602
1 in labs

# ╔═╡ 216ce100-7446-4029-bfea-ef95ab4d3030
begin
xextremes, xdists = na.graph_extremes(graph)
blobs = na.extremes_maxcontents(xextremes, xnd)
end

# ╔═╡ a9e159a0-9aa1-41d0-93e8-f565bbe4565f
begin
mcextremes, mcdists = na.graph_extremes(mcgraph)
mcblobs = na.extremes_maxcontents(mcextremes, mcnd)
end

# ╔═╡ 10f15bcd-6fe3-45ed-a2bb-3e94128acade
imc

# ╔═╡ e5253d55-6f15-4ed5-8fae-ae07da3d3b3a
sum(idf.segclass .== 3)

# ╔═╡ 6930c1ae-d94d-4d53-8ccd-1ba66a289f68
begin
#x = rand(1000,3)
h = SB.fit(SB.Histogram, xcoors, edges)
ucoors = reduce(hcat, xcoors)
bin_ids = map(xi -> SB.binindex(h, Tuple(xi)), eachrow(ucoors))
end

# ╔═╡ 79a78fde-fcff-4e93-b0db-27ff9c4dc9a7
hcat(xcoors)

# ╔═╡ 8dcb907d-2f65-4ec6-81c8-20ada85d1087
xcoors

# ╔═╡ 3b884611-79c2-4e7a-bc1e-3010b02703c4
bin_ids = [SB.binindex(h, Tuple(xi)) for 

# ╔═╡ a46aae36-52bd-40b7-8a82-05e95891466b
Set(idf.segclass)

# ╔═╡ 4adb4ecd-8686-45d1-aedc-247c6a5594f8
histo    = SB.fit(SB.Histogram, xcoors, SB.weights(idf.segclass), edges)

# ╔═╡ 0c9f3947-a399-47ff-8db9-e1fe04427e6e
xcells = findall(x -> x .>= 1, histo.weights)

# ╔═╡ f5c8b2d3-4c1f-40b4-a4d0-df630cd47c2f
CartesianIndex.(xcl.cells)

# ╔═╡ 9fe7aaad-8add-4a80-bb59-aea0831a2f2b
segclass = [histo.weights[index] for index in CartesianIndex.(xcl.cells)]

# ╔═╡ 1f8fe651-2ee4-4b69-a21c-e7773b0f0283
xx = GG.bfs_tree(mcgraph, 1)

# ╔═╡ 3132d75e-84a2-4ebf-9013-5b1274eb4945
GG.eccentricity(xx, 1)

# ╔═╡ 3dac5126-3b9b-4dfb-96f8-5c493cd6ef37
begin
dtrees = Dict()
for i in 1:GG.nv(mcgraph)
	itree = GG.bfs_tree(mcgraph, i)
	idis  = GG.eccentricity(itree, i)
	dtrees[i] = (idis, itree)
end
end

# ╔═╡ 2ac8ee4a-13b3-491e-b077-f724a3a66f4b
dtrees

# ╔═╡ 15496c82-efbc-45f0-acf0-e4c7795aa5f1
itree = dtrees[6][2]

# ╔═╡ e43c2ef5-7c79-4202-a2a2-fc06e30470f5
mc_nlabel

# ╔═╡ e5521919-7178-4d02-ba3d-2f316a68c5cb
GG.distance(itree)

# ╔═╡ aba9805a-7abe-41b9-89fb-7aadc8f5283b
GG.eccentricity(itree)

# ╔═╡ d79eec2a-548c-4d6d-8c0c-2706050e9d4f


# ╔═╡ bf6b5c1b-7302-4d49-8769-59fe00fa0349
function graph_extremes(graph)
	nvertices = GG.nv(graph)
	dists     = zeros(Int64, nvertices, nvertices)
	for i in 1:nvertices
		ds = GG.dijkstra_shortest_paths(graph, i)
		vv = copy(ds.dists)
		vv[vv .> nvertices] .= 0
		dists[i, :] = vv
	end
	extremes = findall(x -> x .== maximum(dists), dists)
	extremes = [(i, j) for (i, j) in Tuple.(extremes) if i < j]
	return extremes, dists
end

# ╔═╡ 16b6403d-18f0-4bd4-9e61-b3da5e3bffb2
graph_extremes(mcgraph)

# ╔═╡ 0af805f0-b6ea-4f1e-b9ae-ab2cb6474f51
begin
xextr, dists = graph_extremes(graph)
end

# ╔═╡ 754c166f-e572-4833-91c0-d61be5b363b4
function extremes_maxcontents(extremes, xnd)
	contents = [xnd.contents[i] + xnd.contents[j] for (i, j) in extremes]
	i = findall(x -> x .== maximum(contents), contents)
	return extremes[i]
end

# ╔═╡ 8b6f177d-1c26-4395-849c-88e3962480b6
extremes_maxcontents(xextr, xnd)

# ╔═╡ 7ef93664-64f1-48b8-9189-211e03526bcf


# ╔═╡ c61ee7b0-6685-4a9e-abd1-591afede39c5
xnd.contents[]

# ╔═╡ a9e47181-d183-4b57-be81-adfcc9052414
dist_gmin

# ╔═╡ 9e6637f1-170d-4efc-b8db-546bfe025072
begin
extremes = findall(x -> x .>= maximum(dist_gmin), dist_gmin)
extremes = [(i, j) for (i, j) in Tuple.(extremes) if i < j]
end

# ╔═╡ 711fcd84-d07f-4bab-abe6-fb93cff817bf
begin
ds = GG.dijkstra_shortest_paths(mcgraph, 1)
vv = copy(ds.dists)
vv[vv .> nvertices] .= 0
vv
end

# ╔═╡ 19b5a414-c674-4735-a007-964916ae19f6
nvertices

# ╔═╡ 13751dfe-80ca-43b2-ae2f-f63dbba26398
a = zeros(Int64, 3, 3)

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

# ╔═╡ a799583a-41df-4d1c-88d0-7442cde49c68
scatter3d(imc.x, imc.y, imc.z, marker_z = Int.(imc.segclass), c = cgrad(:matter, 3, categorical = true), alpha = 0.4,  markersize = 0.5)

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
	plot(sc3, size = (500, 500))
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
	plot(sc3, size = (500, 500))
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
draw_graph(graph, xnd, nlabel)

# ╔═╡ 05848241-42ea-4dab-aff2-a1cdecffaf8d
draw_graph(mcgraph, mcnd, mc_nlabel)

# ╔═╡ 9d6bbe57-73a3-401c-895c-0dcdd8a3d271
draw_graph(mcgraph, mcnd)

# ╔═╡ c1ce0881-420b-4a55-a091-26a5d62f051c
draw_graph(xx, mcnd)

# ╔═╡ 9d9fd9c2-7fa7-4dd3-9a9a-e7da36ae1ef9
draw_graph(dtrees[6][2], mcnd)

# ╔═╡ 9d2ea74b-20c7-4c72-91a6-d136ef1277c5
draw_graph(itree, mcnd, mc_nlabel)

# ╔═╡ 31959beb-7f01-4bb1-a2f5-1103070facf0
nodefillc = [Colors.RGBA(0.0,0.8,0.8, 1) for i in 1:10]

# ╔═╡ 09236d4c-440d-4884-b635-240f420bde23
nodefill = [Colors.RGBA(0.0, 0.8, 0.8, 1.0) for n in GG.nv(graph)]

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
# ╟─f8dbfa77-c88d-42f6-bc2a-600b49f8f98d
# ╟─6e7e284e-88a1-47f9-9da8-4a79797d8ebd
# ╟─03d81ddf-105e-4ea6-a882-e1b40b7ecbfc
# ╟─44b35007-912d-49e3-90bb-09aa1360cbe9
# ╟─14160a98-07a0-4efc-9e71-9aab36ed01b6
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
# ╟─c662109d-8c93-44ce-837e-22c73c64d800
# ╟─c793bb48-6fed-46c6-a68a-0d26c5c57d50
# ╟─3d7b4885-8f3d-4449-aa10-5b481a1f6fba
# ╟─83853cc3-5ec9-4f37-a27e-c7a8563237b8
# ╠═1338cf59-de08-4d7e-a192-b8b2c42efb24
# ╠═c075866b-0b1e-44ef-ac7a-718d1c1d1a4b
# ╠═348bc84f-ad16-4cec-a3df-e1ad7063754e
# ╠═05848241-42ea-4dab-aff2-a1cdecffaf8d
# ╠═5ad5a6a2-2967-4448-8d3e-4bb1c7015a42
# ╠═8e77df59-8a20-4972-a657-f73f25a45d27
# ╠═fb26d9d2-6496-4388-a8fa-8f8157ecf6a3
# ╠═f359135b-5a9c-4640-9026-92f65b1fd525
# ╠═d059f2b3-41a4-4824-90b7-d062e2df92a4
# ╠═4d506bb2-52d8-4eeb-b0d5-04ade6f42d06
# ╟─55d1569e-94f3-43ff-9db4-3531d3964625
# ╠═9d6bbe57-73a3-401c-895c-0dcdd8a3d271
# ╠═4fb921cd-9a6b-4514-aeed-c0d0efa72ab0
# ╠═70d073c5-7186-43a8-81ee-a816eab2c75e
# ╠═2931bf98-f582-44ad-a815-dbb2616d8f17
# ╠═580ed7e6-8633-4241-9a54-24d614238d85
# ╠═5e2aa018-24c0-45f1-8722-48edc915caac
# ╠═4334854e-5e8e-4dab-83e9-5338650b5be1
# ╠═ee14b333-1c4f-4dc2-9cfc-89ed7064b519
# ╠═476e7039-1199-4aa8-a7f3-7e8a14215ade
# ╠═ff69699e-0929-44eb-8fb3-f4ea252a2d91
# ╠═4a469534-0380-41e1-bf6f-805918d1f556
# ╠═944df99e-f941-4454-9ca7-e4c2938d5148
# ╠═333173ae-6bae-4901-988a-289ca89987f9
# ╠═2224f228-be41-427f-8c90-995acf5a4041
# ╠═72ca32c1-d67b-472b-a92d-55a741b0913b
# ╠═ec6d446f-dd80-4832-a833-21e2cfe7ef8e
# ╠═8622dc4d-9f9f-42fe-94a2-9845991ad602
# ╠═216ce100-7446-4029-bfea-ef95ab4d3030
# ╠═a9e159a0-9aa1-41d0-93e8-f565bbe4565f
# ╠═10f15bcd-6fe3-45ed-a2bb-3e94128acade
# ╠═e5253d55-6f15-4ed5-8fae-ae07da3d3b3a
# ╠═6930c1ae-d94d-4d53-8ccd-1ba66a289f68
# ╠═79a78fde-fcff-4e93-b0db-27ff9c4dc9a7
# ╠═8dcb907d-2f65-4ec6-81c8-20ada85d1087
# ╠═3b884611-79c2-4e7a-bc1e-3010b02703c4
# ╠═a46aae36-52bd-40b7-8a82-05e95891466b
# ╠═4adb4ecd-8686-45d1-aedc-247c6a5594f8
# ╠═0c9f3947-a399-47ff-8db9-e1fe04427e6e
# ╠═f5c8b2d3-4c1f-40b4-a4d0-df630cd47c2f
# ╠═9fe7aaad-8add-4a80-bb59-aea0831a2f2b
# ╠═1f8fe651-2ee4-4b69-a21c-e7773b0f0283
# ╠═c1ce0881-420b-4a55-a091-26a5d62f051c
# ╠═3132d75e-84a2-4ebf-9013-5b1274eb4945
# ╠═3dac5126-3b9b-4dfb-96f8-5c493cd6ef37
# ╠═2ac8ee4a-13b3-491e-b077-f724a3a66f4b
# ╠═9d9fd9c2-7fa7-4dd3-9a9a-e7da36ae1ef9
# ╠═15496c82-efbc-45f0-acf0-e4c7795aa5f1
# ╠═9d2ea74b-20c7-4c72-91a6-d136ef1277c5
# ╠═e43c2ef5-7c79-4202-a2a2-fc06e30470f5
# ╠═e5521919-7178-4d02-ba3d-2f316a68c5cb
# ╠═aba9805a-7abe-41b9-89fb-7aadc8f5283b
# ╠═d79eec2a-548c-4d6d-8c0c-2706050e9d4f
# ╠═bf6b5c1b-7302-4d49-8769-59fe00fa0349
# ╠═16b6403d-18f0-4bd4-9e61-b3da5e3bffb2
# ╠═0af805f0-b6ea-4f1e-b9ae-ab2cb6474f51
# ╠═754c166f-e572-4833-91c0-d61be5b363b4
# ╠═8b6f177d-1c26-4395-849c-88e3962480b6
# ╠═7ef93664-64f1-48b8-9189-211e03526bcf
# ╠═c61ee7b0-6685-4a9e-abd1-591afede39c5
# ╠═a9e47181-d183-4b57-be81-adfcc9052414
# ╠═9e6637f1-170d-4efc-b8db-546bfe025072
# ╠═711fcd84-d07f-4bab-abe6-fb93cff817bf
# ╠═19b5a414-c674-4735-a007-964916ae19f6
# ╠═13751dfe-80ca-43b2-ae2f-f63dbba26398
# ╟─a024311d-e35f-472d-94c7-2c894f1d1dca
# ╠═d543b6b1-dcee-4dd5-8613-2105ecb74888
# ╠═a799583a-41df-4d1c-88d0-7442cde49c68
# ╠═e04b9d5c-64d2-4da3-90fb-e79bbf169cc6
# ╠═709da913-34aa-4e87-9fa4-d85eb902cec0
# ╠═8afff426-08e7-49e6-9345-ec727d9110b8
# ╠═d8659fa1-c140-43c9-b2d2-d4ee4c73534e
# ╠═31959beb-7f01-4bb1-a2f5-1103070facf0
# ╠═09236d4c-440d-4884-b635-240f420bde23
# ╟─aa5c46fc-ccf6-49e1-b20b-7714bf2c7630
# ╟─ed0c4b7e-cf34-4b20-9100-43ba5428e096
