#import StatsBase as SB
#import LinearAlgebra as LA
import Graphs as GG

export graph_extremes, extremes_maxcontents

function graph_extremes(graph)
	nvertices = GG.nv(graph)
	dists     = zeros(Int64, nvertices, nvertices)
	for i in 1:nvertices
		ds = GG.dijkstra_shortest_paths(graph, i)
		vv = copy(ds.dists)
		vv[vv .> nvertices] .= 0
		dists[i, :] = vv
	end
	extremes = findall(x -> x .>= maximum(dists), dists)
	extremes = [(i, j) for (i, j) in Tuple.(extremes) if i < j]
	return extremes, dists
end

function extremes_maxcontents(extremes, xnd)
	contents = [xnd.contents[i] + xnd.contents[j] for (i, j) in extremes]
	i = findall(x -> x .== maximum(contents), contents)
	return extremes[i]
end
