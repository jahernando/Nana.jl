using HDF5
using DataFrames
using Clouds
using LinearAlgebra

export _coors, _thekla

function _distance(coors0, coors1)
	x0 = reduce(hcat, coors0)
	x1 = reduce(hcat, coors1)
	vec = x1 - x0
	dis = norm.(eachrow(vec))
	return dis
end

function _dfnodes(xnd, nlabel)

	dd = Dict()
	dd[:contents] = xnd.contents
	dd[:size]     = xnd.size
	for i in 1:length(xnd.coors)
		dd[Symbol(string("x", i, "_mean"))] = xnd.coors[i]
		dd[Symbol(string("x", i, "_std")) ] = xnd.coors_std[i]
	end

	dd[:distance_cell] = _distance(xnd.coors, xnd.coors_cell)

	dd[:maxgrad] = xnd.maxgrad
	dd[:maxlap]  = xnd.maxlap
	dd[:minlap]  = xnd.minlap
	dd[:maxcur]  = xnd.curmax
	dd[:mincur]  = xnd.curmin

	dd[:nedges] = xnd.nedges
	dd[:cloud]  = xnd.cloud
	dd[:ecc]    = xnd.ecc
	dd[:extreme] = Int64.(xnd.extremes)

	dd[:label]   = nlabel

	#return dd
	return DataFrame(dd)

end

function _coors(idf, steps, nsteps = 1.0)

	coors    = Tuple((idf[!, :x], idf[!, :y], idf[!, :z]))
	contents = copy(idf[!, :energy])
	xsteps    = Tuple(nsteps * steps)

	return coors, contents, xsteps
end

function _thekla(idf, imc, steps, event)

	coors, contents, xsteps = _coors(idf, steps)

	cl, nd, graph, edges = clouds(coors, contents, xsteps)
	clabel = label_cell(edges, cl.cells, coors, idf.segclass)
	nlabel = label_node(cl.node, clabel)

	dd = _dfnodes(nd, nlabel)
	dd[!, :event] .= event
	return dd
end

datadir       = "/Users/hernando/work/investigacion/NEXT/data/NEXT100/"
filename      = "bb0nu/v2/beersheba_fixed_label_1_0nubb.h5"
#df, mc, steps = load_data(string(datadir, filename));

#
# function run(filename, config)
#
#     df, mc, steps = load_data(filaname)
#     events        = event_list(df)idf, i
# #    _run          = configure_run(config)
#
#     for event in events:
#         idf  = get_event(df, event)
#         imc  = get_event(mc, event)
#         odfs = _run(idf, imc)
#         #_store(odfs, event)
#     end
#     #_save(odfs, config.output_filename)
# end
#
