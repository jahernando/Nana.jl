import CSV as CSV
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

function _dfnodes(xnd)

	dd = Dict()
	dd[:contents] = xnd.contents
	dd[:size]     = xnd.size
	for i in 1:length(xnd.coors)
		dd[Symbol(string("x", i, "_mean"))] = xnd.coors[i]
		dd[Symbol(string("x", i, "_std")) ] = xnd.coors_std[i]
	end

	dd[:distance_cell] = _distance(xnd.coors, xnd.coors_cell)

	vstd = reduce(hcat, xnd.coors_std)
	dstd = norm.(eachrow(vstd))
	dd[:dst] = dstd

	dd[:maxgrad] = xnd.maxgrad
	dd[:maxlap]  = xnd.maxlap
	dd[:minlap]  = xnd.minlap
	dd[:maxcur]  = xnd.curmax
	dd[:mincur]  = xnd.curmin

	dd[:nedges] = xnd.nedges
	dd[:cloud]  = xnd.cloud
	dd[:ecc]    = xnd.ecc
	dd[:extreme] = Int64.(xnd.extremes)

	#dd[:label]   = nlabel

	#return dd
	return DataFrame(dd)

end

function _coors(idf, steps; nsteps = 1)

	coors    = Tuple((idf[!, :x], idf[!, :y], idf[!, :z]))
	contents = copy(idf[!, :energy])
	xsteps    = Tuple(nsteps * steps)

	return coors, contents, xsteps
end

function _distance_to_blob(nlabel, dist)
	nnodes    = length(nlabel)
	idblobs   = findall(x -> x .== 3, nlabel)
	dist_blob = [minimum([dist[i, k] for k in idblobs]) for i in 1:nnodes]
	return dist_blob
end

function _thekla(idf, steps;
	 			 nsteps = 1, cellnode = false)

	coors, contents, xsteps = _coors(idf, steps; nsteps = nsteps)

	cl, nd, graph, edges = clouds(coors, contents, xsteps;
	  							cellnode = cellnode)
	clabel = label_cell(edges, cl.cells, coors, idf.segclass)
	nlabel = label_node(cl.node, clabel)

	dd         = _dfnodes(nd)
	dd[!, :label] = nlabel

	disttoblob = _distance_to_blob(nlabel, graph.dists)
	dd[!, :disttoblob] = disttoblob

	return dd
end

datadir       = "/Users/hernando/work/investigacion/NEXT/data/NEXT100/"
filename      = "bb0nu/v2/beersheba_fixed_label_1_0nubb.h5"
filenames     = Tuple(string("bb0nu/v2/beersheba_fixed_label_", i, "_0nubb.h5")
				for i in 1:249)
#df, mc, steps = load_data(string(datadir, filename));


function thekla(datadir, filenames;
	 			ofilename = "bb0nu/v2/thekla_nodes_nfiles50",
				reco = true, cellnode = false, nsteps = 1)

	datatype = reco ? "reco" : "mc"
	algoname = cellnode ? "paulina" : "clouds"
	ofile    = string(ofilename,  "_", algoname,
						"_", datatype,
					 	"_nsteps", nsteps,".csv")
	println(ofile)

	nevts = 0
	dfout = DataFrame()
	for (i, filename) in enumerate(filenames)
		df, mc, steps = load_data(string(datadir, filename))
		events        = event_list(df)
		println("events in file ", filename, " : ", length(events))
		for event in events
			idf = get_event(df, event)
			imc = get_event(mc, event)
			xdf = reco ? idf : imc
			odf = _thekla(xdf, steps;
			              nsteps = nsteps, cellnode = cellnode)
			odf[!, :event] .= 1000*i + event
			dfout = nevts == 0 ? odf : vcat(dfout, odf)
			nevts += 1
		end
	end

	ofile = string(datadir, ofile)
	println("write output at  : ", ofile)
	println("processed events : ", nevts)
	CSV.write(ofile, dfout)

	return dfout
end

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
