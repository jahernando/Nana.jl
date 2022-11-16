import CSV as CSV
using DataFrames
using Clouds
using LinearAlgebra

export _coors, _thekla, event_summary

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
	dd[:std] = dstd

	dd[:maxgrad] = xnd.maxgrad
	dd[:maxlap]  = xnd.maxlap
	dd[:minlap]  = xnd.minlap
	dd[:maxcur]  = xnd.curmax
	dd[:mincur]  = xnd.curmin

	dd[:nedges]  = xnd.nedges
	dd[:cloud]   = xnd.cloud
	dd[:ecc]     = xnd.ecc
	dd[:extreme] = Int64.(xnd.extremes)

	#dd[:label]   = nlabel

	#return dd
	return DataFrame(dd)

end

function _connect(contents, emin = 1e-3)
	ene0 = sum(contents)
	sel = contents .== 0.0
	contents[sel] .= emin
	ene1 = sum(contents)
	factor = ene1 > 0 ? ene0/ene1 : 1.0
	contents = factor .* contents
	return contents
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
	if (length(idblobs) <= 0)
		return zeros(Int64, nnodes)
	end
	dist_blob = [minimum([dist[i, k] for k in idblobs]) for i in 1:nnodes]
	return dist_blob
end

function _distance_to_extremes(extrs, dist)
	nnodes = length(extrs)
	idextr = findall(x -> x .== 1, extrs)
	dtoextr1 = zeros(Int64, nnodes)
	dtoextr2 = zeros(Int64, nnodes)
	if ((length(idextr) != 2) || (nnodes <= 1))
		 return dtoextr1, dtoextr2
	end
	ide1, ide2 = idextr[1], idextr[2]
	dtoextr1 = dist[:, ide1]
	dtoextr2 = dist[:, ide2]
	return dtoextr1, dtoextr2
end



function _thekla(idf, steps;
	 			 nsteps = 1, cellnode = false)

	coors, contents, xsteps = _coors(idf, steps; nsteps = nsteps)
	contents = _connect(contents) # fix to connect empty energy voxels

	cl, nd, graph, edges = clouds(coors, contents, xsteps;
	  							  cellnode = cellnode)
	clabel = label_cell(edges, cl.cells, coors, idf.segclass)
	nlabel = label_node(cl.node, clabel)

	dd         = _dfnodes(nd)
	dd[!, :label] = nlabel

	disttoblob = _distance_to_blob(nlabel, graph.dists)
	dd[!, :disttoblob] = disttoblob

	#dtoextr1, dtoextr1 = _distance_to_extremes(dd.extreme, graph.dists)
	#dd[!, :dtoextr1] = dtoextr1
	#dd[!, :dtoextr2] = dtoextr2

	return dd
end

datadir        = "/Users/hernando/work/investigacion/NEXT/data/NEXT100/"
ifiles         = Dict()
ifiles[:bb0nu] = Tuple(string("bb0nu/v2/beersheba_fixed_label_", i, "_0nubb.h5")
				for i in 1:249)
ifiles[:Bi214] = Tuple(string("Bi/Beersheba/fixed_label/beersheba_label_", i,
				"_214Bi.h5") for i in 1:310)
ofiles         = Dict()
ofiles[:bb0nu] = "bb0nu/Thekla/thekla_nodes"
ofiles[:Bi214] = "Bi/Thekla/thekla_nodes"
#df, mc, steps = load_data(string(datadir, filename));

#-----
# Main function

function thekla(; data   = :bb0nu,
	 			  nfiles = -1,
				  reco   = true,
				  cellnode = false,
				  nsteps = 1)

	filenames = nfiles > 0 ? ifiles[data][1:nfiles] : ifiles[data]
	nfiles = length(filenames)
	println("Total number of input files ", nfiles)

	datatype = reco ? "reco" : "mc"
	algoname = cellnode ? "breadth" : "clouds"

	ofilename = ofiles[data]
	ofile     = string(ofilename, "_nfiles", nfiles, "_", algoname,
						"_", datatype,
					 	"_nsteps", nsteps,".csv")
	println("Output filename ", ofile)

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

#---- analysis

function event_summary(df)
	nevents = length(Set(df.event))
	event   = zeros(Int64  , nevents)    # event number
	energy  = zeros(Float64, nevents)    # energy of the event
	ncloud  = zeros(Int64  , nevents)    # number of clouds (tracks)
	nnodes  = zeros(Int64  , nevents)    # number of nodes
	nextrs  = zeros(Int64  , nevents)    # number of extremes
	cextr1  = zeros(Int64  , nevents)    # eccentricity of the extreme 1
	cextr2  = zeros(Int64  , nevents)    # eccentricity of the extreme 2
	eextr1  = zeros(Float64, nevents)    # energy of extreme 1 (most energetic)
	eextr2  = zeros(Float64, nevents)    # energy of extreme 2
	lextr1  = zeros(Int64  , nevents)    # label of the extreme 1
	lextr2  = zeros(Int64  , nevents)    # label of extreme 2
	dextr1  = zeros(Int64  , nevents)    # distance of the extreme 1 to a blob
	dextr2  = zeros(Int64  , nevents)    # distance of the extreme 2 to a blob
	nextrbs = zeros(Int64  , nevents)    # number of extremes labeled as blobs
	nnodebs = zeros(Int64  , nevents)    # number of nodes labeled as blobs
	enodeb1 = zeros(Float64, nevents)    # energy of the node labeled as blob 1-most energetic
	enodeb2 = zeros(Float64, nevents)    # energy of the node labelled as blob 2
	edf     = groupby(df, :event)
	for (i, kdf) in enumerate(edf)
		event[i]   = maximum(kdf.event)
		energy[i]  = sum(kdf.contents)
		ncloud[i]  = maximum(kdf.cloud)
		nnodes[i]  = length(kdf.event)
		idextr     = findall(x -> x == 1, kdf.extreme)
		idblobs    = findall(x -> x == 3, kdf.label)
		nextrs[i]  = length(idextr)
		nnodebs[i] = length(idblobs)
		if (length(idextr) == 2)
			eextr   = kdf.contents[idextr]
			k1, k2  = argmax(eextr), argmin(eextr)
			i1, i2  = idextr[k1], idextr[k2]
			eextr1[i]  = maximum(eextr)
			eextr2[i]  = minimum(eextr)
			lextr1[i]  = kdf.label[i1]
			lextr2[i]  = kdf.label[i2]
			dextr1[i]  = kdf.disttoblob[i1]
			dextr2[i]  = kdf.disttoblob[i2]
			cextr1[i]  = kdf.ecc[i1]
			cextr2[i]  = kdf.ecc[i2]
			nextrbs[i] = (lextr1[i] == 3) + (lextr2[i] == 3)
		end
		if (length(idblobs) >= 2)
			enodebs = sort(kdf.contents[idblobs], rev = true)
			enodeb1[i] = enodebs[1]
			enodeb2[i] = enodebs[2]
		end
	end
	dd = Dict(:event    => event,
			  :energy   => energy,
			  :nclouds  => ncloud,
		      :nnodes   => nnodes,
		      :nextrs   => nextrs,
		      :eextr1   => eextr1,
		      :eextr2   => eextr2,
		      :lextr1   => lextr1,
		      :lextr2   => lextr2,
			  :dextr1   => dextr1,
		      :dextr2   => dextr2,
			  :cextr1   => cextr1,
		      :cextr2   => cextr2,
		      :nextrbs  => nextrbs,
		      :nnodesbs => nnodebs,
		      :enodeb1  => enodeb1,
		      :enodeb2  => enodeb2)
	dd = DataFrame(dd)
return dd
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
