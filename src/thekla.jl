import CSV as CSV
using DataFrames
using Clouds
using LinearAlgebra

export thekla, prod, event_summary, _coors, _thekla

#-----------------
# Production
#-----------------

datadir        = "/Users/hernando/work/investigacion/NEXT/data/NEXT100/"
ifiles         = Dict()
ifiles[:bb0nu] = Tuple(string("bb0nu/v2/beersheba_fixed_label_", i, "_0nubb.h5")
				for i in 1:249)
ifiles[:Bi214] = Tuple(string("Bi/Beersheba/fixed_label/beersheba_label_", i,
				"_214Bi.h5") for i in 1:310)
ofiles         = Dict()
ofiles[:bb0nu] = "bb0nu/Thekla/thekla_nodes"
ofiles[:Bi214] = "Bi/Thekla/thekla_nodes"


function prod(; cellnode = false, nsteps = 1)

	for data in [:bb0nu, :Bi214]
		for reco in (false, true)
			if (reco)
				continue
			end
			thekla(; data = data, evt_min_energy = 2.2, reco = reco,
	         	     cellnode = cellnode,  nsteps = nsteps)
		end
	end
end


#--------------
#  City
#---------------

function thekla(; data           = :bb0nu,
	 			  nfiles         = -1,
				  evt_min_energy = 2.0, # MeV
				  reco           = true,
				  cellnode       = false,
				  nsteps         = 1)

	counters = Dict(:input_hits     => 0,
					:evt_min_energy => 0,
					:evt_min_nodes  => 0,
					:output_nodes   => 0)

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
	dfout = DataFrame() # nodes - frame 
	dfsum = DataFrame() # summary - frame
	for (i, filename) in enumerate(filenames)
		df, mc, steps = load_data(string(datadir, filename))
		events        = event_list(df)
		println("events in file ", filename, " : ", length(events))
		for event in events
			counters[:input_hits] += 1
			idf = get_event(df, event)
			imc = get_event(mc, event)
			xdf = reco ? idf : imc

			ene = sum(xdf.energy)
			if (ene <= evt_min_energy)
				continue
			end
			counters[:evt_min_energy] += 1

			odf = _thekla(xdf, imc, steps;
			              nsteps = nsteps, cellnode = cellnode)
			nnodes = length(odf.contents)
			if (nnodes <= 0)
				continue
			end
			counters[:evt_min_nodes] += 1

			odf[!, :event] .= 1000*i + event

			# event summary of thisfile
			idfsum  = summary(odf)

			dfout = nevts == 0 ? odf     : vcat(dfout, odf)
			dfsum = nevts == 0 ? idfsum  : vcat(dfsum, idfsum) 
			nevts += 1
			counters[:output_nodes] += 1
		end
	end

	# Output
	ofile = string(datadir, ofile)
	println("write node output at  : ", ofile)
	println("processed events : ", nevts)
	CSV.write(ofile, dfout)

	ofile_summary = replace(ofile, "nodes" => "summary")
	println("write summary output at  : ", ofile)
	CSV.write(ofile_summary, dfsum)

	ofile_counters = replace(ofile, "nodes" => "counters")
	println("write counters output at  : ", ofile)
	dcounters = DataFrame(counters)
	CSV.write(ofile_counters, dcounters)


	return dfout, dfsum, dcounters
end


#------------------
#  City per Event
#------------------

function _thekla(idf, imc, steps;
	 			 nsteps = 1,
				 cellnode = false)
	""" 
	main algorithms per event:
		run clouds
		create data-frame of the clouds nodes
		add MC-information into the nodes data frames 
	"""

	# get clouds input
	coors, contents, xsteps = _coors(idf, steps; nsteps = nsteps)
	# fix to connect emtpy energy voxels
	contents = _connect(contents) # fix to connect empty energy voxels

	# call clouds
	cl, nd, graph, edges = clouds(coors, contents, xsteps;
	  							  cellnode = cellnode)

	# convert to DataFrame
	dfcl = _todf(cl, (:node, :contents))
	#dfnd = _todf(nd) 

	# Data Frame with the nodes variables
	dd  = _dfnodes(nd)
	
	# Label the clouds
	mccoors  = Tuple((imc[!, :x], imc[!, :y], imc[!, :z]))
	segclass = idf.segclass
	_, dd = label_clouds!(dfcl, dd, cl.cells, edges, mccoors, coors, segclass)

	# extend the nodes with label information and distances
	dd = _dfnodes_label!(dd, graph.dists)

	return dd
end

#-----
# function to set the nodes Data-Frame with reonstruction (clouds), and MC (labeling)
#-----

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

	return DataFrame(dd)

end

function _dfnodes_label!(dd,      # nodes info DataFrames (to be extended)
						 dists)   # distance matrix between the nodes
				  

	kids = findall(x -> x >= 1, dd[!, :extreme])
	bids = findall(x -> x >= 1, dd[!, :blobindex])
	iids = findall(x -> x >= 1, dd[!, :init]) 
	function _d2(k, ids)
		if (length(ids) <= 0)
			return -1
		end
		return minimum([dists[k, i] for i in ids])
	end

	nnodes = length(dd[!, :contents])
	dd[!, :disttoblob] = [_d2(k, bids) for k in 1:nnodes]
	dd[!, :disttoinit] = [_d2(k, iids) for k in 1:nnodes]

	k1 = kids[argmax(dd[!, :contents][kids])] # extreme with the max energy
	k2 = kids[argmin(dd[!, :contents][kids])] # extreme with the min energy

	dd[!, :disttoextr1] = [_d2(k, (k1,)) for k in 1:nnodes]
	dd[!, :disttoextr2] = [_d2(k, (k2,)) for k in 1:nnodes]
	
	return dd
end

# function _dfnodes_label!(dd,      # nodes info DataFrames (to be extended)
# 	mccoors,      # coordenates of the MC hits
#    coors,    # coordinates of the original hits (Beersehba hits)
# 	 segclass, # label of the original hits (Beersheba hits)
# 	 cl,       # struct with the cells info from the clouds oupout
# 	 graph,    # struct with the graph infrom from the clouds output
# 	 edges)    # edges of the discritization of clouds

# # label cells and nodes
# clabel = label_cell(edges, cl.cells, coors, segclass)
# nlabel = label_node(cl.node, clabel)
# dd[!, :label] = nlabel

# # index the cells as blob 1, 2
# bclabel = label_cell_blobindex(cl.cells, clabel)
# bnlabel = label_node(cl.node, bclabel)
# dd[!, :blobindex] = bnlabel

# # initial cell and node
# cinit  = initial_cell(edges, cl.cells, mccoors)
# ninit  = label_node(cl.node, cinit)
# dd[!, :init]  = ninit

# # distance to blob
# disttoblob = _distance_to_blob(nlabel, graph.dists)
# dd[!, :disttoblob] = disttoblob

# # distance to init
# disttoinit = _distance_to_blob(ninit, graph.dists; label = 1)
# dd[!, :disttoinit] = disttoinit

# #distance to extremes
# dtoextr1, dtoextr2 = _distance_to_extremes(dd.extreme, graph.dists)
# dd[!, :disttoextr1] = dtoextr1
# dd[!, :disttoextr2] = dtoextr2

# return dd
# end




#--------------------------
# Internal Functions
#------------------------

# struct to DF
function _todf(x, names = "")
	names = names == "" ? fieldnames(typeof(x)) : names
	df = DataFrame()
	for name in names
		df[!, name] = getfield(x, name)
	end
	return df
end

function _distance(coors0, coors1)
	x0 = reduce(hcat, coors0)
	x1 = reduce(hcat, coors1)
	vec = x1 - x0
	dis = norm.(eachrow(vec))
	return dis
end

	
function _connect(contents, emin = 1e-3)
	""" fill the null energy beersheba voxels to a minimum energry
	"""
	ene0 = sum(contents)
	sel = contents .== 0.0
	contents[sel] .= emin
	ene1 = sum(contents)
	factor = ene1 > 0 ? ene0/ene1 : 1.0
	contents = factor .* contents
	return contents
end

function _coors(idf, steps; nsteps = 1)
	""" return the hit data (in idf DataFrame) into
	the data required by clouds (coors, contents, steps)
	"""

	coors    = Tuple((idf[!, :x], idf[!, :y], idf[!, :z]))
	contents = copy(idf[!, :energy])
	xsteps   = Tuple(nsteps * steps)

	return coors, contents, xsteps
end

# function _distance_to_blob(nlabel, dist; label = 3)
# 	""" compute the minimun distance of each note to a node labeled as blob
# 	"""
# 	nnodes    = length(nlabel)
# 	idblobs   = findall(x -> x .== label, nlabel)
# 	if (length(idblobs) <= 0)
# 		return zeros(Int64, nnodes)
# 	end
# 	dist_blob = [minimum([dist[i, k] for k in idblobs]) for i in 1:nnodes]
# 	return dist_blob
# end

# function _distance_to_extremes(extrs, dist)
# 	""" return the distance of each note the two extreme nodes, otherwise is 0
# 	"""
# 	nnodes = length(extrs)
# 	idextr = findall(x -> x .== 1, extrs)
# 	dtoextr1 = zeros(Int64, nnodes)
# 	dtoextr2 = zeros(Int64, nnodes)
# 	if ((length(idextr) != 2) || (nnodes <= 1))
# 		 return dtoextr1, dtoextr2
# 	end
# 	ide1, ide2 = idextr[1], idextr[2]
# 	dtoextr1 = dist[:, ide1]
# 	dtoextr2 = dist[:, ide2]
# 	return dtoextr1, dtoextr2
# end


# #---- analysis

# function event_summary(df)
# 	nevents = length(Set(df.event))
# 	event   = zeros(Int64  , nevents)    # event number
# 	energy  = zeros(Float64, nevents)    # energy of the event
# 	ncloud  = zeros(Int64  , nevents)    # number of clouds (tracks)
# 	nnodes  = zeros(Int64  , nevents)    # number of nodes
# 	nextrs  = zeros(Int64  , nevents)    # number of extremes
# 	cextr1  = zeros(Int64  , nevents)    # eccentricity of the extreme 1
# 	cextr2  = zeros(Int64  , nevents)    # eccentricity of the extreme 2
# 	disext  = zeros(Int64  , nevents)    # distance between extreme 1 and 2
# 	eextr1  = zeros(Float64, nevents)    # energy of extreme 1 (most energetic)
# 	eextr2  = zeros(Float64, nevents)    # energy of extreme 2
# 	lextr1  = zeros(Int64  , nevents)    # label of the extreme 1
# 	lextr2  = zeros(Int64  , nevents)    # label of extreme 2
# 	d2bextr1 = zeros(Int64 , nevents)    # distance of the extreme 1 to a blob
# 	d2bextr2 = zeros(Int64 , nevents)    # distance of the extreme 2 to a blob
# 	d2iextr1 = zeros(Int64 , nevents)    # distance of the extreme 1 to the init node
# 	d2iextr2 = zeros(Int64 , nevents)    # distance of the extreme 1 to the init node
# 	nextrbs = zeros(Int64  , nevents)    # number of extremes labeled as blobs
# 	nnodebs = zeros(Int64  , nevents)    # number of nodes labeled as blobs
# 	nnodeis = zeros(Int64  , nevents)    # number of initial nodes
# 	enodeb1 = zeros(Float64, nevents)    # energy of the node labeled as blob 1-most energetic
# 	enodeb2 = zeros(Float64, nevents)    # energy of the node labelled as blob 2
# 	enodei1 = zeros(Float64, nevents)    # energy of the initial node
# 	lnodei1 = zeros(Int64  , nevents)    # label of the initial node
# 	edf     = groupby(df, :event)
# 	for (i, kdf) in enumerate(edf)
# 		event[i]   = maximum(kdf.event)
# 		energy[i]  = sum(kdf.contents)
# 		ncloud[i]  = maximum(kdf.cloud)
# 		nnodes[i]  = length(kdf.event)
# 		idextr     = findall(x -> x == 1, kdf.extreme)
# 		idblobs    = findall(x -> x == 3, kdf.label)
# 		idinit     = findall(x -> x == 1, kdf.init)
# 		nextrs[i]  = length(idextr)
# 		nnodebs[i] = length(idblobs)
# 		nnodeis[i] = length(idinit)
# 		if (length(idextr) >= 2)
# 			eextr   = kdf.contents[idextr]
# 			k1, k2  = argmax(eextr), argmin(eextr)
# 			i1, i2  = idextr[k1], idextr[k2]
# 			eextr1[i]  = maximum(eextr)
# 			eextr2[i]  = minimum(eextr)
# 			lextr1[i]  = kdf.label[i1]
# 			lextr2[i]  = kdf.label[i2]
# 			d2bextr1[i] = kdf.disttoblob[i1]
# 			d2bextr2[i] = kdf.disttoblob[i2]
# 			d2iextr1[i] = kdf.disttoinit[i1]
# 			d2iextr2[i] = kdf.disttoinit[i2]
# 			cextr1[i]  = kdf.ecc[i1]
# 			cextr2[i]  = kdf.ecc[i2]
# 			disext[i]  = maximum(kdf.disttoextr1)
# 			nextrbs[i]  = (lextr1[i] == 3) + (lextr2[i] == 3)
# 		end
# 		if (length(idinit) >= 1)
# 			einit       = kdf.contents[idinit]
# 			k1          = argmax(kdf.contents[idinit])
# 			i1          = idinit[k1]
# 			enodei1[i]  = maximum(einit)
# 			lnodei1[i]  = kdf.label[i1]
# 		end
# 		if (length(idblobs) >= 1)
# 			enodebs = sort(kdf.contents[idblobs], rev = true)
# 			enodeb1[i] = enodebs[1]
# 			if (length(idblobs) >= 2)
# 				enodeb2[i] = enodebs[2]
# 			end
# 		end
# 	end
# 	dd = Dict(:event    => event,
# 			  :energy   => energy,
# 			  :nclouds  => ncloud,
# 		      :nnodes   => nnodes,
# 		      :nextrs   => nextrs,
# 			  :nnodesbs => nnodebs,
# 			  :nnodesis => nnodeis,
# 		      :eextr1   => eextr1,
# 		      :eextr2   => eextr2,
# 		      :lextr1   => lextr1,
# 		      :lextr2   => lextr2,
# 			  :cextr1   => cextr1,
# 		      :cextr2   => cextr2,
# 			  :disext   => disext,
# 			  :d2bextr1 => d2bextr1,
# 			  :d2bextr2 => d2bextr2,
# 			  :d2iextr1 => d2iextr1,
# 			  :d2iextr2 => d2iextr2,
# 		      :nextrbs  => nextrbs,
# 		      :nnodesbs => nnodebs,
# 			  :nnodesis => nnodeis,
# 		      :enodeb1  => enodeb1,
# 		      :enodeb2  => enodeb2,
# 			  :enodei1  => enodei1,
# 			  :lnodei1  => lnodei1
# 			  )
# 	dd = DataFrame(dd)
# return dd
# end

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
