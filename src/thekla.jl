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

