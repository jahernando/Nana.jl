import StatsBase  as SB
using Clouds

export label_cell, label_cell_blobindex, label_node, initial_cell

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
	return 0
end

function label_cell(edges, cells, coors, labels; funct = _label)
	bins = bins_ids(coors, edges)
	xlabel = zeros(Int64, length(cells))
	for (i, icell) in enumerate(cells)
		sel = [b == icell for b in bins]
		labs = Set(labels[sel])
		xlabel[i] = funct(labs)
	end
	return xlabel
end

function label_cell_blobindex(cells, clabel)

	nsize   = length(clabel)
	blabels = zeros(Int64, nsize)

	sel     = clabel .== 3
	if (sum(sel) <= 0)
		return blabels
	end

	bcells  = cells[sel]
	blabel  = clabel[sel]
	bnodes  = cluster_nodes(bcells, blabel)

	blabels[sel] .= bnodes
	return blabels
end

function label_node(cell_node, cell_label)
	nnodes = maximum(cell_node)
	xlabel = zeros(Int64, nnodes)
	for i in 1:nnodes
		sel = cell_node .== i
		labs = Set(cell_label[sel])
		xlabel[i] = _label(labs)
	end
	return xlabel
end

function initial_cell(edges, cells, imc; min_nhits = 8)
	""" set to true the nodes with the first 10 mc-hits,
	(idf dataframe with hits) (imc dataframe with mc hits)
	"""
	mccoors = Tuple((imc[!, :x], imc[!, :y], imc[!, :z]))
	clabel  =  zeros(Int64, length(cells[1]))
	if (length(mccoors[1]) < min_nhits)
		return clabel
	end
	for nhits in 1:min_nhits
		xmccoors  = Tuple((mccoors[i][1:nhits] for i in 1:3))
		inilabel  = ones(Int64, nhits)
		clabel    = label_cell(edges, cells, xmccoors, inilabel)
		if (sum(clabel) >= 1)
			return clabel
		end
		if (nhits == min_nhits)
			return clabel
		end
	end
	return clabel
end
