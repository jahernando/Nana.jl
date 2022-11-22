#---- analysis

using DataFrames
using LinearAlgebra

export summary

function summary(df)
	dd = _summary(df)            # evern summary - reconstruction
	dd = _summary_label!(dd, df) # event summary - mc infor
	return dd
end

function _summary(df)
	
	function _dist(kdf, k1, k2)
		if (k1 <=0) | (k2 <= 0)
			return 0
		end
		return maximum((kdf.disttoextr1[k1], kdf.disttoextr1[k2]))
	end

	edf = groupby(df, :event)
	dd = Dict()
	dd[:event]   = [maximum(kdf.event)   for kdf in edf]
	dd[:energy]  = [sum(kdf.contents)    for kdf in edf]
	dd[:nnodes]  = [length(kdf.contents) for kdf in edf]
	dd[:nclouds] = [maximum(kdf.cloud)   for kdf in edf]
	dd[:nextrs]  = [sum(kdf.extreme)     for kdf in edf]
	ks  = [_idextremes(kdf) for kdf in edf]
	k1s = [ki[1] for ki in ks]
	k2s = [ki[2] for ki in ks]
	dd[:e1ene]  = [_val(kdf.contents, k) for (k, kdf) in zip(k1s, edf)]
	dd[:e2ene]  = [_val(kdf.contents, k) for (k, kdf) in zip(k2s, edf)]
	dd[:e1ecc]  = [_val(kdf.ecc, k)      for (k, kdf) in zip(k1s, edf)]
	dd[:e2ecc]  = [_val(kdf.ecc, k)      for (k, kdf) in zip(k2s, edf)]
	dd[:esdis]  = [_dist(kdf, k...)      for (k, kdf)  in zip(ks , edf)]

	return DataFrame(dd)
end

function _summary_label!(dd, df)

	edf = groupby(df, :event)

	function _d2n(kdf, k, i)
		if (k <= 0) 
			return 0
		end
		if (i <= 0)
			return 0
		end
		dd_ = kdf.disttoextr1[k] == 0 ? kdf.disttoextr1 : kdf.disttoextr2
		return dd_[i]
	end

	function _d2b(kdf, k, bs)
		dds = [_d2n(kdf, k, b) for b in bs]
		return minimum(dds)
	end
		
	ks  = [_idextremes(kdf) for kdf in edf]
	bs  = [_idblobs(kdf)    for kdf in edf]
	i1s = [_idinit(kdf)     for kdf in edf]

	k1s = [ki[1] for ki in ks]
	k2s = [ki[2] for ki in ks]
	b1s = [ki[1] for ki in bs]
	b2s = [ki[2] for ki in bs]

	dd[!, :e1label] = [_val(kdf.label, k) for (k, kdf) in zip(k1s, edf)]
	dd[!, :e2label] = [_val(kdf.label, k) for (k, kdf) in zip(k2s, edf)]
	dd[!, :i1label] = [_val(kdf.label, k) for (k, kdf) in zip(i1s, edf)]
	
	dd[!, :e1isbs]  = [_val(kdf.blobindex, k) > 0 for (k, kdf) in zip(k1s, edf)]
	dd[!, :e2isbs]  = [_val(kdf.blobindex, k) > 0 for (k, kdf) in zip(k2s, edf)]

	dd[!, :b1ene] = [_val(kdf.contents, k) for (k, kdf) in zip(b1s, edf)]
	dd[!, :b2ene] = [_val(kdf.contents, k) for (k, kdf) in zip(b2s, edf)]
	dd[!, :i1ene] = [_val(kdf.contents, k) for (k, kdf) in zip(i1s, edf)]

	dd[!, :e1d2b] = [_d2b(kdf, k, bi) for (k, bi, kdf) in zip(k1s, bs, edf)]
	dd[!, :e2d2b] = [_d2b(kdf, k, bi) for (k, bi, kdf) in zip(k2s, bs, edf)]

	dd[!, :e1d2i] = [_d2n(kdf, k, bi) for (k, bi, kdf) in zip(k1s, i1s, edf)]
	dd[!, :e2d2i] = [_d2n(kdf, k, bi) for (k, bi, kdf) in zip(k2s, i1s, edf)]

	return dd

end


function _idextremes(kdf)
    idextr  = findall(x -> x == 1, kdf.extreme)
    i1, i2 = 0, 0
    if (length(idextr) == 2)
        eextr   = kdf.contents[idextr]
        k1, k2  = argmax(eextr), argmin(eextr)
        i1, i2  = idextr[k1], idextr[k2]
    end
    return (i1, i2)
end

function _idblobs(kdf)
	
	ids  = findall(x -> x >= 1, kdf.blobindex)
	i1 = _idmax(ids, kdf.contents)
	if (i1 == 0)
		return (0, 0)
	end
	bindex =  kdf.blobindex[i1]
	ids  = findall(x -> x != bindex, ids)
	i2   = _idmax(ids, kdf.contents)
	return (i1, i2)
end

function _idinit(kdf)
	ids = findall(x -> x == 1, kdf.init)
	i1  = _idmax(ids, kdf.contents)
	return i1 
end

function _idmax(ids, contents)
	if (length(ids) <= 0)
		return 0
	end
	enes = contents[ids]
	k    = argmax(enes)
	i1   = ids[k]
	return i1
end

function _val(vals, k)
	if k in 1:length(vals)
		return vals[k]
	end
	return 0
end



function event_summary(df)
	nevents = length(Set(df.event))
	event   = zeros(Int64  , nevents)    # event number
	energy  = zeros(Float64, nevents)    # energy of the event
	ncloud  = zeros(Int64  , nevents)    # number of clouds (tracks)
	nnodes  = zeros(Int64  , nevents)    # number of nodes
	nextrs  = zeros(Int64  , nevents)    # number of extremes
	cextr1  = zeros(Int64  , nevents)    # eccentricity of the extreme 1
	cextr2  = zeros(Int64  , nevents)    # eccentricity of the extreme 2
	disext  = zeros(Int64  , nevents)    # distance between extreme 1 and 2
	eextr1  = zeros(Float64, nevents)    # energy of extreme 1 (most energetic)
	eextr2  = zeros(Float64, nevents)    # energy of extreme 2
	lextr1  = zeros(Int64  , nevents)    # label of the extreme 1
	lextr2  = zeros(Int64  , nevents)    # label of extreme 2
	d2bextr1 = zeros(Int64 , nevents)    # distance of the extreme 1 to a blob
	d2bextr2 = zeros(Int64 , nevents)    # distance of the extreme 2 to a blob
	d2iextr1 = zeros(Int64 , nevents)    # distance of the extreme 1 to the init node
	d2iextr2 = zeros(Int64 , nevents)    # distance of the extreme 1 to the init node
	nextrbs = zeros(Int64  , nevents)    # number of extremes labeled as blobs
	nnodebs = zeros(Int64  , nevents)    # number of nodes labeled as blobs
	nnodeis = zeros(Int64  , nevents)    # number of initial nodes
	enodeb1 = zeros(Float64, nevents)    # energy of the node labeled as blob 1-most energetic
	enodeb2 = zeros(Float64, nevents)    # energy of the node labelled as blob 2
	enodei1 = zeros(Float64, nevents)    # energy of the initial node
	lnodei1 = zeros(Int64  , nevents)    # label of the initial node

	edf     = groupby(df, :event)
	for (i, kdf) in enumerate(edf)
		event[i]   = maximum(kdf.event)
		energy[i]  = sum(kdf.contents)
		ncloud[i]  = maximum(kdf.cloud)
		nnodes[i]  = length(kdf.event)
		idextr     = findall(x -> x == 1, kdf.extreme)
		idblobs    = findall(x -> x == 3, kdf.label)
		idinit     = findall(x -> x == 1, kdf.init)
		nextrs[i]  = length(idextr)
		nnodebs[i] = length(idblobs)
		nnodeis[i] = length(idinit)
		if (length(idextr) >= 2)
			eextr   = kdf.contents[idextr]
			k1, k2  = argmax(eextr), argmin(eextr)
			i1, i2  = idextr[k1], idextr[k2]
			eextr1[i]  = maximum(eextr)
			eextr2[i]  = minimum(eextr)
			lextr1[i]  = kdf.label[i1]
			lextr2[i]  = kdf.label[i2]
			d2bextr1[i] = kdf.disttoblob[i1]
			d2bextr2[i] = kdf.disttoblob[i2]
			d2iextr1[i] = kdf.disttoinit[i1]
			d2iextr2[i] = kdf.disttoinit[i2]
			cextr1[i]  = kdf.ecc[i1]
			cextr2[i]  = kdf.ecc[i2]
			disext[i]  = maximum(kdf.disttoextr1)
			nextrbs[i]  = (lextr1[i] == 3) + (lextr2[i] == 3)
		end
		if (length(idinit) >= 1)
			einit       = kdf.contents[idinit]
			k1          = argmax(kdf.contents[idinit])
			i1          = idinit[k1]
			enodei1[i]  = maximum(einit)
			lnodei1[i]  = kdf.label[i1]
		end
		if (length(idblobs) >= 1)
			enodebs = sort(kdf.contents[idblobs], rev = true)
			enodeb1[i] = enodebs[1]
			if (length(idblobs) >= 2)
				enodeb2[i] = enodebs[2]
			end
		end
	end
	dd = Dict(:event    => event,
			  :energy   => energy,
			  :nclouds  => ncloud,
		      :nnodes   => nnodes,
		      :nextrs   => nextrs,
			  :nnodesbs => nnodebs,
			  :nnodesis => nnodeis,
		      :eextr1   => eextr1,
		      :eextr2   => eextr2,
		      :lextr1   => lextr1,
		      :lextr2   => lextr2,
			  :cextr1   => cextr1,
		      :cextr2   => cextr2,
			  :disext   => disext,
			  :d2bextr1 => d2bextr1,
			  :d2bextr2 => d2bextr2,
			  :d2iextr1 => d2iextr1,
			  :d2iextr2 => d2iextr2,
		      :nextrbs  => nextrbs,
		      :nnodesbs => nnodebs,
			  :nnodesis => nnodeis,
		      :enodeb1  => enodeb1,
		      :enodeb2  => enodeb2,
			  :enodei1  => enodei1,
			  :lnodei1  => lnodei1
			  )
	dd = DataFrame(dd)
return dd
end
