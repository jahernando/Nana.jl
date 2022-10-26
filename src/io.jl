using HDF5
using DataFrames

export load_data, get_event

function load_data(filename::String)

	dfs           = _get_dfs(filename)
	x0, x1, steps = _get_dimensions(dfs)

	df            = dfs["BeershebaVoxels"]
	df[!, "x"] = x0[1] .+ (df[!, "xbin"] .+ 0.5) * steps[1]
	df[!, "y"] = x0[2] .+ (df[!, "ybin"] .+ 0.5) * steps[2]
	df[!, "z"] = x0[3] .+ (df[!, "zbin"] .+ 0.5) * steps[3]

	mc            = dfs["MCHits"]
	return df, mc, steps
end;



"""
	get_event()

Returns the DF subset associated to an event
"""
function get_event(df::DataFrame,
	 			   event::Int64)

    sel  = Vector(df.dataset_id .== event)
    #println(sum(sel))
    #xbin = df[!, :xbin]
    idf  = df[sel, :]
	return idf

end;

#--
# Internal
#--

"""
	_get_dfs(filename)

Return a dictionary of DataFrames.

The DFs are `BeershebaVoxels`, `MCHits`, `BinsInfo`

Parameters:
----------
    filename : String, filename

Returns
-------
    dfs      : Dict(String, DF), dictionary of DataFrames

"""
function _get_dfs(filename::String)

	fid = h5open(filename, "r")
	dfs = Dict()

	for key in ["BeershebaVoxels", "MCHits", "BinsInfo"]
		dataset  = fid["DATASET"][key]
	    dfs[key] = _getdf(dataset)
	end

	close(fid)

    return dfs
end


"""
    _getdf(dataset)
"""
function _getdf(dataset)

    dda = read(dataset)

    ddic = Dict()
    for (i, key) in enumerate(keys(dda[1]))
        #println(i, ',', key)
        ddic[key] = [ddi[i] for ddi in dda]
    end

    df = DataFrame(ddic)
end

"""
    returns dimensions: frame origin and vixel size

    Parameters
    ----------
    dfs : Dict{DF}, dictionary with the DataFrame, requires "BinsInfo"

    Returns
    -------
    x0         : Vector{Real}, origin of the frame
	x1         : Vector{Real}, end of the frame
    voxel_size : Tuple{Real}, voxel size


"""
function _get_dimensions(dfs)

	df         = dfs["BinsInfo"]

    voxel_size = [df[!, var][1] for var in ("size_x", "size_y", "size_z")]
    x0         = [df[!, var][1] for var in ("min_x" , "min_y" , "min_z")]
	x1         = [df[!, var][1] for var in ("max_x" , "max_y" , "max_z")]

    return x0, x1, voxel_size
end
