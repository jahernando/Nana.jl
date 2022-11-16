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

# ╔═╡ 13bdcd45-50fd-4ce0-9dba-1f3d4532d589
begin
using Pkg;
Pkg.activate("/Users/hernando/work/investigacion/NEXT/software/julias/Nana")
end

# ╔═╡ 4a04e7a4-5c5f-11ed-3a08-5fcaddecdf9d
begin
using CSV
using DataFrames
using StatsBase
using Plots
import PlutoUI as PUI
import Nana as na
end

# ╔═╡ 68e36c10-a172-449b-b37e-28ecd17ff3a8
md"""

# NEXT-100 thekla  roc

Compute the RoC with the eblob2

J.A Hernando

November 2022

"""

# ╔═╡ a6347485-2dad-4343-bfe0-ed0de5a1e787
plotly();

# ╔═╡ e4f36710-f912-4f6d-a47c-bc7fb4eff6bb
begin
datadir  = "/Users/hernando/work/investigacion/NEXT/data/NEXT100/"
end

# ╔═╡ b16c5925-ad77-41a0-9fec-4f016a4e0b66
md"""
## Blob success rate
"""

# ╔═╡ f7c8c720-7516-4fe9-b0d5-ccf8341870bd
anas = sort(["clouds_mc_nsteps1", "clouds_reco_nsteps1",
	         "clouds_mc_nsteps2", "clouds_reco_nsteps2"])
#	    "breadth_mc_nsteps1", "breadth_mc_nsteps3", 
#	    "breadth_reco_nsteps3"]


# ╔═╡ 73efbbfa-70a8-4e6a-a019-77810f45b2e3
begin
dfbb = Dict()
for afile in anas
	filename = string("bb0nu/Thekla/thekla_nodes_nfiles249_", afile, ".csv")
	println(filename)
	di = DataFrame(CSV.File(string(datadir, filename)))
	dfbb[afile] = na.event_summary(di)
end
dfbi = Dict()
for afile in anas
	filename = string("Bi/Thekla/thekla_nodes_nfiles310_", afile, ".csv")
	println(filename)
	di = DataFrame(CSV.File(string(datadir, filename)))
	dfbi[afile] = na.event_summary(di)
end
end

# ╔═╡ 892239b9-b79a-4118-9ca6-91eb453c8102
begin
bana       = @bind iananames PUI.MultiCheckBox(anas)
bselnames  = @bind iselnames PUI.MultiCheckBox(["RoI", "nclouds", "nextremes"])
md"""

## RoC

### Configuration

Select configuration : $(bana)

Mark selections      : $(bselnames)

"""
end

# ╔═╡ 16e59b98-ef34-4f88-87d5-e0a9a3f8c300
function _select(dfe, selnames)
	xsel = ones(Bool, length(dfe.event))
	if "RoI" in selnames
		xsel = xsel .& (dfe.energy .>= 2.41) .& (dfe.energy .<= 2.5)
	end
	if "nclouds" in selnames
		xsel = xsel .& (dfe.nclouds .== 1)
	end
	if "nclouds" in selnames
		xsel = xsel .& (dfe.nextrs .== 2)
	end
	xdfe = dfe[xsel, :] 
	eff  = round(sum(xsel)/length(dfe.event), digits = 3)
	return eff, xdfe
end

# ╔═╡ 1d10dc9d-84e6-4fe9-a460-55375ae0572b
begin
iname = iananames[1]
effbb, idfbb = _select(dfbb[iname], iselnames)
effbi, idfbi = _select(dfbi[iname], iselnames)
end;

# ╔═╡ bcb2522a-f9f3-4f2f-a091-e909071dc3fb
md"""

### Selection efficiency

| selection    | bb0nu   |   Bi  |
| --:          |  --:    |    --: |
| $(iname) | $(effbb) | $(effbi) |   

"""

# ╔═╡ 41253d0e-82b2-4956-a5d9-dee5ec050edf
function _roc(var1, var2, erange)
	nvar1, nvar2 = length(var1), length(var2)
	eff1 = [sum(var1 .>= ethr)/nvar1 for ethr in erange]
	eff2 = [sum(var2 .>= ethr)/nvar2 for ethr in erange]
	return eff1, eff2
end

# ╔═╡ e174710e-eb52-45c1-a0d8-179cf9818caa
eff1, eff2 = _roc(idfbb.eextr2, idfbi.eextr2, range(0, 1.5, 100));

# ╔═╡ acac3652-344f-4e0b-ad96-f67d14d7c36e
#plot(eff1, eff2, label = iname, lw = 2.)

# ╔═╡ 6ff40eae-44de-4dc8-b257-b110a81e14ee
md"""

## RoC curve

"""

# ╔═╡ 47c070bd-a2e5-4772-b618-bc8344ecf358
begin
xs, ys   = [], []
hs1, hs2 = [], []
ethrs    = range(0., 1.5, 100)
for (i, iname) in enumerate(iananames)
	effbb, idfbb = _select(dfbb[iname], iselnames)
	effbi, idfbi = _select(dfbi[iname], iselnames)
	h1 = histogram(idfbb.eextr2, nbins = ethrs, title = string("bb ", iname))
	push!(hs1, h1)
	h2 = histogram(idfbi.eextr2, nbins = ethrs, title = string("Bi ", iname))
	push!(hs2, h2)
	eff1, eff2 = _roc(idfbb.eextr2, idfbi.eextr2, ethrs)
	push!(xs, eff1)
	push!(ys, eff2)
end
plot(xs, ys, label = Tuple(iananames), lw = 2)
end

# ╔═╡ e615d28f-9323-422d-a671-1a19dbdd9a17
md"""

### eblob bb

"""

# ╔═╡ c889c1f9-02a3-4343-ad5d-c6a1cad3d3b2
plot(hs1..., layout = (2, 2))

# ╔═╡ 233595ab-0662-481f-b93a-b6c1f8e9fb60
md"""
### eblob Bi
"""

# ╔═╡ aa411bd3-8129-4518-a465-9144d6b1b59b
plot(hs2..., layout = (2, 2))

# ╔═╡ 65f166ea-9259-4e07-9baf-a25e568f7acb
md"""
### RoC with blob2 (bb)
"""

# ╔═╡ f8586d5f-323f-4e7a-89f3-d2bed5355941
begin
bxs, bys   = [], []
bhs1, bhs2 = [], []
for (i, iname) in enumerate(iananames)
	effbb, idfbb = _select(dfbb[iname], iselnames)
	effbi, idfbi = _select(dfbi[iname], iselnames)
	h1 = histogram(idfbb.enodeb2, nbins = ethrs, title = string("bb ", iname))
	push!(bhs1, h1)
	h2 = histogram(idfbi.eextr2, nbins = ethrs, title = string("Bi ", iname))
	push!(bhs2, h2)
	eff1, eff2 = _roc(idfbb.enodeb2, idfbi.eextr2, ethrs)
	push!(bxs, eff1)
	push!(bys, eff2)
end
plot(bxs, bys, label = [iananames[1], iananames[2]], lw = 2)
end

# ╔═╡ f40a6b15-5a9a-4596-b910-605257627185
md"""
### eblob2 (bb)
"""

# ╔═╡ 25b565cd-a934-4162-b480-9a423e9ed0b1
plot(bhs1..., layout = (2, 2))

# ╔═╡ 4953f0b4-c268-4d85-8aae-8059b9f1149c
plot(bhs2..., layout = (2, 2))

# ╔═╡ 07f27c84-82de-4590-85c9-c91900469ead
md"""

## Code

"""

# ╔═╡ Cell order:
# ╠═68e36c10-a172-449b-b37e-28ecd17ff3a8
# ╠═13bdcd45-50fd-4ce0-9dba-1f3d4532d589
# ╠═4a04e7a4-5c5f-11ed-3a08-5fcaddecdf9d
# ╟─a6347485-2dad-4343-bfe0-ed0de5a1e787
# ╠═e4f36710-f912-4f6d-a47c-bc7fb4eff6bb
# ╟─b16c5925-ad77-41a0-9fec-4f016a4e0b66
# ╠═f7c8c720-7516-4fe9-b0d5-ccf8341870bd
# ╠═73efbbfa-70a8-4e6a-a019-77810f45b2e3
# ╟─892239b9-b79a-4118-9ca6-91eb453c8102
# ╟─16e59b98-ef34-4f88-87d5-e0a9a3f8c300
# ╠═1d10dc9d-84e6-4fe9-a460-55375ae0572b
# ╟─bcb2522a-f9f3-4f2f-a091-e909071dc3fb
# ╟─41253d0e-82b2-4956-a5d9-dee5ec050edf
# ╟─e174710e-eb52-45c1-a0d8-179cf9818caa
# ╟─acac3652-344f-4e0b-ad96-f67d14d7c36e
# ╟─6ff40eae-44de-4dc8-b257-b110a81e14ee
# ╟─47c070bd-a2e5-4772-b618-bc8344ecf358
# ╟─e615d28f-9323-422d-a671-1a19dbdd9a17
# ╠═c889c1f9-02a3-4343-ad5d-c6a1cad3d3b2
# ╟─233595ab-0662-481f-b93a-b6c1f8e9fb60
# ╠═aa411bd3-8129-4518-a465-9144d6b1b59b
# ╟─65f166ea-9259-4e07-9baf-a25e568f7acb
# ╟─f8586d5f-323f-4e7a-89f3-d2bed5355941
# ╟─f40a6b15-5a9a-4596-b910-605257627185
# ╠═25b565cd-a934-4162-b480-9a423e9ed0b1
# ╠═4953f0b4-c268-4d85-8aae-8059b9f1149c
# ╠═07f27c84-82de-4590-85c9-c91900469ead
