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
end;

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

# NEXT-100 thekla roc

Compute the RoC fo bb0nu and Bi214

Input data: Beersheba mc labeled hits of NEXT-100 Gonzalo production

J.A Hernando

November 2022

"""

# ╔═╡ a9c6da38-0aff-44f9-95b1-e838f451bc96
PUI.TableOfContents()

# ╔═╡ e4f36710-f912-4f6d-a47c-bc7fb4eff6bb
begin
datadir  = "/Users/hernando/work/investigacion/NEXT/data/NEXT100/"
end;

# ╔═╡ b16c5925-ad77-41a0-9fec-4f016a4e0b66
md"""
## Configuration
"""

# ╔═╡ cb1771d3-ec0a-43c8-be19-6db64b6758b5
begin
breco = @bind reco PUI.Select([:mc, :reco])
md"""
Select $(breco)
"""
end

# ╔═╡ f7c8c720-7516-4fe9-b0d5-ccf8341870bd
begin
anas = sort(["clouds_mc_nsteps1" , "clouds_mc_nsteps2", 
		     "breadth_mc_nsteps1", "breadth_mc_nsteps3"])
if reco == :reco
	anas = [replace(ana, "mc" => "reco") for ana in anas]
	anas = sort(["clouds_reco_nsteps1" , "clouds_reco_nsteps2", 
	    	     "breadth_reco_nsteps2", "breadth_reco_nsteps3"])
end
isos = ["bb0nu", "Bi"]
end;

# ╔═╡ 6f846223-6ff2-4fb3-80ea-2a5b1f44d8a7
md"""

Data 

| | |
|:-- | :-- |
| $(isos[1]) | $(isos[2]) |

Configurations

|  |  |  |  |
| :-- | :-- | :-- | :-- |
|$(anas[1]) | $(anas[2]) | $(anas[3]) | $(anas[4]) |

"""

# ╔═╡ 8c79c376-ad14-4a07-8575-6a61eb4b8664
function load_dfs(isos, anas, datatype)
	dfs = Dict{String, Dict{String, DataFrame}}()
	for iso in isos
		idf = Dict{String, DataFrame}()
		for ana in anas
			nfiles = iso == "bb0nu" ? "_nfiles249_" : "_nfiles310_"
			filename = string(iso, "/Thekla/thekla_", datatype, nfiles, ana, ".csv")
			#println("loading ", datatype, " : ", iso, ", ", ana, )
			di = DataFrame(CSV.File(string(datadir, filename)))
			idf[ana] = di
		end
		dfs[iso] = idf
	end
	return dfs
end

# ╔═╡ eda24864-e201-4d21-9dfb-08091c2b4796
begin
dfscounter = load_dfs(isos, anas, "counters")
dfssummary = load_dfs(isos, anas, "summary")

md"""
Data loaded!
"""
end

# ╔═╡ 892239b9-b79a-4118-9ca6-91eb453c8102
begin
#bana       = @bind iananames PUI.MultiCheckBox(anas)
#bana       = @bind sana      PUI.Select(anas)
bselnames  = @bind sselnames PUI.MultiCheckBox(["RoI", "nclouds", "nextremes"])
md"""

## Selection

Mark selections      : $(bselnames)

"""
end

# ╔═╡ 18cfba15-032b-47b9-b43a-1c475fe89647
function _effs_reco()
	effs = Dict{String, Dict{String, Float64}}()
	for iso in isos
		dd = Dict{String, Float64}()
		for ana in sort(anas)
			n0 = dfscounter[iso][ana].input_hits[1]
			n1 = dfscounter[iso][ana].output_nodes[1]
			dd[ana] = round(float(n1/n0), digits = 4)
		end
		effs[iso] = dd
	end
	return effs
end

# ╔═╡ 9a35cc03-d038-4878-a3f6-271ec63b597e
effs_reco = _effs_reco();

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
	eff  = round(sum(xsel)/length(dfe.event), digits = 4)
	return eff, xdfe
end

# ╔═╡ db3c9910-5856-470c-b5a9-3b01bbbe01ba
begin
dfs_sel  = Dict{String, Dict{String, DataFrame}}()
effs_sel = Dict{String, Dict{String, Float64}}()
for iso in isos
	dds = Dict{String, DataFrame}()
	dde = Dict{String, Float64}()
	for ana in anas
		eff, dd = _select(dfssummary[iso][ana], sselnames)
		dde[ana] = eff
		dds[ana] = dd
	end
	effs_sel[iso] = dde
	dfs_sel[iso]  = dds
end
md"""
Selection applied
"""
end

# ╔═╡ e3622633-dfca-4f83-98f5-0d6c13dfbfae
function mk_effs(iso, anas)
	cas = [effs_reco[iso][ana] for ana in anas] 
    cbs = [effs_sel[iso][ana]  for ana in anas]
	ccs = [round(a * b, digits = 4) for (a, b) in zip(cas, cbs)]
	
	n1, n2, n3, n4 = anas
	a1, a2, a3, a4 = cas
	b1, b2, b3, b4 = cbs
	c1, c2, c3, c4 = ccs

	md"""
	**Seleccion efficiency $(iso)**
	
	|type    | nodes | selection | total |
	| :--    |  :--  | :--       | :--   | 
	|$(n1)   | $(a1) | $(b1)     | $(c1) |
	|$(n2)   | $(a2) | $(b2)     | $(c2) |
	|$(n3)   | $(a3) | $(b3)     | $(c3) |
	|$(n4)   | $(a4) | $(b4)     | $(c4) |
	"""
end

# ╔═╡ 1ee5b093-269b-40f9-bc93-4327ca4892c7
mk_effs(isos[1], anas)

# ╔═╡ 70c0e062-1a47-4d4f-b695-ab44cd7db273
mk_effs(isos[2], anas)

# ╔═╡ b37a9175-f7eb-45c1-9770-05946134ce78
md"""

## Label Efficiency

"""

# ╔═╡ 3fa7cc02-a226-4d37-9198-6667cd6b8bed
function _eff_blob_init()
	effs_label_blob = Dict{String, Dict{String, Float64}}()
	effs_label_init = Dict{String, Dict{String, Float64}}()
	for (i, iso) in enumerate(["Bi", "bb0nu"])
		ddb = Dict{String, Float64}()
		ddi = Dict{String, Float64}()
		for ana in anas
			n0 = length(dfs_sel[iso][ana].nnblobs)
			n1 = sum(dfs_sel[iso][ana].nnblobs .== i)
			n2 = sum(dfs_sel[iso][ana].nninit .== 1)
			ddb[ana] = round(float(n1/n0), digits = 4)
			ddi[ana] = round(float(n2/n0), digits = 4)
		end
		effs_label_blob[iso] = ddb
		effs_label_init[iso] = ddi
	end
	return effs_label_blob, effs_label_init
end

# ╔═╡ f4e8ae4c-9307-4bdd-96fe-acc23190b1de
effs_label_blob, effs_label_init = _eff_blob_init();

# ╔═╡ efea7aca-917d-48c9-a42b-54181687b584
function mk_effs_label(iso, anas)
	cas = [effs_label_blob[iso][ana] for ana in anas] 
	cbs = [effs_label_init[iso][ana] for ana in anas] 

	n1, n2, n3, n4 = anas
	a1, a2, a3, a4 = cas
	b1, b2, b3, b4 = cbs
	
	md"""
	**Labeling blob and initial segments efficiency : $(iso)**
	
	|type    | blob  | init  |
	| :--    |  :--  |  :--  |
	|$(n1)   | $(a1) | $(b1) |
	|$(n2)   | $(a2) | $(b2) |
	|$(n3)   | $(a3) | $(b3) |
	|$(n4)   | $(a4) | $(b4) |
	"""
end

# ╔═╡ f63b7904-28d5-4dd3-84fa-4c3c3e72019c
mk_effs_label("bb0nu", anas)

# ╔═╡ 6e2e2c88-6ff1-4fca-aef9-3243b722df46
mk_effs_label("Bi", anas)

# ╔═╡ 6ff40eae-44de-4dc8-b257-b110a81e14ee
md"""

## RoC true

"""

# ╔═╡ 41253d0e-82b2-4956-a5d9-dee5ec050edf
function _roc(sig, bkg, erange)
	nvar1, nvar2 = length(sig), length(bkg)
	effsig = [sum(sig .>= ethr)/nvar1 for ethr in erange]
	effbkg = [sum(bkg .>= ethr)/nvar2 for ethr in erange]
	return 1 .- effbkg, effsig
end

# ╔═╡ 47c070bd-a2e5-4772-b618-bc8344ecf358
function _rocs(signal, bkg)
	xs, ys   = [], []
	hbs, his = [], []
	ethrs    = range(0., 1.5, 100)
	for (i, iname) in enumerate(anas)
		kdfbb = dfs_sel["bb0nu"][iname]
		kdfbi = dfs_sel["Bi"][iname]
		e2sig = kdfbb[!, signal]
		e2bkg = kdfbi[!, bkg]
		e2sig = e2sig #[e2sig .> 0.0]
		e2bkg = e2bkg #[e2bkg .> 0.0]
		emax = maximum(e2sig)
		nbins = range(0, emax, 100)
		hb = histogram(e2sig[e2sig .> 0], nbins = nbins, c = :green,
			title = iname, label = "bb", alpha = 0.5, normed = true, 
			xlabel = "energy (MeV)")
		histogram!(e2bkg[e2bkg .> 0], nbins = nbins, c = :blue,
			title = iname, label = "b", alpha = 0.5, normed = true)
		eff1, eff2 = _roc(e2sig, e2bkg, ethrs)
		push!(xs, eff1)
		push!(ys, eff2)
		push!(hbs, hb)
		#push!(his, hi)
	end
	return (xs, ys), hbs
end

# ╔═╡ c2cdf95f-0469-40cd-94d6-d6fbe6241d72
begin
points, hists = _rocs(:b2ene, :i1ene)
plot(points..., labels = permutedims(anas), lw = 2, legend= (0.12, 0.37),
	 xlabel = "bakground rejection", ylabel = "signal efficiency")
end

# ╔═╡ 43f8f033-6470-4895-98ac-60ccc757d4cf
md"""

Energy distribution of nodes labeled as blob2 and init
"""

# ╔═╡ d212fdd9-5180-4b27-8681-a87a93b48db8
plot(hists..., size= (650, 550), layout = (2, 2), legend = false)

# ╔═╡ d6f521f6-c77c-4c29-a6e9-a1ebd240db22
md"""
## RoC
"""

# ╔═╡ 7eb3d1a1-cd1d-4f1b-8a5a-9d840dc318a6
begin
rpoints, rhists = _rocs(:e2ene, :e2ene)
plot(rpoints..., labels = permutedims(anas), lw = 2, legend= (0.12, 0.37),
    xlabel = "bakground rejection", ylabel = "signal efficiency")

end

# ╔═╡ e7c2650c-13b8-4050-b546-965ab9fe67a1
md"""
Energy distribution of the extreme nodes
"""

# ╔═╡ 906cdbc0-f87f-4d1b-aa5f-a1df3b3eb5f4
begin
plot(rhists..., size= (650, 550), layout = (2, 2), legend = false)
end

# ╔═╡ 45db4fc4-f852-4855-904f-6cf810a4d1ab
md"""
## Blob2  success
"""

# ╔═╡ 75cb0aaf-8afd-4c34-8223-fee0e1536ee4
function _eff_esuccess(;extreme = 2)
	effs = Dict{String, Dict{String, Tuple{Float64, Float64, Float64}}}()
	for iso in isos
		dd = Dict{String, Tuple{Float64, Float64, Float64}}()
		for ana in anas
			df  = dfs_sel[iso][ana]
			l1, l2 = extreme == 1 ? (:e1d2b, :e1d2b) : (:e2d2b, :e2d2i)
			var = iso == "bb0nu" ? df[!, l1] : df[!, l2]
			n0  = length(var)
			n1  = sum(var .== 0)
			n2  = sum(var .<= 1)
			n3  = sum(var .<= 2) 
			ef1 = round(n1/n0, digits = 4)
			ef2 = round(n2/n0, digits = 4)
			ef3 = round(n3/n0, digits = 4)
			dd[ana] = (ef1, ef2, ef3)
		end
		effs[iso] = dd
	end
	return effs
end

# ╔═╡ 879b9f2a-6a9c-49a9-9d1d-06076382f52f
eff_e2success = _eff_esuccess();

# ╔═╡ 9906f064-d2e1-4416-8027-3beb1afe7082
function mk_eff_label(title, labels, effs, iso, anas)
	
	cas = [effs[iso][ana][1] for ana in anas] 
	cbs = [effs[iso][ana][2] for ana in anas] 
	ccs = [effs[iso][ana][3] for ana in anas] 


	l1, l2, l3, l4 = labels
	n1, n2, n3, n4 = anas
	a1, a2, a3, a4 = cas
	b1, b2, b3, b4 = cbs
	c1, c2, c3, c4 = ccs
	
	md"""
	**$(title)**
	
	|$(l1) | $(l2) | $(l3) | $(l4) |
	| :--  | :--   |  :--  |  :--  |
	|$(n1) | $(a1) | $(b1) | $(c1) |
	|$(n2) | $(a2) | $(b2) | $(c2) |
	|$(n3) | $(a3) | $(b3) | $(c3) |
	|$(n4) | $(a4) | $(b4) | $(c4) |
	"""
end

# ╔═╡ a870cecc-8a48-4ba7-ad36-855175709c4c
mk_eff_label(string("extreme 2 as blob ", "bb0nu"),
	         ("type", "dist == 0", "dist <= 1", "dist <= 2"),
			 eff_e2success, "bb0nu", anas)

# ╔═╡ cce1f837-87b9-4303-b8aa-66d2ea9f2eb4
mk_eff_label(string("extreme 2 as init ", "Bi"),
	          ("type", "dist == 0", "dist <= 1", "dist <= 2"),
			  eff_e2success, "Bi", anas)

# ╔═╡ c0547835-8d6d-40d5-97e4-cf116d7f4490
md"""
## Blob1 success
"""

# ╔═╡ 9b666e10-9140-46a2-9bc9-8a6bc3225cca
eff_e1success = _eff_esuccess(extreme = 1);

# ╔═╡ 720f1269-0d15-4aac-8072-95c4e2c4d988
mk_eff_label(string("extreme 1 as blob ", "bb0nu"),
	         ("type", "dist == 0", "dist <= 1", "dist <= 2"),
			 eff_e1success, "bb0nu", anas)

# ╔═╡ 04e3e912-38f4-41e7-a588-f67bf46df868
mk_eff_label(string("extreme 1 as blob ", "Bi"),
	         ("type", "dist == 0", "dist <= 1", "dist <= 2"),
			 eff_e1success, "Bi", anas)

# ╔═╡ 2b5982cd-3233-46d7-8fe9-a5ef908c831b
md"""
## Problematic Events
"""

# ╔═╡ e8c0d32e-e183-4311-99e3-10f090f45ebe
function _evts_condition(condition)
	nevts = Dict{String, Dict{String, NTuple{10, Int64}}}()
	for iso in isos
		dd = Dict{String, NTuple{10, Int64}}()
		for ana in anas
			df = dfs_sel[iso][ana]
			sel = condition(df)
			evts = sum(sel) < 10 ? zeros(Int64, 10) : df.event[sel][1:10]
			dd[ana] = Tuple(evts)
		end
		nevts[iso] = dd
	end
	return nevts
end

# ╔═╡ 6fcad78c-5183-44c1-9f30-1bb4beea213b
function _get_problematic_evts()
	conditions = Dict(:noblobs => x -> x.nnblobs .<= 1,
			          :e1nob   => x -> x.e1d2b   .>= 1,
			          :e2nob   => x -> x.e2d2b   .>= 1,
		    	      :e2noi   => x -> x.e2d2i   .>= 1)

	nevts = Dict()
	for key in keys(conditions)
		con   = conditions[key]
		ievts = _evts_condition(con)
		nevts[key] = ievts
	end
	return nevts
end

# ╔═╡ d68fdc1a-603e-4ef2-9052-307b84dd7ca7
badevts = _get_problematic_evts();

# ╔═╡ 0a46fbc2-fdd0-4f6c-9c13-05d6bb8b9ad5
begin
bbadevts_iso = @bind badevts_iso PUI.Select(isos)
bbadevts_ana = @bind badevts_ana PUI.Select(anas)
bbadevts_con = @bind badevts_con PUI.Select([:noblobs, :e1nob, :e2nob, :e2noi])
md"""

**Select Problematic events**

   problem      $(bbadevts_con) 

   isotope       $(bbadevts_iso)

   configuration $(bbadevts_ana)
"""
end

# ╔═╡ 9f093561-2d32-4f16-ac88-4d08b2fbfb93
println(badevts[badevts_con][badevts_iso][badevts_ana])

# ╔═╡ f3815da6-0bac-46ff-8a6f-5eeb43b298d3
md"""
## Code
"""

# ╔═╡ 6704978f-a68c-4774-acb1-8a10627e7d77
function mk_eff_44(title, labels, effs, iso, anas)

	@assert length(labels) == 4
	@assert length(anas)   == 4
	
	cas = [effs[iso][ana][1] for ana in anas] 
	cbs = [effs[iso][ana][2] for ana in anas] 
	ccs = [effs[iso][ana][3] for ana in anas] 

	l1, l2, l3, l4 = labels
	n1, n2, n3, n4 = anas
	a1, a2, a3, a4 = cas
	b1, b2, b3, b4 = cbs
	c1, c2, c3, c4 = ccs
	
	md"""
	**$(title)**
	
	|$(l1) | $(l2) | $(l3) | $(l4) |
	| :--  | :--   |  :--  |  :--  |
	|$(n1) | $(a1) | $(b1) | $(c1) |
	|$(n2) | $(a2) | $(b2) | $(c2) |
	|$(n3) | $(a3) | $(b3) | $(c3) |
	|$(n4) | $(a4) | $(b4) | $(c4) |
	"""
end

# ╔═╡ Cell order:
# ╟─68e36c10-a172-449b-b37e-28ecd17ff3a8
# ╟─a9c6da38-0aff-44f9-95b1-e838f451bc96
# ╟─13bdcd45-50fd-4ce0-9dba-1f3d4532d589
# ╟─4a04e7a4-5c5f-11ed-3a08-5fcaddecdf9d
# ╟─e4f36710-f912-4f6d-a47c-bc7fb4eff6bb
# ╟─b16c5925-ad77-41a0-9fec-4f016a4e0b66
# ╟─cb1771d3-ec0a-43c8-be19-6db64b6758b5
# ╟─f7c8c720-7516-4fe9-b0d5-ccf8341870bd
# ╟─6f846223-6ff2-4fb3-80ea-2a5b1f44d8a7
# ╟─8c79c376-ad14-4a07-8575-6a61eb4b8664
# ╟─eda24864-e201-4d21-9dfb-08091c2b4796
# ╟─892239b9-b79a-4118-9ca6-91eb453c8102
# ╟─18cfba15-032b-47b9-b43a-1c475fe89647
# ╟─9a35cc03-d038-4878-a3f6-271ec63b597e
# ╟─16e59b98-ef34-4f88-87d5-e0a9a3f8c300
# ╟─db3c9910-5856-470c-b5a9-3b01bbbe01ba
# ╟─e3622633-dfca-4f83-98f5-0d6c13dfbfae
# ╟─1ee5b093-269b-40f9-bc93-4327ca4892c7
# ╟─70c0e062-1a47-4d4f-b695-ab44cd7db273
# ╟─b37a9175-f7eb-45c1-9770-05946134ce78
# ╟─3fa7cc02-a226-4d37-9198-6667cd6b8bed
# ╟─f4e8ae4c-9307-4bdd-96fe-acc23190b1de
# ╟─efea7aca-917d-48c9-a42b-54181687b584
# ╟─f63b7904-28d5-4dd3-84fa-4c3c3e72019c
# ╟─6e2e2c88-6ff1-4fca-aef9-3243b722df46
# ╟─6ff40eae-44de-4dc8-b257-b110a81e14ee
# ╟─41253d0e-82b2-4956-a5d9-dee5ec050edf
# ╟─47c070bd-a2e5-4772-b618-bc8344ecf358
# ╟─c2cdf95f-0469-40cd-94d6-d6fbe6241d72
# ╟─43f8f033-6470-4895-98ac-60ccc757d4cf
# ╟─d212fdd9-5180-4b27-8681-a87a93b48db8
# ╟─d6f521f6-c77c-4c29-a6e9-a1ebd240db22
# ╟─7eb3d1a1-cd1d-4f1b-8a5a-9d840dc318a6
# ╟─e7c2650c-13b8-4050-b546-965ab9fe67a1
# ╟─906cdbc0-f87f-4d1b-aa5f-a1df3b3eb5f4
# ╟─45db4fc4-f852-4855-904f-6cf810a4d1ab
# ╟─75cb0aaf-8afd-4c34-8223-fee0e1536ee4
# ╟─879b9f2a-6a9c-49a9-9d1d-06076382f52f
# ╟─9906f064-d2e1-4416-8027-3beb1afe7082
# ╟─a870cecc-8a48-4ba7-ad36-855175709c4c
# ╟─cce1f837-87b9-4303-b8aa-66d2ea9f2eb4
# ╟─c0547835-8d6d-40d5-97e4-cf116d7f4490
# ╟─9b666e10-9140-46a2-9bc9-8a6bc3225cca
# ╟─720f1269-0d15-4aac-8072-95c4e2c4d988
# ╟─04e3e912-38f4-41e7-a588-f67bf46df868
# ╟─2b5982cd-3233-46d7-8fe9-a5ef908c831b
# ╟─e8c0d32e-e183-4311-99e3-10f090f45ebe
# ╟─6fcad78c-5183-44c1-9f30-1bb4beea213b
# ╟─d68fdc1a-603e-4ef2-9052-307b84dd7ca7
# ╟─0a46fbc2-fdd0-4f6c-9c13-05d6bb8b9ad5
# ╟─9f093561-2d32-4f16-ac88-4d08b2fbfb93
# ╟─f3815da6-0bac-46ff-8a6f-5eeb43b298d3
# ╟─6704978f-a68c-4774-acb1-8a10627e7d77
