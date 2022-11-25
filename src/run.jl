
using Nana


datadir        = "/Users/hernando/work/investigacion/NEXT/data/NEXT100/"
ifiles         = Dict()
ifiles[:bb0nu] = Tuple(string("bb0nu/v2/beersheba_fixed_label_", i, "_0nubb.h5")
				for i in 1:249)
ifiles[:Bi214] = Tuple(string("Bi/Beersheba/fixed_label/beersheba_label_", i,
				"_214Bi.h5") for i in 1:310)
ofiles         = Dict()
ofiles[:bb0nu] = "bb0nu/Thekla/thekla_nodes"
ofiles[:Bi214] = "Bi/Thekla/thekla_nodes"

thekla(; data = :bb0nu, evt_min_energy = 2.2, 
cellnode = true, reco = true, nsteps = 1)
thekla(; data = :Bi214, evt_min_energy = 2.2, 
cellnode = true, reco = true, nsteps = 1)


#-------------

function test()
    thekla(; reco = false, nfiles = 1)
end

function full_prod()

	thekla(; data = :bb0nu, evt_min_energy = 2.2, 
             cellnode = true, reco = false, nsteps = 1)
	thekla(; data = :Bi214, evt_min_energy = 2.2, 
             cellnode = true, reco = false, nsteps = 1)

	prod(; cellnode = false, nsteps = 1)

	prod(; cellnode = false, nsteps = 3)
	prod(; cellnode = true , nsteps = 3)

	prod(; cellnode = false, nsteps = 2)
	prod(; cellnode = true , nsteps = 2)

end


function prod(; cellnode = false, nsteps = 1)

	for data in [:bb0nu, :Bi214]
		for reco in (false, true)
			thekla(; data = data, evt_min_energy = 2.2, reco = reco,
	         	     cellnode = cellnode,  nsteps = nsteps)
		end
	end
end
