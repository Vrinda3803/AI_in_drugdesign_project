reinitialize
load data/external/ftsZ/22HQ.cif, hq_raw
load data/external/ftsZ/22HR.cif, hr_raw
create hq_chainA_gdp, hq_raw and chain A and (polymer.protein or resn GDP)
create hr_chainA_apo, hr_raw and chain A and polymer.protein
save data/external/ftsZ/prepared/22HQ_chainA_GDPbound.pdb, hq_chainA_gdp
save data/external/ftsZ/prepared/22HR_chainA_apo.pdb, hr_chainA_apo
create hr_aligned, hr_chainA_apo
align hr_aligned and polymer.protein and name CA, hq_chainA_gdp and polymer.protein and name CA, cycles=5
hide everything
show cartoon, hq_chainA_gdp and polymer.protein
show cartoon, hr_aligned and polymer.protein
color marine, hq_chainA_gdp and polymer.protein
color gray70, hr_aligned and polymer.protein
set cartoon_transparency, 0.48, hr_aligned and polymer.protein
quit
