# FBAwMC example
### Example to run FBAwMC for all 6 LAB to obtain their growth solution and internal flux
- [f, v, modelIrrev] = FBAwMC_LAB(model,C,biomass,num)
- C, or Cytoplasmic density values for 6 LABs are provided in 
- LAB-probiotic-evaluation-framework/Data/Cytoplasmic density values for 6 LAB.csv
- Number of flux solutions is set at 5000


### Load 6 LAB models
lca=load('lca.mat');<br />
lfe=load('lfe.mat');<br />
lla=load('lla.mat');<br />
lme=load('lme.mat');<br />
lsa=load('lsa.mat');<br />
lpl=load('lpl.mat');<br />


### Set parameters
num = 5000;<br />
biomass = 'biomass';<br />


### Run FBAwMC for all 6 LAB
[lca_f, lca_v, lca_modelIrrev] = FBAwMC_LAB(lca,0.345,biomass,num);<br />
[lfe_f, lfe_v, lfe_modelIrrev] = FBAwMC_LAB(lfe,0.363,biomass,num);<br />
[lla_f, lla_v, lla_modelIrrev] = FBAwMC_LAB(lla,0.436,biomass,num);<br />
[lme_f, lme_v, lme_modelIrrev] = FBAwMC_LAB(lme,0.272,biomass,num);<br />
[lsa_f, lsa_v, lsa_modelIrrev] = FBAwMC_LAB(lsa,0.446,biomass,num);<br />
[lpl_f, lpl_v, lpl_modelIrrev] = FBAwMC_LAB(lpl,0.31,biomass,num);<br />
