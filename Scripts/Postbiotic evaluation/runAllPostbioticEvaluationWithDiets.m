%[ProdRate] = PostbioticEvaluationWithDiets(model, C, biomass, num, targetProds, ex_rxns, lb_values)

% Load 6 LAB models
lca=load('lca.mat');
lfe=load('lfe.mat');
lla=load('lla.mat');
lme=load('lme.mat');
lsa=load('lsa.mat');
lpl=load('lpl.mat');

% Load VMH diet data
load("VMH diets.mat")

% Set parameters
biomass = 'biomass';
num = 5000;

% Load 6 LAB models

% Make cell array consisting list of postbiotics to be evaluated
targetProds = {'lac-L[c]', 'ppa[c]', 'ac[c]', 'diact[c]', 'udpg[c]', 'udpgal[c]','trp-L[c]','thm[c]', 'nac[c]', 'pydx[c]','btn[c]','fol[c]','pnto-R[c]','ribflv[c]', 'mqn7[c]', 'gthrd[c]', 'inost[c]','peptido[c]','CRB[c]','LTA[c]', 'ile-L[c]','val-L[c]','leu-L[c]','phe-L[c]', '3mba[c]', '2mpa[c]', '4abut[c]', 'ptrc[c]', 'spmd[c]', 'hista[c]', 'nh4[c]', 'indole[c]', 'phenol[c]', 'ch4s[c]'};



for i=1:1:length(diet_name)
    
    lca_ProdRate(i).lca_ProdRate = PostbioticEvaluationWithDiets(lca, 0.345, biomass, num, targetProds, ex_rxns, lb_values(:,i))
    lfe_ProdRate(i).lfe_ProdRate = PostbioticEvaluationWithDiets(lfe, 0.363, biomass, num, targetProds, ex_rxns, lb_values(:,i))
    lla_ProdRate(i).lla_ProdRate = PostbioticEvaluationWithDiets(lla, 0.436, biomass, num, targetProds, ex_rxns, lb_values(:,i))
    lme_ProdRate(i).lme_ProdRate = PostbioticEvaluationWithDiets(lme, 0.272, biomass, num, targetProds, ex_rxns, lb_values(:,i))
    lsa_ProdRate(i).lsa_ProdRate = PostbioticEvaluationWithDiets(lsa, 0.446, biomass, num, targetProds, ex_rxns, lb_values(:,i))
    lpl_ProdRate(i).lpl_ProdRate = PostbioticEvaluationWithDiets(lpl, 0.31, biomass, num, targetProds, ex_rxns, lb_values(:,i))
     
end