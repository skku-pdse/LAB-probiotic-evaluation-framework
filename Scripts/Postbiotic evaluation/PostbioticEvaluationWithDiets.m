function [ProdRate] = PostbioticEvaluationWithDiets(model, C, biomass, num, targetProds, ex_rxns, lb_values)
% INPUT
%  model             COBRA model structure with 3 additional vectors of same size as 'rxns':
%                    kcat_f, kcat_b, molwt (if any of the value unknown,
%                    provide '0')
%                    kcat units should be '1/s' and mol wt in 'Dalton'
%  C                 Cytoplasmic density
%  biomass           Name of biomass reaction (to be excluded from enzyme capacity flux constraint)
%  num               Number of flux solutions for computation
%  targetProds       Cell array consisting list of postbiotics to be evaluated
%  ex_rxns           List of exchange reactions for metabolites present in target diet
%  lb_values         Lower bound values corresponding to ex_rxns
% 
% OUTPUT
% ProdRate           Estimated postbiotic production rate  
%
%
% Lokanand Koduru            10/03/18
% Meiyappan Lakshmanan       10/04/18 Generalized the code

temp_model1 = addReaction(model,'EX_TempProd',{'A'},[-1],false);

% Metabolite name conformation
new_mets = {'LTA[c]','LTA[c]','LTA[c]','LTA[c]','LTA[c]','LTA[c]','LTA[c]','peptido[c]','peptido[c]','peptido[c]','peptido[c]','peptido[c]','peptido[c]','CRB[c]','CRB[c]','CRB[c]','CRB[c]','CRB[c]','CRB[c]','CRB[c]'};
original_mets = {'LTA_LCA[c]','LTA_LPL[c]','LTA_LME[c]','LTA_LSA[c]','LTA_LFE[c]','LTAAlaGal_LLA[c]','LTAtotal_LRE[c]','peptido_LCA[c]','peptido_LPL[c]','peptido_LME[c]','peptido_LSA[c]','peptido_LFE[c]','peptido_LRE[c]','CRB_LCA[c]','CPS_LPL2[c]','CRB_LME[c]','CRB_LSA[c]','CRB_LFE[c]','CPS_LLA[c]','CPS_LRE[c]'};
for m=1:length(new_mets)
    met = original_mets(1,m);
    A = find(ismember(temp_model1.mets,met));
    if ~isempty(A)
        temp_model1.mets(A,1) = new_mets(1,m);
    end
end

temp_model = temp_model1;
temp_model = constrainDiets(temp_model,ex_rxns,lb_values);
temp_model = changeRxnBounds(temp_model,'EX_o2(e)',0,'l');
temp_model = changeRxnBounds(temp_model,'ATPH',0,'l');

% Generate crowd positions
% Identify best crowding positions with non-zero growth and lactate production
[f, v, modelIrrev, Best_Crd_positions, EnzConSamples] = FBAwMC_LAB(temp_model,C,biomass,num);

%% Postbiotic evaluation with model constrained with diet and calculated best crowding positions
[ProdRate] = postbioticEvaluation(temp_model, biomass, targetProds, EnzConSamples, Best_Crd_positions);


    
