function [ProdRate] = postbioticEvaluation(model, biomass, targetProds, EnzConSamples, Best_Crd_positions)
% INPUT
%  model                COBRA model structure with 3 additional vectors of same size as 'rxns':
%                       kcat_f, kcat_b, molwt (if any of the value unknown,
%                       provide '0')
%                       kcat units should be '1/s' and mol wt in 'Dalton'
%  biomass              Name of biomass reaction (to be excluded from enzyme capacity flux constraint)
%  targetProds          Cell array consisting list of postbiotics to be evaluated
%  EnzConSamples        Enzyme crowding coefficients
%  Best_Crd_positions   Best crowding positions with non-zero growth and lactate production
% 
% OUTPUT
%  ProdRate             Estimated postbiotic production rate        
%
%
% Lokanand Koduru            10/03/18
% Meiyappan Lakshmanan       10/04/18 Generalized the code

model1 = convertToIrreversible(model);
model1.description = 'model_with_CRD';
revFlag = 'false';
lowerBound = 0;
upperBound = 1;
model1 = addReaction(model1,'Crd_C_Rxn',{'Crd_Met[c]'},-1,revFlag,lowerBound,upperBound,0,'','','','');
Crd_MetInd = find(ismember(model1.mets,'Crd_Met[c]'));
selExc = findExcRxns(model1);
ExRxnInd = find(selExc);
BiomassRxnInd = find(ismember(model1.rxns,biomass));
ATPMRxnInd = find(ismember(model1.rxns,'ATPH'));
ProdRxnInd = find(ismember(model1.rxns,'EX_TempProd'));

Solutions = zeros(length(model1.rxns),length(Best_Crd_positions));
Solutions_GR = zeros(length(Best_Crd_positions),1);
ProdRate = zeros(length(Best_Crd_positions),length(targetProds));
model2 = model1;
for p = 1:length(Best_Crd_positions)
    model1.S(Crd_MetInd,:) = transpose(EnzConSamples(:,Best_Crd_positions(p)));
    model1.S(Crd_MetInd,BiomassRxnInd) = 0;
    model1.S(Crd_MetInd,BiomassRxnInd) = 0;
    model1.S(Crd_MetInd,ATPMRxnInd) = 0;
    model1.S(Crd_MetInd,ExRxnInd) = 0;
    model1.S(Crd_MetInd,ProdRxnInd) = 0;
    model1.S(Crd_MetInd,find(ismember(model1.rxns,'Crd_C_Rxn'))) = -1;
    
    % Removing the target metabolite from biomass objective before maximization
    model1.S(find(model1.S(:,ProdRxnInd)),BiomassRxnInd)=0;
    model1 = changeObjective(model1,biomass);

    sol = optimizeCbModel(model1);
    Solutions(:,p) = sol.x;
    Solutions_GR(p) = sol.f;
    
    % Constrain biomass to 50% of wild type growth
    model1.lb(BiomassRxnInd) = 0.5*Solutions_GR(p);
    model3 = model1;
    for q = 1:length(targetProds)
        Trgt_MetInd = find(ismember(model1.mets,targetProds(q)));    
        model1 = changeObjective(model1,'EX_TempProd');
        model1.S(:,ProdRxnInd) = 0;
        model1.S(Trgt_MetInd,ProdRxnInd) = -1;
        sol = optimizeCbModel(model1);
        ProdRate(p,q) = sol.f;
        model1 = model3;
    end
    model1 = model2;
end
end