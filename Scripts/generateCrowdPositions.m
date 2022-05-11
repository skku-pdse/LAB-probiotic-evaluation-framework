function [modelIrrev, matchRev, irrev2rev, solution] = generateCrowdPositions(model,C,biomass,num)
% INPUT
%  model             COBRA model structure with 3 additional vectors of same size as 'rxns':
%                    kcat_f, kcat_b, molwt (if any of the value unknown,
%                    provide '0')
%                    kcat units should be '1/s' and mol wt in 'Dalton'
%  crd_val           crowding coefficient, or fraction of enzymatic mass in overall dry cell weight 
%  biomass           name of biomass reaction (to be excluded from enzyme
%                    capacity flux constraint)
%  num               number of flux solutions for computation
% 
% OUTPUT
%  modelIrrev        Model in irreversible format with enzyme constraint added as pseudo reaction                
%  matchRev          Matching of forward and backward reactions of a reversible reaction
%  irrev2rev         Matching from irreversible to reversible reactions
%  solution          Flux solutions
%
%
% Lokanand Koduru            10/03/18
% Meiyappan Lakshmanan       10/04/18 Generalized the code

%% Convert to irreversible format
[modelIrrev,matchRev,~,irrev2rev] = convertToIrreversible(model);
modelIrrev.description='model_with_EnzCon';
revFlag='false';

%% Identify reactions with kcat values
kcat_rxns = model.rxns((model.kcat)~=0);

%% Filter out the mol wt of reactions with kcat values
mw = model.mw(find(ismember(model.rxns,kcat_rxns)));

%% Identify reactions in Irreversible model with kcat values
kcat_rxns_irrev = modelIrrev.rxns(ismember(irrev2rev,find(ismember(model.rxns,kcat_rxns))));

%% Assign kcat values to reactions in Irreversible model with kcat values
kcat_irrev = zeros(length(kcat_rxns_irrev),1);

for i=1:1:length(kcat_rxns_irrev)
    if model.kcat(irrev2rev(find(ismember(modelIrrev.rxns,kcat_rxns_irrev(i))))) ~= 0
        kcat_irrev(i) = model.kcat(irrev2rev(find(ismember(modelIrrev.rxns,kcat_rxns_irrev(i)))));
    end
end
for i=1:1:length(kcat_rxns_irrev)
    mw_irrev(i) = model.mw(irrev2rev(find(ismember(modelIrrev.rxns,kcat_rxns_irrev(i)))));
end
kcat_mw = zeros(length(kcat_rxns_irrev),1);

%% Calculate kcat/mw values
for i=1:1:length(kcat_rxns_irrev)
    kcat_mw(i) = (((mw_irrev(i)/1000)*0.73)/(kcat_irrev(i)*3600))*C;
end

kcat_rxns_irrev = kcat_rxns_irrev(~isinf(kcat_mw));
kcat_mw = kcat_mw(~isinf(kcat_mw));
kcat_rxns_irrev = kcat_rxns_irrev(kcat_mw~=0);
kcat_mw = kcat_mw(kcat_mw~=0);
kcat_rxns_irrev = kcat_rxns_irrev(~isnan(kcat_mw));
kcat_mw = kcat_mw(~isnan(kcat_mw));

%% Add dummy reaction to model representing the enzyme constraint
lowerBound=0;
upperBound=1;
modelIrrev=addReaction(modelIrrev,'EnzCon_C_Rxn',{'EnzCon_Met[c]'},-1,revFlag,lowerBound,upperBound,0,'','','','');
EnzCon_MetInd=find(ismember(modelIrrev.mets,'EnzCon_Met[c]'));
Rxn_WithEnzConInd=find(ismember(modelIrrev.rxns,kcat_rxns_irrev));
selExc=findExcRxns(modelIrrev);
ExRxnInd=find(selExc);
BiomassRxnInd=find(ismember(modelIrrev.rxns,biomass));
Solutions_EnzConSamples=zeros(length(modelIrrev.rxns),num);
EnzConSamples=zeros(length(modelIrrev.rxns),num);
RandomEnzConCoeff=zeros(length(modelIrrev.rxns),num);

for i=1:1:num
    RandomEnzConCoeff(:,i)=randsample(kcat_mw,length(modelIrrev.rxns),true);
end

environment=getEnvironment();

parfor j=1:num
    
    restoreEnvironment(environment);
    
    modelIrrev1 = modelIrrev;
    kcat_mw1 = kcat_mw;
    modelIrrev1.S(EnzCon_MetInd,:) = transpose(RandomEnzConCoeff(:,j));
    
    %% Removing EnzCon coeff for biomass equation
    modelIrrev1.S(EnzCon_MetInd,BiomassRxnInd)=0;

    %% Replacing the randomly assigned EnzCon coeff with the original coefficients for those reactions whose EnzCon coeffs are available (calculated list) 
    for k=1:length(Rxn_WithEnzConInd)
         modelIrrev1.S(EnzCon_MetInd,Rxn_WithEnzConInd(k))=kcat_mw1(find(ismember(kcat_rxns_irrev,modelIrrev1.rxns(Rxn_WithEnzConInd(k)))));
    end
    
    %% Removing EnzCon coeff for exchange reactions
    modelIrrev1.S(EnzCon_MetInd,ExRxnInd)=0;
    modelIrrev1.S(EnzCon_MetInd,find(ismember(modelIrrev1.rxns,'EnzCon_C_Rxn')))=-1;
    solution=optimizeCbModel(modelIrrev1,'max');
    if(isempty(solution.x))==0
        Solutions_EnzConSamples(:,j)=solution.x;
        EnzConSamples(:,j) = modelIrrev1.S(EnzCon_MetInd,:)';
    end
    
end
solution = Solutions_EnzConSamples;

end