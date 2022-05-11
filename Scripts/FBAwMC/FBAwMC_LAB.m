function [f, v, modelIrrev] = FBAwMC_LAB(model,C,biomass,num)

% INPUT
%  model             COBRA model structure with 3 additional vectors of same size as 'rxns':
%                    kcat_f, kcat_b, molwt (if any of the value unknown,
%                    provide '0')
%                    kcat units should be '1/s' and mol wt in 'Dalton'
%  crd_val           Cytoplasmic density 
%  biomass           name of biomass reaction (to be excluded from enzyme
%                    capacity flux constraint)
%  num               number of flux solutions for computation
%
% 
% OUTPUT
%  modelIrrev    Model in irreversible format with enzyme constraint added
%                as pseudo reaction
%  f             Objective value
%  v             Reaction rates
%
%
% Lokanand Koduru            10/03/18
% Meiyappan Lakshmanan       10/04/18 Generalized the code

        %% Generate crowd positions
        [modelIrrev, matchRev, irrev2rev, solution] = generateCrowdPositions(model,C,biomass,num);
        
        %% Identify best crowd positions with non-zero growth and lactate production
        NonZero_BiomassCrdInds = find((solution(find(ismember(irrev2rev,find(ismember(model.rxns,biomass)))),:)));
        L_LactateExchInd = find(ismember(modelIrrev.rxns,{'EX_lac-L(e)_f'}));
        D_LactateExchInd = find(ismember(modelIrrev.rxns,{'EX_lac-D(e)_f'}));
        NonZero_L_LactateCrdInds=find((solution(L_LactateExchInd,:)));
        NonZero_D_LactateCrdInds=find((solution(D_LactateExchInd,:)));
        NonZero_Biomass_Lactate_CrdInds = intersect(NonZero_BiomassCrdInds,union(NonZero_L_LactateCrdInds,NonZero_D_LactateCrdInds));
        [Nbins,Edges] = histcounts(log(solution(1,NonZero_Biomass_Lactate_CrdInds)),80);
        minGR=Edges((find(Nbins == max(Nbins))-2));
        maxGR=Edges((find(Nbins == max(Nbins))));
        MinMaxGRIndsOf_NonZero_Biomass_Lactate_CrdInds=[];
        for h=1:length(NonZero_Biomass_Lactate_CrdInds)
            if log(solution(1,NonZero_Biomass_Lactate_CrdInds(h))) >= minGR
                if log(solution(1,NonZero_Biomass_Lactate_CrdInds(h))) <= maxGR
                    MinMaxGRIndsOf_NonZero_Biomass_Lactate_CrdInds=[MinMaxGRIndsOf_NonZero_Biomass_Lactate_CrdInds,NonZero_Biomass_Lactate_CrdInds(h)];
                end
            end
        end
        Best_Crd_positions = MinMaxGRIndsOf_NonZero_Biomass_Lactate_CrdInds;
        
        %% Objective value, averaged from flux solutions
        Solutions_GR = solution(find(ismember(modelIrrev.rxns,biomass)),Best_Crd_positions);
        f = [mean(Solutions_GR,2)];
        
        %% Reaction rates for all reactions in IrrevModel, averaged from flux solutions
        Solutions = solution(:,Best_Crd_positions);
        v = [mean(Solutions,2)];
        
        
end












