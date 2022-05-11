function model_out = constrainDiets(model,ex_rxns,lb_values)
% INPUT
%  model             COBRA model structure with 3 additional vectors of same size as 'rxns':
%                    kcat_f, kcat_b, molwt (if any of the value unknown,
%                    provide '0')
%                    kcat units should be '1/s' and mol wt in 'Dalton'
%  ex_rxns           List of exchange reactions for metabolites present in target diet
%  lb_values         Lower bound values corresponding to ex_rxns
% 
% OUTPUT
%  model_out         Model constrained with target diet        
%
%
% Lokanand Koduru            10/03/18
% Meiyappan Lakshmanan       10/04/18 Generalized the code
    
    model_out = model;
    
    %% Change bounds to constrain diet
    for i=1:1:length(ex_rxns)
        model_out = changeRxnBounds(model_out,ex_rxns(i),lb_values(i),'l');

    end
    
    % Make sure that ions are freely available
    model_out = changeRxnBounds(model_out,{'EX_h2o(e)','EX_ca2(e)','EX_cl(e)','EX_k(e)','EX_mg(e)','EX_na1(e)','EX_pi(e)','EX_cu2(e)','EX_fe2(e)','EX_fe3(e)','EX_mn2(e)','EX_zn(e)'},-1000,'l');
end