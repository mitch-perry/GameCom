% Load Cobra Toolbox and set to compatible solver (CPLEX).
initCobraToolbox(false)
changeCobraSolver('ibm_cplex', 'LP');

% Load 9 species microbiome model from SteadyCom paper.
load('model9.mat');

% Format for use with SteadyCom function from Cobra Toolbox.
model9.indCom = SteadyComSubroutines('infoCom2indCom', model9);

% I"m modifying the following reaction rates;
% don't force any internal reactions for any of the species.
model9.lb(find((model9.lb > 0) & (model9.indCom.rxnSps > 0))) = 0;

[sol, result] = SteadyCom(model9);

biomass = result.BM;
fluxes = result.flux;

lumen_reactions_idx = find(model9.indCom.rxnSps == 0);
lumen_metabolites_idx = find(model9.indCom.metSps == 0);
lumen_reactions = model9.rxns(lumen_reactions_idx);

lumen_uptake_reactions_idx = lumen_reactions_idx(1:558);
lumen_export_reactions_idx = lumen_reactions_idx(559:end);
lumen_uptake_reactions = lumen_reactions(1:558);
lumen_export_reactions = lumen_reactions(559:end);

fluxes(lumen_uptake_reactions_idx) = fluxes(lumen_export_reactions_idx) - fluxes(lumen_uptake_reactions_idx);
fluxes(lumen_export_reactions_idx) = [];

save('biomasses.mat', 'biomass');
save('fluxes.mat', 'fluxes');
