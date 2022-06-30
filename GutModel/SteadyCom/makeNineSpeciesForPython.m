initCobraToolbox(false);
load('model9.mat');
model9.indCom = SteadyComSubroutines('infoCom2indCom', model9);
% I"m modifying the following reaction rates;
% don't force any internal reactions for any of the species.
model9.lb(find((model9.lb > 0) & (model9.indCom.rxnSps > 0))) = 0;
% Now get model parameters to be used in Python code.
S = model9.S;
lb = model9.lb;
ub = model9.ub;

lumen_reactions_idx = find(model9.indCom.rxnSps == 0);
lumen_metabolites_idx = find(model9.indCom.metSps == 0);
lumen_reactions = model9.rxns(lumen_reactions_idx);

lumen_uptake_reactions_idx = lumen_reactions_idx(1:558);
lumen_export_reactions_idx = lumen_reactions_idx(559:end);
lumen_uptake_reactions = lumen_reactions(1:558);
lumen_export_reactions = lumen_reactions(559:end);

lb(lumen_uptake_reactions_idx) = -ub(lumen_uptake_reactions_idx);
ub(lumen_uptake_reactions_idx) = ub(lumen_export_reactions_idx);
lb(lumen_export_reactions_idx) = [];

S(:,lumen_uptake_reactions_idx) = -S(:,lumen_uptake_reactions_idx);
S(:,lumen_export_reactions_idx) = [];

[I,J] = size(S);


save('S.mat', 'S');
save('I.mat', 'I');
save('J.mat', 'J');
save('lb.mat', 'lb');
save('ub.mat', 'ub');

save('lumen_reactions_idx.mat', 'lumen_uptake_reactions_idx');
save('lumen_metabolites_idx.mat', 'lumen_metabolites_idx');
save('lumen_reactions.mat', 'lumen_uptake_reactions');

Bt_reactions_idx = find(model9.indCom.rxnSps == 1);
Bt_metabolites_idx = find(model9.indCom.metSps == 1);
Bt_reactions = model9.rxns(Bt_reactions_idx);
Bt_biomass_idx = find(strcmp(Bt_reactions, model9.infoCom.spBm(1)) == 1);

save('Bt_reactions_idx.mat', 'Bt_reactions_idx');
save('Bt_metabolites_idx.mat', 'Bt_metabolites_idx');
save('Bt_biomass_idx.mat', 'Bt_biomass_idx');
save('Bt_reactions.mat', 'Bt_reactions');

Fp_reactions_idx = find(model9.indCom.rxnSps == 2);
Fp_metabolites_idx = find(model9.indCom.metSps == 2);
Fp_reactions = model9.rxns(Fp_reactions_idx);
Fp_biomass_idx = find(strcmp(Fp_reactions, model9.infoCom.spBm(2)) == 1);

save('Fp_reactions_idx.mat', 'Fp_reactions_idx');
save('Fp_metabolites_idx.mat', 'Fp_metabolites_idx');
save('Fp_biomass_idx.mat', 'Fp_biomass_idx');
save('Fp_reactions.mat', 'Fp_reactions');

Kp_reactions_idx = find(model9.indCom.rxnSps == 3);
Kp_metabolites_idx = find(model9.indCom.metSps == 3);
Kp_reactions = model9.rxns(Kp_reactions_idx);
Kp_biomass_idx = find(strcmp(Kp_reactions, model9.infoCom.spBm(3)) == 1);

save('Kp_reactions_idx.mat', 'Kp_reactions_idx');
save('Kp_metabolites_idx.mat', 'Kp_metabolites_idx');
save('Kp_biomass_idx.mat', 'Kp_biomass_idx');
save('Kp_reactions.mat', 'Kp_reactions');

St_reactions_idx = find(model9.indCom.rxnSps == 4);
St_metabolites_idx = find(model9.indCom.metSps == 4);
St_reactions = model9.rxns(St_reactions_idx);
St_biomass_idx = find(strcmp(St_reactions, model9.infoCom.spBm(4)) == 1);

save('St_reactions_idx.mat', 'St_reactions_idx');
save('St_metabolites_idx.mat', 'St_metabolites_idx');
save('St_biomass_idx.mat', 'St_biomass_idx');
save('St_reactions.mat', 'St_reactions');

Ba_reactions_idx = find(model9.indCom.rxnSps == 5);
Ba_metabolites_idx = find(model9.indCom.metSps == 5);
Ba_reactions = model9.rxns(Ba_reactions_idx);
Ba_biomass_idx = find(strcmp(Ba_reactions, model9.infoCom.spBm(5)) == 1);

save('Ba_reactions_idx.mat', 'Ba_reactions_idx');
save('Ba_metabolites_idx.mat', 'Ba_metabolites_idx');
save('Ba_biomass_idx.mat', 'Ba_biomass_idx');
save('Ba_reactions.mat', 'Ba_reactions');

Ec_reactions_idx = find(model9.indCom.rxnSps == 6);
Ec_metabolites_idx = find(model9.indCom.metSps == 6);
Ec_reactions = model9.rxns(Ec_reactions_idx);
Ec_biomass_idx = find(strcmp(Ec_reactions, model9.infoCom.spBm(6)) == 1);

save('Ec_reactions_idx.mat', 'Ec_reactions_idx');
save('Ec_metabolites_idx.mat', 'Ec_metabolites_idx');
save('Ec_biomass_idx.mat', 'Ec_biomass_idx');
save('Ec_reactions.mat', 'Ec_reactions');

Ef_reactions_idx = find(model9.indCom.rxnSps == 7);
Ef_metabolites_idx = find(model9.indCom.metSps == 7);
Ef_reactions = model9.rxns(Ef_reactions_idx);
Ef_biomass_idx = find(strcmp(Ef_reactions, model9.infoCom.spBm(7)) == 1);

save('Ef_reactions_idx.mat', 'Ef_reactions_idx');
save('Ef_metabolites_idx.mat', 'Ef_metabolites_idx');
save('Ef_biomass_idx.mat', 'Ef_biomass_idx');
save('Ef_reactions.mat', 'Ef_reactions');

Er_reactions_idx = find(model9.indCom.rxnSps == 8);
Er_metabolites_idx = find(model9.indCom.metSps == 8);
Er_reactions = model9.rxns(Er_reactions_idx);
Er_biomass_idx = find(strcmp(Er_reactions, model9.infoCom.spBm(8)) == 1);

save('Er_reactions_idx.mat', 'Er_reactions_idx');
save('Er_metabolites_idx.mat', 'Er_metabolites_idx');
save('Er_biomass_idx.mat', 'Er_biomass_idx');
save('Er_reactions.mat', 'Er_reactions');

Lc_reactions_idx = find(model9.indCom.rxnSps == 9);
Lc_metabolites_idx = find(model9.indCom.metSps == 9);
Lc_reactions = model9.rxns(Lc_reactions_idx);
Lc_biomass_idx = find(strcmp(Lc_reactions, model9.infoCom.spBm(9)) == 1);

save('Lc_reactions_idx.mat', 'Lc_reactions_idx');
save('Lc_metabolites_idx.mat', 'Lc_metabolites_idx');
save('Lc_biomass_idx.mat', 'Lc_biomass_idx');
save('Lc_reactions.mat', 'Lc_reactions');
