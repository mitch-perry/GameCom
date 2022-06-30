initCobraToolbox(false)
changeCobraSolver('ibm_cplex', 'LP');

global CBTDIR
iAF1260 = readCbModel([CBTDIR filesep 'test' filesep 'models' filesep 'mat' filesep 'iAF1260.mat']);

% I'm modifying the following reaction rates.
iAF1260.lb(find(iAF1260.lb > 0)) = 0;

% convert the compartment format from e.g., '_c' to '[c]'
iAF1260.mets = regexprep(iAF1260.mets, '_([^_]+)$', '\[$1\]');
% make all empty cells in cell arrays to be empty string
fieldToBeCellStr = {'metFormulas'; 'genes'; 'grRules'; 'metNames'; 'rxnNames'; 'subSystems'};
for j = 1:numel(fieldToBeCellStr)
iAF1260.(fieldToBeCellStr{j})(cellfun(@isempty, iAF1260.(fieldToBeCellStr{j}))) = {''};
end

iAF1260 = addReaction(iAF1260,{'METt3pp',''},'met__L[c] + h[c] => met__L[p] + h[p]');

argH = {'ARGSL'}; % essential for arginine biosynthesis
lysA = {'DAPDC'}; % essential for lysine biosynthesis
metA = {'HSST'}; % essential for methionine biosynthesis
ilvE = {'PPNDH'}; % essential for phenylalanine biosynthesis

argO = {'ARGt3pp'}; % Evidence for an arginine exporter encoded by yggA (argO) that is regulated by the LysR-type transcriptional regulator ArgP in Escherichia coli.
lysO = {'LYSt3pp'}; % Distinct paths for basic amino acid export in Escherichia coli: YbjE (LysO) mediates export of L-lysine
yjeH = {'METt3pp'}; % YjeH is a novel L-methionine and branched chain amino acids exporter in Escherichia coli
yddG = {'PHEt2rpp'}; % YddG from Escherichia coli promotes export of aromatic amino acids.

% auxotrophic for Lys and Met, not exporting Phe
Ec1 = iAF1260;
Ec1 = changeRxnBounds(Ec1, [lysA; metA; yddG], 0, 'b');
% auxotrophic for Arg and Phe, not exporting Met
Ec2 = iAF1260;
Ec2 = changeRxnBounds(Ec2, [argH; yjeH; ilvE], 0, 'b');
% Auxotrophic for Arg and Phe, not exporting Lys
Ec3 = iAF1260;
Ec3 = changeRxnBounds(Ec3, [argH; lysO; ilvE], 0, 'b');
% Auxotrophic for Lys and Met, not exporting Arg
Ec4 = iAF1260;
Ec4 = changeRxnBounds(Ec4, [argO; lysA; metA], 0, 'b');

% extracellular metabolites (met[e])
metEx = strcmp(getCompartment(iAF1260.mets),'e');
% the corresponding exchange reactions
rxnExAll = find(sum(iAF1260.S ~= 0, 1) == 1);
[rxnEx, ~] = find(iAF1260.S(metEx, rxnExAll)'); % need to be in the same order as metEx
rxnEx = rxnExAll(rxnEx);
% exchange rate
lbEx = iAF1260.lb(rxnEx);

nameTagsModel = {'Ec1'; 'Ec2'; 'Ec3'; 'Ec4'};
cd('/Users/mitchellperry/Dropbox/Microbiome (cleaned up)/matlabHelperFunctions');
EcCom = createMultipleSpeciesModel({Ec1; Ec2; Ec3; Ec4}, nameTagsModel);
EcCom.csense = char('E' * ones(1,numel(EcCom.mets))); % correct the csense
clear Ec1 Ec2 Ec3 Ec4

[EcCom.infoCom, EcCom.indCom] = getMultiSpeciesModelId(EcCom, nameTagsModel);
disp(EcCom.infoCom);

rxnBiomass = strcat(nameTagsModel, 'BIOMASS_Ec_iAF1260_core_59p81M'); % biomass reaction names
rxnBiomassId = findRxnIDs(EcCom, rxnBiomass); % ids
EcCom.infoCom.spBm = rxnBiomass; % .spBm for organism biomass reactions
EcCom.indCom.spBm = rxnBiomassId;

[yn, id] = ismember(strrep(iAF1260.mets(metEx), '[e]', '[u]'), EcCom.infoCom.Mcom); % map the metabolite name
assert(all(yn)); % must be a 1-to-1 mapping
EcCom.lb(EcCom.indCom.EXcom(:,1)) = lbEx(id); % assign community uptake bounds
EcCom.ub(EcCom.indCom.EXcom(:,1)) = 1e5;
EcCom.lb(EcCom.indCom.EXsp) = repmat(lbEx(id), 1, 4); % assign organism-specific uptake bounds

% only allow to take up the amino acids that one is auxotrophic for
exRate = 1; % maximum uptake rate for cross feeding AAs
% Ec1
EcCom = changeRxnBounds(EcCom, {'Ec1IEX_arg__L[u]tr'; 'Ec1IEX_phe__L[u]tr'}, 0, 'l');
EcCom = changeRxnBounds(EcCom, {'Ec1IEX_met__L[u]tr'; 'Ec1IEX_lys__L[u]tr'}, -exRate, 'l');
% Ec2
EcCom = changeRxnBounds(EcCom, {'Ec2IEX_arg__L[u]tr'; 'Ec2IEX_phe__L[u]tr'}, -exRate, 'l');
EcCom = changeRxnBounds(EcCom, {'Ec2IEX_met__L[u]tr'; 'Ec2IEX_lys__L[u]tr'}, 0, 'l');
% Ec3
EcCom = changeRxnBounds(EcCom, {'Ec3IEX_arg__L[u]tr'; 'Ec3IEX_phe__L[u]tr'}, -exRate, 'l');
EcCom = changeRxnBounds(EcCom, {'Ec3IEX_met__L[u]tr'; 'Ec3IEX_lys__L[u]tr'}, 0, 'l');
% Ec4
EcCom = changeRxnBounds(EcCom, {'Ec4IEX_arg__L[u]tr'; 'Ec4IEX_phe__L[u]tr'}, 0, 'l');
EcCom = changeRxnBounds(EcCom, {'Ec4IEX_met__L[u]tr'; 'Ec4IEX_lys__L[u]tr'}, -exRate, 'l');
% allow production of anything for each member
EcCom.ub(EcCom.indCom.EXsp(:)) = 1000;

% Now make 4 E.Coli model to be used in Python code.
[I,J] = size(EcCom.S);
S = EcCom.S;
lb = EcCom.lb;
ub = EcCom.ub;

cd('/Users/mitchellperry/Dropbox/Microbiome (cleaned up)/EColiModelforPython')

save('S.mat', 'S');
save('I.mat', 'I');
save('J.mat', 'J');
save('lb.mat', 'lb');
save('ub.mat', 'ub');

lumen_reactions_idx = find(EcCom.indCom.rxnSps == 0);
lumen_metabolites_idx = find(EcCom.indCom.metSps == 0);
lumen_reactions = EcCom.rxns(lumen_reactions_idx);

save('lumen_reactions_idx.mat', 'lumen_reactions_idx');
save('lumen_metabolites_idx.mat', 'lumen_metabolites_idx');
save('lumen_reactions.mat', 'lumen_reactions');

Ec1_reactions_idx = find(EcCom.indCom.rxnSps == 1);
Ec1_metabolites_idx = find(EcCom.indCom.metSps == 1);
Ec1_reactions = EcCom.rxns(Ec1_reactions_idx);
Ec1_biomass_idx = find(strcmp(Ec1_reactions, EcCom.infoCom.spBm(1)) == 1);

save('Ec1_reactions_idx.mat', 'Ec1_reactions_idx');
save('Ec1_metabolites_idx.mat', 'Ec1_metabolites_idx');
save('Ec1_biomass_idx.mat', 'Ec1_biomass_idx');
save('Ec1_reactions.mat', 'Ec1_reactions');

Ec2_reactions_idx = find(EcCom.indCom.rxnSps == 2);
Ec2_metabolites_idx = find(EcCom.indCom.metSps == 2);
Ec2_reactions = EcCom.rxns(Ec2_reactions_idx);
Ec2_biomass_idx = find(strcmp(Ec2_reactions, EcCom.infoCom.spBm(2)) == 1);

save('Ec2_reactions_idx.mat', 'Ec2_reactions_idx');
save('Ec2_metabolites_idx.mat', 'Ec2_metabolites_idx');
save('Ec2_biomass_idx.mat', 'Ec2_biomass_idx');
save('Ec2_reactions.mat', 'Ec2_reactions');

Ec3_reactions_idx = find(EcCom.indCom.rxnSps == 3);
Ec3_metabolites_idx = find(EcCom.indCom.metSps == 3);
Ec3_reactions = EcCom.rxns(Ec3_reactions_idx);
Ec3_biomass_idx = find(strcmp(Ec3_reactions, EcCom.infoCom.spBm(3)) == 1);

save('Ec3_reactions_idx.mat', 'Ec3_reactions_idx');
save('Ec3_metabolites_idx.mat', 'Ec3_metabolites_idx');
save('Ec3_biomass_idx.mat', 'Ec3_biomass_idx');
save('Ec3_reactions.mat', 'Ec3_reactions');

Ec4_reactions_idx = find(EcCom.indCom.rxnSps == 4);
Ec4_metabolites_idx = find(EcCom.indCom.metSps == 4);
Ec4_reactions = EcCom.rxns(Ec4_reactions_idx);
Ec4_biomass_idx = find(strcmp(Ec4_reactions, EcCom.infoCom.spBm(4)) == 1);

save('Ec4_reactions_idx.mat', 'Ec4_reactions_idx');
save('Ec4_metabolites_idx.mat', 'Ec4_metabolites_idx');
save('Ec4_biomass_idx.mat', 'Ec4_biomass_idx');
save('Ec4_reactions.mat', 'Ec4_reactions');







yada = find(strcmp(EcCom.rxns, 'Ec2METt3pp'));
EcCom.subSystems(yada)

