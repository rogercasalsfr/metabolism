%% 
% _(c) Dr. Maria Pires Pacheco 2016_
% 
% _Example script adapted by Tamara Bintener 2021_
%% Disclaimer
% In order to run FASTCORMICS RNA-seq you will need
% 
% SOFTWARE:
%% 
% * Matlab
% * Compatible IBM ILOG Cplex version, added to the Matlab path
% * Statistics and Machine Learning Toolbox
% * Curve fitting toolbox
% * Cobra Toolbox installed (https://opencobra.github.io/cobratoolbox/latest/installation.html)
%% 
% DATA: 
% 
% The example data was downloaded from <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM1697009 
% https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM1697009> (FPKM values 
% of normal tissue). The first two samples (TCGA06067511A32RA36H07 and TCGA06067811A32RA36H07) 
% will serve as an example was already imported, split into 3 variables, and saved 
% under _example_FPKM.mat ._
%% 
% * colnames:   cell array with the sample names (here TCGA06067511A32RA36H07)
% * rownames:  cell aray with the gene IDs
% * fpkm:           fpkm values for the samples, size(fpkm) = lenghth(geneIDs) 
% length(colnames)
%% 
% The model, dico, and an example medium composition is provided in the exampleData 
% folder.
%% 
% * model:		(consistent) genome-scale metabolic reconstruction in the COBRA 
% format, i.e. Recon 2.04 (from https://vmh.uni.lu/#downloadview )
% * dico:		table which contains corresponding gene identifier information. Needed 
% to map the rownames to the genes in the model. Can be manually assembled in 
% https://www.ensembl.org/biomart/martview and imported into Matlab.
% * medium:		[optional] defines metabolites in the growth medium of cells to 
% constrain the model, see example medium_example.mat
%% Setup
% First, all the needed files should be added to the path, therefore navigate 
% to the _rFASTCORMICS for RNA-seq data_ folder.

addpath(genpath(pwd)) % add all subfolders to path in the current dicrectory.      

feature astheightlimit 2000      


addpath(genpath("/opt/ibm/ILOG/CPLEX_Studio1210/cplex/matlab/x86-64_linux/"))     
addpath("/home/claudia/cobratoolbox/")   
initCobraToolbox(false)  
addpath(genpath("/home/claudia/rFASTCORMICS/"))



data = readtable('DESeq2_TPM_values.txt', "ReadRowNames",true);
fpkm = table2array(data);
rownames = data.Properties.RowNames;
colnames = data.Properties.VariableNames;




%% Data discretization
% The FPKM values will be discretized into expressed, not expressed, and unknown 
% expression status. The discretized data is one input for rFASTCORMICS.

discretized = discretize_FPKM(fpkm, colnames) % no figures
% discretized = discretize_FPKM(fpkm, colnames, 1) % with figures
%% rFASTCORMICS
% Prepare rFASTCORMICS
%% Needed inputs
%% 
% * *Cmodel*:         a consistent model, here: obtained from running fastcc 
% on Recon204 (from <https://vmh.uni.lu/#downloadview) https://vmh.uni.lu/#downloadview)>

load Recon3DModel_301.mat, 
Recon3DModel.description='Recon';
A = fastcc_4_rfastcormics(model, 1e-4,0) % create consistent model by running FASTCC (Vlassis et al., 2014)
Cmodel = removeRxns(model, model.rxns(setdiff(1:numel(model.rxns),A)))
%% 
% * *discretized*:    FPKM values that have been discretized in the previous 
% step
% * *rownames*:        gene names from the FPKM data
% * *dico*:               identifier file, which links the geneIDs from the 
% RNA-seq data to the gene identifiers in the model. Can be created with BioMart 
% and imported to Matlab

load dico_rFASTCORMICS.mat % dictionary to map the rownname identifier to the genes in the model
%% Optional inputs
%% 
% * *biomass_rxn*: name of the biomass reaction in the model if present (check 
% model.rxns). Setting this variable will always enable the biomass to carry a 
% flux

biomass_rxn = {'biomass_maintenance'}; %find the corresponding biomass reaction in your model, here it is for recon 2
%% 
% * *already mapped tag:* 1, if the data was already to the model.rxns in this 
% case data p = n  and  0, if the data has to be mapped  using the GPR rules of 
% the model

already_mapped_tag = 0;
%% 
% * *epsilon*        avoids small number errors

epsilon = 1e-4; %avoid small number errors
%% 
% * *consensus_proportion*

consensus_proportion = 0.9; %gene has to be expressed in 90% of the cases in order to be included. Only relevant if you want to create one generic model from different samples
%% 
% Set *optional_settings* such as:
%% 
% * *unpenalizedSystems*      
% * *func*: reaction(s) forced to be present in the model
% * *not_medium_constrained*
% * *medium* : medium composition, defines metabolites in the growth medium 
% of cells to constrain the model, see example medium_example.mat

load medium_example % need to define medium for the cells used here

unpenalizedSystems = {'Transport, endoplasmic reticular';
    'Transport, extracellular';
    'Transport, golgi apparatus';
    'Transport, mitochondrial';
    'Transport, peroxisomal';
    'Transport, lysosomal';
    'Transport, nuclear'};
Cmodel_subSystems = cellfun(@(x) x{1}, Cmodel.subSystems, 'UniformOutput', false); % unpenalizedSystems cells were char and model.subSystems cells were cells
unpenalized = Cmodel.rxns(ismember(Cmodel_subSystems,unpenalizedSystems));

optional_settings.unpenalized = unpenalized; % biomass rx name = biomass_maintenance
optional_settings.func = {'biomass_maintenance','DM_atp_c_'}; % forced additional reactions into the  model
not_medium_constrained = 'EX_tag_hs[e]';
optional_settings.not_medium_constrained = not_medium_constrained;
optional_settings.medium = medium_example;% remove field if no constraint is provided
%% Create models
% the cobra toolbox is needed for the model creation
% 
% make sure you have a compatible cplex version for your Matlab version
% 
% single models:

for i = 1:numel(colnames) %for each sample
    [ContextModel, A_keep] = fastcormics_RNAseq(Cmodel, discretized(:,i), ...
        rownames, dico, biomass_rxn, already_mapped_tag, consensus_proportion, epsilon, optional_settings);
     filePath = fullfile('results_adrian', strcat('SampleModel_', colnames{i}));
     save(filePath, 'ContextModel');   
        models_keep(A_keep,i) = 1;
end
%% 
% consensus models:

[model_out_consensus, A_keep_consensus] = fastcormics_RNAseq(Cmodel, discretized, rownames, dico, ...
    biomass_rxn, already_mapped_tag, consensus_proportion, epsilon, optional_settings);
%%
delete clone*.log %delete some files generated by cplex
