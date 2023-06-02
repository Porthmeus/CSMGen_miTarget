function [mod_tiss] = extractModels(model_in, model_out, diet, consistent_rxns, core_rxns, sample, cobra_path, zero_cutoff, objective)
    % read the model, extract the genes from the tpm matrix and calculate a
    % transcriptpome specific model per column in the TPM_matrix. Write these
    % models to an SBML file

    % @ model_in - a file which can be read to a CBT model by readCbModel()
    % @ model_out - the file where the results should be stored as a csv file for the list of reactions in the model, which are tissue specific
    % @ diet - a csv file with constrains for the diet - the first column should contain the rxns names, the second column values for the fluxes of these reaction. The function will set the lower bounds to all exchange reaction not in that list. For the reactions in the list the lower bound will be set to -1*value. If such a file is not given, it is assumed, that the model is correctly constrained.
    % @ consistent_rxns - a table of the reaction which make up the consistent model. The column named "rxn" should contain the reaction-IDs of the consistent model. If no file is given the model is assumed to be consistent.
    % @ core_rxns - a csv file which contains the core reactions to be extracted. Can contain several columns with core reactions, @sample will define reactions which are considered core.
    % @ sample - the column in the table of core_rxns which should be used as core reactions. Default is the first column
    % @ cobraPath - the location of the CobraToolbox
    % @ zero_cutoff- the value which is will be considered as zero for model extractions (epsilon in the original paper)
    % @ objective - boolean, whether or not include the objective function(s) in the reaction list of active reactions

       
    cobra_path = char(cobra_path);
    % load the CobraToolbox
    run(fullfile(cobra_path, 'initCobraToolbox(false)'));
    disp('initialized CobraToolbox');

    % initialize the input and output variables
    if ~(exist('model_in','var'))
        model_in = "/home/taube/Work/miTarget/eMed_FUTURE/Pipeline_MM/results/data/CBT_models/colormore22.mat";
    end
    if ~(exist('model_out','var'))
        model_out = 'temp/extractedRxn.csv';
    end
    if ~(exist('objective','var'))
        objective = true;
    end
    if ~(exist('core_rxns','var'))
        core_rxns = 'results/data/coreRxns/coreRxns.GL25-L50-GU75.colormore22.csv';
    end
    if ~(exist('sample','var'))
        sample = 'E02875_L1_S2_L001';
    end
    if ~(exist('zero_cutoff','var'))
        zero_cutoff = getCobraSolverParams('LP','optTol');
    end

    % convert the input into character in case its string, and display the
    % values for the log file
    model_in = char(model_in)
    model_out = char(model_out)
    core_rxns = char(core_rxns)
    sample = char(sample)

    % load the model
    mod = readCbModel(model_in);
    
    % adjust the diet
    if exist('diet','var')
        diet = char(diet)
        % load the diet
        diet = readtable(diet);

        % convert the values to double
        if(~all(class(diet{1:end,2} == 'double')))
            diet.val = str2double(diet{1:end,2});
        else
            diet.val = diet{1:end,2};
        end
        diet.val(isnan(diet.val)) = 0;

        % get the reaction IDs of the model
        diet.rxnID = findRxnIDs(mod,diet{:,1});
        diet = diet(~(diet.rxnID == 0),:);


        % reconstrain model
        ex_rxn = findExcRxns(mod);
        EX_rxn = startsWith(mod.rxns,'EX');
        ex = ex_rxn & EX_rxn;
        % remove possible fluxes on the exchange reactions
        mod = changeRxnBounds(mod, mod.rxns(ex), 0, 'l');
        mod = changeRxnBounds(mod, mod.rxns(diet.rxnID),-1*diet.val,'l');
    end
    modDiet = mod;

    if (exist('consistent_rxns','var'))
        consistent_rxns = char(consistent_rxns)
        % remove inconsistent rxns
        consrxn = readtable(consistent_rxns);
        rxnRm = mod.rxns(~ismember(mod.rxns, consrxn.rxn));
        mod = removeRxns(mod, rxnRm);
    end

    % get the core reactions
    corerxns = readtable(core_rxns,'ReadRowNames',true);
    core_all = corerxns.Row(find(corerxns{:,sample}));
    core = findRxnIDs(mod, core_all);
    core = core(find(~(core == 0)));

    % add the objective functions to the core reactions if needed
    if objective
        obj_rxn = find(mod.c);
        obj_rxn = obj_rxn(~ismember(obj_rxn, core));
        core(length(core) + (1:length(obj_rxn))) = obj_rxn;
    end

    % extract the model with fastcore and write it an sbml-file
    opt= struct('core', core,'solver','fastCore','epsilon',zero_cutoff);
    tic
    try;
        mod_tiss = createTissueSpecificModel(mod, opt);
    catch ME;
        warning('Model extraction failed, trying to generate new consistent model using fastcc');
        consID = fastcc(modDiet);
        rxnRM = mod.rxns(~ismember(mod.rxns, mod.rxns(consID)));
        mod = removeRxns(mod, rxnRM);
        core = findRxnIDs(mod, core_all);
        core = core(find(~(core == 0)));
        opt = struct('core', core,'solver','fastCore','epsilon',zero_cutoff);
        mod_tiss = createTissueSpecificModel(mod, opt);
    end;
    toc

    % write the reactions to a csv file
    consRxn = struct;
    consRxn.rxn = mod_tiss.rxns;
    consRxn.id = findRxnIDs(mod, consRxn.rxn);
    writetable(struct2table(consRxn), model_out);
