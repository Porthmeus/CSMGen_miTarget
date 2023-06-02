function [mod_tiss] = extractModels(model_in, model_out, core_Rxns, diet, consistent_Rxns, cobraPath, objective)
    % read the model, extract the genes from the tpm matrix and calculate a
    % transcriptpome specific model per column in the TPM_matrix. Write these
    % models to an SBML file

    % @ model_in - a file which can be read to a CBT model by readCbModel()
    % @ model_out - the file where the results should be stored as a csv file for the list of reactions in the model, which are tissue specific
    % @ core_Rxns - a table of core reactions, first column is the ID of the reaction, the second column is indicating whether the rxn is part of the core (=1) or not (=0)
    % @ diet - a csv file with constrains for the diet - the first column should contain the rxns names, the second column values for the fluxes of these reaction. The function will set the lower bounds to 0 for all exchange reaction not in that list. For the reactions in the list the lower bound will be set to -1*value. If such a file is not given, it is assumed, that the model is correctly constrained.
    % @ consistent_Rxns - a table of the reaction which make up the consistent model. First column should contain the reaction-IDs of the consistent model. If no file is given the model is assumed to be consistent.
    % @ cobraPath - the location of the CobraToolbox
    % @ biomass - boolean, whether or not include the objective function(s) in the reaction list of active reactions

       
    try
        cobraPath = char(cobraPath);
        % load the CobraToolbox
        run(fullfile(cobraPath, 'initCobraToolbox(false)'));
        disp('initialized CobraToolbox');


        if ~(exist('model_in','var'))
            model_in = "/home/taube/Work/miTarget/eMed_FUTURE/Pipeline_MM/results/data/CBT_models/colormore22.mat";
        end
        if ~(exist('model_out','var'))
            model_out = 'temp/extractedRxn.csv';
        end
        if ~(exist('core_Rxns','var'))
            core_Rxns = "/home/taube/Work/miTarget/eMed_FUTURE/Pipeline_MM/temp/SPLITGL25|L50_colormore22-F02234_L1_S1_L001.csv";
        end
        if ~(exist('objective','var'))
            objective = false;
        end
        
        
        % convert the input into character in case its string, and display the
        % values for the log file
        model_in = char(model_in)
        model_out = char(model_out)
        dietF = char(dietF)
        consistent_Rxns = char(consistent_Rxns)
        core_Rxns = char(core_Rxns)
        cobraPath = char(cobraPath)
        disp(objective)

        % debug
        %consistent_Rxns = 'temp/consistent_rxns.csv';%'results/data/consistentModels/MatjesAbsorption.colormore22_Consistent.csv';
        %dietF = 'resources/diets/MatjesAbsorption.csv';

 
        % load the model
        mod = readCbModel(model_in);
        
        if exist('diet','var')
            % load the diet
            diet = readtable(dietF);
            if(~all(class(diet{1:end,2} == 'double')))
                diet.val = str2double(diet{1:end,2});
            else
                diet.val = diet{1:end,2};
            end
            diet.rxnID = findRxnIDs(mod,diet{:,1});
            diet = diet(~(diet.rxnID == 0),:);
            diet.val(isnan(diet.val)) = 0;


            % reconstrain model
            ex_rxn = findExcRxns(mod);
            EX_rxn = startsWith(mod.rxns,'EX');
            ex = ex_rxn & EX_rxn;
            % remove possible fluxes on the exchange reactions
            mod = changeRxnBounds(mod, mod.rxns(ex), 0, 'l');
            mod = changeRxnBounds(mod, mod.rxns(diet.rxnID),-1*diet.val,'l');
        end

        if (exist('consistent_Rxns','var'))
            % remove inconsistent rxns
            consRxn = readtable(consistent_Rxns);
            rxnRm = mod.rxns(~ismember(mod.rxns, consRxn.rxn));
            mod = removeRxns(mod, rxnRm);
        end
        
       % %%%%%%%%%%% DEPRECATED %%%%%%%%%%%%%%%%%%%
       % % it might happen if the model was created with cobraPy, that gene
       % % naming is screwed after loading the model with matlab. Thus I will
       % % fix, this issue here now, but in a very specific way for the current
       % % pipeline... I hope this will not end up in to much confusion later
       % % on, but otherwise I would have to change the in the source of cobraPy
       % % or cobraToolbox, which are not good options either... sorry for the
       % % ones reading this after encountering problems.
       % mod.geneNames = replace(mod.geneNames, ...
       %     "G_HGNC__58__", "HGNC:");

       % % load the conversion table to match 
       % tbl_ensemble2HGNC = readtable(conversionTable);
       % opt = detectImportOptions(TPM_matrix); % I don't know why, but is required to load the table correctly
       % tbl_FutureReadCounts = readtable(TPM_matrix, opt);

       % % extract relevant information from the transcriptome data
       % tbl_ensemble2HGNC_modOnly = tbl_ensemble2HGNC(...
       %     ismember(tbl_ensemble2HGNC.HGNC_ID, ...
       %     mod.geneNames),:);
       % [dummy, ia, dummy] = unique(tbl_ensemble2HGNC_modOnly.HGNC_ID);
       % tbl_ensemble2HGNC_modOnly = tbl_ensemble2HGNC_modOnly(ia,:);
       % tbl_FutureReadCounts_modOnly = tbl_FutureReadCounts( ismember(tbl_FutureReadCounts{:,1}, ...
       %     tbl_ensemble2HGNC_modOnly.Gene_stable_ID),:);
       % 
       % % set the model gene names to the rownames, delete the original identifier
       % tbl_FutureReadCounts_MONewNames = join(tbl_FutureReadCounts_modOnly, ...
       %     tbl_ensemble2HGNC_modOnly(:,[1,5]), ...
       %     'LeftKey', 1, ...
       %     'RightKey', 'Gene_stable_ID');
       % tbl_FutureReadCounts_MONewNames.Row = tbl_FutureReadCounts_MONewNames.HGNC_ID;
       % tbl_FutureReadCounts_MONewNames(:,1) = [];
       % tbl_FutureReadCounts_MONewNames.HGNC_ID = [];

     
       % % transform the data to reactions of the model
       % data = struct;
       % data.gene = tbl_FutureReadCounts_MONewNames.Row;
       % data.value = tbl_FutureReadCounts_MONewNames{:,1};
        
       % disp("Number of core genes")
       % disp(data)

       %% reactionExpression = mapExpressionToReactions(mod, data);
       %% reaction_IDs = find(reactionExpression >= prctile(reactionExpression, str2double(thrsld)));
       %% another annoying patch, but for some reason fastcore does not run, if
       %% the reaction 'r0364' is in the core set. I will simply remove it from
       %% the core set if it is in it.
       % id2remove = find(ismember(mod.rxns(reaction_IDs), ['r0364']));
       % reaction_IDs(id2remove) = [];

       % %%%%%%%%%%% DEPRECATED END %%%%%%%%%%%%%%%%%%%

        % check for infinity boundaries and correct them
        sel = mod.ub == Inf;
        mod.ub(sel) = 1000;
        sel = mod.lb == Inf;
        mod.lb(sel) = 1000;
        sel = mod.ub == -Inf;
        mod.ub(sel) = -1000;
        sel = mod.lb == -Inf;
        mod.lb(sel) = -1000;

        % just to be sure once more check for inconsistencies
        rm_rxn = 'dummy';
        while length(rm_rxn) >0
            [A,flib,V] = fastcc(mod,1E-4,1);
            cons_rxn = mod.rxns(A);
            rm_rxn = mod.rxns(find(~ismember(mod.rxns,cons_rxn)));
            mod = removeRxns(mod, rm_rxn);
            if length(rm_rxn) == 0
                break
            end
        end


        % load the core rxns
        coreRxns = readtable(core_Rxns);
        reaction_names = coreRxns{find(coreRxns{1:end,2}),1};
        reaction_IDs = findRxnIDs(mod, reaction_names);     
        % this most likely contains reactions which are not present anymore, because of inconsistencies of the model - remove those reactions
        reaction_IDs = reaction_IDs(find(~(reaction_IDs == 0)));


        if objective
            obj_rxn = find(mod.c);
            obj_rxn = obj_rxn(~ismember(obj_rxn, reaction_IDs));
            reaction_IDs(length(reaction_IDs) + (1:length(obj_rxn))) = obj_rxn;
        end


        % extract the model with fastcore and write it an sbml-file
        opt= struct('core', reaction_IDs,'solver','fastCore','epsilon',1E-4);
        tic
        mod_tiss = createTissueSpecificModel(mod, opt);
        toc
        % write the reactions to a csv file

        consRxn = struct;
        consRxn.rxn = mod_tiss.rxns;
        consRxn.id = findRxnIDs(mod, consRxn.rxn);
        writetable(struct2table(consRxn), model_out);

        %writeCbModel(mod_tiss, 'format', 'sbml','fileName', model_out);
        

    catch ME
        disp(getReport(ME));
        %exit;
    end
        
