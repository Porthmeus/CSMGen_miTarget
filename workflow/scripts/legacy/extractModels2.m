function [mod_test] = extractModels(model, TPM_matrix, conversionTable, opts, cobraPath)
    % read the model, extract the genes from the tpm matrix and calculate a
    % transcriptpome specific model per column in the TPM_matrix. Write these
    % models to an SBML file
       
    try
        cobraPath = char(cobraPath);
        % load the CobraToolbox
        cwdir = pwd;
        cd(cobraPath);
        initCobraToolbox(false);
        cd(cwdir);
        disp('initialized CobraToolbox');


        if ~(exist('model','var'))
            model = "~/Work/miTarget/FUTURE/PipelineAnalysis/resources/Recon2.2_edit.sbml";
        end

        if ~(exist('TPM_matrix','var'))
            TPM_matrix = "/home/taube/Work/miTarget/FUTURE/PipelineAnalysis/results/data/SPLIT_TPM_emed_future/SPLIT_F02343_L1_S37_L005.csv"
        end

        if ~(exist('conversionTable','var'))
            conversionTable ='~/Work/miTarget/FUTURE/PipelineAnalysis/resources/ensemble2HGNC.csv';
        end
        
        if ~(exist('opts','var'))
            opts = struct('thrsld', 75 ,... % defines the percentile larger than what is considered a coreset to parse to fastcore
                'outFile', 'extractedModel.xml');
        elseif ~(isfield(opts, 'thrsld'))
            opts.thrsld = 75;
        elseif ~(isfield(opts, 'outFile'))
            opts.outFile ='extractedModel.xml' ;
        end
        
        % convert the input into character in case its string, and display the
        % values for the log file
        model = char(model)
        TPM_matrix = char(TPM_matrix)
        conversionTable = char(conversionTable)
        disp(opts) 
 
        % load the model
        mod_recon2_consistent = readCbModel(model);
        
        % it might happen if the model was created with cobraPy, that gene
        % naming is screwed after loading the model with matlab. Thus I will
        % fix, this issue here now, but in a very specific way for the current
        % pipeline... I hope this will not end up in to much confusion later
        % on, but otherwise I would have to change the in the source of cobraPy
        % or cobraToolbox, which are not good options either... sorry for the
        % ones reading this after encountering problems.
        mod_recon2_consistent.geneNames = replace(mod_recon2_consistent.geneNames, ...
            "G_HGNC__58__", "HGNC:");

        % load the conversion table to match 
        tbl_ensemble2HGNC = readtable(conversionTable);
        opt = detectImportOptions(TPM_matrix); % I don't know why, but is required to load the table correctly
        tbl_FutureReadCounts = readtable(TPM_matrix, opt);

        % extract relevant information from the transcriptome data
        tbl_ensemble2HGNC_modOnly = tbl_ensemble2HGNC(...
            ismember(tbl_ensemble2HGNC.HGNC_ID, ...
            mod_recon2_consistent.geneNames),:);
        [dummy, ia, dummy] = unique(tbl_ensemble2HGNC_modOnly.HGNC_ID);
        tbl_ensemble2HGNC_modOnly = tbl_ensemble2HGNC_modOnly(ia,:);
        tbl_FutureReadCounts_modOnly = tbl_FutureReadCounts( ismember(tbl_FutureReadCounts{:,1}, ...
            tbl_ensemble2HGNC_modOnly.Gene_stable_ID),:);
        
        % set the model gene names to the rownames, delete the original identifier
        tbl_FutureReadCounts_MONewNames = join(tbl_FutureReadCounts_modOnly, ...
            tbl_ensemble2HGNC_modOnly(:,[1,5]), ...
            'LeftKey', 1, ...
            'RightKey', 'Gene_stable_ID');
        tbl_FutureReadCounts_MONewNames.Row = tbl_FutureReadCounts_MONewNames.HGNC_ID;
        tbl_FutureReadCounts_MONewNames(:,1) = [];
        tbl_FutureReadCounts_MONewNames.HGNC_ID = [];

     
        % transform the data to reactions of the model
        data = struct;
        data.gene = tbl_FutureReadCounts_MONewNames.Row;
        data.value = tbl_FutureReadCounts_MONewNames{:,1};
        
        disp("Number of core genes")
        disp(data)

        reactionExpression = mapExpressionToReactions(mod_recon2_consistent, data);
        reaction_IDs = find(reactionExpression >= prctile(reactionExpression, str2num(opts.thrsld)));
        disp(length(reaction_IDs))
        
        % another annoying patch, but for some reason fastcore does not run, if
        % the reaction 'r0364' is in the core set. I will simply remove it from
        % the core set if it is in it.
        id2remove = find(ismember(mod_recon2_consistent.rxns(reaction_IDs), ['r0364']));
        reaction_IDs(id2remove) = [];

        % extract the model with fastcore and write it an sbml-file
        opt= struct('core', reaction_IDs,'solver','fastCore');
        mod_test = createTissueSpecificModel(mod_recon2_consistent, opt);
        
        writeCbModel(mod_test, 'format', 'sbml','fileName', opts.outFile);
        

    catch ME
        disp(getReport(ME));
        exit;
    end
        
