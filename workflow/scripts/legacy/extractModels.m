function [tbl_SamplePath] = extractModels(model, TPM_matrix, conversionTable, opts)
    % read the model, extract the genes from the tpm matrix and calculate a
    % transcriptpome specific model per column in the TPM_matrix. Write these
    % models to an SBML file
    try
        if ~(exist('model','var'))
            model = "~/Work/miTarget/FUTURE/PipelineAnalysis/resources/Recon2.2_edit.sbml";
        end

        if ~(exist('TPM_matrix','var'))
            TPM_matrix = "/home/taube/Work/miTarget/FUTURE/PipelineAnalysis/results/data/TPM_emed_future.csv";
        end

        if ~(exist('conversionTable','var'))
            conversionTable ='~/Work/miTarget/FUTURE/PipelineAnalysis/resources/ensemble2HGNC.csv';
        end
        
        if ~(exist('opts','var'))
            opts = struct('thrsld', 75 ,... % defines the percentile larger than what is considered a coreset to parse to fastcore
                'dirPath', 'extractedModels');
        elseif ~(isfield(opts, 'thrsld'))
            opts.thrsld = 75;
        elseif ~(isfield(opts, 'dirPath'))
            opts.dirPath = 'dirPath';
        end

        model = char(model)
        TPM_matrix = char(TPM_matrix)
        conversionTable = char(conversionTable)
        disp(opts) 
        % load the model
        mod_recon2 = readSBML(model,1000);

        % remove the inconsistent reactions
        [idx_consistentReactions, mod_recon2_flipped, dummy] = fastcc(mod_recon2, 1e-4, 0);
        char_inconsistenReactions = setdiff(mod_recon2.rxns, mod_recon2.rxns(idx_consistentReactions));
        mod_recon2_consistent = removeRxns(mod_recon2, char_inconsistenReactions);

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

     
        % create a table which contains path to the extracted model
        Path = strcat(opts.dirPath, ...
            filesep, ...
            tbl_FutureReadCounts_MONewNames.Properties.VariableNames, ...
            '.sbml');
        Sample = tbl_FutureReadCounts_MONewNames.Properties.VariableNames;
        tbl_SamplePath = cell2table(horzcat(Sample', Path'), 'VariableNames', {'Sample','Path'});
        tbl_MetaData = table(Sample, Path);
        

        for (i = 1:height(tbl_SamplePath))
            disp(tbl_SamplePath(i,:))
            % transform the data to reactions of the model
            data = struct;
            data.gene = tbl_FutureReadCounts_MONewNames.Row;
            data.value = tbl_FutureReadCounts_MONewNames{:,i};
            disp(data)
            reactionExpression = mapExpressionToReactions(mod_recon2_consistent, data);
            reaction_IDs = find(reactionExpression >= prctile(reactionExpression, opts.thrsld));
            disp(length(reaction_IDs))

            % extract the model with fastcore and write it an sbml-file
            opt= struct('core', reaction_IDs,'solver','fastCore');
            mod_test = createTissueSpecificModel(mod_recon2_consistent, opt);
            
            writeCbModel(mod_test, 'format', 'sbml','fileName',char(tbl_SamplePath{i,'Path'}));
        end
        
        % write the table containing the sample names and the file path into the results directory
        csvPath = strcat(opts.dirPath, filesep, "SamplePath.csv");
        writetable(tbl_SamplePath, csvPath);

    catch ME
        disp(getReport(ME));
        exit;
    end
        
