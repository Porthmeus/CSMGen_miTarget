function [consistentModel] = extractConsistentModel(model_in, diet_in, model_out, cobraPath)
        % read an inconsistent model from an sbml file, find incosistent
        % reactions using fastcc, remove the inconsistent reactions and write
        % the consistent model

        
        try
            disp(model_in)
            disp(diet_in)
            disp(model_out)
            disp(cobraPath)
            % load the CobraToolbox
            run(fullfile(cobraPath, 'initCobraToolbox(false)'));
            disp('initialized CobraToolbox');
            
            % load the model
            mod = readCbModel(model_in);
            disp('loaded model')

            % reconstrain diet
            if exist('diet_in', 'var')
                % load the diet
                diet = readtable(diet_in);
                diet.val = str2double(diet{1:end,2});
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
                
            % find inconsistent reactions and remove them
            consistentModel = mod;

            char_inconsistenReactions = 1;
            while length(char_inconsistenReactions) > 0
                [idx_consistentReactions, mod_flipped, dummy] = fastcc(consistentModel, 1e-4, 0);
                char_inconsistenReactions = setdiff(consistentModel.rxns, consistentModel.rxns(idx_consistentReactions));
                catstring = strcat("Found ", num2str(length(char_inconsistenReactions)), " inconsistent reactions");
                disp(catstring)
                if length(char_inconsistenReactions) == 0
                    break
                end
                consistentModel = removeRxns(consistentModel, char_inconsistenReactions);
            end
            disp('extracted model')
            toc


            % write the consistent reactions to disk
            consRxn.rxn = consistentModel.rxns
            consRxn.id = findRxnIDs(mod, consRxn.rxn);
            writetable(struct2table(consRxn), model_out);
            %writeCbModel(consistentModel, 'format', 'sbml', 'fileName',model_out)
            
            disp('wrote consistent rxns')
        catch ME
            disp(ME)
            exit
        end

