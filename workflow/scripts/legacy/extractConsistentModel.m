function [consistentModel] = extractConsistentModel(model_in, model_out, cobraPath)
        % read an inconsistent model from an sbml file, find incosistent
        % reactions using fastcc, remove the inconsistent reactions and write
        % the consistent model

        
        try
            disp(model_in)
            disp(model_out)
            disp(cobraPath)
            % load the CobraToolbox
            cwdir = pwd;
            cd(cobraPath);
            initCobraToolbox(false);
            cd(cwdir);
            disp('initialized CobraToolbox');
            % load the model
            
            mod = readCbModel(model_in);
            disp('loaded model')

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

            % write the model to disk
            writeCbModel(consistentModel, 'format', 'sbml', 'fileName',model_out)
            
            disp('wrote model')
        catch ME
            disp(ME)
            exit
        end

