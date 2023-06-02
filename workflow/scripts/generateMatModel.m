function [consistentModel] = generateMatModel(model_in, model_out, cobraPath)
        % read an inconsistent model from an sbml file, find incosistent
        % reactions using fastcc, remove the inconsistent reactions and write
        % the consistent model

        
        try
            disp(model_in)
            disp(model_out)
            disp(cobraPath)
            % load the CobraToolbox
            run(fullfile(cobraPath, 'initCobraToolbox(false)'));
            disp('initialized CobraToolbox');
            
            % load the model
            mod = readCbModel(model_in);
            disp('loaded model')

            writeCbModel(mod, model_out)
            
            disp('Wrote .mat model')
        catch ME
            disp(ME)
            exit
        end

