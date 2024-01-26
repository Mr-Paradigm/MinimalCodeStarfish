function status=odetpbar(t,~,flag)
    % persistent tf tstart;
    
    if isempty(flag)
        % Integration steps
        ts=mean(t);
        % progress=100*ts/tf;
        status=0;
        
        if mod(ts,5) < 0.1
            disp(ts);
        end
    else
        switch flag
            case 'init'     % Initializing progress bar
                %tstart=tic;
                % tf=max(t);
                disp('Simulation is running!');
            case 'done'     % Finishing status function
                % tf=[];
                % display([ '   Integration time: ' num2str(toc(tstart))]);
                disp('Finished!');
                % tstart=[];
            otherwise
                error('odetpbar:UnknownError',...
                    'Unknown error has occured');
        end
    end
end