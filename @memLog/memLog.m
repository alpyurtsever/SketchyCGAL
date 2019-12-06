%% memLog - by Alp Yurtsever (alp.yurtsever@epfl.ch OR alpy@mit.edu)
% This is a simple class for monitoring the memory usage of MATLAB
% externally in UNIX systems. This simple prototype first gets the PID of
% MATLAB, and then runs a simple bash script on the system level that reads
% the storage usage every second and logs the maximum usage on a hidden tmp
% file. Then, we can read this maximum value at any time using the prompt
% method of the class. Example usage is as follows:
%
%   hmL = memLog('A UNIQUE NAME');  % Create a memLog handle
%   hmL.start();                    % Start memory logging
%   mem1 = hmL.prompt();            % Get the maximum memory usage thus far
%   A = randn(10000,5000);          
%   mem2 = hmL.prompt();            % Get the maximum memory usage thus far
%   B = randn(10000,5000);          
%   mem3 = hmL.prompt();            % Get the maximum memory usage thus far
%   hmL.stop();                     % or (clear hmL): stop logging and delete tmp file
%
% Notes: 
% -> hmL.prompt() pauses 10-20 seconds because it waits for a steady state 
% of memory usage 
% -> Note that the internal process of MATLAB is far more complicated than 
% line by line processing, and this method only gives a rough idea about
% the actual memory usage. It is NOT DESIGNED for accurate measurements.
% The results you will obtain will depend on the system specifications as
% well as the system load, hence will not be reproducible. 
classdef memLog < handle
    
    properties
        ID
        PID
        PIDstr
    end
    
    methods
        
        function obj = memLog(varargin)
            narginchk(0,1);
            if ~isunix, warning('memLog works only in UNIX systems for now!'); end
            if ~isempty(varargin), IDappend = regexprep(varargin{1},'/',''); else, IDappend = ''; end
            obj.ID = [IDappend,datestr(now,30),num2str(randi([10000,99999]))];
            obj.PID = num2str(feature('getpid'));
            obj.PIDstr = num2str(obj.PID);
        end
        
        function start(obj)
            
            if ~isunix, warning('memLog works only in UNIX systems for now!'); return; end
            
            if ~exist('./.memLogTmp','dir'), mkdir('./.memLogTmp'); end
            
            unix(...
                ['logpid() { MAX=0; while sleep 1; do ',...
                'VAR=$(ps -p $1 -o vsz= ); ',...
                'if [ $VAR -gt $MAX ]; then  MAX=$VAR; echo "$MAX"; fi; done; }; ',...
                ...
                'logpid "',obj.PID,'" >> ./.memLogTmp/',obj.ID,'.log.tmp &']);
            
            pause(2);
            
            tic;
            while ~exist(['./.memLogTmp/',obj.ID,'.log.tmp'],'file')
                if toc > 20
                    warning('memLog could not start!');
                    return;
                end
                pause(1);
            end
        end
        
        function stop(obj)
            
            if ~isunix, warning('memLog works only in UNIX systems for now!'); return; end
            
            if exist(['./.memLogTmp/',obj.ID,'.log.tmp'],'file')
                delete(['./.memLogTmp/',obj.ID,'.log.tmp']);
            else
                warning('memLog log file is not found. object might be already stopped!');
            end
        end
        
        function [MAX] = prompt(obj)
            
            if ~isunix, warning('memLog works only in UNIX systems for now!'); MAX = nan; return; end

            if exist(['./.memLogTmp/',obj.ID,'.log.tmp'],'file')
                MAX = -90;
                MAX2 = -89;
                while MAX ~= MAX2
                    [~,VALSTR] = system(['tail -n 1 ','./.memLogTmp/',obj.ID,'.log.tmp']);
                    MAX = str2double(VALSTR);
                    pause(10)
                    [~,VALSTR] = system(['tail -n 1 ','./.memLogTmp/',obj.ID,'.log.tmp']);
                    MAX2 = str2double(VALSTR);
                end
            else
                warning('memLog log file is not found. object might be stopped!');
                MAX = nan;
            end
        end
        
        function delete(obj)
            if exist(['./.memLogTmp/',obj.ID,'.log.tmp'],'file')
                delete(['./.memLogTmp/',obj.ID,'.log.tmp']);
            end
        end
        
    end
    
end
%% Last edit: Alp Yurtsever - November 29, 2019
