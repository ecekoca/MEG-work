%% Maxfilter
% Rik Henson, 2017; Modifications by Ece K 2017

% Script detects and interpolates bad channels, performs signal space
% separation and shifts data to a default space.

addpath /imaging/local/meg_misc/
addpath /neuro/meg_pd_1.2/
addpath /imaging/projects/cbu/ntad/scripts/functions


clear; clc;

% Directories
raw_dir = '/imaging/projects/cbu/ntad/meg_data/'; % Where MEG data lives
ana_dir = '/imaging/projects/cbu/ntad/meg_data/'; % Where you want to put the pre-processed MEG data

subjects = {'P1061','P1063','P1067','P1020','P1045','P3006','P3004','P1042','P1027','P1053'};

all_tasks={
    'avtask'
    'mmn'
    'scenerep'
    'vab'
    'rest_open'
    'rest_closed'
    }; % List of all tasks

do_tasks={
      'avtask'
      'mmn'
      'scenerep'
      'vab'
      'rest_open'
      'rest_closed'
    }; % List of task data you want to maxfilter

session='AF'; %'AF', 'BL' One session at a time.

% Processing options
run_autobad = 1;
run_tsss = 1;
run_trans_default = 1;
OverWrite = 1;

%%=========================================================================
%% Check that file exists and it's named correctly (EK)
%%=========================================================================

if ~exist(ana_dir); mkdir(ana_dir); end

disp('Are these filenames correct? Any typos?') 
for sub=1:length(subjects)
    for ses = 1:length(do_tasks)
        fold=[raw_dir subjects{sub} '/' session '/' do_tasks{ses} '/']; 
        if exist(fold)
        cd(fold)
        disp(subjects{sub})
        ls([do_tasks{ses}(1:2) '*'])
        else
            disp(subjects{sub})
            disp('File does not exist')
        end
            
    end
end

%%=========================================================================
%% Check if maxfiltering has been done already (EK)
%%=========================================================================

if ~OverWrite
    k=0; clear newlist
    for sub=1:length(subjects)
        for ses = 1:length(do_tasks)
            ses_num = find(strcmp(all_tasks,do_tasks{ses}));
            ses_dir = fullfile(ana_dir,subjects{sub},session,do_tasks{ses});
            if exist(ses_dir,'dir')
                cd(ses_dir);
                if ~exist([ses_dir '/transdef.fif'],'file')
                    disp([subjects{sub} ': no Maxfilter']);k=k+1;
                    newlist{k}=subjects{sub};
                end
            else
                disp([subjects{sub} ': no task data']);%k=k+1;
                %newlist{k}=subjects{sub};
            end
        end
    end; disp([num2str(k) '/' num2str(length(subjects)) ' of your subjects are missing Maxfiltering'])
    subjects=newlist;
end

clear newlist k


%%=========================================================================
%% Setting up Matlab Pool
%%=========================================================================
 if length(subjects)>1 %EK
     parpool close force CBU_Cluster
     if parpool('CBU_Cluster')==0
         %MaxNsubs = max([16 length(subjects)]);
         P = cbupool(15);
         P.ResourceTemplate = '-l nodes=^N^,mem=4GB,walltime=140:00:00';
         P.SubmitArguments='-W x="NODESET:ONEOF:FEATURES:MAXFILTER"';
         parpool(P);
     end
 end

%%=========================================================================
%% Begin
%%=========================================================================

basestr = ' -ctc /neuro/databases/ctc/ct_sparse.fif -cal /neuro/databases/sss/sss_cal.dat'; % set autobad settings
basestr = [basestr ' -linefreq 50 -hpisubt amp'];
basestr = [basestr ' -force'];
maxfstr = '!/neuro/bin/util/x86_64-pc-linux-gnu/maxfilter-2.2.12 ';
tSSSstr = ' -st 10 -corr 0.98'; %'tSSSstr = '';
dsstr = ''; % dsstr = ' -ds 4';  % downsample to 250Hz - can cause MF crash?
incEEG = 0;
TransTarget = '';

for sub = 1:length(subjects)
    %parfor sub = 1:length(subjects)
    sub_dir = fullfile(ana_dir,subjects{sub},session);
    
    movfile = fullfile(sub_dir,'trans_move.txt'); % This file will record translations between runs
    rik_eval(sprintf('!touch %s',movfile));
    rik_eval(sprintf('!echo ''%s'' >> %s',datestr(now),movfile));
    
    for ses = 1:length(do_tasks)
        if exist(fullfile(ana_dir,subjects{sub},session,do_tasks{ses}))
        ses_num = find(strcmp(all_tasks,do_tasks{ses}));
        ses_dir = fullfile(ana_dir,subjects{sub},session,do_tasks{ses});
        %date_dir = dir(fullfile(raw_dir,subjects{sub},'*')); date_dir = date_dir(3).name;
        %rawsub_dir = fullfile(raw_dir,subjects{sub},date_dir);
        cd(ses_dir)
        if strcmp(do_tasks{ses},'vab')
            raw_file=dir('*ta*raw*'); % in case there's a typo in the filename EK
        else
            raw_file = dir('*raw*');
        end
        raw_file=raw_file.name;
        
        raw_stem=raw_file(1:end-4);
        if ~exist(raw_file,'file')
            error(raw_file,' does not exist')
        end
        
        % Fit sphere (since better than MaxFilter does)
        [co ki] = hpipoints(raw_file);
        cd(ses_dir)
        tmppoints = co(:,ki>1)';  % don't include the fiducial points
        nloc = find(~(tmppoints(:,3) < 0 & tmppoints(:,2) > 0));  % Remove the nose points if you have any
        headpoints = tmppoints(nloc,:);
        save([ses_dir '/hpi.txt'],'-ASCII','headpoints') % save headpoints
        cmd_fit = ['/imaging/local/software/neuromag/bin/util/fit_sphere_to_points hpi.txt'];
        
        if exist([ses_dir '/fittmp.txt'],'file') %EK
            s=dir('fittmp.txt');
            while s.bytes==0
                try delete('fittmp.txt'), end
                eval(['! ' cmd_fit ' > fittmp.txt']);
                s=dir('fittmp.txt');
            end
        else
            s=[];
            s.bytes=0;
            while s.bytes==0
                try delete('fittmp.txt'), end
                eval(['! ' cmd_fit ' > fittmp.txt']);
                s=dir('fittmp.txt');
            end
        end
        
        spherefit = dlmread('fittmp.txt');
        %orig = str2num(spherefit)*1000;  % m to mm;
        
        if length(spherefit)<5;
            error('Spherefit failed for a different reason!')
        end
        
        orig=spherefit*1000;
        origstr = sprintf(' -origin %d %d %d -frame head',orig(1),orig(2),orig(3))
        
        %%=================================================================
        %% 1. Autobad
        %%=================================================================
        % (this email says important if doing tSSS later 
        % https://www.jiscmail.ac.uk/cgi-bin/webadmin?A2=NEUROMEG;d3f363f3.1205)
        
        if run_autobad
            outfile = fullfile(ses_dir,'bad'); badfile = fullfile(ses_dir,'bad.txt'); logfile = fullfile(ses_dir,'bad.log');
            badstr  = sprintf(' -autobad %d -badlimit %d',1800,7); % 1800s is 30mins - ie enough for all do_sessions EK
            
            if ~exist(logfile,'file') | OverWrite
                filestr = sprintf(' -f %s -o %s.fif',[sub_dir '/' do_tasks{ses} '/' raw_file],outfile);
                posfile = fullfile(ses_dir,'headpos.txt');
                compstr = sprintf(' -headpos -hpistep 10 -hp %s',posfile);
                finstr = [maxfstr filestr origstr basestr badstr compstr sprintf(' -v | tee %s.log',outfile)]
                rik_eval(finstr);
                delete(sprintf('%s.fif',outfile));
            end
            
            % Pull out bad channels from logfile:
            delete(badfile); rik_eval(sprintf('!echo '' '' > %s', badfile));
            rik_eval(sprintf('!cat %s.log | sed -n -e ''/Detected/p'' -e ''/Static/p'' | cut -f 5- -d '' '' >> %s',outfile,badfile));
            tmp=dlmread(badfile,' '); % Nbuf = size(tmp,1)-1; Doesn't work if subset of channels bad for short time
            tmp=reshape(tmp,1,prod(size(tmp)));
            tmp=tmp(tmp>0); % Omit zeros (padded by dlmread):
            
            % Mark bad based on threshold (currently ~5% of buffers (assuming 500 buffers)):
            [frq,allbad] = hist(tmp,unique(tmp));
            badchans = allbad(frq>0.05*max(frq)); %save badchans.mat allbad badchans; %EK assuming that there is always at least one sensor that's bad throughout eg 813
            if isempty(badchans) badstr = '';
            else
                badstr = sprintf(' -bad %s',num2str(badchans))
            end
        end
        
        %%=================================================================
        %% 2. tSSS
        %%=================================================================
        % (ie, align within subject if multiple sessions)
       
        if run_tsss
            if ~exist(sprintf('%s.fif',outfile),'file') | ~exist([outfile '.log'],'file') | OverWrite
                % waitbar(0.66, h, ['Running tSSS and movecomp on Sub:' subjects{sub}(7:10) ' Task:' do_sessions{ses} '...']); %EK
                transfstfile = [outfile '.fif'];
                if ~isempty(TransTarget) & ~strcmp(do_tasks{ses},TransTarget)
                    outfile = fullfile(ses_dir,sprintf('tsss_trans_%s',TransTarget))
                    transtr = sprintf(' -trans %s/tsss.fif',fullfile(sub_dir,TransTarget));
                else
                    transtr = '';
                    outfile = fullfile(ses_dir,sprintf('tsss'))
                end
                
                compstr = sprintf(' -movecomp inter'); %compstr = sprintf(' -movecomp');
                filestr = sprintf(' -f %s -o %s.fif',[ses_dir '/' raw_file],outfile);
                finstr = [maxfstr filestr basestr badstr tSSSstr compstr origstr transtr dsstr sprintf(' -v | tee %s.log',outfile)]
                rik_eval(finstr);
            end
            
            
            if ~isempty(TransTarget) & ~strcmp(do_tasks{ses},TransTarget)
                %fprintf(fp,'\nTransfirst %s to %s: ',do_sessions{ses},TransTarget);
                rik_eval(sprintf('!echo ''Trans %s to %s'' >> %s',do_tasks{ses},TransTarget,movfile));
                rik_eval(sprintf('!cat %s.log | sed -n ''/Position change/p'' | cut -f 7- -d '' '' >> %s',outfile,movfile));
            end
        end
        
        %%=================================================================
        %% 3. Trans default
        %%=================================================================
        
        if run_trans_default
            transdeffile = fullfile(ses_dir,sprintf('transdef'));
                if ~exist([transdeffile '.fif'],'file') | ~exist([transdeffile '.log'],'file') | OverWrite
                    outfile = 'tsss';
                    transtr = [' -trans default -origin ' num2str(orig(1)+0) ' ' num2str(orig(2)-13) ' ' num2str(orig(3)+6) ' -frame head -force'];
                    filestr = sprintf(' -f %s.fif -o %s.fif',outfile,transdeffile);
                    finstr = [maxfstr filestr transtr sprintf(' -v | tee %s.log',transdeffile)]
                    rik_eval(finstr);
                end
                
                rik_eval(sprintf('!echo ''Transdef for %s'' >> %s',do_tasks{ses},movfile));
                rik_eval(sprintf('!cat %s.log | sed -n ''/Position change/p'' | cut -f 7- -d '' '' >> %s',transdeffile,movfile));
            % delete(h); %EK
            %fclose(fp);
        end
        
    end %EK
    end
end

