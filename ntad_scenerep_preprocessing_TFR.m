%% Preprocessing script for Scenerep for evoked TFR analysis
% Rik Henson (2017), Modifications by Ece K (2020).
% Modifications indicated by EK

addpath /imaging/local/meg_misc/
addpath /imaging/local/software/spm_cbu_svn/releases/spm12_latest/;
addpath /hpc-software/matlab/cbu/;
addpath /imaging/projects/cbu/ntad/scripts
addpath (genpath('/imaging/projects/cbu/ntad/scripts/functions'))
% spm eeg;
addpath(genpath('/imaging/local/software/spm_cbu_svn/releases/spm12_latest/external/fieldtrip/')); %EK
addpath(genpath('/imaging/local/software/spm_cbu_svn/releases/spm12_latest/external/mne/')); %EK
addpath(genpath('/imaging/projects/cbu/ntad/scripts/functions/FastICA_25/'))

clear all; clc;

% Directories
ana_dir = '/imaging/projects/cbu/ntad/meg_data/';
task='scenerep';
session = 'TW';
% subjects = {'C1001','C1002','C1003','C1004','C1006','C1007','C1008',...
% 'C1010','C1011','C1012','C1013','C1014','C1015','C1016','C1017',...
% 'P1001','P1002','P1004','P1005','P1007','P1008','P1009','P1010',....
% 'P1011','P1012','P1015','P1016','P1020','P1021','P1022','P1023',...
% 'P1024','P1026','P1027','P1029','P1030','P1031','P1032','P1035',...
% 'P1036','P1037','P1038','P1042','P1043','P1045','P1048','P1049',...
% 'P1053','P1054','P1055','P1056','P1058','P1060','P1061','P1062',...
% 'P1063','P1064','P1065','P1066','P1067','P1070','P3002','P3004',...
% 'P3005','P3006'}; %all subjects-SCENEREP

subjects ={'P1002', 'P1004','P1008','P1009','P1011','P1022','P1035',...
   'P1038','P3002','P3004','P1037','P1049','P1055','P1060'}; %'P1027' patients with two-week data

% subjects = 'C1014','C1016','C1017, P1029, P1035,'P1067' %no EEG!

% Processing switches
% doconvert = 0; % Convert to SPM
% dodown = 0; % Downsample
% dofilter = 0; % Filter
% donotch = 0; % Notch filter
% doartefact = 0; % Artefact rejection,
% doica = 0; % Run ICA
dorename =0; % Copy file and rename with 'tfr' extension, so that we don't replace previous files.
doepoch = 0; % Epoch with padding
doartefact2 = 0; % Second artefact rejection
doaverage = 0; % Average trials and low pass
docontrast = 1; % weight conditions
dotfr= 1; % time frequency decomposition
docrop= 1; % Remove padding
dorescale= 1; % rescale amplitude 
docombine = 1; % Combine planar
dorename2 = 1; % Rename files to something sensible
doconvert2image= 0; % convert to images

% Parameters
old_labels={'NEW','NEW','NEW','NEW','NEW','NEW','NEW','NEW','NEW','NEW','NEW','NEW','NEW','NEW','NEW','NEW','NEW','NEW','NEW','NEW','NEW','NEW','NEW','NEW','NEW','NEW','NEW','OLD','OLD','OLD','OLD','NEW','OLD','OLD','OLD','NEW','NEW','NEW','NEW','OLD','NEW','OLD','OLD','NEW','NEW','OLD','OLD','NEW','OLD','OLD','NEW','OLD','OLD','OLD','OLD','NEW','OLD','NEW','NEW','NEW','NEW','OLD','OLD','OLD','NEW','NEW','OLD','OLD','NEW','NEW','NEW','NEW','OLD','NEW','NEW','OLD','NEW','NEW','NEW','OLD','NEW','OLD','NEW','OLD','OLD','OLD','NEW','OLD','NEW','OLD','OLD','OLD','OLD','NEW','NEW','OLD','NEW','NEW','OLD','NEW','OLD','OLD','OLD','NEW','NEW','OLD','OLD','OLD','NEW','NEW','NEW','NEW','NEW','NEW','OLD','OLD','OLD','NEW','OLD','OLD','OLD','NEW','OLD','NEW','NEW','OLD','OLD','OLD','OLD','NEW','OLD','NEW','OLD','NEW','NEW','OLD','NEW','OLD','NEW','NEW','NEW','NEW','NEW','NEW','OLD','NEW','OLD','NEW','OLD','NEW','OLD','NEW','OLD','OLD','OLD','OLD','OLD','NEW','OLD','NEW','OLD','OLD','NEW','OLD','NEW','NEW','NEW','NEW','NEW','OLD','OLD','NEW','OLD','OLD','OLD','NEW','OLD','OLD','NEW','OLD','OLD','NEW','OLD','NEW','OLD','OLD','NEW','NEW','NEW','NEW','NEW','OLD','NEW','OLD','OLD','OLD','NEW','OLD','OLD','OLD','OLD','OLD','NEW','NEW','NEW','OLD','OLD','OLD','OLD','NEW','NEW','OLD','OLD','NEW','OLD','NEW','NEW','OLD','OLD','NEW','OLD','NEW','OLD','OLD','OLD','NEW','NEW','NEW','OLD','OLD','NEW','OLD','NEW','NEW','NEW','OLD','NEW','OLD','OLD','OLD','OLD','OLD','OLD','OLD','NEW','OLD','NEW','NEW','OLD','NEW','OLD','OLD','NEW','NEW','NEW','OLD','NEW','OLD','NEW','NEW','OLD','OLD','NEW','NEW','NEW','OLD','OLD','NEW','NEW','NEW','NEW','NEW','NEW','NEW','OLD','OLD','NEW','NEW','OLD','OLD','NEW','NEW','OLD','NEW','OLD','NEW','OLD','OLD','NEW','NEW','OLD','NEW','NEW','OLD','OLD','OLD','OLD','OLD','OLD','NEW','OLD','OLD','NEW','OLD','NEW','NEW','OLD','NEW','NEW','OLD','NEW','NEW','OLD','OLD','NEW','OLD','OLD','OLD','OLD','OLD','OLD','NEW','NEW','OLD','NEW','OLD','OLD','NEW','OLD','NEW','OLD','NEW','NEW','NEW','NEW','OLD','NEW','OLD','NEW','NEW','NEW','OLD','NEW','NEW','OLD','OLD','OLD','NEW','OLD','OLD','NEW','NEW','OLD','NEW','OLD','OLD','OLD','OLD','NEW','OLD','OLD','NEW','NEW','NEW','NEW','OLD','OLD','NEW','OLD','NEW','OLD','NEW','OLD','NEW','NEW','OLD','OLD','OLD','NEW','NEW','OLD','NEW','NEW','NEW','NEW','OLD','OLD','NEW','NEW','NEW','OLD','NEW','NEW','NEW','NEW','OLD','NEW','NEW','NEW','NEW','NEW','NEW','NEW','NEW','OLD','NEW','NEW','NEW','OLD','NEW','NEW','NEW','OLD','OLD','NEW','NEW','OLD','NEW','NEW','NEW','OLD','OLD','OLD','NEW','OLD','NEW','OLD','NEW','OLD','OLD','NEW','OLD','OLD','NEW','NEW','NEW','OLD','OLD','OLD','NEW','NEW','NEW','OLD','NEW','NEW','OLD','NEW','OLD','NEW','OLD','NEW','NEW','NEW','OLD','OLD','NEW','OLD','OLD','NEW','NEW','OLD','OLD','OLD','OLD','OLD','NEW','NEW','OLD','OLD','OLD','OLD','NEW','NEW','NEW','OLD','NEW','OLD','OLD','NEW','OLD','OLD','OLD','NEW','NEW','OLD','NEW','OLD','NEW','NEW','OLD','NEW','NEW','NEW','NEW','NEW','OLD','NEW','NEW','OLD','OLD','NEW','NEW','OLD','OLD','OLD','NEW','OLD','NEW','NEW','NEW','OLD','OLD','OLD','OLD','OLD','NEW','OLD','OLD','OLD','OLD','NEW','NEW','NEW','NEW','NEW','OLD','OLD','OLD','OLD','OLD','OLD','OLD','OLD','OLD','NEW'};
new_trl_labels={'NEW','NEW','NEW','NEW','NEW','NEW','NEW','NEW','NEW','NEW','NEW','NEW','NEW','NEW','NEW','NEW','NEW','NEW','NEW','NEW','NEW','NEW','NEW','NEW','NEW','NEW','NEW','OLD-S','OLD-S','OLD-S','OLD-S','NEW','OLD-S','OLD-S','OLD-S','NEW','NEW','NEW','NEW','OLD-S','NEW','OLD-S','OLD-S','NEW','NEW','OLD-S','OLD-S','NEW','OLD-S','OLD-S','NEW','OLD-S','OLD-S','OLD-L','OLD-S','NEW','OLD-S','NEW','NEW','NEW','NEW','OLD-L','OLD-S','OLD-S','NEW','NEW','OLD-S','OLD-S','NEW','NEW','NEW','NEW','OLD-S','NEW','NEW','OLD-S','NEW','NEW','NEW','OLD-L','NEW','OLD-L','NEW','OLD-L','OLD-L','OLD-L','NEW','OLD-L','NEW','OLD-S','OLD-L','OLD-S','OLD-L','NEW','NEW','OLD-S','NEW','NEW','OLD-S','NEW','OLD-L','OLD-L','OLD-S','NEW','NEW','OLD-S','OLD-L','OLD-S','NEW','NEW','NEW','NEW','NEW','NEW','OLD-L','OLD-L','OLD-L','NEW','OLD-L','OLD-L','OLD-L','NEW','OLD-S','NEW','NEW','OLD-S','OLD-L','OLD-L','OLD-S','NEW','OLD-S','NEW','OLD-L','NEW','NEW','OLD-S','NEW','OLD-S','NEW','NEW','NEW','NEW','NEW','NEW','OLD-S','NEW','OLD-S','NEW','OLD-L','NEW','OLD-L','NEW','OLD-S','OLD-L','OLD-S','OLD-L','OLD-L','NEW','OLD-L','NEW','OLD-S','OLD-S','NEW','OLD-L','NEW','NEW','NEW','NEW','NEW','OLD-L','OLD-S','NEW','OLD-S','OLD-S','OLD-S','NEW','OLD-L','OLD-S','NEW','OLD-S','OLD-L','NEW','OLD-L','NEW','OLD-L','OLD-S','NEW','NEW','NEW','NEW','NEW','OLD-L','NEW','OLD-S','OLD-L','OLD-L','NEW','OLD-L','OLD-S','OLD-S','OLD-S','OLD-S','NEW','NEW','NEW','OLD-S','OLD-S','OLD-S','OLD-L','NEW','NEW','OLD-S','OLD-S','NEW','OLD-S','NEW','NEW','OLD-S','OLD-L','NEW','OLD-S','NEW','OLD-L','OLD-S','OLD-S','NEW','NEW','NEW','OLD-S','OLD-S','NEW','OLD-S','NEW','NEW','NEW','OLD-L','NEW','OLD-S','OLD-S','OLD-S','OLD-L','OLD-L','OLD-S','OLD-S','NEW','OLD-S','NEW','NEW','OLD-S','NEW','OLD-S','OLD-S','NEW','NEW','NEW','OLD-S','NEW','OLD-S','NEW','NEW','OLD-S','OLD-S','NEW','NEW','NEW','OLD-L','OLD-S','NEW','NEW','NEW','NEW','NEW','NEW','NEW','OLD-L','OLD-S','NEW','NEW','OLD-S','OLD-S','NEW','NEW','OLD-S','NEW','OLD-S','NEW','OLD-S','OLD-S','NEW','NEW','OLD-S','NEW','NEW','OLD-S','OLD-L','OLD-S','OLD-S','OLD-S','OLD-S','NEW','OLD-S','OLD-S','NEW','OLD-S','NEW','NEW','OLD-S','NEW','NEW','OLD-S','NEW','NEW','OLD-L','OLD-S','NEW','OLD-S','OLD-S','OLD-S','OLD-S','OLD-S','OLD-S','NEW','NEW','OLD-S','NEW','OLD-S','OLD-S','NEW','OLD-S','NEW','OLD-S','NEW','NEW','NEW','NEW','OLD-S','NEW','OLD-L','NEW','NEW','NEW','OLD-S','NEW','NEW','OLD-S','OLD-S','OLD-S','NEW','OLD-S','OLD-L','NEW','NEW','OLD-S','NEW','OLD-S','OLD-S','OLD-L','OLD-S','NEW','OLD-L','OLD-S','NEW','NEW','NEW','NEW','OLD-S','OLD-S','NEW','OLD-S','NEW','OLD-S','NEW','OLD-S','NEW','NEW','OLD-L','OLD-S','OLD-S','NEW','NEW','OLD-L','NEW','NEW','NEW','NEW','OLD-S','OLD-L','NEW','NEW','NEW','OLD-S','NEW','NEW','NEW','NEW','OLD-S','NEW','NEW','NEW','NEW','NEW','NEW','NEW','NEW','OLD-L','NEW','NEW','NEW','OLD-S','NEW','NEW','NEW','OLD-L','OLD-L','NEW','NEW','OLD-L','NEW','NEW','NEW','OLD-S','OLD-L','OLD-L','NEW','OLD-S','NEW','OLD-S','NEW','OLD-L','OLD-L','NEW','OLD-S','OLD-L','NEW','NEW','NEW','OLD-L','OLD-L','OLD-L','NEW','NEW','NEW','OLD-L','NEW','NEW','OLD-L','NEW','OLD-L','NEW','OLD-L','NEW','NEW','NEW','OLD-L','OLD-L','NEW','OLD-L','OLD-L','NEW','NEW','OLD-L','OLD-L','OLD-L','OLD-L','OLD-L','NEW','NEW','OLD-L','OLD-L','OLD-L','OLD-L','NEW','NEW','NEW','OLD-L','NEW','OLD-L','OLD-L','NEW','OLD-L','OLD-L','OLD-L','NEW','NEW','OLD-L','NEW','OLD-L','NEW','NEW','OLD-L','NEW','NEW','NEW','NEW','NEW','OLD-L','NEW','NEW','OLD-L','OLD-L','NEW','NEW','OLD-L','OLD-L','OLD-L','NEW','OLD-L','NEW','NEW','NEW','OLD-L','OLD-L','OLD-L','OLD-L','OLD-L','NEW','OLD-L','OLD-L','OLD-L','OLD-L','NEW','NEW','NEW','NEW','NEW','OLD-L','OLD-L','OLD-L','OLD-L','OLD-L','OLD-L','OLD-L','OLD-L','OLD-L','NEW'};
new_con_labels={'NEW','OLD-S','OLD-L'};

drate = 500; % Downsampling rate
mods = {'MEGMAG', 'MEGPLANAR','EEG'}; % modalities
pads = [-2500 2500]; % extra padding in ms to allow estimation of low frequencies
con_values = [22 33]; % trigger values that correspond to your conditions of interest
con_labels = {'NEW','OLD'}; % condition names in the same order as trigger values
epochWin = [-100 600]; % epoch window in ms
offset = [35 35]; % stimulus delivery delay for each condition 'Auditory','Visual','Motor' delay
ref_chans = {'VEOG', 'HEOG', 'ECG'}; % external references for ICA
basefile = 'transdef'; % 'transdef' Maxfilter output to preprocess
filtwindow = [0.1 100]; % Filter window (decimals will give you filter instability, EK)
PCdim = 60; % Number of PCs for ICA
samprate=1000; %sampling rate in Hz
tfr_fwin=[1 100];
tfr_ncycles= 5;
tfr_twin=[epochWin(1)+pads(1) epochWin(2)+pads(2)];
tfr_rescale_method='LogR';
tfr_rescale_baseline=[epochWin(1) 0];
crop_twin=epochWin; % time window to keep
crop_fwin=[1 100]; % frequency window to keep
% crop_channels={'EEG001', 'EEG003','EEG068','EEG070','EEG071', 'EEG072', 'EEG073','EEG074',...
%     'MEG1741','MEG1931', 'MEG2111','MEG2121', 'MEG2141', 'MEG2131', 'MEG2331', 'MEG2541',...
%     'MEG1742+1743','MEG1932+MEG1933','MEG2112+2113','MEG2122+2123',...
%     'MEG2142+2143', 'MEG2132+2133','MEG2332+2333','MEG2542+2543'}; % channel labels to keep 'all'
ses_nam=session;

%%=========================================================================
%% Setting up Matlab Pool
%%=========================================================================

%parpool(length(subjects));

%%=========================================================================
%% Convert Data to SPM
%%=========================================================================
% 
% if doconvert
%     chanfile = fullfile('chan_select_MEG_EEG_STI101.mat'); %EK
%     % Channels to read: MEG+3EEG+EOG+ECG+STI101 (excludes MISC for subjects with no eye-tracking)
%     chan = load(chanfile);
%     ref_old_names = {'EOG061','EOG062','ECG063'};
%     ref_new_names = {'HEOG';'VEOG';'ECG'};
%     
% %     parfor sub=1:length(subjects)
%     for sub = 1:length(subjects)
%         
%         sub_dir = fullfile(ana_dir,subjects{sub},session,task);
%         all_trl{sub} = {}; all_con{sub} = {};
%         ses_dir = [ana_dir subjects{sub} '/' ses_nam '/' task '/'];
%         outfile = [ses_dir basefile '.mat'];
%         
%         S = [];
%         S.outfile = outfile;
%         S.dataset = fullfile(ses_dir,sprintf('%s.fif',basefile));
%         S.save = 0;  % saved below anywapadsy (if S.save=1, then prompts for new filename!)
%         S.reviewtrials = 0;
%         S.channels = chan.label;
%         S.continuous = 1;
%         S.checkboundary = 0;
%         D = spm_eeg_convert(S);
%         tc = find(strcmp(D.chanlabels,'STI101'));
%         D(tc,:) = D(tc,:)/1e6;  % 1e6 new scaling introduced by SPM12!
%         D = units(D,tc,'V');
%         
%         for c = 1:length(ref_old_names)
%             ch = indchannel(D,ref_old_names{c});
%             D = chanlabels(D,ch,ref_new_names{c});
%         end
%         
%         D.save;
%         
%     end
% end
% 
% %%=========================================================================
% %% Downsample
% %%=========================================================================
% 
% if dodown
%     parfor sub=1:length(subjects)
% %     for sub=1:length(subjects)
% 
%             ses_dir = fullfile(ana_dir,subjects{sub},ses_nam,task,'/')
%             S = [];
%             S.D = [ses_dir '/' basefile '.mat'];
%             S.method = 'downsample';
%             S.fsample_new = drate;
% 
%             D = spm_eeg_downsample(S);
%     end
% end
% %%=========================================================================
% %% Filter
% %%=========================================================================
% 
% if dofilter
%     
%     low=filtwindow(2); high=filtwindow(1);
%     
% %     parfor sub=1:length(subjects)
%     for sub=1:length(subjects)
%         ses_dir = fullfile(ana_dir,subjects{sub},ses_nam,task)
%         S = [];
%         S.D = [ses_dir '/d' basefile '.mat'];
%         S.type = 'butterworth';
%         S.order = 5;
%         S.band = 'low';
%         S.freq = low;
%         D = spm_eeg_filter(S);    
%        
%         S= [];
%         S.D = D;
%         S.type = 'butterworth';
%         S.order = 5;
%         S.band = 'high';
%         S.freq = high;
%         D = spm_eeg_filter(S);
%         
%     end
%     
% end
% 
% %%=========================================================================
% %% Notch filter EK
% %%=========================================================================
% 
% if donotch
%     
% %    parfor sub=1:length(subjects)
%    for sub=1:length(subjects) 
%        
%        ses_dir = fullfile(ana_dir,subjects{sub},ses_nam,task)
%        S = [];
%        S.D = [ses_dir '/ffd'  basefile '.mat'];
%     
%        S.type = 'butterworth';
%        S.order = 5;
%        S.band = 'stop';
%        S.freq = [45 55];
%             
%        D = spm_eeg_filter(S);
%        
%        S = [];
%        S.D = D;
%        S.type = 'butterworth';
%        S.order = 5;
%        S.band = 'stop';
%        S.freq = [95 105];
%             
%        D = spm_eeg_filter(S);
%        
%    end 
% end
% %%=========================================================================
% %% Artefact rejection EK
% %%=========================================================================
% 
% if doartefact % uses Oxford's OSL, works much better than SPM's.
%     for sub=1:length(subjects)
%         
%         ses_dir = [ana_dir subjects{sub} '/' ses_nam '/' task '/']
%         D = spm_eeg_load([ses_dir '/ffffd' basefile '.mat']); %EK
%         D=D.badchannels(1:D.nchannels,0);
%         D2=osl_detect_artefacts(D);D=D2;
%         save([ses_dir '/affffd' basefile '.mat'],'D','-v7.3')
%         
%     end
% end
% 
% %%=========================================================================
% %% ICA
% %%=========================================================================
% 
% if doica
%     
%     ref_old_names = {'EOG061','EOG062','ECG063'};
%     ref_new_names = {'HEOG';'VEOG';'ECG'};
%     ICA = [];
%     ICA.PCA_dim = 60; % Number of PCs for ICA
%     ICA.Nperm = 0; %round((4*sqrt(CorPval)/CorPval)^2); % If want accurate p-values
%     ICA.Rseed = 1; % to make reproducible
%     all_ica_remove = {}; all_ica_TraMat = {}; all_ica_out = {};
%     
%     arttopos = load('/imaging/rh01/Methods/MEGEEG70ArtifactTemplateTopographies');
%     sparefs = {};
%     mind = find(ismember({'MEGMAG', 'MEGPLANAR', 'EEG'},mods));
%     for m=1:length(mods)
%         sparefs{m} = {};
%         for r = 1:length(ref_chans)
%             tmp = getfield(arttopos,upper(ref_chans{r}));  % If doing blinks only
%             sparefs{m}{r} = tmp{mind(m)}';
%         end
%     end
%     
% %     parfor sub=1:length(subjects)
%     for sub=1:length(subjects)
%         ses_dir = fullfile(ana_dir,subjects{sub},ses_nam,task);
%         cd(ses_dir);
%         
%         D=spm_eeg_load(['affffd' basefile '.mat']);
%         disp([subjects{sub} ' +++++++++++++++++++++++++++++++++++++++++++++++++++++++++'])
%         outfile = [ses_dir '/affffd' basefile '.mat'];
%         
%         if size(D,1)>330
%             mod={'MEGMAG', 'MEGPLANAR', 'EEG'};
%         else
%             mod={'MEGMAG', 'MEGPLANAR'};
%         end;
%         
%         icafile = fullfile(ses_dir,sprintf('%s.mat',['Maffffd' basefile]));
%         
%         S = ICA;
%         S.refs.tem = {}; % Reference signal for correlating with ICs
%         if ~isempty(find(ismember(D.chanlabels,'VEOG')))
%             S.refs.tem{1}=D(find(ismember(D.chanlabels,'VEOG')),:); %VEOG
%             S.refs.tem{2}=D(find(ismember(D.chanlabels,'HEOG')),:); %VEOG
%             S.refs.tem{3}=D(find(ismember(D.chanlabels,'ECG')),:); %VEOG
%         else
%             S.refs.tem{1}=D(find(ismember(D.chanlabels,'EOG062')),:); %VEOG
%             S.refs.tem{2}=D(find(ismember(D.chanlabels,'EOG061')),:); %VEOG
%             S.refs.tem{3}=D(find(ismember(D.chanlabels,'ECG063')),:); %VEOG
%         end
%         
%         chans = {}; remove = {}; weights = {}; temcor = {}; spacor = {}; TraMat = {};
%         for m = 1:length(mod)
%             S.refs.spa = sparefs{m};
%             chans = indchantype(D,mod{m});
%             S.d  = D(chans,:);
%             [Out,ICs] = detect_ICA_artefacts_ek(S);
%             Out.chans = chans;
%             ses_nam = session;
%             icafile = fullfile(ses_dir,sprintf('%s_ica.mat',['affffd' basefile '_' mod{m}]));
%             all_subs_ica_out{sub}{m}=Out;
%         end
%     end
%     
%     % Save ICA output
%     ses_nam=session;
%     for sub = 1:length(subjects)
%         ses_dir = fullfile(ana_dir,subjects{sub},ses_nam,task,'/'); cd(ses_dir)
%         outfile = [ses_dir '/affffd' basefile '.mat'];
%         D = spm_eeg_load(outfile);
%         if size(D,1)>330
%             mod={'MEGMAG', 'MEGPLANAR', 'EEG'};
%         else
%             mod={'MEGMAG', 'MEGPLANAR'};
%         end
%         for m = 1:length(mod)
%             icafile = fullfile(ses_dir,sprintf('%s_ica.mat',['affffd' basefile '_' mod{m}]));
%             Out=all_subs_ica_out{sub}{m};
%             parsave(icafile,Out);
%         end
%     end
%     
% %     parfor sub=1:length(subjects)
%     for sub=1:length(subjects)
%         ses_nam = session;
%         ses_dir = fullfile(ana_dir,subjects{sub},ses_nam,task);
%         TraMat = {}; chans = [];
%         D=spm_eeg_load([ses_dir '/affffd' basefile '.mat']); %original file
%         if size(D,1)>330
%             mods={'MEGMAG', 'MEGPLANAR', 'EEG'};
%         else
%             mods={'MEGMAG', 'MEGPLANAR'};
%         end
%         
%         for m = 1:length(mods)
%             icafile = fullfile(ses_dir,sprintf('%s_ica.mat',['affffd' basefile '_' mods{m}]));
%             ica = load(icafile); %ica output
%             TraMat{m} = ica.Out.TraMat;
%             chans = [chans ica.Out.chans];
%         end
%         
%         if sum(TraMat{1,1}(:))~=size(TraMat{1,1},1)
%             S = []; S.D = D;
%             for c = 1:length(chans)
%                 S.montage.labelorg{c} = D.chanlabels{chans(c)};
%             end
%             S.montage.labelnew = S.montage.labelorg;
%             S.montage.tra      = blkdiag(TraMat{:});
%             S.keepothers       = 1;
%             D = spm_eeg_montage(S); % becomes 'Mfftransdef.mat' GOODIES
%         else
%             S=[]; S.D = D;
%             S.outfile = [ses_dir '/Maffffd' basefile '.mat'];
%             D= spm_eeg_copy(S);
%         end
%     end
%     
%     % Plot removed ICs
%     for sub=1:length(subjects)
%         ses_nam = session;
%         ses_dir = fullfile(ana_dir,subjects{sub},ses_nam,task);
%         D=spm_eeg_load([ses_dir '/affffd' basefile '.mat']); %original file
%         if size(D,1)>330; mod={'MEGMAG','MEGPLANAR','EEG'}; else; mod={'MEGMAG','MEGPLANAR'};end
%         for m=1:length(mod)
%             load([ses_dir '/affffd' basefile '_' mod{m} '_ica.mat']);
%             w=pinv(Out.weights)'; tcorrs=Out.temcor;scorrs=Out.spacor;
%             idx=find(strcmp(D.chantype,mod{m}));
%             in.noButtons=1;in.type=mod{m};in.plotpos=0;
%             %w = weights(idx,:)';
%             close all; 
%             rem=unique([Out.temrem, Out.sparem]);
%             for com=1:length(rem)
%                 bc=rem(com);
%                 figure
%                 set(gca,'FontSize',9)
%                 [ZI,f] =spm_eeg_plotScalpData_ek2(w(bc,:)',D.coor2D(idx),D.chanlabels(idx),in);
%                 close(figure(2));
%                 imagesc(flipdim(ZI,1));
%                 colordata=colormap;cc=colordata(1,:);
%                 colordata(4:67,:)=colormap;colordata(2:3,:)=repmat(cc,2,1);colordata(1,:)=[1 1 1];
%                 colormap(colordata); clear colordata cc
%                 axis off;
%                 if tcorrs(1,bc)>=0.3 || scorrs(1,bc)>=0.7
%                     txt=['HEOG_tem_' num2str(tcorrs(1,bc),'%.2f') '_spa_' num2str(scorrs(1,bc),'%.2f')];
%                 elseif tcorrs(2,bc)>=0.3 || scorrs(2,bc)>=0.7
%                     txt=['VEOG_tem_' num2str(tcorrs(2,bc),'%.2f') '_spa_' num2str(scorrs(2,bc),'%.2f')]; 
%                 elseif tcorrs(3,bc)>=0.3 || scorrs(3,bc)>=0.7
%                     txt=['ECG_tem_' num2str(tcorrs(3,bc),'%.2f') '_spa_' num2str(scorrs(3,bc),'%.2f')]; 
%                 end
%                 title({[ 'IC' num2str(bc) ' ' strrep(txt,'_',' ') ' ' mod{m}]});
%                 print(figure(1),'-dtiff', [ses_dir '/' mod{m} '_ic' num2str(bc) '_' subjects{sub} '_' txt '.tiff'],'-r200');
%                 close(figure(1));
%             end
%         end
%     end
% end

%%=========================================================================
%% Rename files
%%=========================================================================

if dorename
    
    for sub=1:length(subjects)
        ses_dir = [ana_dir subjects{sub} '/' ses_nam '/' task '/'];
        S=[]; S.D = spm_eeg_load([ses_dir '/Maffffd' basefile '.mat' ]);
        S.outfile = [ses_dir '/Maffffd' basefile '_tfr.mat'];
        D= spm_eeg_copy(S);
    end
    
end

basefile=[basefile '_tfr'];
%%=========================================================================
%% Read triggers & Epoch
%%=========================================================================

if doepoch
    all_trl = {}; all_con = {};
    padding = pads;
    ses_nam = session;
    
%     parfor sub=1:length(subjects)
    for sub=1:length(subjects)
        if exist(fullfile(ana_dir,subjects{sub},ses_nam,'/',task),'dir')
            
            ses_dir = [ana_dir subjects{sub} '/' ses_nam '/' task '/'];
%             if ~exist([ses_dir 'trial_info.txt'],'file')
%                 filelog=fopen([ses_dir 'trial_info.txt'],'w'); 
%                 fclose(filelog);%EK
%             end
            
%             filelog=fopen([ses_dir 'trial_info.txt'],'a');%EK
            D = spm_eeg_load([ses_dir 'dtransdef.mat']);
            d = D(indchannel(D,'STI101'),:); 
            if max(d(:))>5000; d=d/1000000; end; d=round(d); %EK
            if isequal(unique(d),0); d = D(indchannel(D,'STI101'),:); end
            fprintf('unique trigger codes: %s\n',num2str(unique(d)))  % triggers found
%             fprintf(filelog, 'unique trigger codes: %s\n',num2str(unique(d)));
                   
            dd = [0 diff(d)]; % onset of triggers
            ons = find(ismember(dd,con_values));
            trl = []; conditionlabels = {};
            for o = 1:length(ons)
                con = find(con_values == dd(ons(o)));
                conditionlabels{end+1} = con_labels{con};
                trl(end+1,1) = ons(o) + round((epochWin(1) + padding(1) + offset(con))*(D.fsample/1000));
                trl(end,2)   = ons(o) + round((epochWin(2) + padding(2) + offset(con))*(D.fsample/1000));
                trl(end,3)   = round((epochWin(1) + padding(1))*(D.fsample/1000));
            end
            
            trl(:,2)=trl(:,1)+min(trl(:,2)-trl(:,1));
            
%             fprintf('%s\n', subjects{sub});
%             for con = 1:length(con_values)
%                 fprintf('%d \ttriggers of type %s\n',length(find(strcmp(conditionlabels,con_labels{con}))),con_labels{con});
%                 fprintf(filelog, '%d \ttriggers of type %s\n',length(find(strcmp(conditionlabels,con_labels{con}))),con_labels{con}); %EK
%             end

%             fclose(filelog);%EK

% Change labels EK
labs={}; mis=[];
labs(:,1)=conditionlabels;
labs(:,2)=old_labels(1:length(conditionlabels));
idx=strcmp(labs(:,1),labs(:,2));te=find(idx==0);k=0;
while ~isempty(te)
    te=te(1);
    k=k+1;
    mis=[mis, te+k-1]; % find skips
    l=old_labels;l(mis)=[];
    labs(:,2)=l(1:length(conditionlabels));
    idx=strcmp(labs(:,1),labs(:,2));te=find(idx==0);
end
l=new_trl_labels;l(mis)=[];
conditionlabels=l;

% for con = 1:length(new_con_labels)
%     fprintf('%d \ttriggers of type %s\n',length(find(strcmp(conditionlabels,new_con_labels{con}))),new_con_labels{con});
%     fprintf(filelog, '%d \ttriggers of type %s\n',length(find(strcmp(conditionlabels,new_con_labels{con}))),new_con_labels{con}); %EK
% end




all_trl{sub} = trl; all_con{sub} = conditionlabels;
        end
    end
    
    for sub = 1:length(subjects)
        if exist(fullfile(ana_dir,subjects{sub},ses_nam,task))
            ses_nam = session;
            ses_dir = fullfile(ana_dir,subjects{sub},ses_nam,task);
            trl = all_trl{sub}; conditionlabels = all_con{sub};
            trlfile = fullfile(ses_dir,sprintf('%strl_new_lbl_tfr.mat',''));
            save(trlfile,'trl','conditionlabels');
        end
    end
    
    bline=epochWin(1);
%         parfor sub=1:length(subjects)
    parfor sub=1:length(subjects)
        
        ses_dir = [ana_dir subjects{sub} '/' ses_nam '/' task '/']
        trl = load(fullfile(ses_dir,sprintf('%strl_new_lbl_tfr.mat','')));
        %trl = load(fullfile(ses_dir,sprintf('%s_trl.mat',basefile(1:4))));
        D = spm_eeg_load([ses_dir '/Maffffd' basefile '.mat']); %EK
        S = []; S.D = D;
        S.trl = trl.trl;
        S.conditionlabels = trl.conditionlabels;
        S.eventpadding = 0;  % Already included above!
        S.bc = 0;            % No need to baseline correct - done below...
        D = spm_eeg_epochs(S);
        
        % Baseline correct (ie ignoring padding)
        S = []; S.D = D;
        S.timewin = [bline 0];
        S.save = 1;
        S.prefix='';
        D = spm_eeg_bc(S);
        
    end    
end
%%=========================================================================
%% Artefact rejection and removal 
%%=========================================================================

if doartefact2 % 
    for sub=1:length(subjects)
        
        ses_dir = [ana_dir subjects{sub} '/' ses_nam '/' task '/']
        D = spm_eeg_load([ses_dir '/eMaffffd' basefile '.mat']); %EK
        
        if ~isequal(max(ismember(D.chantype,'MEGMAG')),1)
        load([ses_dir '/eMaffffd' basefile '.mat'])
        chans={'MEG0111';'MEG0121';'MEG0131';'MEG0141';'MEG0211';'MEG0221';'MEG0231';'MEG0241';'MEG0311';'MEG0321';'MEG0331';'MEG0341';'MEG0411';'MEG0421';'MEG0431';'MEG0441';'MEG0511';'MEG0521';'MEG0531';'MEG0541';'MEG0611';'MEG0621';'MEG0631';'MEG0641';'MEG0711';'MEG0721';'MEG0731';'MEG0741';'MEG0811';'MEG0821';'MEG0911';'MEG0921';'MEG0931';'MEG0941';'MEG1011';'MEG1021';'MEG1031';'MEG1041';'MEG1111';'MEG1121';'MEG1131';'MEG1141';'MEG1211';'MEG1221';'MEG1231';'MEG1241';'MEG1311';'MEG1321';'MEG1331';'MEG1341';'MEG1411';'MEG1421';'MEG1431';'MEG1441';'MEG1511';'MEG1521';'MEG1531';'MEG1541';'MEG1611';'MEG1621';'MEG1631';'MEG1641';'MEG1711';'MEG1721';'MEG1731';'MEG1741';'MEG1811';'MEG1821';'MEG1831';'MEG1841';'MEG1911';'MEG1921';'MEG1931';'MEG1941';'MEG2011';'MEG2021';'MEG2031';'MEG2041';'MEG2111';'MEG2121';'MEG2131';'MEG2141';'MEG2211';'MEG2221';'MEG2231';'MEG2241';'MEG2311';'MEG2321';'MEG2331';'MEG2341';'MEG2411';'MEG2421';'MEG2431';'MEG2441';'MEG2511';'MEG2521';'MEG2531';'MEG2541';'MEG2611';'MEG2621';'MEG2631';'MEG2641';'MEG0113';'MEG0112';'MEG0122';'MEG0123';'MEG0132';'MEG0133';'MEG0143';'MEG0142';'MEG0213';'MEG0212';'MEG0222';'MEG0223';'MEG0232';'MEG0233';'MEG0243';'MEG0242';'MEG0313';'MEG0312';'MEG0322';'MEG0323';'MEG0333';'MEG0332';'MEG0343';'MEG0342';'MEG0413';'MEG0412';'MEG0422';'MEG0423';'MEG0432';'MEG0433';'MEG0443';'MEG0442';'MEG0513';'MEG0512';'MEG0523';'MEG0522';'MEG0532';'MEG0533';'MEG0542';'MEG0543';'MEG0613';'MEG0612';'MEG0622';'MEG0623';'MEG0633';'MEG0632';'MEG0642';'MEG0643';'MEG0713';'MEG0712';'MEG0723';'MEG0722';'MEG0733';'MEG0732';'MEG0743';'MEG0742';'MEG0813';'MEG0812';'MEG0822';'MEG0823';'MEG0913';'MEG0912';'MEG0923';'MEG0922';'MEG0932';'MEG0933';'MEG0942';'MEG0943';'MEG1013';'MEG1012';'MEG1023';'MEG1022';'MEG1032';'MEG1033';'MEG1043';'MEG1042';'MEG1112';'MEG1113';'MEG1123';'MEG1122';'MEG1133';'MEG1132';'MEG1142';'MEG1143';'MEG1213';'MEG1212';'MEG1223';'MEG1222';'MEG1232';'MEG1233';'MEG1243';'MEG1242';'MEG1312';'MEG1313';'MEG1323';'MEG1322';'MEG1333';'MEG1332';'MEG1342';'MEG1343';'MEG1412';'MEG1413';'MEG1423';'MEG1422';'MEG1433';'MEG1432';'MEG1442';'MEG1443';'MEG1512';'MEG1513';'MEG1522';'MEG1523';'MEG1533';'MEG1532';'MEG1543';'MEG1542';'MEG1613';'MEG1612';'MEG1622';'MEG1623';'MEG1632';'MEG1633';'MEG1643';'MEG1642';'MEG1713';'MEG1712';'MEG1722';'MEG1723';'MEG1732';'MEG1733';'MEG1743';'MEG1742';'MEG1813';'MEG1812';'MEG1822';'MEG1823';'MEG1832';'MEG1833';'MEG1843';'MEG1842';'MEG1912';'MEG1913';'MEG1923';'MEG1922';'MEG1932';'MEG1933';'MEG1943';'MEG1942';'MEG2013';'MEG2012';'MEG2023';'MEG2022';'MEG2032';'MEG2033';'MEG2042';'MEG2043';'MEG2113';'MEG2112';'MEG2122';'MEG2123';'MEG2133';'MEG2132';'MEG2143';'MEG2142';'MEG2212';'MEG2213';'MEG2223';'MEG2222';'MEG2233';'MEG2232';'MEG2242';'MEG2243';'MEG2312';'MEG2313';'MEG2323';'MEG2322';'MEG2332';'MEG2333';'MEG2343';'MEG2342';'MEG2412';'MEG2413';'MEG2423';'MEG2422';'MEG2433';'MEG2432';'MEG2442';'MEG2443';'MEG2512';'MEG2513';'MEG2522';'MEG2523';'MEG2533';'MEG2532';'MEG2543';'MEG2542';'MEG2612';'MEG2613';'MEG2623';'MEG2622';'MEG2633';'MEG2632';'MEG2642';'MEG2643';'EEG001';'EEG002';'EEG003';'EEG004';'EEG005';'EEG006';'EEG007';'EEG008';'EEG009';'EEG010';'EEG011';'EEG012';'EEG013';'EEG014';'EEG015';'EEG016';'EEG017';'EEG018';'EEG019';'EEG020';'EEG021';'EEG022';'EEG023';'EEG024';'EEG025';'EEG026';'EEG027';'EEG028';'EEG029';'EEG030';'EEG031';'EEG032';'EEG033';'EEG034';'EEG035';'EEG036';'EEG037';'EEG038';'EEG039';'EEG040';'EEG041';'EEG042';'EEG043';'EEG044';'EEG045';'EEG046';'EEG047';'EEG048';'EEG049';'EEG050';'EEG051';'EEG052';'EEG053';'EEG054';'EEG055';'EEG056';'EEG057';'EEG058';'EEG059';'EEG060';'EEG065';'EEG066';'EEG067';'EEG068';'EEG069';'EEG070';'EEG071';'EEG072';'EEG073';'EEG074';'HEOG';'VEOG';'ECG';'STI101'};
        types={'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EOG';'EOG';'ECG';'Other'};
        for ch=1:size(D.data,1)
            D.channels(ch).label=chans{ch};
            D.channels(ch).type=types{ch};
        end; save([ses_dir '/eMaffffd' basefile '.mat'],'D','-v7.3')
        clear D; D = spm_eeg_load([ses_dir '/eMaffffd' basefile '.mat']); %EK
        end
        D=D.badchannels(1:D.nchannels,0);
        D2=osl_detect_artefacts(D);D=D2;
        save([ses_dir '/aeMaffffd' basefile '.mat'],'D','-v7.3')
        D = spm_eeg_load([ses_dir '/aeMaffffd' basefile '.mat']); %EK
        S = []; S.D = D;
        D = spm_eeg_remove_bad_trials(S);
        
    end
end

%%=========================================================================
%% Average
%%=========================================================================

if doaverage
    %parfor
   parfor sub=1:length(subjects)
        %if exist(fullfile(ana_dir,subjects{sub},ses_nam,'/'))
            ses_dir = [ana_dir subjects{sub} '/' ses_nam '/' task '/']
            D = spm_eeg_load([ses_dir '/raeMaffffd' basefile '.mat']);
            
            S = []; S.D = D;
            S.robust.savew = 0;
            S.robust.bycondition = 1;
            S.robust.ks = 3;
            S.robust.removebad = 1;
            
            D = spm_eeg_average(S);
        %end
    end
    
    low=filtwindow(2);
    parfor sub=1:length(subjects)
%     for sub=1:length(subjects)
        ses_dir = fullfile(ana_dir,subjects{sub},ses_nam,task)
        S = [];
        S.D = [ses_dir '/mraeMaffffd' basefile '.mat'];
        S.type = 'butterworth';
        S.order = 5;
        S.band = 'low';
        S.freq = low;
        D = spm_eeg_filter(S);
    end
end

%% Contrast
%%=========================================================================
% 
if docontrast 
       
    parfor sub=1:length(subjects)
  % for sub=1:length(subjects)
        
        ses_dir = [ana_dir subjects{sub} '/' ses_nam '/' task '/']
        D = spm_eeg_load([ses_dir '/fmraeMaffffd' basefile '.mat']);
        
        S = [];
        S.D = D;
        S.label = {'NEW','OLD','OLD-S','OLD-L'};
        S.c = zeros(4,3);
        i=find(strcmp(D.conditions,'NEW'))
        S.c(1,i) = 1;
        j=find(strcmp(D.conditions,'OLD-S'))
        S.c(2,j) = 0.5;
        j=find(strcmp(D.conditions,'OLD-L'))
        S.c(2,j) = 0.5;
        k=find(strcmp(D.conditions,'OLD-S'))
        S.c(3,k) = 1;
        k=find(strcmp(D.conditions,'OLD-L'))
        S.c(4,k) = 1;
        
        S.weighted = 0;
        D = spm_eeg_contrast(S);
        
    end
end

%%=========================================================================
%% Time-frequency decomposition 
%%=========================================================================
if dotfr
    
    f1=tfr_fwin(1); f2=tfr_fwin(2);
    parfor sub=1:length(subjects)
        
        ses_dir = [ana_dir subjects{sub} '/' ses_nam '/' task '/']
        D = spm_eeg_load([ses_dir 'wfmraeMaffffd' basefile '.mat']);
        S   = []; S.D = D;
        S.frequencies=f1:f2;
        S.timewin= tfr_twin;
        S.method='morlet';
        S.settings.ncycles=tfr_ncycles;
        S.settings.subsample=1; % takes samples of the data instead of the whole thing
        S.phase=0;
        D= spm_eeg_tf(S);
        
    end
end

%%=========================================================================
%% Remove padding and select sensor data
%%=========================================================================

if docrop
    
    parfor sub=1:length(subjects)
        
        ses_dir = [ana_dir subjects{sub} '/' ses_nam '/' task '/']
        D = spm_eeg_load([ses_dir 'tf_wfmraeMaffffd' basefile '.mat']);
        S = []; S.D = D;
        
        S.timewin = crop_twin;
        S.freqwin=crop_fwin;
        %S.channels=crop_channels;
        D = spm_eeg_crop(S);
        
    end
    
end

%%=========================================================================
%% Rescale TFRs
%%=========================================================================
if dorescale
    
    parfor sub=1:length(subjects)
        %parfor
        %if exist(fullfile(ana_dir,subjects{sub},ses_nam,'/'))
            ses_dir = [ana_dir subjects{sub} '/' ses_nam '/' task '/']
            D = spm_eeg_load([ses_dir '/ptf_wfmraeMaffffd' basefile '.mat']);
            S=[]; S.D=D;
            S.method='LogR';
            S.timewin=tfr_rescale_baseline;
            
            D=spm_eeg_tf_rescale(S);
        %end
    end
end

%%=========================================================================
%% Combine planar
%%=========================================================================

if docombine
    parfor sub=1:length(subjects)
        ses_dir = [ana_dir subjects{sub} '/' ses_nam '/' task '/'];
        
        S=[];
        S.D =[ses_dir '/rptf_wfmraeMaffffd' basefile '.mat' ];
        S.mode = 'replace';
        D = spm_eeg_combineplanar(S);
        
    end
end

%%=========================================================================
%% Rename files
%%=========================================================================

if dorename2
    
    for sub=1:length(subjects)
        ses_dir = [ana_dir subjects{sub} '/' ses_nam '/' task '/'];
        S=[]; S.D = spm_eeg_load([ses_dir '/Prptf_wfmraeMaffffd' basefile '.mat' ]);
        S.outfile = [ses_dir '/' subjects{sub} '_' ses_nam '_' task '_tfr_new_lbl.mat'];
        D= spm_eeg_copy(S);
    end
    
end

%%=========================================================================
%% Convert to images
%%=========================================================================

if doconvert2image
    for sub=1:length(subjects)
        
        ses_dir = [ana_dir subjects{sub} '/' ses_nam '/' task '/']
        
        D = spm_eeg_load([ses_dir '/' subjects{sub} '_' ses_nam '_' task '_tfr.mat']);
        S=[]; S.D=D;
        S.mode='time x frequency';
        S.timewin=crop_twin;
        S.freqwin=crop_fwin;
        S.channels='MEG';
        S.prefix='MEG_';
        [D,S]=spm_eeg_convert2images(S);
        
        D = spm_eeg_load([ses_dir '/' subjects{sub} '_' ses_nam '_' task '_tfr.mat']);
        S=[]; S.D=D;
        S.mode='time x frequency';
        S.timewin=crop_twin;
        S.freqwin=crop_fwin;
        S.channels='MEGCOMB';
        S.prefix='MEGCOMB_';
        [D,S]=spm_eeg_convert2images(S);
        
        D = spm_eeg_load([ses_dir '/' subjects{sub} '_' ses_nam '_' task '_tfr.mat']);
        S=[]; S.D=D;
        S.mode='time x frequency';
        S.timewin=crop_twin;
        S.freqwin=crop_fwin;
        S.channels='EEG';
        S.prefix='EEG_';
        [D,S]=spm_eeg_convert2images(S);
        
    end
end

