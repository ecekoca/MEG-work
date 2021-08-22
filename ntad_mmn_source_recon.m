%% Automated Source Reconstruction

% (modified by Ece K, inherited from Alex C who inherited from Rik H)

% There is a bug in SPM12 (11/2017), cannot coregister montaged E/MEG files.
% So you need to change the lines in the spm function to put ft_plot_sens in
% a try-catch loop. Also your grads should not be combined.

addpath(genpath('/imaging/local/software/spm_cbu_svn/releases/spm12_latest/external/fieldtrip/'));
addpath /imaging/local/software/spm_cbu_svn/releases/spm12_latest/external/bemcp/;
addpath /imaging/local/software/spm_cbu_svn/releases/spm12_latest/external/mne/;
addpath /imaging/projects/cbu/ntad/scripts/functions 

clear all; clc;
bwd = '/imaging/projects/cbu/ntad/meg_data/';
ana = '/BL/mmn/';

filebegin=''; % 
fileend='_BL_mmn_all_trl.mat';

subs = {'C1001','C1002','C1003','C1004','C1006','C1007',...
'C1008','C1010','C1011','C1012','C1013','C1014','C1015','C1016','C1017',...
'c2','c5','c7','c10','c12','c13','c14',...
'P1002','P1004','P1005','P1007','P1008','P1009',...
'P1010','P1011','P1012','P1015','P1020','P1021','P1022','P1023',...
'P1024','P1026','P1027','P1029','P1030','P1031','P1032','P1035',...
'P1036','P1037','P1038','P1042','P1043','P1045','P1048','P1049',...
'P1053','P1055','P1056','P1058','P1060','P1061','P1062','P1063',...
'P1064','P1065','P1066','P1067','P1070','P3002','P3004','P3005'}; 

% Parameters
epo = [-100 400]; % inversion window in ms, which is your epoch
f_win=[0.1 45];% frequency range to invert 0:all [0 4]:delta [4 8]:theta [8 12]:alpha [12 30]:beta  [30 100]:gamma
inv_typ = 'COH';    % inversion method 'IID' 'GS' 'MSP'
mesh_size = 2;   % [1-3] for coarse-medium-fine
val=1; % arbitrary inversion identifier
hanning=0; % hanning window switch
trialtype='evoked'; %'induced'%depends if it's single trial or average
inv_trls = {'DEV','REP5','DEV-STD'}; % conditions, leave 'Undefined' for resting state
use_headshape = 1; % switch for inclusion of digitisation points
con_win={[100 200], [0 400]}; % focused window,[0 500]

% Processing options
% dodicomconvert = 1; % convert T1's dicom to nifti (already done in ntad_organise.m)
dogetmesh = 1; % process mr
docoregister = 1; % coregistration
autocoregister = 0; % use coregister if you do not have fiducial mat file
doformod = 1; % forward modelling
doinversion= 1; % invert
dowrite= 1; % write niftis of time windows

%% ================================================================

cd(bwd)

for ss = 1:length(subs)
    
    if contains(subs{ss},'c'); sub_dir=[bwd '/TGB_controls/' subs{ss} ana]; else;
    sub_dir=[bwd subs{ss} ana]; end
    
    if exist([sub_dir subs{ss} filebegin fileend])
        cd(sub_dir);
        D = spm_eeg_load([sub_dir subs{ss} filebegin fileend]);
        
        if size(D,1)<350 %EK
            inv_mods={'MEGPLANAR', 'MEG'};
        else
            inv_mods={'MEGPLANAR', 'MEG','EEG'};
        end
        
        % =====================================================================
        %% Dicom conversion EK
        % =====================================================================
        
%         if dodicomconvert %EK
%             
%             if strcmp(submr{ss},'Tem')==0
%                 cd(dicomdir)
%                 l=dir([dicomdir 'CBU' submr{ss} '*']);folname=[dicomdir l.name]; cd(folname)
%                 l=dir([folname '/20*']); folname=[l.folder '/' l.name '/']; cd(folname)
%                 l=dir([folname '*MPRAGE*']);
%                 if length(l)>1
%                     folname=[l(end).folder '/' l(end).name '/'];
%                 else
%                     folname=[l.folder '/' l.name '/'];
%                 end
%                 cd(folname);l=dir(folname);
%                 filenames=[]; flen=[];
%                 for f=1:length(l); flen=[flen, length([l(f).folder '/' l(f).name])]; end
%                 ns=(ones(1,length(l)).*max(flen))-flen;
%                 for f=3:length(l)
%                     filenames=[filenames; [l(f).folder '/' l(f).name blanks(ns(f))]];
%                 end; clear f
%                 hdr=spm_dicom_headers(filenames,0);
%                 mkdir([mridir subs{ss} '/']);
%                 spm_dicom_convert(hdr,'all','flat','nii',[mridir subs{ss} '/']); clear hdr filenames l folname
%             end
%             
%         end
        
        
        % =====================================================================
        %% Get cortical mesh
        % =====================================================================
        
        if dogetmesh
            
            if contains(subs{ss},'c'); submr=[bwd '/TGB_controls/' subs{ss} '/sMRI/']; else
            submr=[bwd subs{ss} '/sMRI/']; end
            
            if exist(submr) && ~isempty(submr) %EK
                cd(submr)
                if contains(subs{ss},'c'); l=dir('c*.nii');else; l=dir('sMR*.nii');end
                mri=l.name;
                D.val = val;
                D.inv{val}.date    = strvcat(date,datestr(now,15));
                D.inv{val}.comment = {sprintf('%s_%s',inv_typ,char(inv_mods)')};    % remember inversion type and sensor configuration
                D.inv{val}.mesh    = [];
                D.inv{val}.datareg = [];
                D.inv{val}.forward = [];
                D.inv{val}.mesh.sMRI = [submr mri];
            else
                D.inv{val}.mesh.sMRI = '/imaging/projects/cbu/ntad/scripts/canonical_brain/single_subj_T1.nii';
            end
            
            if isempty(D.inv{val}.mesh.sMRI)
                error('MRI not found')
            end
            
            if exist(submr) && ~isempty(submr)
                D.inv{val}.mesh.def  = spm_select('FPList',[submr mri],'^y_sM.*\.nii$');% Check whether inverse-normalised surfaces present
            else
                D.inv{val}.mesh.def  = spm_select('FPList','/imaging/projects/cbu/ntad/scripts/canonical_brain/single_subj_T1.nii','^y_sM.*\.nii$');
            end
            
            if isempty(D.inv{val}.mesh.def)
                D.inv{val}.mesh = spm_eeg_inv_mesh(D.inv{val}.mesh.sMRI, mesh_size);
            end
            D.save;
            cd(sub_dir);
            
        elseif ~isequal(val,1)
            D.inv{val}.mesh=D.inv{1}.mesh;
        end
        
        
        % =====================================================================
        %% Co-registration
        % =====================================================================
        
        if docoregister
            
            newmrifid           = [];
            if autocoregister
                newmrifid.fid.pnt   = D.inv{val}.mesh.fid.fid.pnt(1:3,:); % EK force selection
            else
                cd(submr); load('fiducials_num_array.mat');
                newmrifid.fid.pnt = fiducials_num_array;
            end
            
            newmrifid.fid.label = {'Nasion';'LPA';'RPA'};
            newmrifid.pnt       = D.inv{val}.mesh.fid.pnt;
            
            meegfid = D.fiducials;
            meegfid.pnt(find(meegfid.pnt(:,2)>0 & meegfid.pnt(:,3)<0),:) = [];
            
            D = spm_eeg_inv_datareg_ui_ek(D, val, meegfid, newmrifid,use_headshape);
            D.save;
            
        elseif ~isequal(val,1)
            
            D.inv{val}.datareg=D.inv{1}.datareg;
            
        end
        
        % =====================================================================
        %% Forward Modelling
        % =====================================================================
        
        if doformod
            
            fprintf(1, 'Creating forward model\n');
            D.inv{val}.forward = struct([]);
            
            for ind = 1:length(D.inv{val}.datareg) %EK
                if strcmp(D.inv{val}.datareg(ind).modality,'MEG')
                    D.inv{val}.forward(ind).voltype = 'Single Shell';
                elseif strcmp(D.inv{val}.datareg(ind).modality,'EEG')
                    D.inv{val}.forward(ind).voltype = 'EEG BEM';
                end
            end
            
            fprintf(1, 'Computing leadfield\n');
            D = spm_eeg_inv_forward(D,val);
            D.save;
            
        elseif ~isequal(val,1)
            D.inv{val}.forward=D.inv{1}.forward;
        end
        
        
        % =====================================================================
        %% Inversion
        % =====================================================================
        
        if doinversion
            
            D.inv{val}.inverse = [];
            D.inv{val}.inverse.trials = inv_trls ;
            D.inv{val}.inverse.type   = inv_typ;
            D.inv{val}.inverse.lpf    = f_win(1);
            D.inv{val}.inverse.hpf    = f_win(2);
            D.inv{val}.inverse.woi    = epo;% as much as the memory can handle
            D.inv{val}.inverse.Han    = hanning; %no hanning
            D.inv{val}.inverse.modality = inv_mods;
            D.inv{val}.inverse.dplot = 0;
            
            D = spm_eeg_invert_ek(D,val);
            D.save;
            
            if D.inv{val}.inverse.R2<80 % Check if inversion R2 is greater than 80, if not remove EEG and repeat inversion EK
                 
                 D.inv{val}.inverse = [];
                 D.inv{val}.inverse.trials = inv_trls ;
                 D.inv{val}.inverse.type   = inv_typ;
                 D.inv{val}.inverse.lpf    = f_win(1);
                 D.inv{val}.inverse.hpf    = f_win(2);
                 D.inv{val}.inverse.woi    = epo;% as much as the memory can handle
                 D.inv{val}.inverse.Han    = hanning; %no hanning
                 D.inv{val}.inverse.modality = {'MEGPLANAR', 'MEG'};
                 D.inv{val}.inverse.dplot = 0;
                 
                 D = spm_eeg_invert_ek(D,val);
                 
             end
             
            save([sub_dir subs{ss} filebegin fileend],'D','-v7.3'); %EK
            D.save
        end
        
        % =====================================================================
        %% Write niftis
        % =====================================================================
        
        if dowrite
            
            for cont=1:length(con_win)
                
                D.inv{val}.contrast.woi  = con_win{cont};
                D.inv{val}.contrast.fboi = f_win;
                D.inv{val}.contrast.type =  'evoked';
                D = spm_eeg_inv_results(D);
                
                D.con=1;
                spm_eeg_inv_results_display(D);    % Display first condition
                cd(bwd);
                D.inv{val}.contrast.smoothing = 2;
                D = spm_eeg_inv_Mesh2Voxels_ek(D);
                
            end
            D.save;
        end
        
    end
end

%% Check inversion accuracies %EK
k=1;
for ss=1:length(subs)
    sub_dir=[bwd subs{ss} ana];
    if exist([sub_dir subs{ss} filebegin fileend])
        cd(sub_dir);
        sub{k}=subs{ss};
        D=spm_eeg_load([sub_dir subs{ss} filebegin fileend]);
        r2(k)=D.inv{val}.inverse.R2;
        pve(k)=D.inv{val}.inverse.VE;
        k=k+1;
    end
end

%EK
T=table(sub',r2',nonzeros(pve),'VariableNames',{'SS','R2','VE'}) % subject, R squared, variance explained
disp('Do the R2 values look alright for everyone?')


