%% Extracts timecourses from vertices within an ROI mask image
% Alex C (10/2012)
% Modifications Ece K (10/2015)

% General variable setup
%addpath('/work/imaging8/MEG/bananas_2015_EK/scripts/NIfTI_20140122/');
addpath /imaging/local/meg_misc/
addpath /hpc-software/matlab/cbu/;

clear all; clc;
datadir = '/imaging/projects/cbu/ntad/meg_data/';
maskdir='/imaging/ek01/atlases/AAL/names_long/';
outdir='/imaging/ek01/ntad/mmn_dcm/dcm_input/r0_neg_N100/';

subs = {'C1001','C1002','C1003','C1004','C1006','C1007',...
   'C1008','C1010','C1011','C1012','C1013','C1014','C1015',...
   'C1016','C1017','P1002','P1004','P1005','P1007','P1008','P1009',...
    'P1010','P1011','P1012','P1015','P1020','P1021','P1022','P1023',...
    'P1024','P1026','P1027','P1029','P1030','P1031','P1032','P1035',...
    'P1036','P1037','P1038','P1042','P1043','P1048','P1049',...
    'P1053','P1055','P1056','P1058','P1060','P1061','P1062','P1063',...
    'P1064','P1065','P1066','P1067','P1070','P3002','P3004','P3005'}; %MMN 'P1001',

% subs ={'P1002', 'P1004','P1008','P1009','P1011','P1022','P1035',...
% 'P1038','P3002','P3004','P1037','P1049','P1055','P1060'}; % patients with two-week data

midfol='/BL/mmn/';
session='BL';

% Filename of ROI masks and nicname used in naming new files
maskname = {'r1mm_ROI_MNI_V4_Frontal_Inf_L.nii','r1mm_ROI_MNI_V4_Temporal_Sup_L.nii','r1mm_ROI_MNI_V4_Heschl_L.nii',...
    'r1mm_ROI_MNI_V4_Frontal_Inf_R.nii','r1mm_ROI_MNI_V4_Temporal_Sup_R.nii','r1mm_ROI_MNI_V4_Heschl_R.nii'};
masknic = {'LIFG','LSTG','LHG','RIFG','RSTG','RHG'}; % Region names that correspond to the coordinates below

homo = 0; % do RH homologues of LH ROIs?
reslice=0;
newvoxel=[1 1 1];

% Parts of filename with source localised info in
front = '';
lastbit = 'raeMaffffdtransdef.mat';

% Options
wholebrain=0; % normal mesh 8196 vertices. if you want more you need to increase your mesh size
usemask =0;
virtualelectrode=1;% use ROI or single coord
val = 1;  % inversion number
virtual_el_coor=[-46 -20 8;
    -61 -32 8;
    -42 -22 7;
    46 20 8;
    59 -25 8;
    46 -14 8];

%% Get coordinates of vertices in ROIs

% if length(subs)>1 %EK
%     %matlabpool close force CBU_Cluster % this bit works well with matlab 2012b
%     if matlabpool('size')==0
%         MaxNsubs = max([16 length(subs)]);
%         P = cbupool(MaxNsubs);
%         P.ResourceTemplate = '-l nodes=^N^,mem=4GB,walltime=140:00:00';
%         P.SubmitArguments='-W x="NODESET:ONEOF:FEATURES:MAXFILTER"';
%         matlabpool(P);
%     end
% end
if ~exist(outdir)
    mkdir(outdir)
end

XYZ_all=[];
load('/imaging/ek01/mindmaps/scripts/MEG_source_vertices_MNI.mat');
vertices = int16(vertices);
k=1;
for mask = 1:length(masknic)
    
    %disp(['Working on mask: ' masknic{mask}]);
    if usemask
        
        if reslice==0
            P = ([maskdir  maskname{mask}]);
            [M XYZ] = spm_read_vols(spm_vol(P));
            XYZ = XYZ(:,find(M))'; clear P M V
        elseif reslice==1
            nii_fname=[maskdir  maskname{mask}];
            reslice_nii(nii_fname,[maskdir 'r' num2str(newvoxel(1)) 'mm_' maskname{mask}], newvoxel);
            [M XYZ] = spm_read_vols(spm_vol([maskdir 'r' num2str(newvoxel(1)) 'mm_' maskname{mask}]));
            XYZ_all = XYZ(:,find(M))'; clear P M XYZ V
        end
        
        if homo
            x = XYZ_all(:,1).*-1;
            XYZ_all(:,1) = x; clear x
        end
        % Get mask coords for vertices only
        member = ismember(vertices,XYZ, 'rows');
        nchan=length(double(vertices(member,:)));
        maskind(k:(k+nchan-1),1)=mask;
        XYZ_all = [XYZ_all; double(vertices(member,:))]; 
        k=k+nchan; 
    elseif virtualelectrode
        XYZ = virtual_el_coor;
    elseif wholebrain
        load('/imaging/ek01/mindmaps/scripts/MEG_source_vertices_MNI.mat');
        vertices = double(vertices);
        XYZ=vertices;
    end
    
end
% 
% for sub = 1:length(subs)
%     disp(['Extracting data from subject ' subs{sub} ' ...']);
%     theData = [datadir subs{sub} midfol subs{sub} front lastbit];
%     
%     D = spm_eeg_load(theData);
%     S = []; 
%     S.D = D;
%     S.timewin = [-100 0];
%     S.save = 1;
%     S.prefix='';
%     D = spm_eeg_bc(S);
%     
% end

%%
for sub = 1:length(subs)
    disp(['Extracting data from subject ' subs{sub} ' ...']);
    theData = [datadir subs{sub} midfol subs{sub} front lastbit];
    
    D = spm_eeg_load(theData);
    
    D.val = val;
    D.inv{val}.source.fname = [ outdir 'source_r0_neg_N100_' subs{sub} '_' session '_' front lastbit];
    D.inv{val}.source.type  = 'trials'; %'evoked'
    D.inv{val}.source.rad = 0;
    D.inv{val}.source.XYZ = XYZ;
    
    numClusters = size(XYZ,1);
    for r = 1:numClusters
        D.inv{val}.source.label{r} = masknic{r}; %sprintf('ROI %d %d %d chan%d',XYZ(r,1), XYZ(r,2), XYZ(r,3), r);  % Insert your names if you want
    end
    
    [Ds, D] = spm_eeg_inv_extract(D);
    
%     for mask=1:length(masknic)
%         
%         peaktime=85; % around 70 ms
%         ver_ind=find(maskind==mask);
%         
%         
%         ver_means=abs(nanmean(Ds(ver_ind,(peaktime-20):(peaktime+20),1),2));
%         peak_ver=find(ver_means==max(ver_means));
%         peak_ver=peak_ver(1);
%         peaks(mask)=ver_ind(peak_ver);
%         
%         if abs(min(Ds(ver_ind(peak_ver),(peaktime-20):(peaktime+20),1)))>abs(max(Ds(ver_ind(peak_ver),(peaktime-20):(peaktime+20),1)))
%             c(mask)=-1;
%         else
%             c(mask)=1;
%         end  
%     end
%     
%     clear D
%     D = spm_eeg_load(theData);
%     
%     D.val = val;
%     D.inv{val}.source.fname = [ outdir 'source_' subs{sub} front lastbit];
%     D.inv{val}.source.type  = 'evoked'; %'trials'
%     D.inv{val}.source.rad = 0;
%     D.inv{val}.source.XYZ = XYZ_all(peaks,:);
%     
%     for r = 1:length(masknic)
%         D.inv{val}.source.label{r} = masknic{r}; %sprintf('ROI %d %d %d chan%d',XYZ(r,1), XYZ(r,2), XYZ(r,3), r);  % Insert your names if you want
%     end
%     
%     [Ds, D] = spm_eeg_inv_extract(D);
%       
  %  for sub=1:length(subs)
   %Ds = spm_eeg_load([outdir 'source_r0_neg_N100_' subs{sub} '_' session '_' front lastbit]);
    
    for mask=1:length(masknic)
        for tr=1:size(Ds,3)
            if  Ds(mask,75,tr)>0%abs(min(Ds(mask,65:85,tr)))>abs(max(Ds(mask,65:85,tr)))
                Ds(mask,:,tr)=-Ds(mask,:,tr);
            end
        end
    end
    
    D=Ds;
   save([outdir 'source_r0_neg_N100_' subs{sub} '_' session '_' front lastbit],'D','-v7.3')
    %end
end


for sub=1:length(subs)
    
    D = spm_eeg_load([outdir 'source_r0_neg_N100_' subs{sub} '_' session '_'  front lastbit]);
    S = []; 
    S.D = D;
    S.timewin = [-100 0];
    S.save = 1;
    S.prefix='';
    D = spm_eeg_bc(S);
    
end
