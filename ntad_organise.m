%% Organise project imaging space (Ece K, 2019)

% Creates the folder structure, copies the raw data and converts the T1s.
% Copies the T1 niftis under the subject's folder, for ease of access.

clear all; clc;
mrdir='/mridata/cbu_ntad/';
megdir='/megdata/cbu/ntad/';
imgdir='/imaging/projects/cbu/ntad/meg_data/';
megmap='/imaging/projects/cbu/ntad/meg_list.txt'; % To be updated manually
mrmap='/imaging/projects/cbu/ntad/mr_list.txt'; % To be updated manually

%% ================================================================
%% Create folder structure
%% ================================================================

cd(imgdir)
[megsubs megids megses ]=textread(megmap,'%s %s %s');
[mrsubs mrids mrses ]=textread(mrmap,'%s %s %s');

for s=1:length(megsubs)
    if ~exist([imgdir megsubs{s} '/']) % create subject folder
        mkdir([imgdir megsubs{s} '/']);
    end
    ind=find(strcmp(megsubs,megsubs{s}));
    subses=megses(ind);
    for ss=1:length(subses)
       if ~exist([imgdir megsubs{s} '/' subses{ss}]) 
            mkdir([imgdir megsubs{s} '/' subses{ss}]); % create session folder
       end
    end  
end

%% ================================================================
%% Copy raw files
%% ================================================================

for s=1:length(megsubs)
    ind=find(strcmp(megsubs,megsubs{s}));
    subses=megses(ind);
    disp(['Checking data: ' megsubs{s} ' ---------------'])
    for ss=2:length(subses)
        cd([megdir megids{ind(ss)} '/'])
        l=dir([megdir megids{ind(ss)} '/']); fol=[l(3).folder '/' l(3).name '/']; cd(fol);
        tdir=[imgdir megsubs{s} '/' subses{ss} '/avtask/'];
        if ~isempty(dir('av*')) && ~exist(tdir)
            file=dir('av*'); file=file.name;
            mkdir(tdir); cd(tdir); copyfile([fol file],tdir);
            disp('Copied: avtask')
        elseif isempty(dir('av*'))
            disp(['Cant find avtask, skipping'])
        elseif exist(tdir)
            disp(['Avtask already exists'])
        end
        tdir=[imgdir megsubs{s} '/' subses{ss} '/scenerep/'];cd(fol);
        if ~isempty(dir('sc*')) && ~exist(tdir)
            file=dir('sc*'); file=file.name;
            mkdir(tdir); cd(tdir); copyfile([fol file],tdir);
            disp('Copied: scenerep')
        elseif isempty(dir('sc*'))
            disp(['Cant find scenerep, skipping'])
        elseif exist(tdir)
            disp(['Scenerep already exists'])    
        end
        tdir=[imgdir megsubs{s} '/' subses{ss} '/mmn/'];cd(fol);
        if ~isempty(dir('mmn*')) && ~exist(tdir)
            file=dir('mmn*'); file=file.name;
            mkdir(tdir); cd(tdir); copyfile([fol file],tdir);
            disp('Copied: mmn')        
        elseif isempty(dir('mmn*'))
            disp(['Cant find mmn, skipping'])
        elseif exist(tdir)
            disp(['mmn already exists'])     
        end
        tdir=[imgdir megsubs{s} '/' subses{ss} '/vab/'];cd(fol);
        if ~isempty(dir('v*_ta*')) && ~exist(tdir)
            file=dir('v*_ta*'); file=file.name;
            mkdir(tdir); cd(tdir); copyfile([fol file],tdir); cd(fol)
            disp('Copied: vabtask')
        elseif isempty(dir('v*_ta*'))
            disp(['Cant find vabtask, skipping'])
        elseif exist(tdir)
            disp(['vabtask already exists'])      
        end
        tdir=[imgdir megsubs{s} '/' subses{ss} '/vab/'];cd(fol);
        if ~isempty(dir('v*_tr*')) && isempty(dir([tdir '/v*_tr*']))
            file=dir('v*_tr*'); file=file.name;cd(tdir);
            copyfile([fol file],tdir);
            disp('Copied: vabtrain')
        elseif isempty(dir('v*_tr*'))
            disp(['Cant find vabtrain, skipping'])
        elseif exist(tdir)
            disp(['vabtrain already exists'])      
        end
        tdir=[imgdir megsubs{s} '/' subses{ss} '/rest_open/'];cd(fol);
        if ~isempty(dir('*op*')) && ~exist(tdir)
            file=dir('*op*'); file=file.name;
            mkdir(tdir); cd(tdir); copyfile([fol file],tdir);
            disp('Copied: resting eyes open')
        elseif isempty(dir('*op*'))
            disp(['Cant find resting state eyes open, skipping'])
        elseif exist(tdir)
            disp(['Resting state eyes open already exists'])          
        end
        tdir=[imgdir megsubs{s} '/' subses{ss} '/rest_closed/'];cd(fol);
        if ~isempty(dir('*clo*')) && ~exist(tdir)
            file=dir('*clo*'); file=file.name;
            mkdir(tdir); cd(tdir); copyfile([fol file],tdir);
            disp('Copied: resting eyes closed')
        elseif isempty(dir('*clo*'))
            disp(['Cant find resting state eyes closed, skipping'])
        elseif exist(tdir)
            disp(['Resting state eyes closed already exists'])          
        end
    end
end

%% ================================================================
%% Convert dicoms
%% ================================================================

cd(mrdir)
for s=1:length(mrsubs)
    disp(['Checking MR data: ' mrsubs{s} ' ---------------'])
    t1fol=[imgdir mrsubs{s} '/sMRI/'];
    if ~exist(t1fol) % create T1 folder
        mkdir(t1fol);
        l=dir([mrdir mrids{s} '*']);folname=[mrdir l.name]; cd(folname)
        l=dir([folname '/20*']); folname=[l.folder '/' l.name '/']; cd(folname)
        l=dir([folname '*MPRAGE*_iso']);
        if length(l)>1
            folname=[l(end).folder '/' l(end).name '/'];
        else
            folname=[l.folder '/' l.name '/'];
        end
        cd(folname);l=dir(folname);
        filenames=[]; flen=[];
        for f=1:length(l); flen=[flen, length([l(f).folder '/' l(f).name])]; end
        ns=(ones(1,length(l)).*max(flen))-flen;
        for f=3:length(l)
            filenames=[filenames; [l(f).folder '/' l(f).name blanks(ns(f))]];
        end; clear f
        hdr=spm_dicom_headers(filenames,0);
        spm_dicom_convert(hdr,'all','flat','nii',t1fol); clear hdr filenames l folname
        disp(['Converted the MR'])   
    else
        disp(['MR already exists']) 
    end
end
