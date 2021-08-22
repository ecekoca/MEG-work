%% Grand average & plots (Ece K, 2020)

% Calculates the grand average of the averaged data even if the channel
% numbers are different across subjects. Then uses the average to plot
% topoplots.

addpath('/home/ek01/scripts/ntad/')
addpath('/home/ek01/scripts/ntad/functions/')
addpath('/home/ek01/scripts/all_scripts/')
clear all; clc;
indir='/imaging/projects/cbu/ntad/meg_data/';
outdir='/imaging/rowe/users/ek01/prv/ntad/mmn_sensor/roving_mmn_normalised_signal/';
subjects={'C1001','C1002','C1004','C1006',...
    'C1010','C1011','C1012','C1014','C1015','C1017','C2001','C2002','C2008',...
'c2','c5','c7','c10','c12','c13','c14',};
task='mmn';
session='BL';
group='Controls';
times=[50 75 100 125 150 175 200 225 250];% time points to plot, time INDICES not milliseconds, e.g. 50=0ms, 100=100ms, 150=200ms
trials={'DEV','REP5'}; % trials to plot
sensors={'MEGPLANAR'}; %'MEGMAG','EEG', sensors to plot

% subjects = {'P1001','P1005','P1007','P1008','P1009',...
% 'P1010','P1011','P1012','P1015','P1020','P1021','P1022',...
% 'P1024','P1026','P1027','P1029','P1030','P1031','P1032',...
% 'P1036','P1037','P1038','P1042','P1043','P1048','P1049',...
% 'P1053','P1055','P1056','P1058','P1060','P1061','P1062','P1063',...
% 'P1064','P1065','P1066','P1067'}; %MMN

chans={'MEG0111';'MEG0121';'MEG0131';'MEG0141';'MEG0211';'MEG0221';'MEG0231';'MEG0241';'MEG0311';'MEG0321';'MEG0331';'MEG0341';'MEG0411';'MEG0421';'MEG0431';'MEG0441';'MEG0511';'MEG0521';'MEG0531';'MEG0541';'MEG0611';'MEG0621';'MEG0631';'MEG0641';'MEG0711';'MEG0721';'MEG0731';'MEG0741';'MEG0811';'MEG0821';'MEG0911';'MEG0921';'MEG0931';'MEG0941';'MEG1011';'MEG1021';'MEG1031';'MEG1041';'MEG1111';'MEG1121';'MEG1131';'MEG1141';'MEG1211';'MEG1221';'MEG1231';'MEG1241';'MEG1311';'MEG1321';'MEG1331';'MEG1341';'MEG1411';'MEG1421';'MEG1431';'MEG1441';'MEG1511';'MEG1521';'MEG1531';'MEG1541';'MEG1611';'MEG1621';'MEG1631';'MEG1641';'MEG1711';'MEG1721';'MEG1731';'MEG1741';'MEG1811';'MEG1821';'MEG1831';'MEG1841';'MEG1911';'MEG1921';'MEG1931';'MEG1941';'MEG2011';'MEG2021';'MEG2031';'MEG2041';'MEG2111';'MEG2121';'MEG2131';'MEG2141';'MEG2211';'MEG2221';'MEG2231';'MEG2241';'MEG2311';'MEG2321';'MEG2331';'MEG2341';'MEG2411';'MEG2421';'MEG2431';'MEG2441';'MEG2511';'MEG2521';'MEG2531';'MEG2541';'MEG2611';'MEG2621';'MEG2631';'MEG2641';'MEG0113';'MEG0112';'MEG0122';'MEG0123';'MEG0132';'MEG0133';'MEG0143';'MEG0142';'MEG0213';'MEG0212';'MEG0222';'MEG0223';'MEG0232';'MEG0233';'MEG0243';'MEG0242';'MEG0313';'MEG0312';'MEG0322';'MEG0323';'MEG0333';'MEG0332';'MEG0343';'MEG0342';'MEG0413';'MEG0412';'MEG0422';'MEG0423';'MEG0432';'MEG0433';'MEG0443';'MEG0442';'MEG0513';'MEG0512';'MEG0523';'MEG0522';'MEG0532';'MEG0533';'MEG0542';'MEG0543';'MEG0613';'MEG0612';'MEG0622';'MEG0623';'MEG0633';'MEG0632';'MEG0642';'MEG0643';'MEG0713';'MEG0712';'MEG0723';'MEG0722';'MEG0733';'MEG0732';'MEG0743';'MEG0742';'MEG0813';'MEG0812';'MEG0822';'MEG0823';'MEG0913';'MEG0912';'MEG0923';'MEG0922';'MEG0932';'MEG0933';'MEG0942';'MEG0943';'MEG1013';'MEG1012';'MEG1023';'MEG1022';'MEG1032';'MEG1033';'MEG1043';'MEG1042';'MEG1112';'MEG1113';'MEG1123';'MEG1122';'MEG1133';'MEG1132';'MEG1142';'MEG1143';'MEG1213';'MEG1212';'MEG1223';'MEG1222';'MEG1232';'MEG1233';'MEG1243';'MEG1242';'MEG1312';'MEG1313';'MEG1323';'MEG1322';'MEG1333';'MEG1332';'MEG1342';'MEG1343';'MEG1412';'MEG1413';'MEG1423';'MEG1422';'MEG1433';'MEG1432';'MEG1442';'MEG1443';'MEG1512';'MEG1513';'MEG1522';'MEG1523';'MEG1533';'MEG1532';'MEG1543';'MEG1542';'MEG1613';'MEG1612';'MEG1622';'MEG1623';'MEG1632';'MEG1633';'MEG1643';'MEG1642';'MEG1713';'MEG1712';'MEG1722';'MEG1723';'MEG1732';'MEG1733';'MEG1743';'MEG1742';'MEG1813';'MEG1812';'MEG1822';'MEG1823';'MEG1832';'MEG1833';'MEG1843';'MEG1842';'MEG1912';'MEG1913';'MEG1923';'MEG1922';'MEG1932';'MEG1933';'MEG1943';'MEG1942';'MEG2013';'MEG2012';'MEG2023';'MEG2022';'MEG2032';'MEG2033';'MEG2042';'MEG2043';'MEG2113';'MEG2112';'MEG2122';'MEG2123';'MEG2133';'MEG2132';'MEG2143';'MEG2142';'MEG2212';'MEG2213';'MEG2223';'MEG2222';'MEG2233';'MEG2232';'MEG2242';'MEG2243';'MEG2312';'MEG2313';'MEG2323';'MEG2322';'MEG2332';'MEG2333';'MEG2343';'MEG2342';'MEG2412';'MEG2413';'MEG2423';'MEG2422';'MEG2433';'MEG2432';'MEG2442';'MEG2443';'MEG2512';'MEG2513';'MEG2522';'MEG2523';'MEG2533';'MEG2532';'MEG2543';'MEG2542';'MEG2612';'MEG2613';'MEG2623';'MEG2622';'MEG2633';'MEG2632';'MEG2642';'MEG2643';'EEG001';'EEG002';'EEG003';'EEG004';'EEG005';'EEG006';'EEG007';'EEG008';'EEG009';'EEG010';'EEG011';'EEG012';'EEG013';'EEG014';'EEG015';'EEG016';'EEG017';'EEG018';'EEG019';'EEG020';'EEG021';'EEG022';'EEG023';'EEG024';'EEG025';'EEG026';'EEG027';'EEG028';'EEG029';'EEG030';'EEG031';'EEG032';'EEG033';'EEG034';'EEG035';'EEG036';'EEG037';'EEG038';'EEG039';'EEG040';'EEG041';'EEG042';'EEG043';'EEG044';'EEG045';'EEG046';'EEG047';'EEG048';'EEG049';'EEG050';'EEG051';'EEG052';'EEG053';'EEG054';'EEG055';'EEG056';'EEG057';'EEG058';'EEG059';'EEG060';'EEG065';'EEG066';'EEG067';'EEG068';'EEG069';'EEG070';'EEG071';'EEG072';'EEG073';'EEG074'};

%%=========================================================================
%% Grand average
%%=========================================================================
% outname=[group '_' session '_' task '_grand_mean.mat'];
% 
% cd(indir)
% S=[];warning off
% for s=1:length(subjects)
%     if contains(subjects{s},'c')
%         D=spm_eeg_load([indir '/TGB_controls/' subjects{s} '/' session '/' task '/pfmeMaffffdtransdef_all_trl.mat']);
%     else
%         S.D{s,1}=[indir subjects{s} '/' session '/' task '/pfmeMaffffdtransdef_all_trl.mat'];
%     end
% end
% S.outfile= [outdir outname];
% S.weighted=0;
% D=spm_eeg_grandmean_ek(S); % unmarks all bad channels
% 
% D = spm_eeg_load([outdir outname]); %EK
% D=D.badchannels(1:D.nchannels,0);
% D2=osl_detect_artefacts(D);D=D2;
% save([outdir outname],'D','-v7.3') % mark bad channels again so that they won't be included in the plots

%%=========================================================================
%% Manual grand average for data with different sensor orders
%%=========================================================================
D=spm_eeg_load([indir subjects{1} '/' session '/' task '/' subjects{1} '_' session '_' task '_all_trl.mat']);
cn=find(ismember(D.conditions,trials));
cd(indir); 
k=1;
for i=1:length(trials)
    for s=1:length(subjects)
        
        if contains(subjects{s},'c')
            D=spm_eeg_load([indir '/TGB_controls/' subjects{s} '/' session '/' task '/' subjects{s} '_' session '_' task '_all_trl.mat']);
        else
            D=spm_eeg_load([indir subjects{s} '/' session '/' task '/' subjects{s} '_' session '_' task '_all_trl.mat']);
        end
        
        data{1}(s,:,:)=squeeze(D(D.selectchannels(chans(1:102)),:,cn(i))); %MAG
        data{2}(s,:,:)=squeeze(D(D.selectchannels(chans(103:306)),:,cn(i))); %GRAD
        
%         if  size(D,1)>330
%             data{3}(k,:,:)=squeeze(D(D.selectchannels(chans(307:376)),:,cn(i))); %EEG
%             k=k+1;
%         end
    end
    
    means{1}(:,:,i)=squeeze(nanmean(data{1}(:,:,:),1)); %MAG
    means{2}(:,:,i)=squeeze(nanmean(data{2}(:,:,:),1)); %GRAD
%     means{3}(:,:,i)=squeeze(nanmean(data{3}(:,:,:),1)); %EEG
end

%%=========================================================================
%% Topoplots
%%=========================================================================

mkdir([outdir '/topoplots/'])
load('/home/ek01/scripts/ntad/functions/topoplot_coor2D.mat') % 2D coordinates of the sensors, better than what's in the grand mean data

for m=1%length(sensors)
    
    %idx=find(strcmp(D.chantype,sensors{m}));
    in.noButtons=1; % don't put buttons in the image
    in.type=sensors{m}; % what is the type of sensor
    in.plotpos=0; % don't plot channel positions
    
    if strcmp(sensors{m},'MEGMAG')
        coor=megcoor(:,1:102); clims=[-200 200]; cmap='redblue';
        idx=1:102;
    elseif strcmp(sensors{m},'MEGPLANAR')
         coor=megcoor(:,103:306);clims=[-20 20];cmap='jet';
         idx=103:306;
         cmap=customcolormap([0 0.1 0.2 0.3 0.4 0.5 1],...
            [1 1 1;
            1 1 0;
            1 0.5 0;
            1 0 0;
            0.9 0 0; 0.8 0.8 0.8
            0.8 0.8 0.8]); cmap(1:30,:)=repmat([0.8 0.8 0.8],30,1); 
        %clims=[-0.5 1.5];clims2=[-0.1 1];
        clims=[-1.5 2.5];clims2=[-1 1.5];
    else
        coor=eegcoor;clims=[-15 15];cmap='redblue';
        idx=306:376;
    end
    
    for i=1:length(trials)
        %ix=find(ismember(D.conditions,trials{i}));
        for t=times
            d=squeeze(means{2}(:,t,i));
            [ZI,f] =spm_eeg_plotScalpData_ek2(d,coor,chans(idx),in);
            close(gcf);
            figure;pcolor(ZI);shading flat;colormap(cmap);axis square;axis off
            caxis(clims);print(gcf,'-dbmp', '-r0',...
                [outdir '/topoplots/' task '_' session '_' trials{i} '_' sensors{m} '_t'...
                num2str(t) '_' group '.bmp'],'-r200'); close(gcf)
        end
    end
    
    for t=times
        d=squeeze(means{2}(:,t,1))-squeeze(means{2}(:,t,2));
        [ZI,f] =spm_eeg_plotScalpData_ek2(d,coor,D.chanlabels(idx),in);
        close(gcf);
        figure;pcolor(ZI);shading flat;colormap(cmap);axis square;axis off
        caxis(clims2);print(gcf,'-dbmp', '-r0',...
            [outdir '/topoplots/' task '_' session '_REP6-DEV_' sensors{m} '_t'...
            num2str(t) '_' group '.bmp'],'-r200'); close(gcf)
    end
    
end



