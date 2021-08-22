%% NTAD Sensor level RMS analysis for Roving (Ece K, 2019)

addpath /imaging/local/meg_misc/
addpath /imaging/local/software/spm_cbu_svn/releases/spm12_latest/;
addpath /hpc-software/matlab/cbu/;
addpath(genpath('/imaging/local/software/spm_cbu_svn/releases/spm12_latest/external/fieldtrip/')); %EK
addpath /imaging/ek01/ntad/scripts/

close all
clear all; clc
subs = {'C1001','C1002','C1003','C1004','C1006','C1007','C1008',...
'C1010','C1011','C1012','C1013','C1014','C1015','C1016','C1017'...
'P1001','P1002','P1004','P1005','P1007','P1008','P1009','P1010',....
'P1011','P1012','P1015','P1016','P1021','P1022','P1023',...
'P1024','P1026','P1027','P1029','P1030','P1031','P1032','P1035',...
'P1036','P1037','P1042','P1043','P1045','P1048','P1049',...
'P1053','P1054','P1055','P1056','P1058','P1060','P1061','P1062',...
'P1063','P1065','P1066','P1067','P1070','P3002','P3004',...
'P3005','P3006'}; %all subjects-SCENEREP exclude 1020 1064

% subjects ={'P1002', 'P1004','P1008','P1009','P1011','P1022','P1035',...
 %   'P1038','P3002','P3004','P1037','P1049','P1055','P1060'}; %'P1027' patients with two-week data
 
datadir = '/imaging/projects/cbu/ntad/meg_data/';
outdir ='/imaging/projects/cbu/ntad/analyses/scenerep_sensor_TFR/'; %'/imaging/ek01/ntad/scenerep_sensor/scenerep_normalised_signal/';
midfol='/BL/scenerep/';
lastbit ='_scenerep_tfr_new_lbl.mat';
session='BL';
testname='TFR_sensor';
permnum=500;
age=[70.9863013700000;72.7945205500000;53.6438356200000;66.1726027400000;77.1232876700000;64.5616438400000;59.0849315100000;70.0767123300000;67.0712328800000;60.5506849300000;59.7013698600000;64.0684931500000;63.9671232900000;54.9561643800000;73.8547945200000;75.1287671200000;82.3753424700000;84.8657534200000;74.4986301400000;77.5095890400000;72.3452054800000;79.5232876700000;74.1890411000000;72.8246575300000;78.5068493200000;74.9123287700000;72.7753424700000;64.4000000000000;61.2821917800000;80.6520547900000;57.9260274000000;76.0219178100000;69.7616438400000;62.7150684900000;79.9123287700000;69.6684931500000;74.4136986300000;81.4767123300000;51.8082191800000;68.5397260300000;72.4082191800000;77.5890411000000;63.7945205500000;70.2465753400000;73.7726027400000;70.5780821900000;80.6191780800000;74.9232876700000;68.2876712300000;73.2465753400000;63.8027397300000;57.8109589000000;73.1863013700000;74.8000000000000;73.6054794500000;78.8630137000000;76.6986301400000;82.3150684900000;81.6602739700000;81.4191780800000;81.6356164400000;82.7726027400000];

xticks=[1 50  150 250 350 451];
xlabels={'-100' '0' '200'  '400' '600' };
conditions={'NEW', 'OLD-S','OLD-L'};

%% Channel selection

fr_mag_L={'MEG0121','MEG0341','MEG0321','MEG0331','MEG0641','MEG0621','MEG0611','MEG0541','MEG0311','MEG0511','MEG0531','MEG0821','MEG0521'};
fr_mag_R={'MEG1031','MEG1241','MEG1231','MEG1221','MEG1411','MEG1011','MEG1021','MEG0931','MEG1211','MEG0921','MEG0911','MEG0811','MEG0941'};
tem_mag_L={'MEG0111','MEG0221','MEG0211','MEG0131','MEG0141','MEG1511','MEG1521','MEG1531','MEG1541','MEG0231','MEG0241','MEG1621','MEG1611'};
tem_mag_R={'MEG1421','MEG1441','MEG1321','MEG1311','MEG1431','MEG2611','MEG1331','MEG1341','MEG2621','MEG2641','MEG2421','MEG2411','MEG2631'};
par_mag_L={'MEG0411','MEG0421','MEG0631','MEG0441','MEG0431','MEG0711','MEG1811','MEG1821','MEG0741','MEG1631','MEG1841','MEG1831','MEG2011'};
par_mag_R={'MEG1041','MEG1111','MEG1121','MEG0721','MEG1141','MEG1131','MEG0731','MEG2211','MEG2221','MEG2241','MEG2231','MEG2021','MEG2441'};
occ_mag_L={'MEG1641','MEG1721','MEG1711','MEG1731','MEG1941','MEG1911','MEG2041','MEG1921','MEG1931','MEG1741','MEG2141','MEG2111'};
occ_mag_R={'MEG2121','MEG2131','MEG2031','MEG2341','MEG2331','MEG2541','MEG2311','MEG2321','MEG2511','MEG2531','MEG2521','MEG2431'};

fr_grad_L={'MEG0123+0122','MEG0342+0343','MEG0323+0322','MEG0332+0333','MEG0643+0642','MEG0623+0622','MEG0312+0313','MEG0543+0542','MEG0612+0613','MEG0512+0513','MEG0533+0532','MEG0823+0822','MEG0522+0523'};
fr_grad_R={'MEG0812+0813','MEG0912+0913','MEG0922+0923','MEG0943+0942','MEG0933+0932','MEG1212+1213','MEG1012+1013','MEG1022+1023','MEG1033+1032','MEG1242+1243','MEG1233+1232','MEG1222+1223','MEG1413+1412'};
tem_grad_L={'MEG0112+0113','MEG0142+0143','MEG0133+0132','MEG0212+0213','MEG0223+0222','MEG0233+0232','MEG0242+0243','MEG1512+1513','MEG1542+1543','MEG1623+1622','MEG1612+1613','MEG1523+1522','MEG1532+1533'};
tem_grad_R={'MEG1313+1312','MEG1322+1323','MEG1443+1442','MEG1422+1423','MEG1432+1433','MEG2613+2612','MEG1332+1333','MEG1343+1342','MEG2622+2623','MEG2643+2642','MEG2423+2422','MEG2413+2412','MEG2632+2633'};
par_grad_L={'MEG0412+0413','MEG0422+0423','MEG0632+0633','MEG0442+0443','MEG0432+0433','MEG0712+0713','MEG1812+1813','MEG1822+1823','MEG0742+0743','MEG1632+1633','MEG1842+1843','MEG1832+1833','MEG2012+2013'};
par_grad_R={'MEG1042+1043','MEG1112+1113','MEG1122+1123','MEG0722+0723','MEG1142+1143','MEG1132+1133','MEG0732+0733','MEG2212+2213','MEG2222+2223','MEG2242+2243','MEG2232+2233','MEG2022+2023'};
occ_grad_L={'MEG1642+1643','MEG1722+1723','MEG1712+1713','MEG1912+1913','MEG1942+1943','MEG1732+1733','MEG1742+1743','MEG1932+1933','MEG1922+1923','MEG2042+2043','MEG2112+2113','MEG2142+2143'};
occ_grad_R={'MEG2032+2033','MEG2122+2123','MEG2132+2133','MEG2342+2343','MEG2332+2333','MEG2542+2543','MEG2312+2313','MEG2322+2323','MEG2512+2513','MEG2532+2533','MEG2522+2523','MEG2432+2433'};

%% Gather data

c_ind=find(contains(subs,'C'));
p_ind=find(contains(subs,'P'));

cd(datadir);
mkdir(outdir);


for sub=1:length(subs)
    filen=[datadir subs{sub} midfol subs{sub} '_' session lastbit];
    D=spm_eeg_load(filen);
    
    %origdata(sub,1:size(D,1),:,:)=D(:,:,:);
%     if size(D,1)>320
%         coormeg(sub,:,1:380)=D.coor2D(1:380);
%     else
%         coormeg(sub,:,1:306)=D.coor2D;
%     end
    data{1}(sub,:,:,:,:)=D(1:102,:,:,:);
    data{2}(sub,:,:,:,:)=D(D.selectchannels(fr_mag_L),:,:,:);
    data{3}(sub,:,:,:,:)=D(D.selectchannels(fr_mag_R),:,:,:);
    data{4}(sub,:,:,:,:)=D(D.selectchannels(tem_mag_L),:,:,:);
    data{5}(sub,:,:,:,:)=D(D.selectchannels(tem_mag_R),:,:,:);
    data{6}(sub,:,:,:,:)=D(D.selectchannels(par_mag_L),:,:,:);
    data{7}(sub,:,:,:,:)=D(D.selectchannels(par_mag_R),:,:,:);
    data{8}(sub,:,:,:,:)=D(D.selectchannels(occ_mag_L),:,:,:);
    data{9}(sub,:,:,:,:)=D(D.selectchannels(occ_mag_R),:,:,:);
    
    data{10}(sub,:,:,:,:)=D([103:204],:,:,:);
    data{11}(sub,:,:,:,:)=D(D.selectchannels(fr_grad_L),:,:,:);
    data{12}(sub,:,:,:,:)=D(D.selectchannels(fr_grad_R),:,:,:);
    data{13}(sub,:,:,:,:)=D(D.selectchannels(tem_grad_L),:,:,:);
    data{14}(sub,:,:,:,:)=D(D.selectchannels(tem_grad_R),:,:,:);
    data{15}(sub,:,:,:,:)=D(D.selectchannels(par_grad_L),:,:,:);
    data{16}(sub,:,:,:,:)=D(D.selectchannels(par_grad_R),:,:,:);
    data{17}(sub,:,:,:,:)=D(D.selectchannels(occ_grad_L),:,:,:);
    data{18}(sub,:,:,:,:)=D(D.selectchannels(occ_grad_R),:,:,:);
    
%     if size(D,1)>350
%         data{11}(sub,:,:,:)=D([306:376],:,:);
%         data{12}(sub,:,:,:)=D(D.selectchannels(fr_eeg),:,:);
% %         data{9}(sub,:,:,:)=D(D.selectchannels(tem_eeg),:,:);
%     else
%         data{11}(sub,:,:,:)=NaN;
%         data{12}(sub,:,:,:)=NaN;
% %         data{9}(sub,:,:,:)=NaN;
%     end
end

% for i=1:18
%     for s=1:size(data{i},1)
%         for tr=1:size(data{i},4)
%             for ch=1:size(data{i},2)
%                 if data{i}(s,ch,100,tr)<0
%                     data{i}(s,ch,:,tr)=-(data{i}(s,ch,:,tr)); % P100 is positive
%                     %data{i}(s,ch,:,tr)=NaN(1,1,501);
%                 end
%             end
%         end
%     end
% end

for i=1:18
    %datarms{i}(:,:,:,:)=permute(data{i}(:,:,:,:),[2 1 4 3]);
    %datarms{i}=squeeze(rms_ek(datarms{i},1,[1 100]));
    d=squeeze(nanmean(data{i}(:,:,:,:,:),2));
    d=permute(d,[1 4 2 3]);
    datarms{i}=d;
    
    for j=1:size(datarms{i},1)
        for k=1:size(datarms{i},2)
            d=squeeze(datarms{i}(j,k,:,:));
            d=smooth2a(d,2,10);
%             d=squeeze(datarms{i}(j,k,:,:));
%             bl=nanmean(d(1:51)); d=d-bl; % baseline correct again just in case
            datarms{i}(j,k,:,:)=d;
        end
    end
    
%     ind=find(isnan(datarms{i}));
%     for s=1:size(datarms{i},1)
%         datarms{i}(s,:,:)=mat2gray(datarms{i}(s,:,:));
%     end
%     
%     datarms{i}(ind)=NaN;
%     datarms{i}=datarms{i}(:,:,1:351); % up to 600 ms
%     
end

for i=1:18
    d=datarms{i};
    d2(:,1,:,:)=d(:,1,:,:)-d(:,2,:,:); %NEW-OLD
    d2(:,2,:,:)=d(:,1,:,:)-d(:,3,:,:); %NEW-OLD-S
    d2(:,3,:,:)=d(:,1,:,:)-d(:,4,:,:); %NEW-OLD-L
    d=d2;
    
    for j=1:size(d,1)
        for k=1:size(d,2)
            ts=squeeze(d(j,k,:,:));
%             bl=nanmean(ts(1:51)); ts=ts-bl; % baseline correct again just in case
            diffrms{i}(j,k,:,:)=smooth2a(ts,2,5);
        end
    end
    
end

%% Plot TFRs
clear data
varnm={'MAGs','L frontal MAGs','R frontal MAGs','L temporal MAGs','R temporal MAGs','L central MAGs','R central MAGs',...
    'L occipital MAGs','R occipital MAGs','GRADs','L frontal GRADs','R frontal GRADs',...
    'L temporal GRADs','R temporal GRADs','L central GRADs','R central GRADs','L occipital GRADs','R occipital GRADs'};

for i=[10 12 14] %1:length(varnm)
    
    data=datarms{i}(c_ind,:,:,:);
    m1=squeeze(nanmean(data,1));
    data=diffrms{i}(c_ind,:,:,:);
    m2=squeeze(nanmean(data,1));
    
    figure('color','w') 
    imagesc(flip(squeeze(m1(1,:,:)))); colorbar
    caxis([0 15]); colormap('parula')
    set(gca, 'XTick',xticks, 'XTickLabel',xlabels,'Fontsize',14); 
    set(gca, 'YTick',[1 55 70 88 92 96 100], 'YTickLabel',{'100','45','30','12', '8', '4', '0'},'Fontsize',14); 
    xlabel('Time (ms)'); ylabel('Frequency (Hz)'); ylim([55 100])
    title(['NEW controls mean in ' varnm{i}])
    print(gcf,'-dbmp', '-r0',[outdir  strrep(varnm{i},' ','_') '_' session  '_mean_controls_NEW_' testname '.bmp'],'-r100'); close(gcf)
    
    figure('color','w') 
    imagesc(flip(squeeze(m1(2,:,:)))); colorbar
    caxis([0 15]); colormap('parula')
    set(gca, 'XTick',xticks, 'XTickLabel',xlabels,'Fontsize',14); 
    set(gca, 'YTick',[1 55 70 88 92 96 100], 'YTickLabel',{'100','45','30','12', '8', '4', '0'},'Fontsize',14); 
    xlabel('Time (ms)'); ylabel('Frequency (Hz)'); ylim([55 100])
    title(['OLD controls mean in ' varnm{i}])
    print(gcf,'-dbmp', '-r0',[outdir  strrep(varnm{i},' ','_') '_' session  '_mean_controls_OLD_' testname '.bmp'],'-r100'); close(gcf)
    
    figure('color','w') 
    imagesc(flip(squeeze(m1(3,:,:)))); colorbar
    caxis([0 15]); colormap('parula')
    set(gca, 'XTick',xticks, 'XTickLabel',xlabels,'Fontsize',14); 
    set(gca, 'YTick',[1 55 70 88 92 96 100], 'YTickLabel',{'100','45','30','12', '8', '4', '0'},'Fontsize',14); 
    xlabel('Time (ms)'); ylabel('Frequency (Hz)'); ylim([55 100])
    title(['OLD-S controls mean in ' varnm{i}])
    print(gcf,'-dbmp', '-r0',[outdir  strrep(varnm{i},' ','_') '_' session  '_mean_controls_OLD-S_' testname '.bmp'],'-r100'); close(gcf)
    
    figure('color','w') 
    imagesc(flip(squeeze(m1(4,:,:)))); colorbar
    caxis([0 15]); colormap('parula')
    set(gca, 'XTick',xticks, 'XTickLabel',xlabels,'Fontsize',14); 
    set(gca, 'YTick',[1 55 70 88 92 96 100], 'YTickLabel',{'100','45','30','12', '8', '4', '0'},'Fontsize',14); 
    xlabel('Time (ms)'); ylabel('Frequency (Hz)'); ylim([55 100])
    title(['OLD-L controls mean in ' varnm{i}])
    print(gcf,'-dbmp', '-r0',[outdir  strrep(varnm{i},' ','_') '_' session  '_mean_controls_OLD-L_' testname '.bmp'],'-r100'); close(gcf)
    
    figure('color','w') 
    imagesc(flip(squeeze(m2(1,:,:)))); colorbar
    caxis([-3 3]); colormap('parula')
    set(gca, 'XTick',xticks, 'XTickLabel',xlabels,'Fontsize',14); 
    set(gca, 'YTick',[1 55 70 88 92 96 100], 'YTickLabel',{'100','45','30','12', '8', '4', '0'},'Fontsize',14); 
    xlabel('Time (ms)'); ylabel('Frequency (Hz)'); ylim([55 100])
    title(['NEW-OLD controls mean in ' varnm{i}])
    print(gcf,'-dbmp', '-r0',[outdir  strrep(varnm{i},' ','_') '_' session  '_mean_controls_NEW-OLD_' testname '.bmp'],'-r100'); close(gcf)
    
    figure('color','w') 
    imagesc(flip(squeeze(m2(2,:,:)))); colorbar
    caxis([-3 3]); colormap('parula')
    set(gca, 'XTick',xticks, 'XTickLabel',xlabels,'Fontsize',14); 
    set(gca, 'YTick',[1 55 70 88 92 96 100], 'YTickLabel',{'100','45','30','12', '8', '4', '0'},'Fontsize',14); 
    xlabel('Time (ms)'); ylabel('Frequency (Hz)'); ylim([55 100])
    title(['NEW-OLD-S controls mean in ' varnm{i}])
    print(gcf,'-dbmp', '-r0',[outdir  strrep(varnm{i},' ','_') '_' session  '_mean_controls_NEW-OLD-S_' testname '.bmp'],'-r100'); close(gcf)
    
    figure('color','w') 
    imagesc(flip(squeeze(m2(3,:,:)))); colorbar
    caxis([-3 3]); colormap('parula')
    set(gca, 'XTick',xticks, 'XTickLabel',xlabels,'Fontsize',14); 
    set(gca, 'YTick',[1 55 70 88 92 96 100], 'YTickLabel',{'100','45','30','12', '8', '4', '0'},'Fontsize',14); 
    xlabel('Time (ms)'); ylabel('Frequency (Hz)'); ylim([55 100])
    title(['NEW-OLD-L controls mean in ' varnm{i}])
    print(gcf,'-dbmp', '-r0',[outdir  strrep(varnm{i},' ','_') '_' session  '_mean_controls_NEW-OLD-L_' testname '.bmp'],'-r100'); close(gcf)

end


for i=[10 12 14]% 1:length(varnm)
    
    data=datarms{i}(p_ind,:,:,:);
    m1=squeeze(nanmean(data,1));
    data=diffrms{i}(p_ind,:,:,:);
    m2=squeeze(nanmean(data,1));
    
    figure('color','w') 
    imagesc(flip(squeeze(m1(1,:,:)))); colorbar
    caxis([0 15]); colormap('parula')
    set(gca, 'XTick',xticks, 'XTickLabel',xlabels,'Fontsize',14); 
    set(gca, 'YTick',[1 55 70 88 92 96 100], 'YTickLabel',{'100','45','30','12', '8', '4', '0'},'Fontsize',14); 
    xlabel('Time (ms)'); ylabel('Frequency (Hz)'); ylim([55 100])
    title(['NEW patients mean in ' varnm{i}])
    print(gcf,'-dbmp', '-r0',[outdir  strrep(varnm{i},' ','_') '_' session  '_mean_patients_NEW_' testname '.bmp'],'-r100'); close(gcf)
    
    figure('color','w') 
    imagesc(flip(squeeze(m1(2,:,:)))); colorbar
    caxis([0 15]); colormap('parula')
    set(gca, 'XTick',xticks, 'XTickLabel',xlabels,'Fontsize',14); 
    set(gca, 'YTick',[1 55 70 88 92 96 100], 'YTickLabel',{'100','45','30','12', '8', '4', '0'},'Fontsize',14); 
    xlabel('Time (ms)'); ylabel('Frequency (Hz)'); ylim([55 100])
    title(['OLD patients mean in ' varnm{i}])
    print(gcf,'-dbmp', '-r0',[outdir  strrep(varnm{i},' ','_') '_' session  '_mean_patients_OLD_' testname '.bmp'],'-r100'); close(gcf)
    
    figure('color','w') 
    imagesc(flip(squeeze(m1(3,:,:)))); colorbar
    caxis([0 15]); colormap('parula')
    set(gca, 'XTick',xticks, 'XTickLabel',xlabels,'Fontsize',14); 
    set(gca, 'YTick',[1 55 70 88 92 96 100], 'YTickLabel',{'100','45','30','12', '8', '4', '0'},'Fontsize',14); 
    xlabel('Time (ms)'); ylabel('Frequency (Hz)'); ylim([55 100])
    title(['OLD-S patients mean in ' varnm{i}])
    print(gcf,'-dbmp', '-r0',[outdir  strrep(varnm{i},' ','_') '_' session  '_mean_patients_OLD-S_' testname '.bmp'],'-r100'); close(gcf)
    
    figure('color','w') 
    imagesc(flip(squeeze(m1(4,:,:)))); colorbar
    caxis([0 15]); colormap('parula')
    set(gca, 'XTick',xticks, 'XTickLabel',xlabels,'Fontsize',14); 
    set(gca, 'YTick',[1 55 70 88 92 96 100], 'YTickLabel',{'100','45','30','12', '8', '4', '0'},'Fontsize',14); 
    xlabel('Time (ms)'); ylabel('Frequency (Hz)'); ylim([55 100])
    title(['OLD-L patients mean in ' varnm{i}])
    print(gcf,'-dbmp', '-r0',[outdir  strrep(varnm{i},' ','_') '_' session  '_mean_patients_OLD-L_' testname '.bmp'],'-r100'); close(gcf)
    
    figure('color','w') 
    imagesc(flip(squeeze(m2(1,:,:)))); colorbar
    caxis([-3 3]); colormap('parula')
    set(gca, 'XTick',xticks, 'XTickLabel',xlabels,'Fontsize',14); 
    set(gca, 'YTick',[1 55 70 88 92 96 100], 'YTickLabel',{'100','45','30','12', '8', '4', '0'},'Fontsize',14); 
    xlabel('Time (ms)'); ylabel('Frequency (Hz)'); ylim([55 100])
    title(['NEW-OLD patients mean in ' varnm{i}])
    print(gcf,'-dbmp', '-r0',[outdir  strrep(varnm{i},' ','_') '_' session  '_mean_patients_NEW-OLD_' testname '.bmp'],'-r100'); close(gcf)
    
    figure('color','w') 
    imagesc(flip(squeeze(m2(2,:,:)))); colorbar
    caxis([-3 3]); colormap('parula')
    set(gca, 'XTick',xticks, 'XTickLabel',xlabels,'Fontsize',14); 
    set(gca, 'YTick',[1 55 70 88 92 96 100], 'YTickLabel',{'100','45','30','12', '8', '4', '0'},'Fontsize',14); 
    xlabel('Time (ms)'); ylabel('Frequency (Hz)'); ylim([55 100])
    title(['NEW-OLD-S patients mean in ' varnm{i}])
    print(gcf,'-dbmp', '-r0',[outdir  strrep(varnm{i},' ','_') '_' session  '_mean_patients_NEW-OLD-S_' testname '.bmp'],'-r100'); close(gcf)
    
    figure('color','w') 
    imagesc(flip(squeeze(m2(3,:,:)))); colorbar
    caxis([-3 3]); colormap('parula')
    set(gca, 'XTick',xticks, 'XTickLabel',xlabels,'Fontsize',14); 
    set(gca, 'YTick',[1 55 70 88 92 96 100], 'YTickLabel',{'100','45','30','12', '8', '4', '0'},'Fontsize',14); 
    xlabel('Time (ms)'); ylabel('Frequency (Hz)'); ylim([55 100])
    title(['NEW-OLD-L patients mean in ' varnm{i}])
    print(gcf,'-dbmp', '-r0',[outdir  strrep(varnm{i},' ','_') '_' session  '_mean_patients_NEW-OLD-L_' testname '.bmp'],'-r100'); close(gcf)
    
end


for i=1:18
    
    data=diffrms{i}(c_ind,:,:,:);
    m1=squeeze(nanmean(data,1));
    data=diffrms{i}(p_ind,:,:,:);
    m2=squeeze(nanmean(data,1));
    m2=m1-m2; % controls- patients
    
    figure('color','w') 
    imagesc(flip(squeeze(m2(1,:,:)))); colorbar
    caxis([-3 3]); colormap('parula')
    set(gca, 'XTick',xticks, 'XTickLabel',xlabels,'Fontsize',14); 
    set(gca, 'YTick',[1 55 70 88 92 96 100], 'YTickLabel',{'100','45','30','12', '8', '4', '0'},'Fontsize',14); 
    xlabel('Time (ms)'); ylabel('Frequency (Hz)'); ylim([55 100])
    title(['NEW-OLD interaction in ' varnm{i}])
    print(gcf,'-dbmp', '-r0',[outdir  strrep(varnm{i},' ','_') '_' session  '_Controls-patients_NEW-OLD_' testname '.bmp'],'-r100'); close(gcf)
    
    figure('color','w') 
    imagesc(flip(squeeze(m2(2,:,:)))); colorbar
    caxis([-3 3]); colormap('parula')
    set(gca, 'XTick',xticks, 'XTickLabel',xlabels,'Fontsize',14); 
    set(gca, 'YTick',[1 55 70 88 92 96 100], 'YTickLabel',{'100','45','30','12', '8', '4', '0'},'Fontsize',14); 
    xlabel('Time (ms)'); ylabel('Frequency (Hz)'); ylim([55 100])
    title(['NEW-OLD-S interaction in ' varnm{i}])
    print(gcf,'-dbmp', '-r0',[outdir  strrep(varnm{i},' ','_') '_' session  '_Controls-patients_NEW-OLD-S_' testname '.bmp'],'-r100'); close(gcf)
    
    figure('color','w') 
    imagesc(flip(squeeze(m2(3,:,:)))); colorbar
    caxis([-3 3]); colormap('parula')
    set(gca, 'XTick',xticks, 'XTickLabel',xlabels,'Fontsize',14); 
    set(gca, 'YTick',[1 55 70 88 92 96 100], 'YTickLabel',{'100','45','30','12', '8', '4', '0'},'Fontsize',14); 
    xlabel('Time (ms)'); ylabel('Frequency (Hz)'); ylim([55 100])
    title(['NEW-OLD-L interaction in ' varnm{i}])
    print(gcf,'-dbmp', '-r0',[outdir  strrep(varnm{i},' ','_') '_' session  '_Controls-patients_NEW-OLD-L_' testname '.bmp'],'-r100'); close(gcf)
    
end


%% Task effects in controls

group{1}=c_ind;group{2}=p_ind;
gnm={'Controls','Patients'};
kl=1; testn=''; roi=''; gr=''; tsum=[]; pmin=[]; time=[]; freq=[];

for i=[10 12 14] % NEW-OLD in controls
    disp(['NEW-OLD in controls and in ' varnm{i}])
    
    data1 = squeeze(datarms{i}(c_ind,1,1:45,:));
    data2 = squeeze(datarms{i}(c_ind,2,1:45,:));
    
    [results, perm_max_cl,tmap]=perm_stats_TFR(data1,data2,permnum,0.05,1);
    
    figure('color','w'); tmap=tmap.*(tmap>abs(tinv(0.05,size(data1,1))));
    imagesc(flip(tmap)); colorbar
    caxis([-5 5]); colormap('redblue')
    set(gca, 'XTick',xticks, 'XTickLabel',xlabels,'Fontsize',14); 
    set(gca, 'YTick',[1 16 33 37 41 45], ...
        'YTickLabel',{'45','30','12', '8', '4', '0'},'Fontsize',14); 
    xlabel('Time (ms)'); ylabel('Frequency (Hz)'); 
    title(['NEW-OLD tmap in controls in ' varnm{i}])
    print(gcf,'-dbmp', '-r0',[outdir  strrep(varnm{i},' ','_') '_' session  '_Controls_NEW-OLD_tmap_' testname '.bmp'],'-r100'); 
    close(gcf)
    
    
    if min(results(:,2))>0.07 % nothing significant
        testn{kl,1}='NEW-OLD task';
        roi{kl,1}=varnm{i};
        gr{kl,1}='Controls';
        tsum(kl,1)=NaN; pmin(kl,1)=NaN; time(kl,:)=[NaN NaN]; freq(kl,:)=[NaN NaN];
        kl=kl+1;
    else %significant clusters
        idx=find(results(:,2)<0.07);
        for j=1:length(idx)
            time(kl,1)=(results(idx(j),5)+51)*2-100;
            time(kl,2)=(results(idx(j),6)+51)*2-100;
            freq(kl,1)=results(idx(j),3);
            freq(kl,2)=results(idx(j),4);
            testn{kl,1}='NEW-OLD task';
            roi{kl,1}=varnm{i};
            gr{kl,1}='Controls';
            tsum(kl,1)=results(idx(j),1);
            pmin(kl,1)=results(idx(j),2);
            kl=kl+1; 
        end
    end; clear d ix cl_m tmap idx
    
end

tbl=table(testn,roi,gr,tsum,pmin,time,freq,'VariableNames',{'Test','Sensor','Group','tsum','pcor','time','freq'}); disp(tbl)
save([outdir '/scenerep_TFR_sensor_stats_controls_NEW-OLD_task_effects.mat'],'tbl');

kl=1; testn=''; roi=''; gr=''; tsum=[]; pmin=[]; time=[]; freq=[];

for i=[10 12 14] % NEW-OLD in controls
    disp(['NEW-OLD-S in controls and in ' varnm{i}])
    
    data1 = squeeze(datarms{i}(c_ind,1,1:45,:));
    data2 = squeeze(datarms{i}(c_ind,3,1:45,:));
    
    [results, perm_max_cl,tmap]=perm_stats_TFR(data1,data2,permnum,0.05,1);
    
    figure('color','w'); tmap=tmap.*(tmap>abs(tinv(0.05,size(data1,1))));
    imagesc(flip(tmap)); colorbar
    caxis([-5 5]); colormap('redblue')
    set(gca, 'XTick',xticks, 'XTickLabel',xlabels,'Fontsize',14); 
    set(gca, 'YTick',[1 16 33 37 41 45], ...
        'YTickLabel',{'45','30','12', '8', '4', '0'},'Fontsize',14); 
    xlabel('Time (ms)'); ylabel('Frequency (Hz)'); 
    title(['NEW-OLD-S tmap in controls in ' varnm{i}])
    print(gcf,'-dbmp', '-r0',[outdir  strrep(varnm{i},' ','_') '_' session  '_Controls_NEW-OLD-S_tmap_' testname '.bmp'],'-r100'); 
    %close(gcf)
    
    if min(results(:,2))>0.07 % nothing significant
        testn{kl,1}='NEW-OLD-S task';
        roi{kl,1}=varnm{i};
        gr{kl,1}='Controls';
        tsum(kl,1)=NaN; pmin(kl,1)=NaN; time(kl,:)=[NaN NaN]; freq(kl,:)=[NaN NaN];
        kl=kl+1;
    else %significant clusters
        idx=find(results(:,2)<0.07);
        for j=1:length(idx)
            time(kl,1)=(results(idx(j),5)+51)*2-100;
            time(kl,2)=(results(idx(j),6)+51)*2-100;
            freq(kl,1)=results(idx(j),3);
            freq(kl,2)=results(idx(j),4);
            testn{kl,1}='NEW-OLD-S task';
            roi{kl,1}=varnm{i};
            gr{kl,1}='Controls';
            tsum(kl,1)=results(idx(j),1);
            pmin(kl,1)=results(idx(j),2);
            kl=kl+1; 
        end
    end; clear d ix cl_m tmap idx
    
end

tbl=table(testn,roi,gr,tsum,pmin,time,freq,'VariableNames',{'Test','Sensor','Group','tsum','pcor','time','freq'}); disp(tbl)
save([outdir '/scenerep_TFR_sensor_stats_controls_NEW-OLD-S_task_effects.mat'],'tbl');

kl=1; testn=''; roi=''; gr=''; tsum=[]; pmin=[]; time=[]; freq=[];

for i=[10 12 14] % NEW-OLD in controls
    disp(['NEW-OLD-L in controls and in ' varnm{i}])
    
    data1 = squeeze(datarms{i}(c_ind,1,1:45,:));
    data2 = squeeze(datarms{i}(c_ind,4,1:45,:));
    
    [results, perm_max_cl,tmap]=perm_stats_TFR(data1,data2,permnum,0.05,1);
    
    figure('color','w'); tmap=tmap.*(tmap>abs(tinv(0.05,size(data1,1))));
    imagesc(flip(tmap)); colorbar
    caxis([-5 5]); colormap('redblue')
    set(gca, 'XTick',xticks, 'XTickLabel',xlabels,'Fontsize',14); 
    set(gca, 'YTick',[1 16 33 37 41 45], ...
        'YTickLabel',{'45','30','12', '8', '4', '0'},'Fontsize',14); 
    xlabel('Time (ms)'); ylabel('Frequency (Hz)'); 
    title(['NEW-OLD-L tmap in controls in ' varnm{i}])
    print(gcf,'-dbmp', '-r0',[outdir  strrep(varnm{i},' ','_') '_' session  '_Controls_NEW-OLD-L_tmap_' testname '.bmp'],'-r100'); 
    %close(gcf)
    
    if min(results(:,2))>0.07 % nothing significant
        testn{kl,1}='NEW-OLD-L task';
        roi{kl,1}=varnm{i};
        gr{kl,1}='Controls';
        tsum(kl,1)=NaN; pmin(kl,1)=NaN; time(kl,:)=[NaN NaN]; freq(kl,:)=[NaN NaN];
        kl=kl+1;
    else %significant clusters
        idx=find(results(:,2)<0.07);
        for j=1:length(idx)
            time(kl,1)=(results(idx(j),5)+51)*2-100;
            time(kl,2)=(results(idx(j),6)+51)*2-100;
            freq(kl,1)=results(idx(j),3);
            freq(kl,2)=results(idx(j),4);
            testn{kl,1}='NEW-OLD-L task';
            roi{kl,1}=varnm{i};
            gr{kl,1}='Controls';
            tsum(kl,1)=results(idx(j),1);
            pmin(kl,1)=results(idx(j),2);
            kl=kl+1; 
        end
    end; clear d ix cl_m tmap idx
    
end

tbl=table(testn,roi,gr,tsum,pmin,time,freq,'VariableNames',{'Test','Sensor','Group','tsum','pcor','time','freq'}); disp(tbl)
save([outdir '/scenerep_TFR_sensor_stats_controls_NEW-OLD-L_task_effects.mat'],'tbl');

%% Task effects in patients

group{1}=c_ind;group{2}=p_ind;
gnm={'Controls','Patients'};
kl=1; testn=''; roi=''; gr=''; tsum=[]; pmin=[]; time=[]; freq=[];

for i=[10 12 14] % NEW-OLD in patients
    disp(['NEW-OLD in patients and in ' varnm{i}])
    
    data1 = squeeze(datarms{i}(p_ind,1,1:45,:));
    data2 = squeeze(datarms{i}(p_ind,2,1:45,:));
     
    [results, perm_max_cl,tmap]=perm_stats_TFR(data1,data2,permnum,0.05,1);
    
    figure('color','w'); tmap=tmap.*(tmap>abs(tinv(0.05,size(data1,1))));
    imagesc(flip(tmap)); colorbar
    caxis([-5 5]); colormap('redblue')
    set(gca, 'XTick',xticks, 'XTickLabel',xlabels,'Fontsize',14); 
    set(gca, 'YTick',[1 16 33 37 41 45], ...
        'YTickLabel',{'45','30','12', '8', '4', '0'},'Fontsize',14); 
    xlabel('Time (ms)'); ylabel('Frequency (Hz)'); 
    title(['NEW-OLD tmap in patients in ' varnm{i}])
    print(gcf,'-dbmp', '-r0',[outdir  strrep(varnm{i},' ','_') '_' session  '_Patients_NEW-OLD_tmap_' testname '.bmp'],'-r100'); 
    %close(gcf)
    
    if min(results(:,2))>0.07 % nothing significant
        testn{kl,1}='NEW-OLD task';
        roi{kl,1}=varnm{i};
        gr{kl,1}='Patients';
        tsum(kl,1)=NaN; pmin(kl,1)=NaN; time(kl,:)=[NaN NaN]; freq(kl,:)=[NaN NaN];
        kl=kl+1;
    else %significant clusters
        idx=find(results(:,2)<0.07);
        for j=1:length(idx)
            time(kl,1)=(results(idx(j),5)+51)*2-100;
            time(kl,2)=(results(idx(j),6)+51)*2-100;
            freq(kl,1)=results(idx(j),3);
            freq(kl,2)=results(idx(j),4);
            testn{kl,1}='NEW-OLD task';
            roi{kl,1}=varnm{i};
            gr{kl,1}='Patients';
            tsum(kl,1)=results(idx(j),1);
            pmin(kl,1)=results(idx(j),2);
            kl=kl+1; 
        end
    end; clear d ix cl_m tmap idx
   
end

tbl=table(testn,roi,gr,tsum,pmin,time,freq,'VariableNames',{'Test','Sensor','Group','tsum','pcor','time','freq'}); disp(tbl)
save([outdir '/scenerep_TFR_sensor_stats_patients_NEW-OLD_task_effects.mat'],'tbl');

kl=1; testn=''; roi=''; gr=''; tsum=[]; pmin=[]; time=[]; freq=[];

for i=[10 12 14] % NEW-OLD in patients
    disp(['NEW-OLD-S in patients and in ' varnm{i}])
    
    data1 = squeeze(datarms{i}(p_ind,1,1:45,:));
    data2 = squeeze(datarms{i}(p_ind,3,1:45,:));
     
    [results, perm_max_cl,tmap]=perm_stats_TFR(data1,data2,permnum,0.05,1);
    
    figure('color','w'); tmap=tmap.*(tmap>abs(tinv(0.05,size(data1,1))));
    imagesc(flip(tmap)); colorbar
    caxis([-5 5]); colormap('redblue')
    set(gca, 'XTick',xticks, 'XTickLabel',xlabels,'Fontsize',14); 
    set(gca, 'YTick',[1 16 33 37 41 45], ...
        'YTickLabel',{'45','30','12', '8', '4', '0'},'Fontsize',14); 
    xlabel('Time (ms)'); ylabel('Frequency (Hz)'); 
    title(['NEW-OLD-S tmap in patients in ' varnm{i}])
    print(gcf,'-dbmp', '-r0',[outdir  strrep(varnm{i},' ','_') '_' session  '_Patients_NEW-OLD-S_tmap_' testname '.bmp'],'-r100'); 
    %close(gcf)
    
    if min(results(:,2))>0.07 % nothing significant
        testn{kl,1}='NEW-OLD-S task';
        roi{kl,1}=varnm{i};
        gr{kl,1}='Patients';
        tsum(kl,1)=NaN; pmin(kl,1)=NaN; time(kl,:)=[NaN NaN]; freq(kl,:)=[NaN NaN];
        kl=kl+1;
    else %significant clusters
        idx=find(results(:,2)<0.07);
        for j=1:length(idx)
            time(kl,1)=(results(idx(j),5))*2-100;
            time(kl,2)=(results(idx(j),6))*2-100;
            freq(kl,1)=results(idx(j),3);
            freq(kl,2)=results(idx(j),4);
            testn{kl,1}='NEW-OLD-S task';
            roi{kl,1}=varnm{i};
            gr{kl,1}='Patients';
            tsum(kl,1)=results(idx(j),1);
            pmin(kl,1)=results(idx(j),2);
            kl=kl+1; 
        end
    end; clear d ix cl_m tmap idx
   
end

tbl=table(testn,roi,gr,tsum,pmin,time,freq,'VariableNames',{'Test','Sensor','Group','tsum','pcor','time','freq'}); disp(tbl)
save([outdir '/scenerep_TFR_sensor_stats_patients_NEW-OLD-S_task_effects.mat'],'tbl');

kl=1; testn=''; roi=''; gr=''; tsum=[]; pmin=[]; time=[]; freq=[];

for i=[10 12 14] % NEW-OLD in patients
    disp(['NEW-OLD-L in patients and in ' varnm{i}])
    
    data1 = squeeze(datarms{i}(p_ind,1,1:45,:));
    data2 = squeeze(datarms{i}(p_ind,4,1:45,:));
     
    [results, perm_max_cl,tmap]=perm_stats_TFR(data1,data2,permnum,0.05,1);
    
    figure('color','w'); tmap=tmap.*(tmap>abs(tinv(0.05,size(data1,1))));
    imagesc(flip(tmap)); colorbar
    caxis([-5 5]); colormap('redblue')
    set(gca, 'XTick',xticks, 'XTickLabel',xlabels,'Fontsize',14); 
    set(gca, 'YTick',[1 16 33 37 41 45], ...
        'YTickLabel',{'45','30','12', '8', '4', '0'},'Fontsize',14); 
    xlabel('Time (ms)'); ylabel('Frequency (Hz)'); 
    title(['NEW-OLD-L tmap in patients in ' varnm{i}])
    print(gcf,'-dbmp', '-r0',[outdir  strrep(varnm{i},' ','_') '_' session  '_Patients_NEW-OLD-L_tmap_' testname '.bmp'],'-r100'); 
    %close(gcf)
    
    if min(results(:,2))>0.07 % nothing significant
        testn{kl,1}='NEW-OLD-L task';
        roi{kl,1}=varnm{i};
        gr{kl,1}='Patients';
        tsum(kl,1)=NaN; pmin(kl,1)=NaN; time(kl,:)=[NaN NaN]; freq(kl,:)=[NaN NaN];
        kl=kl+1;
    else %significant clusters
        idx=find(results(:,2)<0.07);
        for j=1:length(idx)
            time(kl,1)=(results(idx(j),5))*2-102;
            time(kl,2)=(results(idx(j),6))*2-102;
            freq(kl,1)=results(idx(j),3);
            freq(kl,2)=results(idx(j),4);
            testn{kl,1}='NEW-OLD-L task';
            roi{kl,1}=varnm{i};
            gr{kl,1}='Patients';
            tsum(kl,1)=results(idx(j),1);
            pmin(kl,1)=results(idx(j),2);
            kl=kl+1; 
        end
    end; clear d ix cl_m tmap idx
   
end

tbl=table(testn,roi,gr,tsum,pmin,time,freq,'VariableNames',{'Test','Sensor','Group','tsum','pcor','time','freq'}); disp(tbl)
save([outdir '/scenerep_TFR_sensor_stats_patients_NEW-OLD-L_task_effects.mat'],'tbl');

%% Interaction effects in controls-patients

group{1}=c_ind;group{2}=p_ind;
gnm={'Controls','Patients'};
kl=1; testn=''; roi=''; gr=''; tsum=[]; pmin=[]; time=[]; freq=[];

for i=[10 12 14] % NEW-OLD 
    disp(['NEW-OLD interaction effects in ' varnm{i}])
    
    data1 = squeeze(diffrms{i}(c_ind,1,1:45,:));
    data2 = squeeze(diffrms{i}(p_ind,1,1:45,:));
     
    [results, perm_max_cl,tmap]=perm_stats_TFR(data1,data2,permnum,0.05,2);
    
    figure('color','w'); tmap=tmap.*(tmap>abs(tinv(0.05,length(subs)-2)));
    imagesc(flip(tmap)); colorbar
    caxis([-5 5]); colormap('redblue')
    set(gca, 'XTick',xticks, 'XTickLabel',xlabels,'Fontsize',14); 
    set(gca, 'YTick',[1 16 33 37 41 45], ...
        'YTickLabel',{'45','30','12', '8', '4', '0'},'Fontsize',14); 
    xlabel('Time (ms)'); ylabel('Frequency (Hz)'); 
    title(['NEW-OLD tmap interaction in ' varnm{i}])
    print(gcf,'-dbmp', '-r0',[outdir  strrep(varnm{i},' ','_') '_' session  '_Controls-patients_NEW-OLD_tmap_' testname '.bmp'],'-r100'); 
    %close(gcf)
   
    if min(results(:,2))>0.07 % nothing significant
        testn{kl,1}='NEW-OLD interaction';
        roi{kl,1}=varnm{i};
        gr{kl,1}='Control-Patients';
        tsum(kl,1)=NaN; pmin(kl,1)=NaN; time(kl,:)=[NaN NaN]; freq(kl,:)=[NaN NaN];
        kl=kl+1;
    else %significant clusters
        idx=find(results(:,2)<0.07);
        for j=1:length(idx)
            time(kl,1)=(results(idx(j),5)+51)*2-100;
            time(kl,2)=(results(idx(j),6)+51)*2-100;
            freq(kl,1)=results(idx(j),3);
            freq(kl,2)=results(idx(j),4);
            testn{kl,1}='NEW-OLD interaction';
            roi{kl,1}=varnm{i};
            gr{kl,1}='Control-Patients';
            tsum(kl,1)=results(idx(j),1);
            pmin(kl,1)=results(idx(j),2);
            kl=kl+1; 
        end
    end; clear d ix cl_m tmap idx
   
end

tbl=table(testn,roi,gr,tsum,pmin,time,freq,'VariableNames',{'Test','Sensor','Group','tsum','pcor','time','freq'}); disp(tbl)
save([outdir '/scenerep_TFR_sensor_stats_control-patients_NEW-OLD_interaction_effects.mat'],'tbl');

kl=1; testn=''; roi=''; gr=''; tsum=[]; pmin=[]; time=[]; freq=[];

for i=[10 12 14] % NEW-OLD in patients
    disp(['NEW-OLD-S interaction effects in ' varnm{i}])
    
    data1 = squeeze(diffrms{i}(c_ind,2,1:45,:));
    data2 = squeeze(diffrms{i}(p_ind,2,1:45,:));
     
    [results, perm_max_cl,tmap]=perm_stats_TFR(data1,data2,permnum,0.05,2);
    
    figure('color','w'); tmap=tmap.*(tmap>abs(tinv(0.05,length(subs)-2)));
    imagesc(flip(tmap)); colorbar
    caxis([-5 5]); colormap('redblue')
    set(gca, 'XTick',xticks, 'XTickLabel',xlabels,'Fontsize',14); 
    set(gca, 'YTick',[1 16 33 37 41 45], ...
        'YTickLabel',{'45','30','12', '8', '4', '0'},'Fontsize',14); 
    xlabel('Time (ms)'); ylabel('Frequency (Hz)'); 
    title(['NEW-OLD-S tmap interaction in ' varnm{i}])
    print(gcf,'-dbmp', '-r0',[outdir  strrep(varnm{i},' ','_') '_' session  '_Controls-patients_NEW-OLD-S_tmap_' testname '.bmp'],'-r100'); 
    %close(gcf)

    
    if min(results(:,2))>0.07 % nothing significant
        testn{kl,1}='NEW-OLD-S interaction';
        roi{kl,1}=varnm{i};
        gr{kl,1}='Control-Patients';
        tsum(kl,1)=NaN; pmin(kl,1)=NaN; time(kl,:)=[NaN NaN]; freq(kl,:)=[NaN NaN];
        kl=kl+1;
    else %significant clusters
        idx=find(results(:,2)<0.07);
        for j=1:length(idx)
            time(kl,1)=(results(idx(j),5))*2-102;
            time(kl,2)=(results(idx(j),6))*2-102;
            freq(kl,1)=results(idx(j),3);
            freq(kl,2)=results(idx(j),4);
            testn{kl,1}='NEW-OLD-S interaction';
            roi{kl,1}=varnm{i};
            gr{kl,1}='Control-Patients';
            tsum(kl,1)=results(idx(j),1);
            pmin(kl,1)=results(idx(j),2);
            kl=kl+1; 
        end
    end; clear d ix cl_m tmap idx
   
end

tbl=table(testn,roi,gr,tsum,pmin,time,freq,'VariableNames',{'Test','Sensor','Group','tsum','pcor','time','freq'}); disp(tbl)
save([outdir '/scenerep_TFR_sensor_stats_control-patients_NEW-OLD-S_interaction_effects.mat'],'tbl');

kl=1; testn=''; roi=''; gr=''; tsum=[]; pmin=[]; time=[]; freq=[];

for i=[10 12 14] % NEW-OLD in patients
    disp(['NEW-OLD-L interaction effects in ' varnm{i}])
    
    data1 = squeeze(diffrms{i}(c_ind,3,1:45,:));
    data2 = squeeze(diffrms{i}(p_ind,3,1:45,:));
     
    [results, perm_max_cl,tmap]=perm_stats_TFR(data1,data2,permnum,0.05,2);
    
    figure('color','w'); tmap=tmap.*(tmap>abs(tinv(0.05,length(subs)-2)));
    imagesc(flip(tmap)); colorbar
    caxis([-5 5]); colormap('redblue')
    set(gca, 'XTick',xticks, 'XTickLabel',xlabels,'Fontsize',14); 
    set(gca, 'YTick',[1 16 33 37 41 45], ...
        'YTickLabel',{'45','30','12', '8', '4', '0'},'Fontsize',14); 
    xlabel('Time (ms)'); ylabel('Frequency (Hz)'); 
    title(['NEW-OLD-L tmap interaction in ' varnm{i}])
    print(gcf,'-dbmp', '-r0',[outdir  strrep(varnm{i},' ','_') '_' session  '_Controls-patients_NEW-OLD-L_tmap_' testname '.bmp'],'-r100'); 
    %close(gcf)
  
    if min(results(:,2))>0.07 % nothing significant
        testn{kl,1}='NEW-OLD-L interaction';
        roi{kl,1}=varnm{i};
        gr{kl,1}='Control-Patients';
        tsum(kl,1)=NaN; pmin(kl,1)=NaN; time(kl,:)=[NaN NaN]; freq(kl,:)=[NaN NaN];
        kl=kl+1;
    else %significant clusters
        idx=find(results(:,2)<0.07);
        for j=1:length(idx)
            time(kl,1)=(results(idx(j),5))*2-100;
            time(kl,2)=(results(idx(j),6))*2-100;
            freq(kl,1)=results(idx(j),3);
            freq(kl,2)=results(idx(j),4);
            testn{kl,1}='NEW-OLD-L interaction';
            roi{kl,1}=varnm{i};
            gr{kl,1}='Control-Patients';
            tsum(kl,1)=results(idx(j),1);
            pmin(kl,1)=results(idx(j),2);
            kl=kl+1; 
        end
    end; clear d ix cl_m tmap idx
   
end

tbl=table(testn,roi,gr,tsum,pmin,time,freq,'VariableNames',{'Test','Sensor','Group','tsum','pcor','time','freq'}); disp(tbl)
save([outdir '/scenerep_TFR_sensor_stats_control-patients_NEW-OLD-L_interaction_effects.mat'],'tbl');

%% Calculate metrics
band=[1 4; 4 8; 8 12; 12 30; 30 45];

for i=1:18
    for s=1:size(datarms{i},1)
        for j=1:4
            for b=1:size(band,1)
                tw1_con_mean(i,s,b,j)=nanmean(nanmean(datarms{i}(s,j,band(b,1):band(b,2),75:125)));
                tw2_con_mean(i,s,b,j)=nanmean(nanmean(datarms{i}(s,j,band(b,1):band(b,2),150:200)));
                tw3_con_mean(i,s,b,j)=nanmean(nanmean(datarms{i}(s,j,band(b,1):band(b,2),200:300)));
            end
        end
    end
end


for i=1:18
    for s=1:size(diffrms{i},1)
        for j=1:3
            for b=1:size(band,1)
                tw1_dif_mean(i,s,b,j)=nanmean(nanmean(diffrms{i}(s,j,band(b,1):band(b,2),75:125)));
                tw2_dif_mean(i,s,b,j)=nanmean(nanmean(diffrms{i}(s,j,band(b,1):band(b,2),150:200)));
                tw3_dif_mean(i,s,b,j)=nanmean(nanmean(diffrms{i}(s,j,band(b,1):band(b,2),200:300)));
            end
        end
    end
end


%% Test metrics

for i=1:18
    for c=1:4
        for b=1:size(band,1)
            [h p ci stats]=ttest2(squeeze(tw1_con_mean(i,c_ind,b,c)),squeeze(tw1_con_mean(i,p_ind,b,c)),'tail','both','vartype','unequal');
            tmat_tw1_con(i,b,c)=stats.tstat; pmat_tw1_con(i,b,c)=p;
            [h p ci stats]=ttest2(squeeze(tw2_con_mean(i,c_ind,b,c)),squeeze(tw2_con_mean(i,p_ind,b,c)),'tail','both','vartype','unequal');
            tmat_tw2_con(i,b,c)=stats.tstat; pmat_tw2_con(i,b,c)=p;
            [h p ci stats]=ttest2(squeeze(tw3_con_mean(i,c_ind,b,c)),squeeze(tw3_con_mean(i,p_ind,b,c)),'tail','both','vartype','unequal');
            tmat_tw3_con(i,b,c)=stats.tstat; pmat_tw3_con(i,b,c)=p;
        end
    end
end


for i=1:18
    for c=1:3
        for b=1:size(band,1)
            [h p ci stats]=ttest2(squeeze(tw1_dif_mean(i,c_ind,b,c)),squeeze(tw1_dif_mean(i,p_ind,b,c)),'tail','both','vartype','unequal');
            tmat_tw1_dif(i,b,c)=stats.tstat; pmat_tw1_dif(i,b,c)=p;
            [h p ci stats]=ttest2(squeeze(tw2_dif_mean(i,c_ind,b,c)),squeeze(tw2_dif_mean(i,p_ind,b,c)),'tail','both','vartype','unequal');
            tmat_tw2_dif(i,b,c)=stats.tstat; pmat_tw2_dif(i,b,c)=p;
            [h p ci stats]=ttest2(squeeze(tw3_dif_mean(i,c_ind,b,c)),squeeze(tw3_dif_mean(i,p_ind,b,c)),'tail','both','vartype','unequal');
            tmat_tw3_dif(i,b,c)=stats.tstat; pmat_tw3_dif(i,b,c)=p;
        end
    end
end

%% Plot results
band_nm={'Delta','Theta','Alpha','Beta','Gamma'};

for x=1:size(band,1)
    
    t=[squeeze(tmat_tw1_con(:,x,:)), squeeze(tmat_tw2_con(:,x,:)), squeeze(tmat_tw3_con(:,x,:))]; 
    p=[squeeze(pmat_tw1_con(:,x,:)), squeeze(pmat_tw2_con(:,x,:)), squeeze(pmat_tw3_con(:,x,:))];
    
    figure('Color','w'); imagesc(t); axis equal; colormap('redblue'); caxis([-3 3]); xlim([0.5 12.5]); colorbar
    set(gca,'XTick',1:12,'XTickLabel',{'NEW','OLD','OLD-S','OLD-L'}); xtickangle(40)
    set(gca,'YTick',1:18,'YTickLabel',varnm);
    [a,b]=ind2sub(size(p),find(p<=0.05)); hold on;
    for d=1:length(a)
        plot(b(d),a(d),'-s','markersize',8,'markeredgecolor','k','markerfacecolor','k')
    end
    [a,b]=ind2sub(size(p),find(p<=0.01)); hold on;
    for d=1:length(a)
        plot(b(d),a(d),'-s','markersize',12,'markeredgecolor','k','markerfacecolor','k')
    end
    title(['Group differences in ' band_nm{x} ' single condition'])
    print(gcf,'-dbmp', '-r0',[outdir 'scenerep_group_dif_tmap_TFR_' band_nm{x} '_single_condition.bmp'],'-r300'); %close(gcf)
    
end

for x=1:size(band,1)
    
    t=[squeeze(tmat_tw1_dif(:,x,:)), squeeze(tmat_tw2_dif(:,x,:)), squeeze(tmat_tw3_dif(:,x,:))]; 
    p=[squeeze(pmat_tw1_dif(:,x,:)), squeeze(pmat_tw2_dif(:,x,:)), squeeze(pmat_tw3_dif(:,x,:))];
    
    figure('Color','w'); imagesc(t); axis equal; colormap('redblue'); caxis([-3 3]); xlim([0.5 9.5]); colorbar
    set(gca,'XTick',1:9,'XTickLabel',{'NEW-OLD','NEW-OLD-S','NEW-OLD-L'}); xtickangle(40)
    set(gca,'YTick',1:18,'YTickLabel',varnm);
    [a,b]=ind2sub(size(p),find(p<=0.05)); hold on;
    for d=1:length(a)
        plot(b(d),a(d),'-s','markersize',8,'markeredgecolor','k','markerfacecolor','k')
    end
    [a,b]=ind2sub(size(p),find(p<=0.01)); hold on;
    for d=1:length(a)
        plot(b(d),a(d),'-s','markersize',12,'markeredgecolor','k','markerfacecolor','k')
    end
    title(['Group differences in ' band_nm{x} ' contrasts'])
    print(gcf,'-dbmp', '-r0',[outdir 'scenerep_group_dif_tmap_TFR_' band_nm{x} '_contrasts.bmp'],'-r300'); %close(gcf)
    %close(gcf)
    
end

%% Calculate AUC

for i=1:18
    for c=1:4
        for b=1:size(band,1)
            param=roc_curve(squeeze(tw1_con_mean(i,p_ind,b,c)),squeeze(tw1_con_mean(i,c_ind,b,c)));
            auc_tw1_con(i,b,c)=param.param.AROC;
            param=roc_curve(squeeze(tw2_con_mean(i,p_ind,b,c)),squeeze(tw2_con_mean(i,c_ind,b,c)));
            auc_tw2_con(i,b,c)=param.param.AROC;
            param=roc_curve(squeeze(tw3_con_mean(i,p_ind,b,c)),squeeze(tw3_con_mean(i,c_ind,b,c)));
            auc_tw3_con(i,b,c)=param.param.AROC;
        end
    end
end

for i=1:18
    for c=1:3
        for b=1:size(band,1)
            param=roc_curve(squeeze(tw1_dif_mean(i,p_ind,b,c)),squeeze(tw1_dif_mean(i,c_ind,b,c)));
            auc_tw1_dif(i,b,c)=param.param.AROC;
            param=roc_curve(squeeze(tw2_dif_mean(i,p_ind,b,c)),squeeze(tw2_dif_mean(i,c_ind,b,c)));
            auc_tw2_dif(i,b,c)=param.param.AROC;
            param=roc_curve(squeeze(tw3_dif_mean(i,p_ind,b,c)),squeeze(tw3_dif_mean(i,c_ind,b,c)));
            auc_tw3_dif(i,b,c)=param.param.AROC;
        end
    end
end

%% Plot AUC results
cmap=customcolormap([0 0.9 1],[29/255 180/255 191/255 ;1 1 1;1 1 1 ]);

for x=1:size(band,1)
    
    auc=[squeeze(auc_tw1_con(:,x,:)), squeeze(auc_tw2_con(:,x,:)), squeeze(auc_tw3_con(:,x,:))];
    idx=find(auc<0.5);
    auc(idx)=abs(1-auc(idx));
    
    figure('Color','w'); imagesc(auc); axis equal; caxis([0.5 1]); xlim([0.5 12.5]); colorbar
    set(gca,'XTick',1:12,'XTickLabel',{'NEW','OLD','OLD-S','OLD-L'}); xtickangle(40)
    set(gca,'YTick',1:18,'YTickLabel',varnm); colormap(cmap)
    
    [a,b]=ind2sub(size(auc),find(auc>=0.7)); hold on;
    for d=1:length(a)
        plot(b(d),a(d),'-s','markersize',8,'markeredgecolor','k','markerfacecolor','k')
    end
    title(['AUC in ' band_nm{x} ' single condition'])
    print(gcf,'-dbmp', '-r0',[outdir 'scenerep_AUC_' band_nm{x} '_single_condition.bmp'],'-r300'); %close(gcf)
    
end

for x=1:size(band,1)
    
    auc=[squeeze(auc_tw1_dif(:,x,:)), squeeze(auc_tw2_dif(:,x,:)), squeeze(auc_tw3_dif(:,x,:))];
    idx=find(auc<0.5);
    auc(idx)=abs(1-auc(idx));
    
    figure('Color','w'); imagesc(auc); axis equal; colormap(cmap); caxis([0.5 1]); xlim([0.5 9.5]); colorbar
    set(gca,'XTick',1:9,'XTickLabel',{'NEW-OLD','NEW-OLD-S','NEW-OLD-L'}); xtickangle(40)
    set(gca,'YTick',1:18,'YTickLabel',varnm);
    
    [a,b]=ind2sub(size(auc),find(auc>=0.7)); hold on;
    for d=1:length(a)
        plot(b(d),a(d),'-s','markersize',8,'markeredgecolor','k','markerfacecolor','k')
    end
    title(['AUC in ' band_nm{x} ' contrasts'])
    print(gcf,'-dbmp', '-r0',[outdir 'scenerep_AUC_' band_nm{x} '_contrasts.bmp'],'-r300');
    %close(gcf)
    
end