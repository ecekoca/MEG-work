%% NTAD Sensor level RMS analysis for Roving (Ece K, 2019)

addpath /imaging/local/meg_misc/
addpath /imaging/local/software/spm_cbu_svn/releases/spm12_latest/;
addpath /hpc-software/matlab/cbu/;
addpath(genpath('/imaging/local/software/spm_cbu_svn/releases/spm12_latest/external/fieldtrip/')); %EK
addpath('/home/ek01/scripts/ntad/functions/')
addpath('/home/ek01/scripts/all_scripts/')
addpath('/imaging/projects/cbu/ntad/scripts/functions/')
addpath('/imaging/rowe/users/ek01/collab/MKL_Delshad_Rik/BioFIND-data-paper-master/')

close all
clear all; clc
subs = {'C1001','C1002','C1004','C1006','C1010','C1011','C1012','C1014',...
'C1015','C1017','C2001','C2002','C2008',...
'c2','c5','c7','c10','c12','c13','c14',...
'P1001','P1005','P1007','P1008','P1009','P1010','P1011',...
'P1012','P1015','P1020','P1021','P1022','P1023','P1024','P1026',...
'P1027','P1029','P1030','P1031','P1032','P1036','P1037',...
'P1038','P1042','P1043','P1048','P1049','P1053','P1055',...
'P1056','P1058','P1060','P1061','P1062','P1063',...
'P1064','P1065','P1066','P1067','P3004'}; % c* are TGB controls, C2* are Oxford NTAD controls

% No MMN: 'P1016' BL exclude: 'C1007', 'c8'

datadir = '/imaging/projects/cbu/ntad/meg_data/';
datadir2= '/imaging/projects/cbu/ntad/meg_data/TGB_controls/';
outdir = '/imaging/rowe/users/ek01/prv/ntad/mmn_sensor/roving_mmn_normalised_signal/';
midfol='/BL/mmn/';
lastbit ='_mmn_all_trl.mat';
session='BL';
testname='norm_signal';

xticks=[1 50 100 150 200 251];
xlabels={'-100' '0' '100' '200' '300' '400'};
conditions={'DEV', 'REP1','REP2','REP3','REP4','REP5','REP6','REP7','REP8','REP9','REP10'};

load([outdir 'input_values.mat']) % should be in the same order as the subs above
age=T.AGE;acer=T.ACER; mmse=T.MMSE; mem=T.ACER_MEM;
edu=T.YOE; gm=T.GM; wm=T.WM; hip=T.HIP;

green=[0 201 121]./255;%[0 130 144]./255
blue=[18 27 194]./255;
navy=[29 41 91]./255;
yellow=[254 210 1]./255;
grey=[94 94 94]./255;
red=[198 4 4]./255;


%% Channel selection
fr_mag={'MEG0111','MEG0121','MEG0341','MEG0131','MEG1411','MEG1421','MEG1221','MEG1441',...
    'MEG0211','MEG0221','MEG0321','MEG1431','MEG1321','MEG1311','MEG1231'};
fr_grad={'MEG0112','MEG0113','MEG0123','MEG0122','MEG0342','MEG0343','MEG0133','MEG0132',...
    'MEG1412','MEG1413','MEG1422','MEG1423','MEG1222','MEG1223','MEG1442','MEG1443',...
    'MEG0212','MEG0213','MEG0222','MEG0223','MEG0322','MEG0323',...
    'MEG1432','MEG1433','MEG1322','MEG1323','MEG1312','MEG1313','MEG1232','MEG1233'};

% fr_mag={'MEG0111','MEG0121','MEG0341','MEG0131','MEG1411','MEG1421','MEG1221','MEG1441'};
% fr_grad={'MEG0112','MEG0113','MEG0123','MEG0122','MEG0342','MEG0343','MEG0133','MEG0132',...
%     'MEG1412','MEG1413','MEG1422','MEG1423','MEG1222','MEG1223','MEG1442','MEG1443'};
fr_eeg={'EEG002','EEG004','EEG005','EEG006','EEG007','EEG008','EEG009','EEG010','EEG011','EEG012','EEG013','EEG014','EEG015',...
    'EEG016','EEG017','EEG019','EEG020','EEG021','EEG022','EEG023','EEG024','EEG025','EEG026','EEG027'};

tem_mag={'MEG0241','MEG0231','MEG0441','MEG1611','MEG1621','MEG1811','MEG1821','MEG1641','MEG1631','MEG1841',...
    'MEG1131','MEG1341','MEG1331','MEG2211','MEG2221','MEG2411','MEG2431','MEG2441'};
tem_grad={'MEG0242','MEG0243','MEG0233','MEG0232','MEG0442','MEG0443','MEG1612','MEG1613','MEG1623','MEG1622','MEG1812',...
    'MEG1813','MEG1823','MEG1822','MEG1642','MEG1643','MEG1633','MEG1632','MEG1842','MEG1843',...
    'MEG1132','MEG1133','MEG1342','MEG1343','MEG1332','MEG1333','MEG2212','MEG2213','MEG2222','MEG2223','MEG2412','MEG2413',...
    'MEG2432','MEG2433','MEG2442','MEG2443'};
tem_eeg={'EEG030','EEG031','EEG032','EEG033','EEG034','EEG035','EEG036','EEG037','EEG041','EEG042','EEG043','EEG044','EEG045',...
    'EEG046','EEG047','EEG048'};

%% Gather data

c_ind=[find(contains(subs,'C')) find(contains(subs,'c'))];
p_ind=find(contains(subs,'P'));

cd(datadir);
mkdir(outdir);
mkdir([outdir '/GLMs/']); mkdir([outdir '/Corr/']); mkdir([outdir '/Par_corr/']); mkdir([outdir '/T_tests/']);
mkdir([outdir '/AUC/']); mkdir([outdir '/ERFs/']);

cn=1:9;

for sub=1:length(subs)
    
    if contains(subs{sub},'c')
        filen=[datadir2 subs{sub} midfol subs{sub} '_' session lastbit];
    else
        filen=[datadir subs{sub} midfol subs{sub} '_' session lastbit];
    end
    D=spm_eeg_load(filen);
    bad=D.badchannels;
    %origdata(sub,1:size(D,1),:,:)=D(:,:,1:10);
%     if size(D,1)>320
%         coormeg(sub,:,1:380)=D.coor2D(1:380);
%     else
%         coormeg(sub,:,1:306)=D.coor2D;
%     end
    x=find(ismember(D.selectchannels('MEGMAG'),bad));
    data{1}(sub,:,:,:)=D(D.selectchannels('MEGMAG'),:,cn);
    if ~isempty(x); data{1}(sub,x,:,:)=NaN; end
    x=find(ismember(D.selectchannels(fr_mag),bad));
    data{2}(sub,:,:,:)=D(D.selectchannels(fr_mag),:,cn);
    if ~isempty(x); data{2}(sub,x,:,:)=NaN; end
    x=find(ismember(D.selectchannels(tem_mag),bad));
    data{3}(sub,:,:,:)=D(D.selectchannels(tem_mag),:,cn);
    if ~isempty(x); data{3}(sub,x,:,:)=NaN; end
    x=find(ismember(D.selectchannels('MEGPLANAR'),bad));
    data{4}(sub,:,:,:)=D(D.selectchannels('MEGPLANAR'),:,cn);
    if ~isempty(x); data{4}(sub,x,:,:)=NaN; end
    x=find(ismember(D.selectchannels(fr_grad),bad));
    data{5}(sub,:,:,:)=D(D.selectchannels(fr_grad),:,cn);
    if ~isempty(x); data{5}(sub,x,:,:)=NaN; end
    x=find(ismember(D.selectchannels(tem_grad),bad));
    data{6}(sub,:,:,:)=D(D.selectchannels(tem_grad),:,cn);
    if ~isempty(x); data{6}(sub,x,:,:)=NaN; end
    
    if size(D,1)>350 && ~contains(subs{sub},'C2')
        x=find(ismember(D.selectchannels('EEG'),bad));
        data{7}(sub,:,:,:)=D(D.selectchannels('EEG'),:,cn);
        if ~isempty(x); data{7}(sub,x,:,:)=NaN; end
        x=find(ismember(D.selectchannels(fr_eeg),bad));
        data{8}(sub,:,:,:)=D(D.selectchannels(fr_eeg),:,cn);
        if ~isempty(x); data{8}(sub,x,:,:)=NaN; end
        x=find(ismember(D.selectchannels(tem_eeg),bad));
        data{9}(sub,:,:,:)=D(D.selectchannels(tem_eeg),:,cn);
        if ~isempty(x); data{9}(sub,x,:,:)=NaN; end
    else
        data{7}(sub,:,:,:)=NaN;
        data{8}(sub,:,:,:)=NaN;
        data{9}(sub,:,:,:)=NaN;
    end
end

for i=1:6
    for s=1:size(data{i},1)
        for tr=1:size(data{i},4)
            for ch=1:size(data{i},2)
                if data{i}(s,ch,100,tr)>0  
                    data{i}(s,ch,:,tr)=-(data{i}(s,ch,:,tr)); % N100 is negative, P300 is positive
                end
            end
        end
    end
end

for i=7:9
    for s=1:size(data{i},1)
        for tr=1:size(data{i},4)
            for ch=1:size(data{i},2)
                d=smooth(data{i}(s,ch,:,tr),20,'moving');
                %if data{i}(s,ch,85,tr)>0
                if d(100)>0    
                    %data{i}(s,ch,:,tr)=-(data{i}(s,ch,:,tr)); % N100 is negative
                    data{i}(s,ch,:,tr)=-d; % N100 is negative
                    %data{i}(s,ch,:,tr)=NaN(1,1,501);
                else
                    data{i}(s,ch,:,tr)=d; % N100 is negative
                end             
            end
        end
    end
end

for i=1:9
    %datarms{i}(:,:,:,:)=permute(data{i}(:,:,:,:),[2 1 4 3]);
    %datarms{i}=squeeze(rms_ek(datarms{i},1,[1 100]));
    d=squeeze(nanmean(data{i}(:,:,:,:),2));
    d=permute(d,[1 3 2]);
    datarms{i}=d;
    
    for j=1:size(datarms{i},1)
        for k=1:size(datarms{i},2)
            d=squeeze(datarms{i}(j,k,:));
            bl=nanmean(d(1:51)); d=d-bl; % baseline correct again just in case
            if i<7
                datarms{i}(j,k,:)=smooth(d,20,'moving');
            end
        end
    end
    
    %ind1=find(isnan(datarms{i}(:,[1 4 7],:)));
    ind=find(isnan(datarms{i}));
    for s=1:size(datarms{i},1)
        d=mat2gray(datarms{i}(s,:,:));
        bl=nanmean(d(1:51)); d=d-bl;
        datarms2{i}(s,:,:)=d;
    end

    datarms2{i}(ind)=NaN;
    datarms2{i}=datarms2{i}(:,:,1:251);
%     datarms2{i}(ind2)=NaN;
%     datarms2{i}=datarms2{i}(:,:,1:251);
end; datarms=datarms2;

%% Plot Means

clear data
varnm={'MAG','Frontal MAG','Temporal MAG','GRAD','Frontal GRAD','Temporal GRAD','EEG','Frontal EEG','Temporal EEG'};

for i=1:9
    
    std_c=6;
    
    data=datarms{i}(c_ind,:,:);
    
    m1=squeeze(nanmean(data,1));
    s1=squeeze(nanstd(data)./sqrt(size(data,1)));
    
    [h,p,ci,stats]=ttest(squeeze(data(:,1,:)),squeeze(data(:,std_c,:)));
    %[pfdr ~]=fdr(p,0.05);
    k=max(max([m1(std_c,:); m1(1,:)]))+max(max([s1(std_c,:);s1(1,:)]));k=k*1.1;
    p=0.3.*(p<0.05); p(find(p==0))=NaN; p(1:51)=NaN;
    
    figure('Color','w');set(gcf, 'Units','normalized', 'Position', [0 0 .25 .4]);
    plot(m1(std_c,:),'LineWidth',3,'Color',grey); hold on
    plot(m1(1,:),'LineWidth',3,'Color',green); legendflex({'REP6','DEV'},'fontsize',14)
    boundedline([1:size(m1,2)],m1(std_c,:),s1(std_c,:),'cmap',grey, 'transparency', 0.3, 'orientation', 'vert', 'alpha');
    boundedline([1:size(m1,2)],m1(1,:),s1(1,:),'cmap',green, 'transparency', 0.3, 'orientation', 'vert', 'alpha');
    plot(p,'-o','LineWidth',2,'MarkerEdgeColor',[0.6 0.6 0.6],'MarkerFaceColor',[0.6 0.6 0.6],'MarkerSize',8)
    set(gca, 'XTick',xticks, 'XTickLabel',xlabels,'Fontsize',16); xlabel('Time (ms)'); ylabel('au')
    title([varnm{i} ' mean in controls']); xlim([0 251]); ylim([-0.6 0.4]);
    %if i>6; ylim([-0.3 0.3]); else; ylim([-0.6 0.4]); end 
    print(gcf,'-dbmp', '-r0',[outdir '/ERFs/' strrep(varnm{i},' ','_') '_' session  '_mean_controls_' testname '.bmp'],'-r200'); close(gcf)
    
    
    data=datarms{i}(p_ind,:,:);
    
    m2=squeeze(nanmean(data,1));
    s2=squeeze(nanstd(data)./sqrt(size(data,1)));
    
    [h,p,ci,stats]=ttest(squeeze(data(:,1,:)),squeeze(data(:,std_c,:)));
    %[pfdr ~]=fdr(p,0.05);
    k=max(max([m1(std_c,:); m1(1,:)]))+max(max([s1(std_c,:);s1(1,:)]));k=k*1.1;
    p=0.3.*(p<0.05); p(find(p==0))=NaN; p(1:51)=NaN;
    
    figure('Color','w');set(gcf, 'Units','normalized', 'Position', [0 0 .25 .4]);
    plot(m2(std_c,:),'LineWidth',3,'Color',grey); hold on
    plot(m2(1,:),'LineWidth',3,'Color',green); legendflex({'REP6','DEV'},'fontsize',14)
    boundedline([1:size(m2,2)],m2(std_c,:),s2(std_c,:),'cmap',grey, 'transparency', 0.3, 'orientation', 'vert', 'alpha');
    boundedline([1:size(m2,2)],m2(1,:),s2(1,:),'cmap',green, 'transparency', 0.3, 'orientation', 'vert', 'alpha');
    plot(p,'-o','LineWidth',2,'MarkerEdgeColor',[0.6 0.6 0.6],'MarkerFaceColor',[0.6 0.6 0.6],'MarkerSize',8)
    set(gca, 'XTick',xticks, 'XTickLabel',xlabels,'Fontsize',16); xlabel('Time (ms)'); ylabel('au')
    title([varnm{i} ' mean in patients']); xlim([0 251]); ylim([-0.6 0.4]);
    %if i>6; ylim([-0.3 0.3]); else; ylim([-0.6 0.4]); end 
    print(gcf,'-dbmp', '-r0',[outdir '/ERFs/' strrep(varnm{i},' ','_') '_' session   '_mean_patients_' testname '.bmp'],'-r200'); close(gcf)
  
    c_data=squeeze(datarms{i}(c_ind,std_c,:)-datarms{i}(c_ind,1,:));
    m1=squeeze(nanmean(c_data,1));
    s1=squeeze(nanstd(c_data)./sqrt(size(c_data,1)));
%     
    p_data=squeeze(datarms{i}(p_ind,std_c,:)-datarms{i}(p_ind,1,:));
    m2=squeeze(nanmean(p_data,1));
    s2=squeeze(nanstd(p_data)./sqrt(size(p_data,1)));
    
    [h,p,ci,stats]=ttest2(c_data,p_data,'tail','right');
    %[pfdr ~]=fdr(p,0.05);
    %k=max(max([m1(1,:); m2(1,:)]))+max(max([s1(1,:);s2(1,:)]));k=k*1.2;
    p=0.25.*(p<0.05); p(find(p==0))=NaN; p(1:51)=NaN;

    figure('Color','w');set(gcf, 'Units','normalized', 'Position', [0 0 .25 .4]);
    plot(m1(1,:),'LineWidth',3,'Color','b'); hold on
    plot(m2(1,:),'LineWidth',3,'Color','r');legendflex({'C','P'},'fontsize',14)
    boundedline([1:size(m1,2)],m1(1,:),s1(1,:),'cmap','b', 'transparency', 0.3, 'orientation', 'vert', 'alpha');
    boundedline([1:size(m2,2)],m2(1,:),s2(1,:),'cmap','r', 'transparency', 0.3, 'orientation', 'vert', 'alpha');
    plot(p,'-o','LineWidth',2,'MarkerEdgeColor',[0.6 0.6 0.6],'MarkerFaceColor',[0.6 0.6 0.6],'MarkerSize',8)
    set(gca, 'XTick',xticks, 'XTickLabel',xlabels,'Fontsize',16); xlabel('Time (ms)'); ylabel('au')
    title([varnm{i} ' MMN']); xlim([0 251]); 
    ylim([-0.3 0.3]); %axis square
    print(gcf,'-dbmp', '-r0',[outdir '/ERFs/' strrep(varnm{i},' ','_') '_' session '_MMN_all_' testname '.bmp'],'-r200'); close(gcf)
    
end

%% Develop metrics

clear data
window{1}=100:150; window{2}=175:250; window_nm={'100-200ms','250-400ms'};

for w=1%:2
for i=1:size(datarms,2)
    
    for s=1:size(datarms{1},1)
        
        data=squeeze(datarms{i}(s,:,:));
        
        dev_m(s)=nanmean(data(1,window{w})); % Mean deviant
        rep2_m(s)=nanmean(data(2,window{w})); % Mean 1st rep
        rep6_m(s)=nanmean(data(6,window{w})); % Mean 5th repetition
        amp_m(s,:)=nanmean(data(cn,window{w}),2);
        
        if ~isnan(dev_m(s))
                coefs=fit([1:2]',[amp_m(s,1:2)]','poly1','normalize','on');
                dev_rep2_dif(s)=coefs.p1;
                rep6_dev_m(s)=rep6_m(s)-dev_m(s);   
                t=80:170; %0 - 300 ms 
                a=data(6,t)-data(1,t);
                rep6_dev_lat(s)= (t(find(a==max(a)))*2-100)/1000;
                coefs=fit([1:6]',[amp_m(s,1:6)]','poly1','normalize','on');
                slope_d36(s)=coefs.p1;
                coefs=fit([1:6]',[amp_m(s,1:6)]','poly3','normalize','on'); % cubic polynomial
                slope_d36_q(s)=coefs.p1; %
                %cs(s,:)=[coefs.p1,coefs.p2,coefs.p3];
        else
            rep6_dev_m(s)=NaN;
            rep6_dev_lat(s)=NaN;
            slope_d36(s)=NaN; 
            slope_d36_q(s)=NaN;
            dev_rep2_dif(s)=NaN;
        end
    end
    
    metrics=[rep6_dev_m' rep6_dev_lat' dev_rep2_dif' slope_d36'  slope_d36_q'];
    
    met_by_sens{1}(:,i)=rep6_dev_m';
    met_by_sens{2}(:,i)=rep6_dev_lat'; 
    met_by_sens{3}(:,i)=dev_rep2_dif';
    met_by_sens{4}(:,i)=slope_d36'; 
    met_by_sens{5}(:,i)=slope_d36_q';
    
    met_names={'REP6_DEV','REP6-DEV_LAT','SLOPE_D2','SLOPE_D6','CSLOPE_D6'};
    varnames={'REP6-DEV','REP6-DEV LAT','SLOPE D2','SLOPE D6','CSLOPE D6'};
    T=table(subs',rep6_dev_m', rep6_dev_lat',dev_rep2_dif', slope_d36', slope_d36_q',...
        'VariableNames',{'SS','REP6_DEV','REP6_DEV_LAT','SLOPE_D2','SLOPE_D6','CSLOPE_D6'}) 
    save([outdir 'roving_metrics_' session '_' strrep(varnm{i},' ','_') '_' testname '_' window_nm{w} '.mat'],'T')
    writetable(T,[outdir 'roving_metrics_' session '_' strrep(varnm{i},' ','_') '_' testname '_' window_nm{w} '.xlsx'])
    
    %% Amplitude change plots
    
    if ~exist([outdir '/amp_change/'],'dir'); mkdir([outdir '/amp_change/']); end
    
    grey2=[0.2 0.2 0.2];
    c_d=amp_m(c_ind,:); p_d=amp_m(p_ind,:);
    m1=nanmean(c_d);m2=nanmean(p_d);
    s1=squeeze(nanstd(c_d)./sqrt(size(c_d,1)));
    s2=squeeze(nanstd(p_d)./sqrt(size(p_d,1)));
    
    figure; errorbar(1:length(m1),m1,s1,'-s','MarkerSize',9,'MarkerEdgeColor',navy,...
        'MarkerFaceColor',navy,'Color',navy, 'LineWidth',1.5);
    hold on; xlim([0 10])
    errorbar(1:length(m2),m2,s2,'-v','MarkerSize',9,'MarkerEdgeColor',red,...
        'MarkerFaceColor',red,'Color',red, 'LineWidth',1.5);
    set(gca, 'XTick',1:9, 'XTickLabel',{'DEV','R2','R3','R4','R5','STD','R7','R8','R9'},'Fontsize',11);
    ylabel('Mean (au)'); title([varnm{i} ' in ' window_nm{w}]); legend({'C','P'},'fontsize',11); axis square
    print(gcf,'-dbmp', '-r0',[outdir '/amp_change/' strrep(varnm{i},' ','_') '_' session '_amplitude_change_' window_nm{w} '.bmp'],'-r200'); close(gcf)
    
    %% Scatterplots
    
    var=[mmse acer mem hip gm wm];
    met=metrics;
    cov_nm={'MMSE','ACER','ACER-mem','HIP','GMV','WMV'};
    c_col=[0 192 163]./255;
    p_col=[0.3 0.3 0.3];
    
    
    for m=1:length(met_names)
        for v=1:length(cov_nm)
            figure('color','w');
            scatter(var(:,v),met(:,m),80,'w','filled');h=lsline; h.LineWidth=3; h.Color='k';hold on;
            scatter(var(c_ind,v),met(c_ind,m),120,c_col,'filled');
            scatter(var(p_ind,v),met(p_ind,m),120,p_col,'filled');
            axis square; box on;
            xlabel(cov_nm{v}); ylabel(varnames{m});
            %title([varnm{m} ' x ptau181']);
            set(gca,'FontSize',14); grid on;
            print(gcf,[outdir '/scatterplots/scatterplot_' strrep(cov_nm{v},' ','_') '_' met_names{m} '_' varnm{i} '.bmp'],'-dbmp','-r300'); close(gcf)
        end
    end
    
    for m=1:length(met_names)
        for v=1:length(cov_nm)
            a=met(:,m); b=var(:,v);
            idx=unique([find(isnan(a)) find(isnan(b))]);
            a(idx)=[]; b(idx)=[];
            [r p]=corrcoef(a,b);
            disp([varnm{i} ': ' met_names{m} 'x' cov_nm{v} ' r= ' num2str(r(1,2)) ' p=' num2str(p(1,2)/2)])
        end
    end
    
    %% Boxplots
    if ~exist([outdir '/boxplots/'],'dir'); mkdir([outdir '/boxplots/']); end
    gr = [0.1*ones(length(c_ind),1); 0.5*ones(length(p_ind),1)];
    for m=1:length(met_names)
        met=metrics;
        figure('color','w');
        H=notBoxPlot([met(c_ind,m); met(p_ind,m)],gr,'jitter',0.3) % plot the two groups data points using notboxplot function
        set([H.data],'MarkerSize',5,...
            'markerFaceColor',[1,1,1]*0.25,...
            'markerEdgeColor', 'none')
        set([H([1:2:end]).semPtch],...
            'FaceColor',[100 222 187 ]./255,...
            'EdgeColor','none')
        set([H([2:2:end]).semPtch],...
            'FaceColor',[170 170 170]./255,...
            'EdgeColor','none')
        
        set([H([1:2:end]).sdPtch],...
            'FaceColor',[163 255 229]./255,...
            'EdgeColor','none')
        set([H([2:2:end]).sdPtch],...
            'FaceColor',[200 200 200]./255,...
            'EdgeColor','none')
        set([H.mu],...
            'Color','k')
        set(gca, 'XtickLabel', {'C','MCI/AD'})
        
        title([varnames{m}]); set(gca,'fontsize',13); box on
        print(gcf,'-dbmp', '-r0',[outdir '/boxplots/Boxplots_' strrep(varnm{i},' ','_') '_' met_names{m} '.bmp'],'-r200'); close(gcf)    
    end
   
    %% Intercorrelations
    
        group(c_ind,1)=1; group(p_ind,1)=0;
        [r, p]=corrcoef([group,age,acer,mmse,mem,gm,wm,hip,metrics],'rows','complete');
        r=tril(r);p=p/2;
        figure('color','w'); imagesc(r);
        axis square; colormap('redblue');
        caxis([-0.4 0.4]);
        yticks(1:13); 
        vnm={'Group','Age','ACER','MMSE','ACER-mem','GMV','WMV','HIP','REP6-DEV','REP6-DEV LAT','SLOPE D2','SLOPE D6','CSLOPE D6'};
        set(gca,'YTickLabel',vnm,...
            'XTick',1:length(vnm),'XTickLabel',vnm)
        xtickangle(30);colorbar; box off
        title(['Intercorrelations in ' varnm{i} ' in ' window_nm{w}])
        [a,b]=ind2sub(size(squeeze(p')),find(squeeze(p'<=0.05))); hold on;
        for d=1:length(a)
            plot(b(d),a(d),'-s','markersize',8,'markeredgecolor','w','markerfacecolor','w')
        end
        [a,b]=ind2sub(size(squeeze(p')),find(squeeze(p'<=0.01))); hold on;
        for d=1:length(a)
            plot(b(d),a(d),'-s','markersize',12,'markeredgecolor','w','markerfacecolor','w')
        end
        print(gcf,'-dbmp', '-r0',[outdir '/Corr/' varnm{i} '_metric_intercorrelations_' window_nm{w} '.bmp'],'-r200'); close(gcf)
        
        [r,p]=partialcorr([group,age,acer,mmse,mem,gm,wm,hip,metrics],age,'rows','complete');
        r=tril(r);p=p/2;
        figure('color','w'); imagesc(r);
        axis square; colormap('redblue');
        caxis([-0.4 0.4]);
        yticks(1:13); 
        vnm={'Group','Age','ACER','MMSE','ACER-mem','GMV','WMV','HIP','REP6-DEV','REP6-DEV LAT','SLOPE D2','SLOPE D6','CSLOPE D6'};
        set(gca,'YTickLabel',vnm,...
            'XTick',1:length(vnm),'XTickLabel',vnm)
        xtickangle(30);colorbar; box off
        title(['Partial intercorrelations in ' varnm{i} ' in ' window_nm{w}])
        [a,b]=ind2sub(size(squeeze(p')),find(squeeze(p'<=0.05))); hold on;
        for d=1:length(a)
            plot(b(d),a(d),'-s','markersize',8,'markeredgecolor','w','markerfacecolor','w')
        end
        [a,b]=ind2sub(size(squeeze(p')),find(squeeze(p'<=0.01))); hold on;
        for d=1:length(a)
            plot(b(d),a(d),'-s','markersize',12,'markeredgecolor','w','markerfacecolor','w')
        end
        print(gcf,'-dbmp', '-r0',[outdir '/Corr/' varnm{i} '_metric_partial_intercorrelations_' window_nm{w} '.bmp'],'-r200'); close(gcf)

    %% Basic comparisons
        
    for m=1:size(metrics,2)
        
        met=metrics(:,m);
        cmet=metrics(c_ind,m);
        idx=find(isoutlier(cmet));cmet(idx)=[];
        pmet=metrics(p_ind,m);
        idx=find(isoutlier(pmet));pmet(idx)=[];
        
        [h p ci stats]=ttest2(cmet,pmet,'vartype','unequal');
        hyp(m)=h;
        pval(m)=p;
        degf(m)=stats.df;
        tval(m)=stats.tstat;
        tot_t(i,m)=stats.tstat;
        
        [h p ci stats]=ttest2(cmet,pmet,'vartype','unequal','tail','right'); % C>P
        hyp_r(m)=h;
        pval_r(m)=p;
        tot_p1(i,m)=p;
        
        [h p ci stats]=ttest2(cmet,pmet,'vartype','unequal','tail','left'); % P>C
        hyp_l(m)=h;
        pval_l(m)=p;
        tot_p2(i,m)=p;
        
        X=[];
        X(c_ind,1)=1;
        X(p_ind,2)=1;
        %X(:,3)=zscore(age);
        %X(:,4)=zscore(edu);
        met=metrics(:,m); idx=find(isnan(met));
        met(idx)=[]; X(idx,:)=[]; 
        %if m==3; met=zscore(met); end
        idx=find(isoutlier(met));met(idx)=[];X(idx,:)=[];
        [t,F,p,df,R2,cR2,B] = glm_rik(met,X,[1 -1],1);
        tot_t_glm1(i,m)=t; tot_p_glm1(i,m)=p;
%         
        met=metrics(:,m);
        [r,p]=corrcoef(met(p_ind),acer(p_ind),'rows','complete');
        acer_r(i,m)=r(1,2); acer_p(i,m)=p(1,2)/2;
        [r,p]=corrcoef(met(p_ind),mmse(p_ind),'rows','complete');
        mmse_r(i,m)=r(1,2); mmse_p(i,m)=p(1,2)/2;
        [r,p]=corrcoef(met(p_ind),gm(p_ind),'rows','complete');
        gm_r(i,m)=r(1,2); gm_p(i,m)=p(1,2)/2;
        [r,p]=corrcoef(met(p_ind),wm(p_ind),'rows','complete');
        wm_r(i,m)=r(1,2); wm_p(i,m)=p(1,2)/2;
        [r,p]=corrcoef(met(p_ind),hip(p_ind),'rows','complete');
        hip_r(i,m)=r(1,2); hip_p(i,m)=p(1,2)/2;
        
        met=metrics(:,m);
        [r,p]=corrcoef(met,acer,'rows','complete');
        acer_r2(i,m)=r(1,2); acer_p2(i,m)=p(1,2)/2;
        [r,p]=corrcoef(met,mmse,'rows','complete');
        mmse_r2(i,m)=r(1,2); mmse_p2(i,m)=p(1,2)/2;
        [r,p]=corrcoef(met,gm,'rows','complete');
        gm_r2(i,m)=r(1,2); gm_p2(i,m)=p(1,2)/2;
        [r,p]=corrcoef(met,wm,'rows','complete');
        wm_r2(i,m)=r(1,2); wm_p2(i,m)=p(1,2)/2;
        [r,p]=corrcoef(met,hip,'rows','complete');
        hip_r2(i,m)=r(1,2); hip_p2(i,m)=p(1,2)/2;
        [r,p]=corrcoef(met,mem,'rows','complete');
        mem_r2(i,m)=r(1,2); mem_p2(i,m)=p(1,2)/2;
        
        met=metrics(:,m);
        [r,p]=partialcorr([met,acer],age,'rows','complete');
        acer_par_r(i,m)=r(1,2); acer_par_p(i,m)=p(1,2)/2; % 1 tailed
        [r,p]=partialcorr([met,mmse],age,'rows','complete');
        mmse_par_r(i,m)=r(1,2); mmse_par_p(i,m)=p(1,2)/2;
        [r,p]=partialcorr([met,gm],age,'rows','complete');
        gm_par_r(i,m)=r(1,2); gm_par_p(i,m)=p(1,2)/2;
        [r,p]=partialcorr([met,wm],age,'rows','complete');
        wm_par_r(i,m)=r(1,2); wm_par_p(i,m)=p(1,2)/2;
        [r,p]=partialcorr([met,hip],age,'rows','complete');
        hip_par_r(i,m)=r(1,2); hip_par_p(i,m)=p(1,2)/2;
        [r,p]=partialcorr([met,mem],age,'rows','complete');
        mem_par_r(i,m)=r(1,2); mem_par_p(i,m)=p(1,2)/2;
  
    end
            
%   
%     T=table(met_names',c_m',c_sd',p_m',p_sd',hyp',tval',degf',pval',hyp_r',pval_r',hyp_l',pval_l',...
%         'VariableNames',{'TEST','C_M','C_SD','P_M','P_SD','H','T','DF','P','H_CP','P_CP','H_PC','P_PC'}) 
%     save([outdir 'roving_stats_' session '_' strrep(varnm{i},' ','_') '_' testname '.mat'],'T')
%     writetable(T,[outdir 'roving_stats_' session '_' strrep(varnm{i},' ','_') '_' testname '.xlsx'])
    
%% GLMs to correct for age 

for m=1:size(metrics,2)
    
    X=[];
    X(c_ind,1)=1;
    X(p_ind,2)=1;
    X(:,3)=zscore(age);
    %X(:,4)=zscore(edu);
    met=metrics(:,m); idx=find(isnan(met)); 
    met(idx)=[]; X(idx,:)=[]; %met=zscore(met); 
    met=zscore(met); 
    idx=find(isoutlier(met));met(idx)=[];X(idx,:)=[];
    [t,F,p,df,R2,cR2,B] = glm_rik(met,X,[1 -1 0],1);
    tot_t_glm2(i,m)=t; tot_p_glm2(i,m)=p;
    
end

%% AUCs & SVMs

for m=1:size(metrics,2)
    
    cmet=metrics(c_ind,m);
    idx=find(isoutlier(cmet));cmet(idx)=[];
    pmet=metrics(p_ind,m);
    idx=find(isoutlier(pmet));pmet(idx)=[];

    param=roc_curve(cmet,pmet);
    auc(i,m)=param.param.AROC;
    
    clear A
    cmet=metrics(c_ind,m);pmet=metrics(p_ind,m); cl(c_ind,1)=1;cl(p_ind,1)=0;
    A(:,1)=[cmet;pmet]; A(:,2)=cl; idx=find(isnan(A(:,1))); A(idx,:)=[];
    acc=repeated_CV(A,5,250);
    svmcl(i,m)=mean(acc(:));
end

end

%% Plot summary t values

figure('color','w');
imagesc(squeeze(tot_t)'); axis equal off
caxis([-2 2]); axis on;
xlim([0.5 9.5]);yticks(1:length(met_names));ylim([0.5 length(met_names)+0.5])
set(gca,'YTickLabel',varnames)
set(gca,'XTick',1:9,'XTickLabel',{'MAG','Fr-MAG','Tem-MAG','GRAD','Fr-GRAD','Tem-GRAD','EEG','Fr-EEG','Tem-EEG'})
xtickangle(30); colorbar; colormap(redblue)
title(['Controls-Patients ' window_nm{w}])

[a,b]=ind2sub(size(squeeze(tot_p1')),find(squeeze(tot_p1'<=0.05))); hold on;
for d=1:length(a)
    plot(b(d),a(d),'-s','markersize',8,'markeredgecolor','w','markerfacecolor','w')
end

[a,b]=ind2sub(size(squeeze(tot_p1')),find(squeeze(tot_p1'<=0.01))); hold on;
for d=1:length(a)
    plot(b(d),a(d),'-s','markersize',12,'markeredgecolor','w','markerfacecolor','w')
end

[a,b]=ind2sub(size(squeeze(tot_p2')),find(squeeze(tot_p2'<=0.05))); hold on;
for d=1:length(a)
    plot(b(d),a(d),'-s','markersize',8,'markeredgecolor','w','markerfacecolor','w')
end

[a,b]=ind2sub(size(squeeze(tot_p2')),find(squeeze(tot_p2'<=0.01))); hold on;
for d=1:length(a)
    plot(b(d),a(d),'-s','markersize',12,'markeredgecolor','w','markerfacecolor','w')
end
print(gcf,'-dbmp', '-r0',[outdir '/T_tests/' 't_map_' session  '_controls-patients_' testname '_' window_nm{w} '.bmp'],'-r200');% close(gcf)

%% Plot GLM results

figure('color','w');
imagesc(squeeze(tot_t_glm1)'); axis equal off
caxis([-2 2]); axis on;
xlim([0.5 9.5]);yticks(1:length(met_names)); ylim([0.5 length(met_names)+0.5])
set(gca,'YTickLabel',varnames)
set(gca,'XTick',1:9,'XTickLabel',{'MAG','Fr-MAG','Tem-MAG','GRAD','Fr-GRAD','Tem-GRAD','EEG','Fr-EEG','Tem-EEG'})
xtickangle(30); colorbar; colormap(redblue)
title(['Controls-Patients ' window_nm{w} ' GLM'])

[a,b]=ind2sub(size(squeeze(tot_p_glm1')),find(squeeze(tot_p_glm1'<=0.05))); hold on;
for d=1:length(a)
    plot(b(d),a(d),'-s','markersize',8,'markeredgecolor','w','markerfacecolor','w')
end

[a,b]=ind2sub(size(squeeze(tot_p_glm1')),find(squeeze(tot_p_glm1'<=0.01))); hold on;
for d=1:length(a)
    plot(b(d),a(d),'-s','markersize',12,'markeredgecolor','w','markerfacecolor','w')
end
print(gcf,'-dbmp', '-r0',[outdir '/GLMs/' 't_map_' session  '_controls-patients_GLM_' testname '_' window_nm{w} '.bmp'],'-r200');% close(gcf)


figure('color','w');
imagesc(squeeze(tot_t_glm2)'); axis equal off
caxis([-2 2]); axis on;
xlim([0.5 9.5]);yticks(1:length(met_names)); ylim([0.5 length(met_names)+0.5])
set(gca,'YTickLabel',varnames)
set(gca,'XTick',1:9,'XTickLabel',{'MAG','Fr-MAG','Tem-MAG','GRAD','Fr-GRAD','Tem-GRAD','EEG','Fr-EEG','Tem-EEG'})
xtickangle(30); colorbar; colormap(redblue)
title(['Controls-Patients ' window_nm{w} ' GLM age cor'])

[a,b]=ind2sub(size(squeeze(tot_p_glm2')),find(squeeze(tot_p_glm2'<=0.05))); hold on;
for d=1:length(a)
    plot(b(d),a(d),'-s','markersize',8,'markeredgecolor','w','markerfacecolor','w')
end

[a,b]=ind2sub(size(squeeze(tot_p_glm2')),find(squeeze(tot_p_glm2'<=0.01))); hold on;
for d=1:length(a)
    plot(b(d),a(d),'-s','markersize',12,'markeredgecolor','w','markerfacecolor','w')
end
print(gcf,'-dbmp', '-r0',[outdir '/GLMs/' 't_map_' session  '_controls-patients_GLM_age_cor_' testname '_' window_nm{w} '.bmp'],'-r200'); %close(gcf)

%% Plot AUCs
cmap=customcolormap([0 0.9 1],[29/255 180/255 191/255 ;1 1 1;1 1 1 ]);
idx=find(auc<0.5);
auc(idx)=abs(1-auc(idx));

figure('color','w');
imagesc(auc'); axis equal off
caxis([0.5 0.8]); axis on;
xlim([0.5 9.5]);yticks(1:length(met_names));ylim([0.5 length(met_names)+0.5])
set(gca,'YTickLabel',varnames)
set(gca,'XTick',1:9,'XTickLabel',{'MAG','Fr-MAG','Tem-MAG','GRAD','Fr-GRAD','Tem-GRAD','EEG','Fr-EEG','Tem-EEG'})
xtickangle(30); colorbar; colormap(cmap)
title(['AUC ' window_nm{w}])

[a,b]=ind2sub(size(auc),find(auc>=0.7)); hold on;
for d=1:length(a)
    plot(a(d),b(d),'-s','markersize',8,'markeredgecolor','k','markerfacecolor','k')
end

print(gcf,'-dbmp', '-r0',[outdir '/AUC/' 'AUC_' session  '_baseline_' testname '_' window_nm{w} '.bmp'],'-r200'); close(gcf)

%% Plot correlations: Patients only

figure('color','w');
imagesc(squeeze(mmse_r)'); axis equal off
caxis([-0.2 0.2]); axis on;
xlim([0.5 9.5]);yticks(1:length(met_names));ylim([0.5 length(met_names)+0.5])
set(gca,'YTickLabel',varnames)
set(gca,'XTick',1:9,'XTickLabel',{'MAG','Fr-MAG','Tem-MAG','GRAD','Fr-GRAD','Tem-GRAD','EEG','Fr-EEG','Tem-EEG'})
xtickangle(30); colorbar; colormap(redblue)
title(['Correlations with MMSE in ' window_nm{w}])

[a,b]=ind2sub(size(squeeze(mmse_p')),find(squeeze(mmse_p'<=0.05))); hold on;
for d=1:length(a)
    plot(b(d),a(d),'-s','markersize',8,'markeredgecolor','w','markerfacecolor','w')
end

[a,b]=ind2sub(size(squeeze(mmse_p')),find(squeeze(mmse_p'<=0.01))); hold on;
for d=1:length(a)
    plot(b(d),a(d),'-s','markersize',12,'markeredgecolor','w','markerfacecolor','w')
end

print(gcf,'-dbmp', '-r0',[outdir '/Corr/' 'mmse_corr_' session  '_patients_' testname '_' window_nm{w} '.bmp'],'-r200'); close(gcf)

figure('color','w');
imagesc(squeeze(acer_r)'); axis equal off
caxis([-0.2 0.2]); axis on;
xlim([0.5 9.5]);yticks(1:length(met_names));ylim([0.5 length(met_names)+0.5])
set(gca,'YTickLabel',varnames)
set(gca,'XTick',1:9,'XTickLabel',{'MAG','Fr-MAG','Tem-MAG','GRAD','Fr-GRAD','Tem-GRAD','EEG','Fr-EEG','Tem-EEG'})
xtickangle(30); colorbar; colormap(redblue)
title(['Correlations with ACER in ' window_nm{w}])

[a,b]=ind2sub(size(squeeze(acer_p')),find(squeeze(acer_p'<=0.05))); hold on;
for d=1:length(a)
    plot(b(d),a(d),'-s','markersize',8,'markeredgecolor','w','markerfacecolor','w')
end

[a,b]=ind2sub(size(squeeze(acer_p')),find(squeeze(acer_p'<=0.01))); hold on;
for d=1:length(a)
    plot(b(d),a(d),'-s','markersize',12,'markeredgecolor','w','markerfacecolor','w')
end

print(gcf,'-dbmp', '-r0',[outdir '/Corr/' 'acer_corr_' session  '_patients_' testname '_' window_nm{w} '.bmp'],'-r200'); close(gcf)

figure('color','w');
imagesc(squeeze(gm_r)'); axis equal off
caxis([-0.2 0.2]); axis on;
xlim([0.5 9.5]);yticks(1:length(met_names));ylim([0.5 length(met_names)+0.5])
set(gca,'YTickLabel',varnames)
set(gca,'XTick',1:9,'XTickLabel',{'MAG','Fr-MAG','Tem-MAG','GRAD','Fr-GRAD','Tem-GRAD','EEG','Fr-EEG','Tem-EEG'})
xtickangle(30); colorbar; colormap(redblue)
title(['Correlations with GMV in ' window_nm{w}])

[a,b]=ind2sub(size(squeeze(gm_p')),find(squeeze(gm_p'<=0.05))); hold on;
for d=1:length(a)
    plot(b(d),a(d),'-s','markersize',8,'markeredgecolor','w','markerfacecolor','w')
end

[a,b]=ind2sub(size(squeeze(gm_p')),find(squeeze(gm_p'<=0.01))); hold on;
for d=1:length(a)
    plot(b(d),a(d),'-s','markersize',12,'markeredgecolor','w','markerfacecolor','w')
end

print(gcf,'-dbmp', '-r0',[outdir '/Corr/' 'gmv_corr_' session  '_patients_' testname '_' window_nm{w} '.bmp'],'-r200'); close(gcf)

figure('color','w');
imagesc(squeeze(wm_r)'); axis equal off
caxis([-0.2 0.2]); axis on;
xlim([0.5 9.5]);yticks(1:length(met_names));ylim([0.5 length(met_names)+0.5])
set(gca,'YTickLabel',varnames)
set(gca,'XTick',1:9,'XTickLabel',{'MAG','Fr-MAG','Tem-MAG','GRAD','Fr-GRAD','Tem-GRAD','EEG','Fr-EEG','Tem-EEG'})
xtickangle(30); colorbar; colormap(redblue)
title(['Correlations with WMV in ' window_nm{w}])

[a,b]=ind2sub(size(squeeze(wm_p')),find(squeeze(wm_p'<=0.05))); hold on;
for d=1:length(a)
    plot(b(d),a(d),'-s','markersize',8,'markeredgecolor','w','markerfacecolor','w')
end

[a,b]=ind2sub(size(squeeze(wm_p')),find(squeeze(wm_p'<=0.01))); hold on;
for d=1:length(a)
    plot(b(d),a(d),'-s','markersize',12,'markeredgecolor','w','markerfacecolor','w')
end

print(gcf,'-dbmp', '-r0',[outdir '/Corr/' 'wmv_corr_' session  '_patients_' testname '_' window_nm{w} '.bmp'],'-r200'); close(gcf)

figure('color','w');
imagesc(squeeze(hip_r)'); axis equal off
caxis([-0.2 0.2]); axis on;
xlim([0.5 9.5]);yticks(1:length(met_names));ylim([0.5 length(met_names)+0.5])
set(gca,'YTickLabel',varnames)
set(gca,'XTick',1:9,'XTickLabel',{'MAG','Fr-MAG','Tem-MAG','GRAD','Fr-GRAD','Tem-GRAD','EEG','Fr-EEG','Tem-EEG'})
xtickangle(30); colorbar; colormap(redblue)
title(['Correlations with HIP in ' window_nm{w}])

[a,b]=ind2sub(size(squeeze(hip_p')),find(squeeze(hip_p'<=0.05))); hold on;
for d=1:length(a)
    plot(b(d),a(d),'-s','markersize',8,'markeredgecolor','w','markerfacecolor','w')
end

[a,b]=ind2sub(size(squeeze(hip_p')),find(squeeze(hip_p'<=0.01))); hold on;
for d=1:length(a)
    plot(b(d),a(d),'-s','markersize',12,'markeredgecolor','w','markerfacecolor','w')
end

print(gcf,'-dbmp', '-r0',[outdir '/Corr/' 'hip_corr_' session  '_patients_' testname '_' window_nm{w} '.bmp'],'-r200'); close(gcf)

%% Plot correlations: Whole sample

figure('color','w');
imagesc(squeeze(mmse_r2)'); axis equal off
caxis([-0.4 0.4]); axis on;
xlim([0.5 9.5]);yticks(1:length(met_names));ylim([0.5 length(met_names)+0.5])
set(gca,'YTickLabel',varnames)
set(gca,'XTick',1:9,'XTickLabel',{'MAG','Fr-MAG','Tem-MAG','GRAD','Fr-GRAD','Tem-GRAD','EEG','Fr-EEG','Tem-EEG'})
xtickangle(30); colorbar; colormap(redblue)
title(['Correlations with MMSE in ' window_nm{w}])

[a,b]=ind2sub(size(squeeze(mmse_p2')),find(squeeze(mmse_p2'<=0.05))); hold on;
for d=1:length(a)
    plot(b(d),a(d),'-s','markersize',8,'markeredgecolor','w','markerfacecolor','w')
end

[a,b]=ind2sub(size(squeeze(mmse_p2')),find(squeeze(mmse_p2'<=0.01))); hold on;
for d=1:length(a)
    plot(b(d),a(d),'-s','markersize',12,'markeredgecolor','w','markerfacecolor','w')
end

print(gcf,'-dbmp', '-r0',[outdir '/Corr/' 'mmse_corr_' session  '_whole_sample_' testname '_' window_nm{w} '.bmp'],'-r200'); close(gcf)

figure('color','w');
imagesc(squeeze(acer_r2)'); axis equal off
caxis([-0.4 0.4]); axis on;
xlim([0.5 9.5]);yticks(1:length(met_names));ylim([0.5 length(met_names)+0.5])
set(gca,'YTickLabel',varnames)
set(gca,'XTick',1:9,'XTickLabel',{'MAG','Fr-MAG','Tem-MAG','GRAD','Fr-GRAD','Tem-GRAD','EEG','Fr-EEG','Tem-EEG'})
xtickangle(30); colorbar; colormap(redblue)
title(['Correlations with ACER in ' window_nm{w}])

[a,b]=ind2sub(size(squeeze(acer_p2')),find(squeeze(acer_p2'<=0.05))); hold on;
for d=1:length(a)
    plot(b(d),a(d),'-s','markersize',8,'markeredgecolor','w','markerfacecolor','w')
end

[a,b]=ind2sub(size(squeeze(acer_p2')),find(squeeze(acer_p2'<=0.01))); hold on;
for d=1:length(a)
    plot(b(d),a(d),'-s','markersize',12,'markeredgecolor','w','markerfacecolor','w')
end

print(gcf,'-dbmp', '-r0',[outdir '/Corr/' 'acer_corr_' session  '_whole_sample_' testname '_' window_nm{w} '.bmp'],'-r200'); close(gcf)

figure('color','w');
imagesc(squeeze(gm_r2)'); axis equal off
caxis([-0.4 0.4]); axis on;
xlim([0.5 9.5]);yticks(1:length(met_names));ylim([0.5 length(met_names)+0.5])
set(gca,'YTickLabel',varnames)
set(gca,'XTick',1:9,'XTickLabel',{'MAG','Fr-MAG','Tem-MAG','GRAD','Fr-GRAD','Tem-GRAD','EEG','Fr-EEG','Tem-EEG'})
xtickangle(30); colorbar; colormap(redblue)
title(['Correlations with GMV in ' window_nm{w}])

[a,b]=ind2sub(size(squeeze(gm_p2')),find(squeeze(gm_p2'<=0.05))); hold on;
for d=1:length(a)
    plot(b(d),a(d),'-s','markersize',8,'markeredgecolor','w','markerfacecolor','w')
end

[a,b]=ind2sub(size(squeeze(gm_p2')),find(squeeze(gm_p2'<=0.01))); hold on;
for d=1:length(a)
    plot(b(d),a(d),'-s','markersize',12,'markeredgecolor','w','markerfacecolor','w')
end

print(gcf,'-dbmp', '-r0',[outdir '/Corr/' 'gmv_corr_' session  '_whole_sample_' testname '_' window_nm{w} '.bmp'],'-r200'); close(gcf)

figure('color','w');
imagesc(squeeze(wm_r2)'); axis equal off
caxis([-0.4 0.4]); axis on;
xlim([0.5 9.5]);yticks(1:length(met_names));ylim([0.5 length(met_names)+0.5])
set(gca,'YTickLabel',varnames)
set(gca,'XTick',1:9,'XTickLabel',{'MAG','Fr-MAG','Tem-MAG','GRAD','Fr-GRAD','Tem-GRAD','EEG','Fr-EEG','Tem-EEG'})
xtickangle(30); colorbar; colormap(redblue)
title(['Correlations with WMV in ' window_nm{w}])

[a,b]=ind2sub(size(squeeze(wm_p2')),find(squeeze(wm_p2'<=0.05))); hold on;
for d=1:length(a)
    plot(b(d),a(d),'-s','markersize',8,'markeredgecolor','w','markerfacecolor','w')
end

[a,b]=ind2sub(size(squeeze(wm_p2')),find(squeeze(wm_p2'<=0.01))); hold on;
for d=1:length(a)
    plot(b(d),a(d),'-s','markersize',12,'markeredgecolor','w','markerfacecolor','w')
end

print(gcf,'-dbmp', '-r0',[outdir '/Corr/' 'wmv_corr_' session  '_whole_sample_' testname '_' window_nm{w} '.bmp'],'-r200'); close(gcf)

figure('color','w');
imagesc(squeeze(hip_r2)'); axis equal off
caxis([-0.4 0.4]); axis on;
xlim([0.5 9.5]);yticks(1:length(met_names));ylim([0.5 length(met_names)+0.5])
set(gca,'YTickLabel',varnames)
set(gca,'XTick',1:9,'XTickLabel',{'MAG','Fr-MAG','Tem-MAG','GRAD','Fr-GRAD','Tem-GRAD','EEG','Fr-EEG','Tem-EEG'})
xtickangle(30); colorbar; colormap(redblue)
title(['Correlations with HIP in ' window_nm{w}])

[a,b]=ind2sub(size(squeeze(hip_p2')),find(squeeze(hip_p2'<=0.05))); hold on;
for d=1:length(a)
    plot(b(d),a(d),'-s','markersize',8,'markeredgecolor','w','markerfacecolor','w')
end

[a,b]=ind2sub(size(squeeze(hip_p2')),find(squeeze(hip_p2'<=0.01))); hold on;
for d=1:length(a)
    plot(b(d),a(d),'-s','markersize',12,'markeredgecolor','w','markerfacecolor','w')
end

print(gcf,'-dbmp', '-r0',[outdir '/Corr/' 'hip_corr_' session  '_whole_sample_' testname '_' window_nm{w} '.bmp'],'-r200'); close(gcf)

figure('color','w');
imagesc(squeeze(mem_r2)'); axis equal off
caxis([-0.4 0.4]); axis on;
xlim([0.5 9.5]);yticks(1:length(met_names));ylim([0.5 length(met_names)+0.5])
set(gca,'YTickLabel',varnames)
set(gca,'XTick',1:9,'XTickLabel',{'MAG','Fr-MAG','Tem-MAG','GRAD','Fr-GRAD','Tem-GRAD','EEG','Fr-EEG','Tem-EEG'})
xtickangle(30); colorbar; colormap(redblue)
title(['Correlations with ACER MEM in ' window_nm{w}])

[a,b]=ind2sub(size(squeeze(mem_p2')),find(squeeze(mem_p2'<=0.05))); hold on;
for d=1:length(a)
    plot(b(d),a(d),'-s','markersize',8,'markeredgecolor','w','markerfacecolor','w')
end

[a,b]=ind2sub(size(squeeze(mem_p2')),find(squeeze(mem_p2'<=0.01))); hold on;
for d=1:length(a)
    plot(b(d),a(d),'-s','markersize',12,'markeredgecolor','w','markerfacecolor','w')
end

print(gcf,'-dbmp', '-r0',[outdir '/Corr/' 'acer_mem_corr_' session  '_whole_sample_' testname '_' window_nm{w} '.bmp'],'-r200'); close(gcf)

%% Plot partial correlations

figure('color','w');
imagesc(squeeze(mmse_par_r)'); axis equal off
caxis([-0.4 0.4]); axis on;
xlim([0.5 9.5]);yticks(1:length(met_names));ylim([0.5 length(met_names)+0.5])
set(gca,'YTickLabel',varnames)
set(gca,'XTick',1:9,'XTickLabel',{'MAG','Fr-MAG','Tem-MAG','GRAD','Fr-GRAD','Tem-GRAD','EEG','Fr-EEG','Tem-EEG'})
xtickangle(30); colorbar; colormap(redblue)
title(['Partial correlations with MMSE in ' window_nm{w}])

[a,b]=ind2sub(size(squeeze(mmse_par_p')),find(squeeze(mmse_par_p'<=0.05))); hold on;
for d=1:length(a)
    plot(b(d),a(d),'-s','markersize',8,'markeredgecolor','w','markerfacecolor','w')
end

[a,b]=ind2sub(size(squeeze(mmse_par_p')),find(squeeze(mmse_par_p'<=0.01))); hold on;
for d=1:length(a)
    plot(b(d),a(d),'-s','markersize',12,'markeredgecolor','w','markerfacecolor','w')
end

print(gcf,'-dbmp', '-r0',[outdir '/Par_corr/' 'mmse_par_corr_' session  '_controls-patients_' testname '_' window_nm{w} '.bmp'],'-r200'); close(gcf)

figure('color','w');
imagesc(squeeze(acer_par_r)'); axis equal off
caxis([-0.4 0.4]); axis on;
xlim([0.5 9.5]);yticks(1:length(met_names));ylim([0.5 length(met_names)+0.5])
set(gca,'YTickLabel',varnames)
set(gca,'XTick',1:9,'XTickLabel',{'MAG','Fr-MAG','Tem-MAG','GRAD','Fr-GRAD','Tem-GRAD','EEG','Fr-EEG','Tem-EEG'})
xtickangle(30); colorbar; colormap(redblue)
title(['Partial correlations with ACER in ' window_nm{w}])

[a,b]=ind2sub(size(squeeze(acer_par_p')),find(squeeze(acer_par_p'<=0.05))); hold on;
for d=1:length(a)
    plot(b(d),a(d),'-s','markersize',8,'markeredgecolor','w','markerfacecolor','w')
end

[a,b]=ind2sub(size(squeeze(acer_par_p')),find(squeeze(acer_par_p'<=0.01))); hold on;
for d=1:length(a)
    plot(b(d),a(d),'-s','markersize',12,'markeredgecolor','w','markerfacecolor','w')
end

print(gcf,'-dbmp', '-r0',[outdir '/Par_corr/' 'acer_par_corr_' session  '_controls-patients_' testname '_' window_nm{w} '.bmp'],'-r200'); close(gcf)

figure('color','w');
imagesc(squeeze(gm_par_r)'); axis equal off
caxis([-0.4 0.4]); axis on;
xlim([0.5 9.5]);yticks(1:length(met_names));ylim([0.5 length(met_names)+0.5])
set(gca,'YTickLabel',varnames)
set(gca,'XTick',1:9,'XTickLabel',{'MAG','Fr-MAG','Tem-MAG','GRAD','Fr-GRAD','Tem-GRAD','EEG','Fr-EEG','Tem-EEG'})
xtickangle(30); colorbar; colormap(redblue)
title(['Partial correlations with GMV in ' window_nm{w}])

[a,b]=ind2sub(size(squeeze(gm_par_p')),find(squeeze(gm_par_p'<=0.05))); hold on;
for d=1:length(a)
    plot(b(d),a(d),'-s','markersize',8,'markeredgecolor','w','markerfacecolor','w')
end

[a,b]=ind2sub(size(squeeze(gm_par_p')),find(squeeze(gm_par_p'<=0.01))); hold on;
for d=1:length(a)
    plot(b(d),a(d),'-s','markersize',12,'markeredgecolor','w','markerfacecolor','w')
end

print(gcf,'-dbmp', '-r0',[outdir '/Par_corr/' 'gmv_par_corr_' session  '_controls-patients_' testname '_' window_nm{w} '.bmp'],'-r200'); close(gcf)

figure('color','w');
imagesc(squeeze(wm_par_r)'); axis equal off
caxis([-0.4 0.4]); axis on;
xlim([0.5 9.5]);yticks(1:length(met_names));ylim([0.5 length(met_names)+0.5])
set(gca,'YTickLabel',varnames)
set(gca,'XTick',1:9,'XTickLabel',{'MAG','Fr-MAG','Tem-MAG','GRAD','Fr-GRAD','Tem-GRAD','EEG','Fr-EEG','Tem-EEG'})
xtickangle(30); colorbar; colormap(redblue)
title(['Partial correlations with WMV in ' window_nm{w}])

[a,b]=ind2sub(size(squeeze(wm_par_p')),find(squeeze(wm_par_p'<=0.05))); hold on;
for d=1:length(a)
    plot(b(d),a(d),'-s','markersize',8,'markeredgecolor','w','markerfacecolor','w')
end

[a,b]=ind2sub(size(squeeze(wm_par_p')),find(squeeze(wm_par_p'<=0.01))); hold on;
for d=1:length(a)
    plot(b(d),a(d),'-s','markersize',12,'markeredgecolor','w','markerfacecolor','w')
end

print(gcf,'-dbmp', '-r0',[outdir '/Par_corr/' 'wmv_par_corr_' session  '_controls-patients_' testname '_' window_nm{w} '.bmp'],'-r200'); close(gcf)

figure('color','w');
imagesc(squeeze(hip_par_r)'); axis equal off
caxis([-0.4 0.4]); axis on;
xlim([0.5 9.5]);yticks(1:length(met_names));ylim([0.5 length(met_names)+0.5])
set(gca,'YTickLabel',varnames)
set(gca,'XTick',1:9,'XTickLabel',{'MAG','Fr-MAG','Tem-MAG','GRAD','Fr-GRAD','Tem-GRAD','EEG','Fr-EEG','Tem-EEG'})
xtickangle(30); colorbar; colormap(redblue)
title(['Partial correlations with HIP in ' window_nm{w}])

[a,b]=ind2sub(size(squeeze(hip_par_p')),find(squeeze(hip_par_p'<=0.05))); hold on;
for d=1:length(a)
    plot(b(d),a(d),'-s','markersize',8,'markeredgecolor','w','markerfacecolor','w')
end

[a,b]=ind2sub(size(squeeze(hip_par_p')),find(squeeze(hip_par_p'<=0.01))); hold on;
for d=1:length(a)
    plot(b(d),a(d),'-s','markersize',12,'markeredgecolor','w','markerfacecolor','w')
end

print(gcf,'-dbmp', '-r0',[outdir '/Par_corr/' 'hip_par_corr_' session  '_controls-patients_' testname '_' window_nm{w} '.bmp'],'-r200'); close(gcf)

figure('color','w');
imagesc(squeeze(mem_par_r)'); axis equal off
caxis([-0.4 0.4]); axis on;
xlim([0.5 9.5]);yticks(1:length(met_names));ylim([0.5 length(met_names)+0.5])
set(gca,'YTickLabel',varnames)
set(gca,'XTick',1:9,'XTickLabel',{'MAG','Fr-MAG','Tem-MAG','GRAD','Fr-GRAD','Tem-GRAD','EEG','Fr-EEG','Tem-EEG'})
xtickangle(30); colorbar; colormap(redblue)
title(['Partial correlations with ACER mem in ' window_nm{w}])

[a,b]=ind2sub(size(squeeze(mem_par_p')),find(squeeze(mem_par_p'<=0.05))); hold on;
for d=1:length(a)
    plot(b(d),a(d),'-s','markersize',8,'markeredgecolor','w','markerfacecolor','w')
end

[a,b]=ind2sub(size(squeeze(mem_par_p')),find(squeeze(mem_par_p'<=0.01))); hold on;
for d=1:length(a)
    plot(b(d),a(d),'-s','markersize',12,'markeredgecolor','w','markerfacecolor','w')
end

print(gcf,'-dbmp', '-r0',[outdir '/Par_corr/' 'acer_mem_par_corr_' session  '_controls-patients_' testname '_' window_nm{w} '.bmp'],'-r200'); close(gcf)

%% Plot kernel density estimations
sens={'MAGs','Fr MAGs', 'Tem MAGs','GRADs','Fr GRADs', 'Tem GRADs','EEGs','Fr EEGs','Tem EEGs'};
% 
% for m=1:size(met_by_sens,2) %tem grads
%     mets=met_by_sens{m};
%     
%     figure('Color','w');set(gcf, 'Units','normalized', 'Position', [0 0 .3 .8]);
%     for j=1:size(mets,2)
%         
%         c=mets(c_ind,j);
%         p=mets(p_ind,j);
%         
%         if ismember(m,[7 8 9]); bw=0.25; else; bw=0.04; end
%         
%         %subplot(size(mets,2),1,j)
%         
%         [a,b,bw]=ksdensity(p,'Bandwidth',bw);
%         h2=area(b,a);
%         set(h2, 'FaceColor', [170 170 170]./255,'EdgeColor',[140 140 140]./255 , 'FaceAlpha',0.6,'LineWidth',2);
%         hold on;
%         pm=nanmean(p); ix=find(b>pm);ix=ix(1);
%         line([pm pm],[0 a(ix)],'Color',[120 120 120]./255,'LineWidth',2.1,'LineStyle','--')
%         [a,b]=ksdensity(c,'Bandwidth',bw);
%         h1=area(b,a);
%         set(h1, 'FaceColor', [247 187 5]./255, 'EdgeColor', [247 187 5]./255, 'FaceAlpha',0.6,'LineWidth',2);
%         cm=nanmean(c); ix=find(b>cm);ix=ix(1);
%         line([cm cm],[0 a(ix)],'Color',[219 146 0]./255,'LineWidth',2.1,'LineStyle','--'), box off
%         ylabel(sens{j})
%         
%         if j==1; title(varnames{m},'fontsize',12); end
%         
%     end
%     print(gcf,'-dbmp', '-r0',[outdir 'Group_kernel_density_' session  '_' strrep(varnames{m},'/','_by_') '_' testname '_' window_nm{w} '.bmp'],'-r200');% close(gcf)
% 
% end
%end

%% Plot ACER MMSE GM distributions

% c=acer(c_ind);
% p=acer(p_ind);

% 
% figure('Color','w');
% [a,b,bw]=ksdensity(p,'Bandwidth',3);
% h2=area(b,a);
% set(h2, 'FaceColor', [170 170 170]./255,'EdgeColor',[140 140 140]./255 , 'FaceAlpha',0.6,'LineWidth',2);
% hold on;
% pm=nanmean(p); ix=find(b>pm);ix=ix(1);
% line([pm pm],[0 a(ix)],'Color',[120 120 120]./255,'LineWidth',2.1,'LineStyle','--')
% [a,b]=ksdensity(c,'Bandwidth',3);
% h1=area(b,a);
% set(h1, 'FaceColor', [247 187 5]./255, 'EdgeColor', [247 187 5]./255, 'FaceAlpha',0.6,'LineWidth',2);
% cm=nanmean(c); ix=find(b>cm);ix=ix(1);
% line([cm cm],[0 a(ix)],'Color',[219 146 0]./255,'LineWidth',2.1,'LineStyle','--'), box off
% title('ACER','fontsize',12);
% print(gcf,'-dbmp', '-r0',[outdir 'Group_kernel_density_' session  '_ACER_' testname '_' window_nm{w} '.bmp'],'-r200'); close(gcf)
% 
% c=mmse(c_ind);
% p=mmse(p_ind);
% 
% 
% figure('Color','w');
% [a,b,bw]=ksdensity(p,'Bandwidth',1);
% h2=area(b,a);
% set(h2, 'FaceColor', [170 170 170]./255,'EdgeColor',[140 140 140]./255 , 'FaceAlpha',0.6,'LineWidth',2);
% hold on;
% pm=nanmean(p); ix=find(b>pm);ix=ix(1);
% line([pm pm],[0 a(ix)],'Color',[120 120 120]./255,'LineWidth',2.1,'LineStyle','--')
% [a,b]=ksdensity(c,'Bandwidth',1);
% h1=area(b,a);
% set(h1, 'FaceColor', [247 187 5]./255, 'EdgeColor', [247 187 5]./255, 'FaceAlpha',0.6,'LineWidth',2);
% cm=nanmean(c); ix=find(b>cm);ix=ix(1);
% line([cm cm],[0 a(ix)],'Color',[219 146 0]./255,'LineWidth',2.1,'LineStyle','--'), box off
% title('MMSE','fontsize',12);
% print(gcf,'-dbmp', '-r0',[outdir 'Group_kernel_density_' session  '_MMSE_' testname '_' window_nm{w} '.bmp'],'-r200'); close(gcf)
% 
% c=gm(c_ind);
% p=gm(p_ind);
% 
% figure('Color','w');
% [a,b,bw]=ksdensity(p,'Bandwidth',0.6);
% h2=area(b,a);
% set(h2, 'FaceColor', [170 170 170]./255,'EdgeColor',[140 140 140]./255 , 'FaceAlpha',0.6,'LineWidth',2);
% hold on;
% pm=nanmean(p); ix=find(b>pm);ix=ix(1);
% line([pm pm],[0 a(ix)],'Color',[120 120 120]./255,'LineWidth',2.1,'LineStyle','--')
% [a,b]=ksdensity(c,'Bandwidth',0.6);
% h1=area(b,a);
% set(h1, 'FaceColor', [247 187 5]./255, 'EdgeColor', [247 187 5]./255, 'FaceAlpha',0.6,'LineWidth',2);
% cm=nanmean(c); ix=find(b>cm);ix=ix(1);
% line([cm cm],[0 a(ix)],'Color',[219 146 0]./255,'LineWidth',2.1,'LineStyle','--'), box off
% title('%GMV','fontsize',12);
% print(gcf,'-dbmp', '-r0',[outdir 'Group_kernel_density_' session  '_GMV_' testname '_' window_nm{w} '.bmp'],'-r200'); close(gcf)
% 
% c=wm(c_ind);
% p=wm(p_ind);
% 
% figure('Color','w');
% [a,b,bw]=ksdensity(p,'Bandwidth',0.6);
% h2=area(b,a);
% set(h2, 'FaceColor', [170 170 170]./255,'EdgeColor',[140 140 140]./255 , 'FaceAlpha',0.6,'LineWidth',2);
% hold on;
% pm=nanmean(p); ix=find(b>pm);ix=ix(1);
% line([pm pm],[0 a(ix)],'Color',[120 120 120]./255,'LineWidth',2.1,'LineStyle','--')
% [a,b]=ksdensity(c,'Bandwidth',0.6);
% h1=area(b,a);
% set(h1, 'FaceColor', [247 187 5]./255, 'EdgeColor', [247 187 5]./255, 'FaceAlpha',0.6,'LineWidth',2);
% cm=nanmean(c); ix=find(b>cm);ix=ix(1);
% line([cm cm],[0 a(ix)],'Color',[219 146 0]./255,'LineWidth',2.1,'LineStyle','--'), box off
% title('%WMV','fontsize',12);
% print(gcf,'-dbmp', '-r0',[outdir 'Group_kernel_density_' session  '_WMV_' testname '_' window_nm{w} '.bmp'],'-r200'); close(gcf)
% 
% c=hip(c_ind);
% p=hip(p_ind);
% 
% figure('Color','w');
% [a,b,bw]=ksdensity(p,'Bandwidth',0.02);
% h2=area(b,a);
% set(h2, 'FaceColor', [170 170 170]./255,'EdgeColor',[140 140 140]./255 , 'FaceAlpha',0.6,'LineWidth',2);
% hold on;
% pm=nanmean(p); ix=find(b>pm);ix=ix(1);
% line([pm pm],[0 a(ix)],'Color',[120 120 120]./255,'LineWidth',2.1,'LineStyle','--')
% [a,b]=ksdensity(c,'Bandwidth',0.02);
% h1=area(b,a);
% set(h1, 'FaceColor', [247 187 5]./255, 'EdgeColor', [247 187 5]./255, 'FaceAlpha',0.6,'LineWidth',2);
% cm=nanmean(c); ix=find(b>cm);ix=ix(1);
% line([cm cm],[0 a(ix)],'Color',[219 146 0]./255,'LineWidth',2.1,'LineStyle','--'), box off
% title('%HIP GMV','fontsize',12);
% print(gcf,'-dbmp', '-r0',[outdir 'Group_kernel_density_' session  '_HIP_GMV_' testname '_' window_nm{w} '.bmp'],'-r200'); close(gcf)

%% Plot Cognition - MMN metric GLMs

for m=1:size(met_by_sens,2) % metric
    
    mets=met_by_sens{m};
    
    for j=1:size(mets,2) %sensor
        
        mmn=mets(:,j);
        ace=acer;
        minim=mmse;
        p_age=age;
        
        X=[];
        X(:,1)=ace;
        X(:,2)=age;
        idx=unique([find(isnan(mmn)); find(isnan(ace)); find(isnan(age))]);
        mmn(idx)=[]; X(idx,:)=[]; %mmn=zscore(mmn); 
        met=zscore(mmn); X(:,1)=zscore(X(:,1)); X(:,2)=zscore(X(:,2));
        idx=find(isoutlier(mmn));mmn(idx)=[];X(idx,:)=[];
        [t,F,p,df,R2,cR2,B] = glm(mmn,X,[1 0],0);
        acer_mmn_cor_t_glm(m,j)=t; acer_mmn_cor_p_glm(m,j)=p;
        
        X=[];mmn=mets(:,j);
        X(:,1)=ace;
        %X(:,2)=age;
        idx=unique([find(isnan(mmn)); find(isnan(ace))]);
        mmn(idx)=[]; X(idx,:)=[]; %mmn=zscore(mmn); 
        met=zscore(mmn); X(:,1)=zscore(X(:,1));
        idx=find(isoutlier(mmn));mmn(idx)=[];X(idx,:)=[];
        [t,F,p,df,R2,cR2,B] = glm(mmn,X,[1],0);
        acer_mmn_t_glm(m,j)=t; acer_mmn_p_glm(m,j)=p;
    
        mmn=mets(:,j);p_age=age;
        X=[];
        X(:,1)=minim;
        X(:,2)=age;
        idx=unique([find(isnan(mmn)); find(isnan(minim)); find(isnan(age))]);
        mmn(idx)=[]; X(idx,:)=[]; %mmn=zscore(mmn); 
         met=zscore(mmn); X(:,1)=zscore(X(:,1)); X(:,2)=zscore(X(:,2));
        idx=find(isoutlier(mmn));mmn(idx)=[];X(idx,:)=[];
        [t,F,p,df,R2,cR2,B] = glm(mmn,X,[1 0],0);
        mmse_mmn_cor_t_glm(m,j)=t; mmse_mmn_cor_p_glm(m,j)=p;
        
        mmn=mets(:,j);p_age=age;
        X=[];
        X(:,1)=minim;
        %X(:,2)=age;
        idx=unique([find(isnan(mmn)) ; find(isnan(minim))]);
        mmn(idx)=[]; X(idx,:)=[]; %mmn=zscore(mmn); 
        met=zscore(mmn); X(:,1)=zscore(X(:,1));
        idx=find(isoutlier(mmn));mmn(idx)=[];X(idx,:)=[];
        [t,F,p,df,R2,cR2,B] = glm(mmn,X,[1],0);
        mmse_mmn_t_glm(m,j)=t; mmse_mmn_p_glm(m,j)=p;
        
    end
end

figure('color','w');
imagesc(squeeze(acer_mmn_cor_t_glm)); axis equal off
caxis([-3 3]); axis on;
xlim([0.5 9.5]);yticks(1:length(met_names)); ylim([0.5 length(met_names)+0.5])
set(gca,'YTickLabel',varnames)
set(gca,'XTick',1:9,'XTickLabel',{'MAG','Fr-MAG','Tem-MAG','GRAD','Fr-GRAD','Tem-GRAD','EEG','Fr-EEG','Tem-EEG'})
xtickangle(30); colorbar; colormap(redblue)
title(['ACER x MMN -Age ' window_nm{w}])

[a,b]=ind2sub(size(squeeze(acer_mmn_cor_p_glm)),find(squeeze(acer_mmn_cor_p_glm<=0.05))); hold on;
for d=1:length(a)
    plot(b(d),a(d),'-s','markersize',8,'markeredgecolor','w','markerfacecolor','w')
end

[a,b]=ind2sub(size(squeeze(acer_mmn_cor_p_glm)),find(squeeze(acer_mmn_cor_p_glm<=0.01))); hold on;
for d=1:length(a)
    plot(b(d),a(d),'-s','markersize',12,'markeredgecolor','w','markerfacecolor','w')
end
print(gcf,'-dbmp', '-r0',[outdir '/GLMs/' 't_map_' session  '_ACER_MMN_GLM_age-cor_' testname '_' window_nm{w} '.bmp'],'-r200'); close(gcf)

figure('color','w');
imagesc(squeeze(acer_mmn_t_glm)); axis equal off
caxis([-3 3]); axis on;
xlim([0.5 9.5]);yticks(1:length(met_names));ylim([0.5 length(met_names)+0.5])
set(gca,'YTickLabel',varnames)
set(gca,'XTick',1:9,'XTickLabel',{'MAG','Fr-MAG','Tem-MAG','GRAD','Fr-GRAD','Tem-GRAD','EEG','Fr-EEG','Tem-EEG'})
xtickangle(30); colorbar; colormap(redblue)
title(['ACER x MMN ' window_nm{w}])

[a,b]=ind2sub(size(squeeze(acer_mmn_p_glm)),find(squeeze(acer_mmn_p_glm<=0.05))); hold on;
for d=1:length(a)
    plot(b(d),a(d),'-s','markersize',8,'markeredgecolor','w','markerfacecolor','w')
end

[a,b]=ind2sub(size(squeeze(acer_mmn_p_glm)),find(squeeze(acer_mmn_p_glm<=0.01))); hold on;
for d=1:length(a)
    plot(b(d),a(d),'-s','markersize',12,'markeredgecolor','w','markerfacecolor','w')
end
print(gcf,'-dbmp', '-r0',[outdir '/GLMs/' 't_map_' session  '_ACER_MMN_GLM_' testname '_' window_nm{w} '.bmp'],'-r200'); close(gcf)

figure('color','w');
imagesc(squeeze(mmse_mmn_cor_t_glm)); axis equal off
caxis([-3 3]); axis on;
xlim([0.5 9.5]);yticks(1:length(met_names)); ylim([0.5 length(met_names)+0.5])
set(gca,'YTickLabel',varnames)
set(gca,'XTick',1:9,'XTickLabel',{'MAG','Fr-MAG','Tem-MAG','GRAD','Fr-GRAD','Tem-GRAD','EEG','Fr-EEG','Tem-EEG'})
xtickangle(30); colorbar; colormap(redblue)
title(['MMSE x MMN -Age ' window_nm{w}])

[a,b]=ind2sub(size(squeeze(mmse_mmn_cor_p_glm)),find(squeeze(mmse_mmn_cor_p_glm<=0.05))); hold on;
for d=1:length(a)
    plot(b(d),a(d),'-s','markersize',8,'markeredgecolor','w','markerfacecolor','w')
end

[a,b]=ind2sub(size(squeeze(mmse_mmn_cor_p_glm)),find(squeeze(mmse_mmn_cor_p_glm<=0.01))); hold on;
for d=1:length(a)
    plot(b(d),a(d),'-s','markersize',12,'markeredgecolor','w','markerfacecolor','w')
end
print(gcf,'-dbmp', '-r0',[outdir '/GLMs/' 't_map_' session  '_MMSE_MMN_GLM_age-cor_' testname '_' window_nm{w} '.bmp'],'-r200'); close(gcf)

figure('color','w');
imagesc(squeeze(mmse_mmn_t_glm)); axis equal off
caxis([-3 3 ]); axis on;
xlim([0.5 9.5]);yticks(1:length(met_names));ylim([0.5 length(met_names)+0.5])
set(gca,'YTickLabel',varnames)
set(gca,'XTick',1:9,'XTickLabel',{'MAG','Fr-MAG','Tem-MAG','GRAD','Fr-GRAD','Tem-GRAD','EEG','Fr-EEG','Tem-EEG'})
xtickangle(30); colorbar; colormap(redblue)
title(['MMSE x MMN ' window_nm{w}])

[a,b]=ind2sub(size(squeeze(mmse_mmn_p_glm)),find(squeeze(mmse_mmn_p_glm<=0.05))); hold on;
for d=1:length(a)
    plot(b(d),a(d),'-s','markersize',8,'markeredgecolor','w','markerfacecolor','w')
end

[a,b]=ind2sub(size(squeeze(mmse_mmn_p_glm)),find(squeeze(mmse_mmn_p_glm<=0.01))); hold on;
for d=1:length(a)
    plot(b(d),a(d),'-s','markersize',12,'markeredgecolor','w','markerfacecolor','w')
end
print(gcf,'-dbmp', '-r0',[outdir '/GLMs/' 't_map_' session  '_MMSE_MMN_GLM_' testname '_' window_nm{w} '.bmp'],'-r200'); close(gcf)

%% Plot MR x MMN relationship 

for m=1:size(met_by_sens,2) % metric
    
    mets=met_by_sens{m};
    
    for j=1:size(mets,2) %sensor
        
        mmn=mets(:,j);
        X=[];
        X(:,1)=gm;
        X(:,2)=age;
        idx=unique([find(isnan(mmn)); find(isnan(gm))]);
        mmn(idx)=[];X(idx,:)=[]; 
        %mmn=zscore(mmn); %X(:,1)=zscore(X(:,1)); X(:,2)=zscore(X(:,2));
        met=zscore(mmn);  X(:,1)=zscore(X(:,1)); X(:,2)=zscore(X(:,2)); 
        idx=find(isoutlier(mmn));mmn(idx)=[];X(idx,:)=[];
        [t,F,p,df,R2,cR2,B] = glm(mmn,X,[1 0],0);
        gm_mmn_cor_t_glm(m,j)=t; gm_mmn_cor_p_glm(m,j)=p;
        
        X=[];mmn=mets(:,j);
        X(:,1)=gm;
        idx=unique([find(isnan(mmn)); find(isnan(gm))]);
        mmn(idx)=[]; X(idx,:)=[]; 
        %X(:,1)=zscore(X(:,1));mmn=zscore(mmn); 
        met=zscore(mmn); X(:,1)=zscore(X(:,1)); 
        idx=find(isoutlier(mmn));mmn(idx)=[];X(idx,:)=[];
        [t,F,p,df,R2,cR2,B] = glm(mmn,X,[1],0);
        gm_mmn_t_glm(m,j)=t; gm_mmn_p_glm(m,j)=p;
    
        mmn=mets(:,j);
        X=[];
        X(:,1)=wm;
        X(:,2)=age;
        idx=unique([find(isnan(mmn)); find(isnan(wm))]);
        mmn(idx)=[]; X(idx,:)=[];
        %X(:,1)=zscore(X(:,1)); X(:,2)=zscore(X(:,2));mmn=zscore(mmn);
        met=zscore(mmn); X(:,1)=zscore(X(:,1)); X(:,2)=zscore(X(:,2)); 
        idx=find(isoutlier(mmn));mmn(idx)=[];X(idx,:)=[];
        [t,F,p,df,R2,cR2,B] = glm(mmn,X,[1 0],0);
        wm_mmn_cor_t_glm(m,j)=t; wm_mmn_cor_p_glm(m,j)=p;
        
        mmn=mets(:,j);
        X=[];
        X(:,1)=wm;
        idx=unique([find(isnan(mmn)); find(isnan(wm))]);
        mmn(idx)=[]; X(idx,:)=[];
        %X(:,1)=zscore(X(:,1));mmn=zscore(mmn);
        met=zscore(mmn); X(:,1)=zscore(X(:,1));
        idx=find(isoutlier(mmn));mmn(idx)=[];X(idx,:)=[];
        [t,F,p,df,R2,cR2,B] = glm(mmn,X,[1],0);
        wm_mmn_t_glm(m,j)=t; wm_mmn_p_glm(m,j)=p;
        
        mmn=mets(:,j);
        X=[];
        X(:,1)=hip;
        X(:,2)=age;
        idx=unique([find(isnan(mmn)) ;find(isnan(hip))]);
        mmn(idx)=[];X(idx,:)=[];% mmn=zscore(mmn); 
        %X(:,1)=zscore(X(:,1)); X(:,2)=zscore(X(:,2));
        met=zscore(mmn); X(:,1)=zscore(X(:,1)); X(:,2)=zscore(X(:,2)); 
        idx=find(isoutlier(mmn));mmn(idx)=[];X(idx,:)=[];
        [t,F,p,df,R2,cR2,B] = glm(mmn,X,[1 0],0);
        hip_mmn_cor_t_glm(m,j)=t; hip_mmn_cor_p_glm(m,j)=p;
        
        mmn=mets(:,j);
        X=[];
        X(:,1)=hip;
        idx=unique([find(isnan(mmn)); find(isnan(hip))]);
        mmn(idx)=[]; X(idx,:)=[]; 
        %X(:,1)=zscore(X(:,1));mmn=zscore(mmn); 
        met=zscore(mmn); X(:,1)=zscore(X(:,1));
        idx=find(isoutlier(mmn));mmn(idx)=[];X(idx,:)=[];
        [t,F,p,df,R2,cR2,B] = glm(mmn,X,[1],0);
        hip_mmn_t_glm(m,j)=t; hip_mmn_p_glm(m,j)=p;
        
    end
end

figure('color','w');
imagesc(squeeze(gm_mmn_cor_t_glm)); axis equal off
caxis([-3 3]); axis on;
xlim([0.5 9.5]);yticks(1:length(met_names)); ylim([0.5 length(met_names)+0.5])
set(gca,'YTickLabel',varnames)
set(gca,'XTick',1:9,'XTickLabel',{'MAG','Fr-MAG','Tem-MAG','GRAD','Fr-GRAD','Tem-GRAD','EEG','Fr-EEG','Tem-EEG'})
xtickangle(30); colorbar; colormap(redblue)
title(['%GMV x MMN -Age ' window_nm{w}])

[a,b]=ind2sub(size(squeeze(gm_mmn_cor_p_glm)),find(squeeze(gm_mmn_cor_p_glm<=0.05))); hold on;
for d=1:length(a)
    plot(b(d),a(d),'-s','markersize',8,'markeredgecolor','w','markerfacecolor','w')
end

[a,b]=ind2sub(size(squeeze(gm_mmn_cor_p_glm)),find(squeeze(gm_mmn_cor_p_glm<=0.01))); hold on;
for d=1:length(a)
    plot(b(d),a(d),'-s','markersize',12,'markeredgecolor','w','markerfacecolor','w')
end
print(gcf,'-dbmp', '-r0',[outdir '/GLMs/' 't_map_' session  '_GMV_MMN_GLM_age-cor_' testname '_' window_nm{w} '.bmp'],'-r200');close(gcf)

figure('color','w');
imagesc(squeeze(gm_mmn_t_glm)); axis equal off
caxis([-3 3]); axis on;
xlim([0.5 9.5]);yticks(1:length(met_names));ylim([0.5 length(met_names)+0.5])
set(gca,'YTickLabel',varnames)
set(gca,'XTick',1:9,'XTickLabel',{'MAG','Fr-MAG','Tem-MAG','GRAD','Fr-GRAD','Tem-GRAD','EEG','Fr-EEG','Tem-EEG'})
xtickangle(30); colorbar; colormap(redblue)
title(['%GMV x MMN ' window_nm{w}])

[a,b]=ind2sub(size(squeeze(gm_mmn_p_glm)),find(squeeze(gm_mmn_p_glm<=0.05))); hold on;
for d=1:length(a)
    plot(b(d),a(d),'-s','markersize',8,'markeredgecolor','w','markerfacecolor','w')
end

[a,b]=ind2sub(size(squeeze(gm_mmn_p_glm)),find(squeeze(gm_mmn_p_glm<=0.01))); hold on;
for d=1:length(a)
    plot(b(d),a(d),'-s','markersize',12,'markeredgecolor','w','markerfacecolor','w')
end
print(gcf,'-dbmp', '-r0',[outdir '/GLMs/' 't_map_' session  '_GMV_MMN_GLM_' testname '_' window_nm{w} '.bmp'],'-r200'); close(gcf)


figure('color','w');
imagesc(squeeze(wm_mmn_cor_t_glm)); axis equal off
caxis([-3 3]); axis on;
xlim([0.5 9.5]);yticks(1:length(met_names)); ylim([0.5 length(met_names)+0.5])
set(gca,'YTickLabel',varnames)
set(gca,'XTick',1:9,'XTickLabel',{'MAG','Fr-MAG','Tem-MAG','GRAD','Fr-GRAD','Tem-GRAD','EEG','Fr-EEG','Tem-EEG'})
xtickangle(30); colorbar; colormap(redblue)
title(['%WMV x MMN -Age ' window_nm{w}])

[a,b]=ind2sub(size(squeeze(wm_mmn_cor_p_glm)),find(squeeze(wm_mmn_cor_p_glm<=0.05))); hold on;
for d=1:length(a)
    plot(b(d),a(d),'-s','markersize',8,'markeredgecolor','w','markerfacecolor','w')
end

[a,b]=ind2sub(size(squeeze(wm_mmn_cor_p_glm)),find(squeeze(wm_mmn_cor_p_glm<=0.01))); hold on;
for d=1:length(a)
    plot(b(d),a(d),'-s','markersize',12,'markeredgecolor','w','markerfacecolor','w')
end
print(gcf,'-dbmp', '-r0',[outdir '/GLMs/' 't_map_' session  '_WMV_MMN_GLM_age-cor_' testname '_' window_nm{w} '.bmp'],'-r200'); close(gcf)

figure('color','w');
imagesc(squeeze(wm_mmn_t_glm)); axis equal off
caxis([-3 3 ]); axis on;
xlim([0.5 9.5]);yticks(1:length(met_names));ylim([0.5 length(met_names)+0.5])
set(gca,'YTickLabel',varnames)
set(gca,'XTick',1:9,'XTickLabel',{'MAG','Fr-MAG','Tem-MAG','GRAD','Fr-GRAD','Tem-GRAD','EEG','Fr-EEG','Tem-EEG'})
xtickangle(30); colorbar; colormap(redblue)
title(['%WMV x MMN ' window_nm{w}])

[a,b]=ind2sub(size(squeeze(wm_mmn_p_glm)),find(squeeze(wm_mmn_p_glm<=0.05))); hold on;
for d=1:length(a)
    plot(b(d),a(d),'-s','markersize',8,'markeredgecolor','w','markerfacecolor','w')
end

[a,b]=ind2sub(size(squeeze(wm_mmn_p_glm)),find(squeeze(wm_mmn_p_glm<=0.01))); hold on;
for d=1:length(a)
    plot(b(d),a(d),'-s','markersize',12,'markeredgecolor','w','markerfacecolor','w')
end
print(gcf,'-dbmp', '-r0',[outdir '/GLMs/' 't_map_' session  '_WMV_MMN_GLM_' testname '_' window_nm{w} '.bmp'],'-r200'); close(gcf)

figure('color','w');
imagesc(squeeze(hip_mmn_cor_t_glm)); axis equal off
caxis([-3 3]); axis on;
xlim([0.5 9.5]);yticks(1:length(met_names));ylim([0.5 length(met_names)+0.5])
set(gca,'YTickLabel',varnames)
set(gca,'XTick',1:9,'XTickLabel',{'MAG','Fr-MAG','Tem-MAG','GRAD','Fr-GRAD','Tem-GRAD','EEG','Fr-EEG','Tem-EEG'})
xtickangle(30); colorbar; colormap(redblue)
title(['%HIP x MMN -Age ' window_nm{w}])

[a,b]=ind2sub(size(squeeze(hip_mmn_cor_p_glm)),find(squeeze(hip_mmn_cor_p_glm<=0.05))); hold on;
for d=1:length(a)
    plot(b(d),a(d),'-s','markersize',8,'markeredgecolor','w','markerfacecolor','w')
end

[a,b]=ind2sub(size(squeeze(hip_mmn_cor_p_glm)),find(squeeze(hip_mmn_cor_p_glm<=0.01))); hold on;
for d=1:length(a)
    plot(b(d),a(d),'-s','markersize',12,'markeredgecolor','w','markerfacecolor','w')
end
print(gcf,'-dbmp', '-r0',[outdir '/GLMs/' 't_map_' session  '_HIP_MMN_GLM_age-cor_' testname '_' window_nm{w} '.bmp'],'-r200'); close(gcf)

figure('color','w');
imagesc(squeeze(hip_mmn_t_glm)); axis equal off
caxis([-3 3 ]); axis on;
xlim([0.5 9.5]);yticks(1:length(met_names));ylim([0.5 length(met_names)+0.5])
set(gca,'YTickLabel',varnames)
set(gca,'XTick',1:9,'XTickLabel',{'MAG','Fr-MAG','Tem-MAG','GRAD','Fr-GRAD','Tem-GRAD','EEG','Fr-EEG','Tem-EEG'})
xtickangle(30); colorbar; colormap(redblue)
title(['%HIP x MMN ' window_nm{w}])

[a,b]=ind2sub(size(squeeze(hip_mmn_p_glm)),find(squeeze(hip_mmn_p_glm<=0.05))); hold on;
for d=1:length(a)
    plot(b(d),a(d),'-s','markersize',8,'markeredgecolor','w','markerfacecolor','w')
end

[a,b]=ind2sub(size(squeeze(hip_mmn_p_glm)),find(squeeze(hip_mmn_p_glm<=0.01))); hold on;
for d=1:length(a)
    plot(b(d),a(d),'-s','markersize',12,'markeredgecolor','w','markerfacecolor','w')
end
print(gcf,'-dbmp', '-r0',[outdir '/GLMs/' 't_map_' session  '_HIP_MMN_GLM_' testname '_' window_nm{w} '.bmp'],'-r200'); close(gcf)

end


