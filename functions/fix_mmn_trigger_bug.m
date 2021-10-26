function [dd, ons]=fix_mmn_triggers(fname)

% Function to fix the MMN trigger bug introduced with the Elekta Triux
% scanner upgrade. It reads the triggers from the channel, removes inserted
% triggers, corrects the value of the incorrect triggers.

% dd : corrected time series of the trigger values
% ons : onset of each trigger value in the time series

% Ece K 2021

D=spm_eeg_load(fname);
tr=D.events;
clear mat m

k=1;
for i=1:length(tr)
    if strcmp(tr(i).type,'STI101_up')
        mat(k,1)=tr(i).value;
        mat(k,2)=tr(i).time;
        k=k+1;
    end
end

mat(1,3)=0;
mat(2:end,3)=mat(2:end,2)-mat(1:end-1,2);

mat2=mat(find(mat(:,3)>0.4),:); %triggers split by at least 400ms remove extra triggers
mat2(find(mat2(:,1)>=200),1)=200;% make sure end of tone trains have the same trigger value

while max(mat2(:,3))>1 % insert missing 200 triggers to mark the end of the tone trains
    idx=find(mat2(:,3)>1);
    if idx(1)==1
        m(1,:)=[200 NaN 0.5];
        m=[m;mat2];
        m(2,3)=0.5;
    else
        m(1:idx(1)-1,:)=mat2(1:idx(1)-1,:); 
        m(idx(1),:)=[200 NaN 0.5];
        m=[m;mat2(idx(1):end,:)];
        m(idx(1)+1,3)=0.5;
    end
    mat2=m; m=[];
end; mat2=[mat2;200 NaN 0.5]; % don't forget the last train

idx=find(mat2(:,1)==200); %find the tone trains
tones=100:150; % correct trigger values
for i=1:length(idx) % replace incorrect triggers with correct ones
    if idx(i)~=idx(end)
        mat2((idx(i)+1):(idx(i+1)-1),4)=tones(1:length((idx(i)+1):(idx(i+1)-1)));
    end
end

mat2(:,1)=mat2(:,4);
mat2(:,3:4)=[]; % tidy up
mat2(:,2)=round(mat2(:,2)*1000/2);
mat2(find(mat2(:,1)==0),:)=[];
ons=mat2(:,2)';
dd=zeros(1,ons(end));
dd(ons)=mat2(:,1);




