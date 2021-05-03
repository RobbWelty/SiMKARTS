%% File Selection Step

clear all;
close all;

t=0;
done=0;
while done ==0
    MPathName=[];
    MFileName=[];
    storemat=[];
    storemat2=[];
    storemat3=[];
    t=t+1;
    folder_name{t,1} = uigetdir;
    Tokens = {};
    A=folder_name{t,1};
    while ~isempty(A)
        [Tokens{end+1} A] = strtok(A,'\');
    end
    condition_store{t,1}=Tokens{end};

%% Lets user select appropriate files for  %%
fprintf(1,'Please select all your path files \n')
[MFileName, MPathName] = uigetfile(strcat(folder_name{t,1},'\','*.dat',';','*.datR'),'Please select all your PATH files','Multiselect','on');
% This is just in case the user selects only one file (converts array ->
% cell array), because the rest of the code expects a cell array.
if ~iscell(MFileName)
    G = cell(1);
    G{1} = MFileName;
    MFileName = G;
end
store_filenames{t,1}=MFileName;
done=input('are we done selecting data? [0] - No [1] - Yes');
    if isempty(done)
        done=0;
    else
    end

end


%% Works for 2 state data (an unbound state and bound state)
for j=1:size(folder_name,1)
    nfiles=[];
    nfiles = size(store_filenames{j,1},2);
    TestFileName = store_filenames{j,1};
    TestPathName = strcat(folder_name{j,1},'\');
    Listraw2=cell(1);
        for n=1:nfiles
            Listraw2{n,1}=dlmread(strcat(TestPathName,TestFileName{1,n}));
        end        
% Reassign top states to 0.75 bottom state to 0.00
newrawout=cell(1);
newraw=cell(1);
for p=1:size(Listraw2,1)
    raw=Listraw2{p,1}(:,1);  %the (:,1) section refers to which column the data is taken from where 1 is the first column
    states=unique(raw);      %determines the number of unique value in the raw data
    for n=1:size(states,1)
        if n==1
            raw(raw==states(n),1)=0.00;  % sets the lowest unique value from raw to zero
        else
            raw(raw==states(n),1)=0.75;  % sets every value that isn't the lowest unique value from raw to 0.75
        end
    end
    newrawout{p,1}=raw;  
end
Listraw2store{j,1}=newrawout;  % this is a vector of all of the raw idealized traces

end 


%% Quantify Bursts
set(0, 'DefaultFigureVisible', 'off');

% global analysis and rank score
storetn=[];
positionstore=[];
startpos=[];
finalpos=[];
for j=1:size(Listraw2store,1)    %lets the subsequent loops work on one file at a time    
    for n=1:size(Listraw2store{j,1},1) %       
        tn=[];
        Listraw2temp=Listraw2store{j,1}{n,1}(:,1);
        tn=find(diff(Listraw2temp')>0);
        startpos=size(storetn,2)+1;
        finalpos=size(storetn,2)+size(tn,2);
        if isempty(storetn)
            storetn=horzcat(storetn,tn);
        else            
            storetn=horzcat(storetn,1000+storetn(end)+tn);
        end       
        positionstore{j,1}{n,1}=[startpos finalpos];
    end    
end

%%%%%%%%  Robb's Notes %%%
% position store is just a collection vectors of spikes starting positions.  
%%%%%%%%%%
[a,b,c]=burst(storetn,200,2); %last two parameters 1) max ISI to be included in burst and 2) burst rank score
temp=cell(1);
for p=1:size(c,2)
    if p==1
        temp{p,1}=c(p)+1:(c(p)+b(p)-1);
    else
        temp{p,1}=c(p)+1:(c(p)+b(p)-1);
    end 
end


%% Extract bursts corresponding to segment
storelimits=cell(1);
for j=1:size(Listraw2store,1)
     for n=1:size(Listraw2store{j,1},1)
    where=(c>=positionstore{j,1}{n,1}(1) & c<=positionstore{j,1}{n,1}(2));
    storelimits{j,1}{n,1}=c(where);        
    storelimits{j,1}{n,2}=b(where);
     end
end

%%
positionstorerenum=[];
for j=1:size(positionstore,1)
    for n=1:size(positionstore{j,1},1)
       if isempty(positionstorerenum)
           positionstorerenum{j,1}{n,1}=positionstore{j,1}{n,1};
       else
           if n==1
           positionstorerenum{j,1}{n,1}=positionstore{j,1}{n,1}-positionstore{j-1,1}{end,1};
           else           
           positionstorerenum{j,1}{n,1}=positionstore{j,1}{n,1}-positionstore{j,1}{n-1,1}(2);
           end
       end        
    end
end

%% Extract bursts corresponding to segment
storelimitsrenum=[];
for j=1:size(storelimits,1)
    for n=1:size(storelimits{j,1},1)
        if isempty(storelimitsrenum)
            storelimitsrenum{j,1}{n,1}=storelimits{j,1}{n,1};
            storelimitsrenum{j,1}{n,2}=storelimits{j,1}{n,2};
            storelimitsrenumnormalized{j,1}{n,1}=storelimits{j,1}{n,2}./size(Listraw2store{j,1}{n,1}(:,1),1);

        else
            
            if n==1
                storelimitsrenum{j,1}{n,1}=storelimits{j,1}{n,1}-positionstore{j-1,1}{end,1}(2);
                storelimitsrenum{j,1}{n,2}=storelimits{j,1}{n,2};
                storelimitsrenumnormalized{j,1}{n,1}=storelimits{j,1}{n,2}./size(Listraw2store{j,1}{n,1}(:,1),1);
            else                
                storelimitsrenum{j,1}{n,1}=storelimits{j,1}{n,1}-positionstore{j,1}{n-1,1}(2);
                storelimitsrenum{j,1}{n,2}=storelimits{j,1}{n,2};
                storelimitsrenumnormalized{j,1}{n,1}=storelimits{j,1}{n,2}./size(Listraw2store{j,1}{n,1}(:,1),1);
            end
                        
        end
    end
end


%%
for j=1:size(Listraw2store,1)
    for n=1:size(Listraw2store{j,1},1)
        tn=[];
        tn=find(diff(Listraw2store{j,1}{n,1}(:,1)')>0);
        dwells=diff(horzcat(0,tn));
        tn2=[];
        tn2=find(diff(Listraw2store{j,1}{n,1}(:,1)')<0);
        
        if ~isempty(storelimitsrenum{j,1}{n,1})
        store=cell(1);
        c=[];
        b=[];
        c=storelimitsrenum{j,1}{n,1};
        b=storelimitsrenum{j,1}{n,2};        
        
        for p=1:size(c,2)
            lengthofburst=[];
            if p==1
                store{p,1}=dwells(c(p)+1:(c(p)+b(p)-1));                
                
                start=tn2(c(p));               
                if (c(p)+b(p)-1)>= size(tn2,2)                                                    
                endspot=tn2(end);
                else                
                endspot=tn2(c(p)+b(p)-1);
                end                                                
            else                
                store{p,1}=dwells(c(p)+1:(c(p)+b(p)-1));                
                start=tn2(c(p));                
                if (c(p)+b(p)-1)>= size(tn2,2)                                    
                endspot=tn2(end);
                else
                endspot=tn2(c(p)+b(p)-1);
                end                
            end
            lengthofburst=endspot-start;
            burstlengthstore{j,1}{n,p}=lengthofburst;
            burstlengthstore2{j,1}{n,p}=lengthofburst./size(Listraw2store{j,1}{n,1}(:,1),1);
        end        
        
        dwellsinburst=cat(2,store{:});
        
        dwellsinburststore{j,1}{n,1}=dwellsinburst;
        
        dwellsoutofburst=[];
        dwellstemp=dwells;
        for h=1:size(dwellsinburst,2)
            where=[];
            where=find(dwellstemp==dwellsinburst(1,h));
            dwellstemp(where)=[];
        end
        dwellsoutofburst=dwellstemp;                
        
        dwellsoutofburststore{j,1}{n,1}=dwellsoutofburst;
        else
            
        tn=[];
        tn=find(diff(Listraw2store{j,1}{n,1}(:,1)')>0);
        dwells=diff(horzcat(0,tn));
        dwellsinburst=[];
        dwellsoutofburst=[];
        dwellsoutofburst=dwells;
        dwellsinburststore{j,1}{n,1}=dwellsinburst;
        dwellsoutofburststore{j,1}{n,1}=dwellsoutofburst;            
        end
    end
end

%%
for n=1:size(dwellsinburststore,1)
    [y,x]=hist(cat(2,dwellsinburststore{n,1}{:}),1:20:1000);   
    [y2,x2]=hist(cat(2,dwellsoutofburststore{n,1}{:}),1:20:1000);
    figure;
    bar(x./10,y./sum(y+y2))    
    title(strcat('Condition',num2str(n),'ISI In Burst'));   
    dlmwrite(strcat('Condition',num2str(n),'ISI In Burst','.txt'),cat(2,dwellsinburststore{n,1}{:})','delimiter','\t')
    xlim([0 75])
    xlabel('Time (frames)') 
    ylabel('Count')
    h=gcf;
    print(h,'-djpeg',strcat('Condition',num2str(n),'ISI In Burst'))
    
    figure;
    bar(x2./10,y2./sum(y+y2))
    title(strcat('Condition',num2str(n),'ISI Out of Burst'));  
    xlim([0 75])
    xlabel('Time (frames)') 
    ylabel('Count')
    h=gcf;
    print(h,'-djpeg',strcat('Condition',num2str(n),'ISI Out of Burst'))
    dlmwrite(strcat('Condition',num2str(n),'ISI Out of Burst','.txt'),cat(2,dwellsoutofburststore{n,1}{:})','delimiter','\t')
    
    figure;
    stairs(x/10,cumsum(y)./sum(y+y2),'r','LineWidth',2)   
    hold on;
    stairs(x/10,cumsum(y2)./sum(y+y2),'k','LineWidth',2)
    title(strcat('Condition',num2str(n),'CDF'));  
    legend('Dwells in busrt','Dwells out of burst','Location','southeast')
    xlabel('Time (frames)') 
    ylabel('Count')
    h=gcf;
    print(h,'-djpeg',strcat('Condition',num2str(n),'CDF'))
    
end


%% Number of events per burst
storehistnumeventsperburst=cell(1);
for n=1:size(storelimitsrenum,1)
    figure;    
    [y,x]=hist(cat(2,storelimitsrenum{n,1}{:,2}),1:1:50);
    bar(x,y./sum(y))
    xlabel('# of Binding Events')
    ylabel('Count')
    title(strcat('Condition',num2str(n),'EventsPerBurst'))
    axis tight
    ylim([0 0.5])    
    h=gcf;
    print(h,'-djpeg',strcat('Condition',num2str(n),'EventsPerBurst'))
    storehistnumeventsperburst{n,1}=[x' y'];
    dlmwrite(strcat('Condition',num2str(n),'NumEventsBurst','.txt'),cat(2,storelimitsrenum{n,1}{:,2})','delimiter','\t')   
end
%close all

%% Total length of burst
storehistburstlength=cell(1);
for n=1:size(burstlengthstore,1)    
    figure;
    [y,x]=hist(cat(2,burstlengthstore{n,1}{:,:}),1:20:2000);
    bar(x/10,y./sum(y))
    xlabel('Total Length of Burst (s)')
    ylabel('Count')
    title(strcat('Condition',num2str(n),'Total Length of Burst'))
    axis tight
    ylim([0 0.35])
    h=gcf;
    print(h,'-djpeg',strcat('Condition',num2str(n),'TotalLengthofBurst'))    
    storehistburstlength{n,1}=[x' y'];             
    dlmwrite(strcat('Condition',num2str(n),'BurstLength','.txt'),cat(2,burstlengthstore{n,1}{:,:})','delimiter','\t')       
end

function [archive_burst_RS,archive_burst_length,archive_burst_start]=burst(tn,limit,RSalpha)
%"Rank Surprise" method for burst detection. Requires statistical toolbox
%of Matlab 6.5 or any program computing the Gaussian CDF.
%
%USE : [archive_burst_RS,archive_burst_length,archive_burst_start]=burst(tn,limit,RSalpha)
%
%INPUT :  tn - spike times
%         limit - ISI value not to include in a burst
%         RSalpha - minimum surprise value to consider
%
%OUTPUT : archive_burst_RS - "rank surprise" values for each burst detected
%         archive_burst_length - burst length for each burst detected (in spikes)
%         archive_burst_start – Spike number of burst start for each burst detected
%
% Base code provided in Supplementary Information associated with 
% Gourévitch, B. & Eggermont, J. J. A nonparametric approach for detection of bursts in spike trains. 
% J. Neurosci. Meth. 160, 349-358 (2007). 
% Minor changes made by MRB 2015 to adapt for appication to SiM-KARTS data
%% Checking input

%risk level
if nargin<3,
    RSalpha=-log(0.01);
end;

%% General parameters

%limit for using the real distribution
q_lim=30;
%minimum length of a burst in spikes
l_min=2;

%% General vectors

alternate=ones(400,1);
alternate(2:2:end)=-1;
%log factorials
log_fac=cumsum(log(1:q_lim));
%to make tn an horizontal vector
tn=tn(:)';

%% Ranks computation
%compute the ISI
ISI=diff(tn);
%ISI=diff([0 tn]); %MRB edit
N=length(ISI);
%ISI value not to include in a burst
if nargin<2,
    %percentile 75% (default)
    limit=prctile(ISI,75);
end;
%compute ranks
R=val2rk(ISI);

%% Find sequences of ISI under 'limit'
ISI_limit=diff(ISI<limit);
%first time stamp of these intervals
begin_int=find(ISI_limit==1)+1;
%manage the first ISI
if ISI(1)<limit,
    begin_int=[1 begin_int];%the first IS is under limit
end;
%last time stamp of these intervals
end_int=find(ISI_limit==-1);
%manage the last ISI
if length(end_int)<length(begin_int),
    end_int=[end_int N];
end;
%length of intervals of interest
length_int=end_int-begin_int+1;

%% Initializations
archive_burst_RS=[];
archive_burst_length=[];
archive_burst_start=[];

%% Going through the intervals of interest
indic=0;
for n_j=begin_int
    indic=indic+1;
    p_j=length_int(indic);
    
    subseq_RS=[];    
        %test each set of spikes
        for i=0:p_j-(l_min-1)
        %for i=0:p_j-(l_min)    %MRB edit
            %length of burst tested
            q=l_min-2;
            while (q<p_j-i)
                q=q+1;
                %statistic
                u=sum(R(n_j+i:n_j+i+q-1));
                u=floor(u);
                if q<q_lim,
                    %exact discrete distribution
                    k=0:(u-q)/N;
                    length_k=length(k);
                    prob=exp((sum(log(u-repmat(k,q,1)*N-repmat((0:q-1)',1,length_k)))...
                        -log_fac([1 k(2:end)])-log_fac(q-k))-q*log(N))*alternate(1:length_k);
                else
                    %approximate Gaussian distribution
                    prob=normcdf((u-q*(N+1)/2)/sqrt(q*(N^2-1)/12));
                end;
                RS=-log(prob);
                %archive results for each subsequence [RSstatistic beginning length]
                if RS>RSalpha,
                    subseq_RS(end+1,:)=[RS i q];                                    
                end;                
            end;
        end;
  
    %vet results archive to extract most significant bursts
    if ~isempty(subseq_RS),
        %sort RS for all subsequences
        subseq_RS=-sortrows(-subseq_RS,1);        
        while ~isempty(subseq_RS),
            %extract most surprising burst
            current_burst=subseq_RS(1,:);            
%             %keep only other bursts non-overlapping with this burst
%           subseq_RS=subseq_RS(subseq_RS(:,2)+subseq_RS(:,3)-1<current_burst(2)|subseq_RS(:,2)>current_burst(2)+current_burst(3)-1,:);
%             
            archive_burst_RS(end+1)=current_burst(1);
            archive_burst_length(end+1)=current_burst(3)+1;%number of ISI involved + 1
            archive_burst_start(end+1)=n_j+current_burst(2);
            %remove most surprising burst from the set
            %subseq_RS=subseq_RS(2:end,:);
            
            %keep only other bursts non-overlapping with this burst
            subseq_RS=subseq_RS(subseq_RS(:,2)+subseq_RS(:,3)-1<current_burst(2)|subseq_RS(:,2)>current_burst(2)+current_burst(3)-1,:);                       
        end;
    end;
    
end;

%sort bursts by ascending time
[archive_burst_start,ind_sort]=sort(archive_burst_start);
archive_burst_RS=archive_burst_RS(ind_sort);
archive_burst_length=archive_burst_length(ind_sort);

%% Utility - Rank computation

    function [ranks]=val2rk(values)
        %Convert values to ranks, with mean of ranks for tied values. Alternative
        %and faster version of "tiedrank" in statistical toolbox of Matlab 6-7.
        lp=length(values);
        [y,cl]=sort(values);
        rk(cl)=(1:lp);
        [y,cl2]=sort(-values);
        rk2(cl2)=(1:lp);
        ranks=(lp+1-rk2+rk)/2;
    end
end