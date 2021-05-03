
%% Description
%Calculates Fano factor values for a given set of time windows (Ln 59). 
%Uses a defined random stream to select a time
%window of speicifed length from each trace in the data set. Re-running the
%script with the same set of traces and number of time windows will select
%the same winodws from those traces and thus return the same fano factor
%values for each time window.
%
%% Mario Blanco, 06/2015
% Put the path files (5 column data with time/donor/acceptor/FRET/ideal) into a folder called ideal (Ln 23) 
% Example: C:\Users\mrblanc\Dropbox\ArlieCluster (1)\For_Clustering\100uMQ\ideal 
% Add the path files into ideal, then when the dialog box comes up select the 100uMQ folder (one up from ideal). 
% The code will search for a folder named ideal, then create a new folder called remainder and spit out the 95% 
% confidence interval and Fano factors in the command window in matlab. 
%% Paul Lund, 06/2015
% Add option to calculate Fano values a user-specified number of times for
% a given set of time windows using randomly selected window, then
% calculate average value and std dev of the fano factor for each time window in the set.
% Remove remainder calculation. 
%% Import Data
clear all;
close all;

% folder_name = uigetdir('Z:\Robb\Papers\Sim-Karts protocol\test');
%find subdirectories; requires subdir.m function from Elmar Tarajan [Elmar.Tarajan@Mathworks.de]
% [subdirect, files]=subdir(folder_name);

% for n=1:size(subdirect,1)
%     Tokens = {};
%     A=subdirect{1,n};
%     while ~isempty(A)
%         [Tokens{end+1} A] = strtok(A,'\');
%     end
%     
%     if strcmp(Tokens{end},'ideal')
%         %         mkdir(subdirect{1,n},'remainder');
%         fprintf('Time to Analyze\n');
%         sprintf([Tokens{1,end-3:end}],'\n')
%         listingall=dir(fullfile(subdirect{1,n},'*.dat'));
%         for j=1:size(listingall,1)
%             %             if exist(strcat(subdirect{1,n},'\remainder\',strtok(listingall(j,1).name,'.'),'_remain.dat'), 'file')
%             %                 fprintf('Already analyzed\n')
%             %             else
%             RAW=[];
%             RAW=importdata(strcat(subdirect{1,n},'\',listingall(j).name));
%             RAWnormstore{j,n}=RAW;
%         end     
%     end
% end
%% New file loading scheme
[file_name, path] = uigetfile('*.dat','MultiSelect','on');
folder_name = path;
for n=1:length(file_name)
    filename=strcat(path,file_name{n});
    RAW=[];
    RAW=importdata(strcat(filename));
    RAWnormstore{1,n}=RAW;
end   

%% Setup options for Fano calculations
Count=[];
lengthvect=[50 100 150 250 500 1000 1500 2000]; %List of time windows in frames for which to calculate Fano factor 
trial = input ('Number of trials if not using defined random stream? [none, use defined stream]');
if isempty(trial)
    myStream=RandStream('mt19937ar');
    RandStream.setGlobalStream(myStream);
    trial = 1;
    definedStream=1;
else 
    definedStream=0;
    Fanotrials = cell(length(lengthvect),1);
    mkdir(folder_name,'ideal\trials');
    verbose = input ('Save result from each trial? [n]','s');
    if  strcmp(verbose,'y')
        verbose = 1;
    else
        verbose = 0;
    end
end
%% Calculate Fano factors
for  q=1:trial
    for j=1:size(RAWnormstore,2)  % number of experiment files
        for n=1:size(RAWnormstore,1) % obsolete as there is only 1 vector of experimental file data
            if ~isempty(RAWnormstore{n,j})
                for p=1:length(lengthvect)
                    states=unique(RAWnormstore{n,j}(:,1));
                    states=sort(states);
                    if (size(RAWnormstore{n,j},1)-lengthvect(p)-1)>0
                        Int = randi([1 size(RAWnormstore{n,j},1)-lengthvect(p)-1],1,1);
                        if numel(states)<2
                            Count = 0;
                        else
                            %Count{n,p} = numel(find(RAWnormstore{n,j}(Int:Int+lengthvect(p),5)==states(2)));
                            %figure;
                            %plot(RAWnormstore{n,j}(Int:Int+lengthvect(p),5))
                            %hold on;
                            Count{j,p}=round(numel(find(diff(RAWnormstore{n,j}(Int:Int+lengthvect(p)))~=0))/2); % one half of the number of transitions 
                            %title(num2str(Count{n,p}));
                            %hold off;
                        end
                    else
                    end
                end
            else
            end
        end
    end
    
    Fano=[];
    for n=1:size(Count,2)
        Fano(n,1)=var(cat(1,Count{:,n}))./mean(cat(1,Count{:,n}));
        Fano(n,2)=mean(cat(1,Count{:,n}));
        Fano(n,3)=var(cat(1,Count{:,n}));
    end
    
  figure;
    plot(lengthvect,Fano(:,1),'o');
    hold on;
    bounds=gaminv([.025,.975],(size(Count,1)-1)/2,2/(size(Count,1)-1));
%     display(bounds)  
    Fanoout=horzcat(lengthvect'./10,Fano(:,1));
%     display(Fanoout);
    axis([min(lengthvect) max(lengthvect) 0 2.5]);
    plot(lengthvect,repmat(bounds(1),1,length(lengthvect)),'r--');
    plot(lengthvect,repmat(bounds(2),1,length(lengthvect)),'r--');
    ylabel('Fano Factor')
    xlabel('Window length (frames)')
    title(['Fano vs Poisson: lengthvect = ' char(num2str(lengthvect))]);
        
    if definedStream==1
        dlmwrite(strcat(folder_name,'Fano.txt'),Fano,'delimiter','\t');
        % save picture of graph
        h =gcf;
        fname = [folder_name strrep(char(num2str(lengthvect)),'  ','-') '.jpg'];
        print(h,'-r150','-djpeg',fname)
    else
        Fanotrials{q} = Fanoout;
        if verbose == 1
            % save result
            output = vertcat(bounds, [Fanoout(:,1) Fanoout(:,2)]);
            save(strcat(folder_name,'\ideal\trials\','Fano', '_trial',num2str(q),'.txt'),'output','-ascii');
            
            %save picture of plot
            h =gcf;
            fname = [folder_name '\ideal\trials\FanoPlot' '_trial',num2str(q),'.jpg'];
            print(h,'-r150','-djpeg',fname)
        else
        end
        close(gcf);
    end
end

% calculate avergage fano value from trials and std deviation
if definedStream ==0
    avgFanoValues = zeros(length(lengthvect),3);
    for k=1:length(lengthvect)
        b=[];
        for m=1:length(Fanotrials) %collect all fano values for time window k
            b = [b Fanotrials{m}(k,2)];
        end
        avgFanoValues(k,1)=lengthvect(k); % store time window
        avgFanoValues(k,2)=mean(b); % store mean Fano value from trials
        avgFanoValues(k,3)=std(b); % store std of Fano values from trials
    end
    close all;
    %plot avg Fano values with std dev
    figure;
    errorbar(avgFanoValues(:,1),avgFanoValues(:,2),avgFanoValues(:,3),'d');
    hold on;
    display(bounds)
    axis([min(lengthvect)*0.95 max(lengthvect)*1.05 0 2.5]);
    plot(lengthvect,repmat(bounds(1),1,length(lengthvect)),'r--');
    plot(lengthvect,repmat(bounds(2),1,length(lengthvect)),'r--');
    ylabel('Fano factor')
    xlabel('Window length (frames)')
    title(strcat('Average Fano factor:', num2str(trial),' trials'));
    
    %save plot
    h =gcf;
    fname2 = [folder_name '\ideal\trials\avg_FanoPlot_',num2str(trial), '_trials.jpg'];
    print(h,'-r150','-djpeg',fname2)
    
    %save result
    boundsOut = [bounds 0];
    output = vertcat([avgFanoValues(:,1)./10 avgFanoValues(:,2) avgFanoValues(:,3)], boundsOut);
    save(strcat(folder_name,'\ideal\trials\','avgFano_', num2str(trial),'_trials.txt'),'output','-ascii');
end