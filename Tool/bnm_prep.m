function out = bnm_prep(doc,layerfolder,show,threshold)
tic
if nargin<4
    threshold=0.8;
end

if istable(doc)
    T = doc;
else
    T = readtable(doc);
end

if isstring(layerfolder) || ischar(layerfolder)
    temp=dir(layerfolder);
    layers={temp.name};
    comp=strncmp('bio',layers,3);
    layers=layers(comp);
    disp('----Reading layers----')

    [Z,R] = arcgridread(strcat(layerfolder,layers{1}));
    N=length(layers);
    Z=zeros(size(Z,1),size(Z,2),N);

    for i =1:N
       progressbar(i/N)
       Z(:,:,i)=arcgridread(strcat(layerfolder,layers{i}));
       eval(strcat('T.bio',num2str(i),'=ltln2val(Z(:,:,i),R,[T.LAT],[T.LONG]);'))
    end
    toc
else
    Z = layerfolder.Z;
    R = layerfolder.R;
    N = size(Z,3);
    for i =1:N
       progressbar(i/N)
       eval(strcat('T.bio',num2str(i),'=ltln2val(Z(:,:,i),R,[T.LAT],[T.LONG]);'))
    end
end


disp('----Finding correlation----')

if show
    correlationCircles(T{:,4:end},'varnames',T.Properties.VariableNames(4:end))
end
[Corr,~] = corr(T{:,4:end},'rows','complete');
%Construction of new variables
Corr=abs(Corr);
T2=T;
T=T(:,4:end);
Tdata=table();
ex=[];
k = 5;
vars=[];
for i=1:N
    if sum(i==ex)==0
        if sum(Corr(i,:)>threshold)>1
            f = find(Corr(i,:)>threshold);
            f = setdiff(f,i);
            ex=[ex,f];
            model=ridge([T{:,i}],[T{:,f}],k,0);
            vars=[vars,i];
            Tdata=addprop(Tdata,{strcat('m',num2str(i))},{'table'});
            eval(strcat("Tdata.Properties.CustomProperties.m",num2str(i),"=@(X) X{:,i}-(X{:,f}*model(2:end)+model(1));"))
        end
    end
end
toc
disp('----Creating predictors----')
figure('Name','Histogram')
%indicators=setdiff(1:size(T,2),union(ex,vars));
indicators=setdiff(1:size(T,2),ex);
siz=size(indicators,2);
D1 = floor(sqrt(siz)); % Number of rows of subplot
D2 = D1+ceil((siz-D1^2)/D1);
count=1;
%outlier remotion must starts before saving Tdata as predictors.
%It is highly recommended to start outlier remotion from normalized PCA
%table
for i=indicators
    subplot(D1,D2,count)
    h=histfit(T{:,i},10,'kernel');
    title(strcat('bio',num2str(i)))
    eval(strcat("Tdata.bio",num2str(i),"= [get(h(2),'xdata');get(h(2),'ydata')]';"))
    count=count+1;
end
figure('Name','Histogram2')
siz=size(vars,2);
D1 = floor(sqrt(siz)); % Number of rows of subplot
D2 = D1+ceil((siz-D1^2)/D1);
count=1;
for i=vars
    subplot(D1,D2,count)
    h=histfit(eval(strcat("Tdata.Properties.CustomProperties.m",num2str(i),"(T)")),10,'kernel');
    title(strcat('var',num2str(i)))
    eval(strcat("Tdata.var",num2str(i),"= [get(h(2),'xdata');get(h(2),'ydata')]';"))
    count=count+1; 
end
for i=1:size(Tdata,2)
    maxim=max(Tdata{:,i}(:,2));
    minim=min(Tdata{:,i}(:,2));
    Tdata{:,i}(:,2)=(Tdata{:,i}(:,2)-minim)./(maxim-minim);
end
toc
disp('¡All done!')
out=struct();
out.Tdata=Tdata;
out.Indicators=indicators;
out.Vars=vars;
out.T2=T2;
out.Z=Z;
out.R=R;
end