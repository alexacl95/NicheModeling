function in = bnm_modeling(in,layerfolder,show,method,parallel)
tic
    if nargin<4
        method=4;
    end

    if nargin<5
        parallel=false;
    end

    if parallel
        parforArg=inf;
    else
        parforArg=0;
    end

if ~isempty(layerfolder)
    
    disp('----Reading layers----')
    temp=dir(layerfolder);
    layers={temp.name};
    comp=strncmp('bio',layers,3);
    layers=layers(comp);

    [Z,R] = arcgridread(strcat(layerfolder,layers{1}));
    %[Z,R] = readgeoraster(strcat(layerfolder,layers{1}));
    
    N=length(layers);
    Z=zeros(size(Z,1),size(Z,2),N);
    
    parfor (i =1:N,parforArg)
%         progressbar(i/N)
       Z(:,:,i)=arcgridread(strcat(layerfolder,layers{i}));
    end
    toc
else
    R=in.R;
    Z=in.Z;
end

Tdata=in.Tdata;
indicators=in.Indicators;
vars=in.Vars;
T2=in.T2;

disp('----Modeling----')
[map,response,minimize,roc] = predictor2(Tdata,Z,R,indicators,vars,T2,show,method,[]);
in.Map=map;
in.Response=response;
in.Minimize=minimize;
in.Method=method;
in.Roc=roc;
in.Z=Z;
in.R=R;
toc

if show
    D1 = 2; % Number of rows of subplot
    D2 = 2;
    figure('Name','Map gradient')
    clf
    count=1;
    for i=0.2:0.2:0.8
        map2=map;
        map2(map2<i)=0;
        cont=map2(:);
        cont=cont(~isnan(cont));
        indice=sum(cont>0)/length(cont)*100;
        subplot(D1,D2,count)
        geoshow(map2,R,'DisplayType','surface')
        contourcmap('jet',0:0.05:1,'colorbar','on','location','vertical')
        title(strcat('Score: ',num2str(round(indice,2)),'%'))
        count=count+1;
    end
end
end