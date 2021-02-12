function [map,response,minimize]=predictor2(Tdata,Z,R,indicators,vars,T2,show,method,coeff)
stot=length(indicators);
tot=stot+length(vars);
if isempty(coeff)
    coeff=ones(1,tot);
end
reps=size(Z);
map=ones(reps(1),reps(2));
caps=reps(3);
template=Z(:,:,1);
data=NaN(length(template(:)),caps);
for i=1:caps
    template=Z(:,:,i);
    data(:,i)=template(:);
end
response=NaN(length(template(:)),length(tot));
i=1;
for j=indicators
    pointer=~isnan(data(:,j));
    xdata=Tdata{:,i}(:,1);
    ydata=Tdata{:,i}(:,2);
    response(pointer,i)=interp1(xdata,ydata.^(coeff(i)),data(pointer,j),'pchip',0);
    i=i+1;
end
for j=vars
    pointer=~isnan(data(:,i));
    xdata=Tdata{:,i}(:,1);
    ydata=Tdata{:,i}(:,2);
    response(pointer,i)=interp1(xdata,ydata.^(coeff(i)),eval(strcat("Tdata.Properties.CustomProperties.m",num2str(j),"(array2table(data(pointer,:)));")),'pchip',0);
    i=i+1;
end
final=(prod(response,2)).^(1/sum(coeff(1:stot)));
map(:)=final(:);
minimize=curverock(map,R,T2,show,method);

if show
    figure(3)
    clf
    geoshow(map,R,'DisplayType','surface');
    contourcmap('jet',0:0.05:1,'colorbar','on','location','vertical')
    geoshow(T2.LAT, T2.LONG, 'DisplayType', 'Point', 'Marker', 'o',...
        'MarkerSize',5,'MarkerFaceColor',[.95 .9 .8],'MarkerEdgeColor',...
        'black','Zdata', 2*ones(length(T2.LONG),1));
end
end
    

