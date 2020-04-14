function irr=irradiance_extrapol(data,cells)

for i=1:length(cells)
    try
  tmp(:,i)=data(cells(i)).step_red.neg_irr_red(2,:) 
    catch
        tmp(:,i)=ones(11,1)*NaN;
    end
end
red2=find(max(tmp)>1)
templ=nanmean(tmp(:,red2),2);
%res=interp1(1:5,templ(3:7),6:9,'linear','extrap')
[t c]=polyfit([1:5]',templ(3:7),1)
ou=[6 7 8 9];
for m=1:4
out(m)=t(1)*ou(m)+t(2)
end

irr=[templ(1:7)' out]
end