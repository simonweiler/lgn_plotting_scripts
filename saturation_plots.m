function saturation_plots(r_idx,i_idx,data,ampa)
if ampa==1

for i=1:length(r_idx)
%red contra
red_sat(:,i)=abs(data(r_idx(i)).step_red.neg_peak1(2,:));
max_idx(i)=find(red_sat(:,i)==max(red_sat(:,i)));

end
for i=1:length(i_idx)
%green ipsi
green_sat(:,i)=abs(data(i_idx(i)).step_blue.neg_peak2(1,:));
max_idx_ipsi(i)=find(green_sat(:,i)==max(green_sat(:,i)));
end
norm_red=red_sat./max(red_sat);
peak_sub=red_sat-max(red_sat);

norm_green=green_sat./max(green_sat);
peak_sub_green=green_sat-max(green_sat);


% for i=1:length(i_idx)
% %red ipsi
% red_sat_i(:,i)=abs(data(i_idx(i)).step_red.neg_peak1(2,:));
% max_idx_i(i)=find(red_sat_i(:,i)==max(red_sat_i(:,i)));
% %green contra
% green_sat_i(:,i)=abs(data(i_idx(i)).step_blue.neg_peak2(1,:));
% max_idx_ipsi_i(i)=find(green_sat_i(:,i)==max(green_sat_i(:,i)));
% end
% 
% norm_red_i=red_sat_i./max(red_sat_i);
% peak_sub_i=red_sat_i-max(red_sat_i);
% 
% norm_green=green_sat_i./max(green_sat_i);
% peak_sub_green=green_sat_i-max(green_sat_i);

else
end

set(0,'DefaultLegendAutoUpdate','off')
fig6= figure;set(fig6, 'Name', 'Step protocol');set(fig6, 'Position', [200, 100, 400, 400]);set(gcf,'color','w');
subplot(2,1,1);
hold on;
me=nanmean(norm_red,2);
st=nanstd(norm_red,[],2)/sqrt(length(norm_red))
for i=1:length(norm_red)
exp=plot(norm_red(:,i),'-r')
exp.Color(4) = 0.05;
end
hold on;
errorbar(me,st,'--or','MarkerSize',5,...
    'MarkerEdgeColor','k','MarkerFaceColor','r');hold on;
box off;ylabel('Peak normalized');xticks([1:1:11]);
exp=[];
set(gca,'FontSize',12);

hold on;
me=nanmean(norm_green,2);
st=nanstd(norm_green,[],2)/sqrt(length(norm_green))
for i=1:size(norm_green,2)
exp=plot(norm_green(:,i),'-g')
exp.Color(4) = 0.05;
end
hold on;
errorbar(me,st,'--og','MarkerSize',5,...
    'MarkerEdgeColor','k','MarkerFaceColor','g');hold on;

subplot(2,1,2);
hold on;
me=nanmean(peak_sub/1000,2);
st=nanstd(peak_sub/1000,[],2)/sqrt(length(peak_sub))

errorbar(me,st,'--or','MarkerSize',5,...
    'MarkerEdgeColor','k','MarkerFaceColor','r');hold on;
box off;ylabel('Peak subtracted (nA)');xticks([1:1:11]);
exp=[];
set(gca,'FontSize',12);
hold on;
me=nanmean(peak_sub_green/1000,2);
st=nanstd(peak_sub_green/1000,[],2)/sqrt(length(peak_sub_green))
errorbar(me,st,'--og','MarkerSize',5,...
    'MarkerEdgeColor','k','MarkerFaceColor','g');hold on;
legend('Chrimson','Chronos','Location','southeast');
legend boxoff;set(gca,'FontSize',10);
for i=1:size(peak_sub_green,2)
exp=plot(peak_sub_green(:,i)/1000,'-g')
exp.Color(4) = 0.05;
end
hold on
for i=1:length(peak_sub)
exp=plot(peak_sub(:,i)/1000,'-r')
exp.Color(4) = 0.05;
end
hold on;
xlabel('Step');
set(gca,'FontSize',12);
end