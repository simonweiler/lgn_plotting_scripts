%SW191127
%PLEASE READ
%Script that needs the LGN structure with all extracted features: updated
%structure are located on I

%Script displays one by one polished ephys figures for the manuscript

%For supplementary figures: it also needs the Chrimson only injection animals which are
%in a separete data strcuture called :  

% Data_Chrimson_only_MF.mat


%Additional functions needed:
    %uipickfiles
    
% %% Load LGN structure into workspace, loaded strcuture should have the name "data"
% str_LGN     = 'C:\Users\Simon-localadmin\Documents\MargrieLab\PhDprojects\LGN';
% folder_list = uipickfiles('FilterSpec',str_LGN);
% load(char(folder_list));
% %check if the structure is called "data" otherwise rename to data
% if exist('data')==1
% data=data;
% else
% wto=whos;wt={wto.class};wt_idx=find(contains(wt,'struct')==1);
% data=eval(wto(wt_idx).name);
% end    
%% Load relevant data from structure

[data] = load_StrArray('Full_Data_26-Feb-2020','C:\Users\Simon-localadmin\Documents\MargrieLab\PhDprojects\LGN','IncludeField',{'step_red','step_blue','ODI_AMPA_step','ODI_NMDA_step','red_failure_AMPA','red_failure_NMDA','blue_failure_AMPA','blue_failure_NMDA'},[])
%% Read out general features
for i=1:length(data)
    %Ocular category based on ODI_step for AMPA/NMDA
    %contra only
    if data(i).ODI_AMPA_step==1 & data(i).ODI_NMDA_step==1;
    ocular_cat(i)=1;
    %ipsi only
    elseif data(i).ODI_AMPA_step==-1 & data(i).ODI_NMDA_step==-1;
    ocular_cat(i)=2;
    %bino
    elseif data(i).ODI_AMPA_step ~=1 & data(i).ODI_NMDA_step~=1 &  data(i).ODI_AMPA_step ~=-1 & data(i).ODI_NMDA_step~=-1;
    ocular_cat(i)=3;
    %contra silent
    elseif data(i).ODI_AMPA_step==1 & data(i).ODI_NMDA_step<1 & data(i).ODI_NMDA_step>0; 
    ocular_cat(i)=4;   
    %ipsi silent
    elseif data(i).ODI_AMPA_step==-1 & data(i).ODI_NMDA_step>-1 & data(i).ODI_NMDA_step<0; 
    ocular_cat(i)=5;
    elseif data(i).ODI_AMPA_step<1 & data(i).ODI_AMPA_step>0 & data(i).ODI_NMDA_step==1;
    ocular_cat(i)=6;
    elseif data(i).ODI_AMPA_step==-1 & data(i).ODI_NMDA_step>-1;
    ocular_cat(i)=7;
    %not recorded or sth else
    else
    ocular_cat(i)=NaN;    
    end
    %Contra red or green injected
    contra_r(i)=data(i).brain_contra_ipsi;
    %Eye injection order
    eye(i)=data(i).eye_inj_ord;
    %Hemipshere
    hemi{i}=data(i).hemisphere;
    %MD 
    md(i)=data(i).MD;
end
%% idx cells based on their ODI, injection order, and md; intersectional idx
ocat1=find(ocular_cat==1);ocat2=find(ocular_cat==2);ocat3=find(ocular_cat==3);ocat4=find(ocular_cat==4);ocat5=find(ocular_cat==5);
c_red=find(contra_r==1);
i_red=find(contra_r==0);
md_n=find(md==0);
%% 

[r1 ~]=intersect(ocat1,c_red);[i1 ~]=intersect(ocat1,i_red);
[r2 ~]=intersect(ocat2,c_red);[i2 ~]=intersect(ocat2,i_red);
[r3 ~]=intersect(ocat3,c_red);[i3 ~]=intersect(ocat3,i_red);
[r4 ~]=intersect(ocat4,c_red);[i4 ~]=intersect(ocat4,i_red);
[r5 ~]=intersect(ocat5,c_red);[i5 ~]=intersect(ocat5,i_red);
[no_md1 ~]=intersect(ocat1,md_n);[no_md2 ~]=intersect(ocat2,md_n);[no_md3 ~]=intersect(ocat3,md_n);
[no_md4 ~]=intersect(ocat4,md_n);[no_md5 ~]=intersect(ocat5,md_n);
%% Get the intersectional cell idx for ocular category and non md contra red or ipsi red
[r_nm1 ~]=intersect(no_md1,c_red);[i_nm1 ~]=intersect(no_md1,i_red);
[r_nm2 ~]=intersect(no_md2,c_red);[i_nm2 ~]=intersect(no_md2,i_red);
[r_nm3 ~]=intersect(no_md3,c_red);[i_nm3 ~]=intersect(no_md3,i_red);
[r_nm4 ~]=intersect(no_md4,c_red);[i_nm4 ~]=intersect(no_md4,i_red);
[r_nm5 ~]=intersect(no_md5,c_red);[i_nm5 ~]=intersect(no_md5,i_red);
%% Observed response types
%Choose different example cells for the different category 
%cat3=11
cat1=14; cat2=3; cat3=19; cat4=2 ; cat5=8;
%linewidth of traces
tr_l=1.2;
%linewidth of red blue indicator
i_l=5;
%which trace from the 11 steps 
trace_nr=7;
%Choose whether red should be contra or ipsi
contra_red=1;
base_start          =   1;
base_end            =   99;
redpeak_start       =   100;
redpeak_end         =   350;%looking for the peak shoud be the same as for blue even though the window is longer?
bluepeak_start      =   351;
bluepeak_end        =   400;
scale_x= 100;
scale_y= 400;
%Find the min/max for the traces for each category
sb_min=[min(data(r1(cat1)).step_red.ephys_traces_70(:,trace_nr,2)) min(data(r2(cat2)).step_red.ephys_traces_70(:,trace_nr,2))...
        min(data(r3(cat3)).step_red.ephys_traces_70(:,trace_nr,2)) min(data(r4(cat4)).step_red.ephys_traces_70(:,trace_nr,2))...
        min(data(ocat5(cat5)).step_red.ephys_traces_70(:,trace_nr,2))];
sb_max=[max(data(r1(cat1)).step_red.ephys_traces_40(:,trace_nr,2)) max(data(r2(cat2)).step_red.ephys_traces_40(:,trace_nr,2))...
        max(data(r3(cat3)).step_red.ephys_traces_40(:,trace_nr,2)) max(data(r4(cat4)).step_red.ephys_traces_40(:,trace_nr,2))...
        max(data(ocat5(cat5)).step_red.ephys_traces_40(:,trace_nr,2))];
ov_min=min(sb_min);
ov_max=max(sb_max);
%figure;set(gcf,'color','w');plot(sb_min,'--or');hold on;plot(sb_max,'--ob');ylabel('Synaptic Input peak (pA)');xlabel('Category');box off
%The offset between the 70 and 40 mV
yoffset=ov_max/8;
ov_maxo=ov_max+ov_max/4;

%FIGURE for paper
%RED=Contra; GREEN=ipsi
%set up figure
fig1= figure;set(fig1, 'Name', 'Observed response types');set(fig1, 'Position', [400, 600, 600, 200]);set(gcf,'color','w');

%Contra only
sb1=subplot(1,5,1);
if data(r1(cat1)).experimentator=='SW'
    srF=1;
else data(r1(cat1)).experimentator=='MF'  
    srF=2;
end
plot(data(r1(cat1)).step_red.ephys_traces_70(:,trace_nr,2),'Color','k','LineWidth',tr_l);box off;axis off
hold on;plot(data(r1(cat1)).step_red.ephys_traces_40(:,trace_nr,2)+yoffset,'Color',[0.5 0.5 0.5],'LineWidth',tr_l);box off;axis off
set(sb1,'ylim',[ov_min ov_maxo]);
hold on; x1= redpeak_start*srF;x2=redpeak_end*srF;p1=plot([x1 x2],[ov_maxo ov_maxo],'-','Color','r','LineWidth',i_l);
hold on; x1= bluepeak_start*srF;x2=bluepeak_end*srF; hold on;p1=plot([x1 x2],[ov_maxo ov_maxo],'-','Color','b','LineWidth',i_l);
hold on;title('contra only');
%Ipsi only
sb2=subplot(1,5,2);
if data(r2(cat2)).experimentator=='SW'
    srF=1;
else data(r2(cat2)).experimentator=='MF'  
    srF=2;
end
plot(data(r2(cat2)).step_red.ephys_traces_70(:,trace_nr,2),'Color','k','LineWidth',tr_l);box off;axis off
hold on;plot(data(r2(cat2)).step_red.ephys_traces_40(:,trace_nr,2)+yoffset,'Color',[0.5 0.5 0.5],'LineWidth',tr_l);box off;axis off
set(sb2,'ylim',[ov_min ov_maxo]);
hold on; x1= redpeak_start*srF;x2=redpeak_end*srF;p1=plot([x1 x2],[ov_maxo ov_maxo],'-','Color','r','LineWidth',i_l);
hold on; x1= bluepeak_start*srF;x2=bluepeak_end*srF; hold on;p1=plot([x1 x2],[ov_maxo ov_maxo],'-','Color','b','LineWidth',i_l);
hold on;title('ipsi only');

%Bino at 70 and 40 mV
sb3=subplot(1,5,3);  
if data(r3(cat3)).experimentator=='SW'
    srF=1;
else data(r3(cat3)).experimentator=='MF'  
    srF=2;
end
plot(data(r3(cat3)).step_red.ephys_traces_70(:,trace_nr,2),'Color','k','LineWidth',tr_l);box off;axis off
hold on;plot(data(r3(cat3)).step_red.ephys_traces_40(:,trace_nr,2)+yoffset,'Color',[0.5 0.5 0.5],'LineWidth',tr_l);box off;axis off
set(sb3,'ylim',[ov_min ov_maxo]);
hold on; x1= redpeak_start*srF;x2=redpeak_end*srF;p1=plot([x1 x2],[ov_maxo ov_maxo],'-','Color','r','LineWidth',i_l);
hold on; x1= bluepeak_start*srF;x2=bluepeak_end*srF; hold on;p1=plot([x1 x2],[ov_maxo ov_maxo],'-','Color','b','LineWidth',i_l);
hold on;title('binocular');

%contra silent
sb4=subplot(1,5,4);  
if data(r4(cat4)).experimentator=='SW'
    srF=1;
else data(r4(cat4)).experimentator=='MF'  
    srF=2;
end
plot(data(r4(cat4)).step_red.ephys_traces_70(:,trace_nr,2),'Color','k','LineWidth',tr_l);box off;axis off
hold on;plot(data(r4(cat4)).step_red.ephys_traces_40(:,trace_nr,2)+yoffset,'Color',[0.5 0.5 0.5],'LineWidth',tr_l);box off;axis off
set(sb4,'ylim',[ov_min ov_maxo]);
hold on; x1= redpeak_start*srF;x2=redpeak_end*srF;p1=plot([x1 x2],[ov_maxo ov_maxo],'-','Color','r','LineWidth',i_l);
hold on; x1= bluepeak_start*srF;x2=bluepeak_end*srF; hold on;p1=plot([x1 x2],[ov_maxo ov_maxo],'-','Color','b','LineWidth',i_l);
hold on;title('contra silent'); 

%ipsi silent
sb5=subplot(1,5,5);  
if data(ocat5(cat5)).experimentator=='SW'
    srF=1;
else data(ocat5(cat5)).experimentator=='MF'  
    srF=2;
end
plot(data(ocat5(cat5)).step_red.ephys_traces_70(:,trace_nr,2),'Color','k','LineWidth',tr_l);box off;axis off
hold on;plot(data(ocat5(cat5)).step_red.ephys_traces_40(:,trace_nr,2)+yoffset,'Color',[0.5 0.5 0.5],'LineWidth',tr_l);box off;axis off
set(sb5,'ylim',[ov_min ov_maxo]);
hold on; x1= redpeak_start*srF;x2=redpeak_end*srF;p1=plot([x1 x2],[ov_maxo ov_maxo],'-','Color','r','LineWidth',i_l);
hold on; x1= bluepeak_start*srF;x2=bluepeak_end*srF; hold on;p1=plot([x1 x2],[ov_maxo ov_maxo],'-','Color','b','LineWidth',i_l);
hold on;title('ipsi silent'); 

%scale barx for all
hold on;x1= 900*srF;x2=1000*srF;p1=plot([x1 x2],[ov_min ov_min],'-','Color','k','LineWidth',2);
%scale bary for all
hold on;y2= ov_min+scale_y;y1=ov_min;p1=plot([x2 x2],[y1 y2],'-','Color','k','LineWidth',2); 



%% Example ramp for display purpose, atm binocular cells: read out data 
%Blue AMPA/NMDA
%choose cell %potrntial: 19,27;--30--,--33--,--39--
cellnr=2;
type_c=r_nm2;
%Concatenate all 11 traces AMPA blue only
temp_cat=data(type_c(cellnr)).step_red.ephys_traces_70(:,:,1);
cat_bluesteps_a=temp_cat(:);
%Concatenate all 11 traces AMPA red with blue
temp_cat=data(type_c(cellnr)).step_red.ephys_traces_70(:,:,2);
cat_redsteps_a=temp_cat(:);
%Concatenate all 11 traces NMDA blue only
temp_cat=data(type_c(cellnr)).step_red.ephys_traces_40(:,:,1);
cat_bluesteps_n=temp_cat(:);
%Concatenate all 11 traces NMDA blue only
temp_cat=data(type_c(cellnr)).step_red.ephys_traces_40(:,:,2);
cat_redsteps_n=temp_cat(:);
%% Read out irradiance and extrapolate for SW
irr_red2=irradiance_extrapol(data,[1:135])
irr_blue1=data(type_c(cellnr)).step_blue.neg_irr_blue(1,:);
irr_blue2=data(type_c(cellnr)).step_blue.neg_irr_blue(2,:);
irr_red1=data(type_c(cellnr)).step_red.neg_irr_red(1,:);
%% 
% Display figure for example ramp
%FIGURE for paper
if data(type_c(cellnr)).experimentator=='SW'
    srF=1;
else data(type_c(cellnr)).experimentator=='MF'  
    srF=2;
end
fig7= figure;set(fig7, 'Name', 'Approach step display');set(fig7, 'Position', [200, 100, 600, 400]);set(gcf,'color','w');
%Plot blue only 
subplot(4,1,1);
plot(cat_bluesteps_a,'Color','k','LineWidth',1);box off;axis off;ylim([min(cat_bluesteps_a) abs(min(cat_bluesteps_a))/2]);
hold on;co=[0:1:10];
for i=1:11
x1(i)=bluepeak_start*srF+length(cat_bluesteps_a)/11*co(i);
x2(i)=bluepeak_end*srF+length(cat_bluesteps_a)/11*co(i);
p1=plot([x1(i) x2(i)],[abs(min(cat_bluesteps_a))/2-(abs(min(cat_bluesteps_a))/2)/3 abs(min(cat_bluesteps_a))/2-(abs(min(cat_bluesteps_a))/2)/3],'-','Color','b','LineWidth',4);
end
subplot(4,1,2);
for i=1:11
  plot(x1(i),irr_blue1(i),'sw','MarkerFaceColor','w','MarkerSize',5);box off
hold on;
plot([x1(i) x1(i)],[0 irr_blue1(i)],'Color','b','LineWidth',1.5);axis off
end

%Plot red/blue 
subplot(4,1,3);
plot(cat_redsteps_a,'Color','k','LineWidth',1);box off;axis off;ylim([min(cat_redsteps_a) abs(min(cat_redsteps_a))/2]);
hold on;co=[0:1:10];
for i=1:11
x1(i)=bluepeak_start*srF+length(cat_bluesteps_a)/11*co(i);
x2(i)=bluepeak_end*srF+length(cat_bluesteps_a)/11*co(i);

p1=plot([x1(i) x2(i)],[abs(min(cat_redsteps_a))/2-(abs(min(cat_redsteps_a))/2)/3 abs(min(cat_redsteps_a))/2-(abs(min(cat_redsteps_a))/2)/3],'-','Color','b','LineWidth',4);
xr1(i)=redpeak_start*srF+length(cat_redsteps_a)/11*co(i);
xr2(i)=redpeak_end*srF+length(cat_redsteps_a)/11*co(i);
p2=plot([xr1(i) xr2(i)],[abs(min(cat_redsteps_a))/2-(abs(min(cat_redsteps_a))/2)/3 abs(min(cat_redsteps_a))/2-(abs(min(cat_redsteps_a))/2)/3],'-','Color','r','LineWidth',4);
end
text(0,abs(min(cat_redsteps_a))/1.5,'637 nm','Color','r');text(1500*srF,abs(min(cat_redsteps_a))/1.5,'472 nm','Color','b');

%Global Scale bars
scale_x= 100;
scale_y= 400;
ov_min=min([min(cat_redsteps_a) min(cat_bluesteps_a)]);
%scale barx
hold on;xm1= 10900*srF;xm2=11000*srF;p1=plot([xm1 xm2],[ov_min ov_min],'-','Color','k','LineWidth',1.5);
%scale bary
hold on;ym2= ov_min+scale_y;ym1=ov_min;p1=plot([xm2 xm2],[ym1 ym2],'-','Color','k','LineWidth',1.5); 

subplot(4,1,4);
for i=1:11
  plot(x1(i),irr_blue2(i),'sw','MarkerFaceColor','w','MarkerSize',5);box off
hold on;
plot([x1(i) x1(i)],[0 irr_blue2(i)],'Color','b','LineWidth',1.5);axis off
hold on;
  plot(xr1(i),irr_red2(i),'sw','MarkerFaceColor','w','MarkerSize',5);box off
  hold on;
plot([xr1(i) xr1(i)],[0 irr_red2(i)],'Color','r','LineWidth',1.5);axis off
end
plot(xr1(1:2),irr_red2(1:2),'sw','MarkerFaceColor','r','MarkerSize',5);box off
%% ODI distributions and category plot
%FIGURE for paper
%Check the categories again: false positive rate is high,
%overrepresentation of bino cells I guess
%NO MD only atm
no_md=length([no_md1 no_md2 no_md3 no_md4 no_md5]);
fig2= figure;set(fig2, 'Name', 'ODI/category histogram');set(fig2, 'Position', [200, 600, 1000, 300]);set(gcf,'color','w');
%Step AMPA/NMDA
subplot(1,3,1);
b=bar([1 2 3 4 5],[length(no_md1)/no_md length(no_md2)/no_md length(no_md3)/no_md length(no_md4)/no_md length(no_md5)/no_md]*100);
box off;ylabel('Fraction of cells (%)');
b.FaceColor = 'flat';ylim([0 60]);
b.CData(2,:) = [1 0 0];b.CData(3,:) = [1 1 1]; b.CData(5,:) = [1 0 0];
yticks([0:10:60]); xticklabels({'contra','ipsi','bino','ipsi silent','contra silent'});set(gca,'FontSize',10)
xtickangle(45);set(gca,'FontSize',10);
%Histograms of ODIs AMPA
subplot(1,3,2);
histogram([data([no_md1 no_md2 no_md3 no_md4 no_md5]).ODI_AMPA_step], 7, 'Normalization','probability','FaceColor',[0.7 0.7 0.7],'EdgeColor',[0 0 0]);
%legend('AMPA');legend boxoff
hold on;title('AMPA');set(gca,'FontSize',10);
ylim([0 0.6]);yticks([0:0.2:0.6]); 
yticklabels({'0','20','40','60','80'});ylabel('Fraction of cells (%)');xlabel('ODI');box off;
%Histograms of ODIs NMDA
subplot(1,3,3);
histogram([data([no_md1 no_md2 no_md3 no_md4 no_md5]).ODI_NMDA_step], 7, 'Normalization','probability','FaceColor',[0.7 0.7 0.7],'EdgeColor',[0 0 0]);
%legend('NMDA');legend boxoff
ylim([0 0.6]);yticks([0:0.2:0.6]); 
yticklabels({'0','20','40','60','80'});ylabel('Fraction of cells (%)');xlabel('ODI');box off;
hold on;title('NMDA');set(gca,'FontSize',10);
%% Step protocol visualization prep
 rm={r_nm1 r_nm2 r_nm3 r_nm4 r_nm5};
 im={i_nm1 i_nm2 i_nm3 i_nm4 i_nm5};
 %% 
 r_idx=[r_nm1 r_nm2 r_nm3 r_nm4 r_nm5];
 
 %% 
 
 saturation_plots([i_nm1 i_nm3 i_nm4],[i_nm2 i_nm3 i_nm5],data,1);
%% Step protocol plots AMPA
tit={'contra only AMPA','ipsi only','bino','ipsi silent','contra silent'};
 fig4= figure;set(fig4, 'Name', 'Step protocol');set(fig4, 'Position', [200, 100, 1400, 200]);set(gcf,'color','w');
 for j=1:5
 subplot(1,5,j);
 title(tit{j});
for i=1:length(rm{j})
    tem=rm{j};
    if j==1 | j==4
    red_contra=data(tem(i)).step_red.neg_peak1(2,:);
    red_contra_norm=red_contra/max(abs(red_contra));
    r_c_n(:,i)=red_contra_norm;
    red_c_norm{:,j,i}=red_contra_norm;
    r_c_mean=nanmean(r_c_n,2);
    r_c_sem=nanstd(r_c_n,1,2)/sqrt(length(tem));
    all_rc(:,j)=r_c_mean;
    all_rc_sem(:,j)=r_c_sem;
    %Plot individual curves that are contra and red
    plot(abs(red_contra_norm),'--r.');hold on;  box off;
    elseif j==2 | j==5
    green_ipsi=data(tem(i)).step_blue.neg_peak2(1,:); 
    green_ipsi_norm=green_ipsi/max(abs(green_ipsi));
    g_i_n(:,i)=green_ipsi_norm;
    green_i_norm{:,j,i}=green_ipsi_norm;
    g_i_mean=nanmean(g_i_n,2);
    g_i_sem=nanstd(g_i_n,1,2)/sqrt(length(tem));
    all_gi(:,j)=g_i_mean;
    all_gi_sem(:,j)=g_i_sem;
   %Plot individual curves that are ipsi and green
     plot(abs(green_ipsi_norm),'--g.');hold on; box off;
    else j==3
    red_contra=data(tem(i)).step_red.neg_peak1(2,:);
    red_contra_norm=red_contra/max(abs(red_contra));
    red_c_norm{:,j,i}=red_contra_norm;
    r_c_n(:,i)=red_contra_norm;
    r_c_mean=nanmean(r_c_n,2);
    r_c_sem=nanstd(r_c_n,1,2)/sqrt(length(tem));
    all_rc(:,j)=r_c_mean;
    all_rc_sem(:,j)=r_c_sem;
    plot(abs(red_contra_norm),'--r.');hold on;  box off;
    green_ipsi=data(tem(i)).step_blue.neg_peak2(1,:); 
    green_ipsi_norm=green_ipsi/max(abs(green_ipsi));
    g_i_n(:,i)=green_ipsi_norm;
    green_i_norm{:,j,i}=green_ipsi_norm;
    g_i_mean=nanmean(g_i_n,2);
    g_i_sem=nanstd(g_i_n,1,2)/sqrt(length(tem));
    all_gi(:,j)=g_i_mean;
    all_gi_sem(:,j)=g_i_sem;
    plot(abs(green_ipsi_norm),'--g.');hold on; box off;
    end
    ylabel('Peak normalized');xlabel('Step nummber');xticks([1:1:11])        
end
hold on;
 title(tit{j});
for i=1:length(im{j})
    tem_g=im{j};
    if j==1 | j==4
    green_contra=data(tem_g(i)).step_blue.neg_peak2(1,:);
    green_contra_norm=green_contra/max(abs(green_contra));
    g_c_n(:,i)=green_contra_norm;
    green_c_norm{:,j,i}=green_contra_norm;
    g_c_mean=nanmean(g_c_n,2);
    g_c_sem=nanstd(g_c_n,1,2)/sqrt(length(tem));
    all_gc(:,j)=g_c_mean;
    all_gc_sem(:,j)=g_c_sem;
    %Plot individual curves that are contra and green
    plot(abs(green_contra_norm),'--g.');hold on; box off;
    elseif j==2 | j==5
    red_ipsi=data(tem_g(i)).step_red.neg_peak1(2,:);
    red_ipsi_norm=red_ipsi/max(abs(red_ipsi));
    r_i_n(:,i)=red_ipsi_norm;
     red_i_norm{:,j,i}=green_ipsi_norm;
    r_i_mean=nanmean(r_i_n,2);
    r_i_sem=nanstd(r_i_n,1,2)/sqrt(length(tem));
    all_ri(:,j)=r_i_mean;
    all_ri_sem(:,j)=r_i_sem;
    %Plot individual curves that are ipsi and red
    plot(abs(red_ipsi_norm),'--r.');hold on;  box off;
    else j==3
    green_contra=data(tem_g(i)).step_blue.neg_peak2(1,:);
    green_contra_norm=green_contra/max(abs(green_contra));
    g_c_n(:,i)=green_contra_norm;
     green_c_norm{:,j,i}=green_contra_norm;
    g_c_mean=nanmean(g_c_n,2);
    g_c_sem=nanstd(g_c_n,1,2)/sqrt(length(tem));
    all_gc(:,j)=g_c_mean;
    all_gc_sem(:,j)=g_c_sem;
    %Plot individual curves that are contra and green
    plot(abs(green_contra_norm),'--g.');hold on; box off;
    red_ipsi=data(tem_g(i)).step_red.neg_peak1(2,:);
    red_ipsi_norm=red_ipsi/max(abs(red_ipsi));
    r_i_n(:,i)=red_ipsi_norm;
    red_i_norm{:,j,i}=red_ipsi_norm;
    r_i_mean=nanmean(r_c_n,2);
    r_i_sem=nanstd(r_i_n,1,2)/sqrt(length(tem));
    all_ri(:,j)=r_i_mean;
    all_ri_sem(:,j)=r_i_sem;
    %Plot individual curves that are ipsi and red
    plot(abs(red_ipsi_norm),'--r.');hold on;  box off;
    end
    ylabel('Peak normalized');xlabel('Step nummber');xticks([1:1:11]);
end
 end

 %% Plot Contra red (Chrimson) and Contra green (Chronos) irrespective of category AMPA
 %FIGURE for paper
 
 set(0,'DefaultLegendAutoUpdate','off')
%Contra red
fig6= figure;set(fig6, 'Name', 'Step protocol');set(fig6, 'Position', [200, 100, 400, 400]);set(gcf,'color','w');
subplot(2,1,1);
hold on;title('Contra','Color','b');
me=abs(nanmean(reshape([red_c_norm{:}],11,length([red_c_norm{:}])/11),2));
st=nanstd(reshape([red_c_norm{:}],11,length([red_c_norm{:}])/11),1,2)/sqrt(length([red_c_norm{:}])/11);
itr=abs(reshape([red_c_norm{:}],11,length([red_c_norm{:}])/11))
for i=1:length(itr)
exp=plot(itr(:,i),'-r')
exp.Color(4) = 0.05;
hold on;
end
hold on;
errorbar(me,st,'--or','MarkerSize',5,...
    'MarkerEdgeColor','k','MarkerFaceColor','r');hold on;
box off;ylabel('Peak normalized');xticks([1:1:11]);

% Contra green
hold on;
me=abs(nanmean(reshape([green_c_norm{:}],11,length([green_c_norm{:}])/11),2));
st=nanstd(reshape([green_c_norm{:}],11,length([green_c_norm{:}])/11),1,2)/sqrt(length([green_c_norm{:}])/11);
itg=abs(reshape([green_c_norm{:}],11,length([green_c_norm{:}])/11))
for i=1:length(itg)
exp=plot(itg(:,i),'-g')
exp.Color(4) = 0.2;
hold on;
end
errorbar(me,st,'--og','MarkerSize',5,...
    'MarkerEdgeColor','k','MarkerFaceColor','g');hold on;
box off;ylabel('Peak normalized');xticks([1:1:11]);
set(gca,'FontSize',10);

subplot(2,1,2);
hold on;title('Ipsi','Color','r');
%Ipsi red
me=abs(nanmean(reshape([red_i_norm{:}],11,length([red_i_norm{:}])/11),2));
st=nanstd(reshape([red_i_norm{:}],11,length([red_i_norm{:}])/11),1,2)/sqrt(length([red_i_norm{:}])/11);
hold on;
errorbar(me,st,'--or','MarkerSize',5,...
    'MarkerEdgeColor','k','MarkerFaceColor','r');hold on;
box off;ylabel('Peak normalized');xlabel('Step number');xticks([1:1:11]);

% Ipsi green
me=abs(nanmean(reshape([green_i_norm{:}],11,length([green_i_norm{:}])/11),2));
st=nanstd(reshape([green_i_norm{:}],11,length([green_i_norm{:}])/11),1,2)/sqrt(length([green_i_norm{:}])/11);
hold on;
errorbar(me,st,'--og','MarkerSize',5,...
    'MarkerEdgeColor','k','MarkerFaceColor','g');hold on;
box off;ylabel('Peak normalized');xlabel('Step number');xticks([1:1:11]);
legend('Chrimson','Chronos','Location','southeast');
legend boxoff;set(gca,'FontSize',10);

%Ipsi red individuals
itr_i=abs(reshape([red_i_norm{:}],11,length([red_i_norm{:}])/11))
for i=1:length(itr_i)
exp=plot(itr_i(:,i),'-r')
exp.Color(4) = 0.1;
hold on;
end
% Ipsi green individuals
itg_i=abs(reshape([green_i_norm{:}],11,length([green_i_norm{:}])/11))
for i=1:length(itg_i)
exp=plot(itg_i(:,i),'-g')
exp.Color(4) = 0.05;
hold on;
end


%% 
% for i=1:length(itr)
%   c_itr(i)=find(itr(:,i)==1)  
% end
% for i=1:length(itg)
%   c_itg(i)=find(itg(:,i)==1)  
% end
% for i=1:length(itr_i)
%   c_itr_i(i)=find(itr_i(:,i)==1)  
% end
% for i=1:length(itg_i)
%   c_itg_i(i)=find(itg_i(:,i)==1)  
% end
% figure;histogram(c_itr);hold on;histogram(c_itg);hold on;histogram(c_itr_i);hold on;histogram(c_itg_i)

%% %% Step protocol plots NMDA
red_contra=[];
red_contra_norm=[];
r_c_mean=[];
r_c_sem=[];
green_ipsi=[];
green_contra=[];
red_ipsi=[];

tit={'contra only NMDA','ipsi only','bino','ipsi silent','contra silent'};
fig5= figure;set(fig5, 'Name', 'Step protocol');set(fig5, 'Position', [200, 100, 1400, 200]);set(gcf,'color','w');
 for j=1:5
 subplot(1,5,j);
 title(tit{j});
for i=1:length(rm{j})
    tem=rm{j};
    if j==1 | j==4
    if isempty(data(tem(i)).step_red.pos_peak1(2,:))==0   
    red_contra=data(tem(i)).step_red.pos_peak1(2,:);
     else 
    red_contra=ones(1,11)*NaN;
    end
    red_contra_norm=red_contra/max(abs(red_contra));
    r_c_n(:,i)=red_contra_norm;
    red_c_norm{:,j,i}=red_contra_norm;
    r_c_mean=nanmean(r_c_n,2);
    r_c_sem=nanstd(r_c_n,1,2)/sqrt(length(tem));
    all_rc(:,j)=r_c_mean;
    all_rc_sem(:,j)=r_c_sem;
   
    %Plot individual curves that are contra and red
    plot(abs(red_contra_norm),'--r.');hold on;  box off;
    elseif j==2 | j==5
     if isempty(data(tem(i)).step_blue.pos_peak2(1,:))==0      
    green_ipsi=data(tem(i)).step_blue.pos_peak2(1,:); 
     else
     green_ipsi=ones(1,11)*NaN;
     end
    green_ipsi_norm=green_ipsi/max(abs(green_ipsi));
    g_i_n(:,i)=green_ipsi_norm;
    green_i_norm{:,j,i}=green_ipsi_norm;
    g_i_mean=nanmean(g_i_n,2);
    g_i_sem=nanstd(g_i_n,1,2)/sqrt(length(tem));
    all_gi(:,j)=g_i_mean;
    all_gi_sem(:,j)=g_i_sem;
   %Plot individual curves that are ipsi and green
     plot(abs(green_ipsi_norm),'--g.');hold on; box off;
    else j==3
        if isempty(data(tem(i)).step_red.pos_peak1)==0
      if isempty(data(tem(i)).step_red.pos_peak1(2,:))==0   
    red_contra=data(tem(i)).step_red.pos_peak1(2,:);
     else 
    red_contra=ones(1,11)*NaN;
      end 
        else
            red_contra=ones(1,11)*NaN;
        end
   
    red_contra_norm=red_contra/max(abs(red_contra));
    red_c_norm{:,j,i}=red_contra_norm;
    r_c_n(:,i)=red_contra_norm;
    r_c_mean=nanmean(r_c_n,2);
    r_c_sem=nanstd(r_c_n,1,2)/sqrt(length(tem));
    all_rc(:,j)=r_c_mean;
    all_rc_sem(:,j)=r_c_sem;
    plot(abs(red_contra_norm),'--r.');hold on;  box off;
    
     if isempty(data(tem(i)).step_blue.pos_peak2)==0
      if isempty(data(tem(i)).step_blue.pos_peak2(1,:))==0   
      green_ipsi=data(tem(i)).step_blue.pos_peak2(1,:); 
     else 
    green_ipsi=ones(1,11)*NaN;
      end 
        else
            green_ipsi=ones(1,11)*NaN;
        end
  
    green_ipsi_norm=green_ipsi/max(abs(green_ipsi));
    g_i_n(:,i)=green_ipsi_norm;
    green_i_norm{:,j,i}=green_ipsi_norm;
    g_i_mean=nanmean(g_i_n,2);
    g_i_sem=nanstd(g_i_n,1,2)/sqrt(length(tem));
    all_gi(:,j)=g_i_mean;
    all_gi_sem(:,j)=g_i_sem;
    plot(abs(green_ipsi_norm),'--g.');hold on; box off;
        
    end
    ylabel('Peak normalized');xlabel('Step nummber');xticks([1:1:11])        
end


hold on;
 title(tit{j});
for i=1:length(im{j})
    tem_g=im{j};
    if j==1 | j==4
         if isempty(data(tem_g(i)).step_blue.pos_peak2(1,:))==0   
    green_contra=data(tem_g(i)).step_blue.pos_peak2(1,:);
     else 
    green_contra=ones(1,11)*NaN;
         end
    green_contra_norm=green_contra/max(abs(green_contra));
    g_c_n(:,i)=green_contra_norm;
    green_c_norm{:,j,i}=green_contra_norm;
    g_c_mean=nanmean(g_c_n,2);
    g_c_sem=nanstd(g_c_n,1,2)/sqrt(length(tem));
    all_gc(:,j)=g_c_mean;
    all_gc_sem(:,j)=g_c_sem;
    %Plot individual curves that are contra and green
    plot(abs(green_contra_norm),'--g.');hold on; box off;
    elseif j==2 | j==5
       if isempty(data(tem_g(i)).step_red.pos_peak1(2,:))==0   
    red_ipsi=data(tem_g(i)).step_red.pos_peak1(2,:);
     else 
    red_ipsi=ones(1,11)*NaN;
    end    
   
    red_ipsi_norm=red_ipsi/max(abs(red_ipsi));
    r_i_n(:,i)=red_ipsi_norm;
     red_i_norm{:,j,i}=green_ipsi_norm;
    r_i_mean=nanmean(r_i_n,2);
    r_i_sem=nanstd(r_i_n,1,2)/sqrt(length(tem));
    all_ri(:,j)=r_i_mean;
    all_ri_sem(:,j)=r_i_sem;
    %Plot individual curves that are ipsi and red
    plot(abs(red_ipsi_norm),'--r.');hold on;  box off;
    else j==3
        if isempty(data(tem_g(i)).step_blue.pos_peak2)==0
         if isempty(data(tem_g(i)).step_blue.pos_peak2(1,:))==0   
    green_contra=data(tem_g(i)).step_blue.pos_peak2(1,:);
     else 
    green_contra=ones(1,11)*NaN;
    end 
        else
            green_contra=ones(1,11)*NaN;
        end
    green_contra_norm=green_contra/max(abs(green_contra));
    g_c_n(:,i)=green_contra_norm;
     green_c_norm{:,j,i}=green_contra_norm;
    g_c_mean=nanmean(g_c_n,2);
    g_c_sem=nanstd(g_c_n,1,2)/sqrt(length(tem));
    all_gc(:,j)=g_c_mean;
    all_gc_sem(:,j)=g_c_sem;
    %Plot individual curves that are contra and green
    plot(abs(green_contra_norm),'--g.');hold on; box off;
    
    if isempty(data(tem_g(i)).step_red.pos_peak1)==0
         if isempty(data(tem_g(i)).step_red.pos_peak1(2,:))==0   
     red_ipsi=data(tem_g(i)).step_red.pos_peak1(2,:);
     else 
    red_ipsi=ones(1,11)*NaN;
    end 
        else
            red_ipsi=ones(1,11)*NaN;
        end
   
    red_ipsi_norm=red_ipsi/max(abs(red_ipsi));
    r_i_n(:,i)=red_ipsi_norm;
    red_i_norm{:,j,i}=red_ipsi_norm;
    r_i_mean=nanmean(r_c_n,2);
    r_i_sem=nanstd(r_i_n,1,2)/sqrt(length(tem));
    all_ri(:,j)=r_i_mean;
    all_ri_sem(:,j)=r_i_sem;
    %Plot individual curves that are ipsi and red
    plot(abs(red_ipsi_norm),'--r.');hold on;  box off;
      
        end
   
    ylabel('Peak normalized');xlabel('Step nummber');xticks([1:1:11]);
end
 end

%% Plot Contra red (Chrimson) and Contra green (Chronos) irrespective of category NMDA
 %FIGURE for paper
 itr=[];
 set(0,'DefaultLegendAutoUpdate','off')
%Contra red
fig6= figure;set(fig6, 'Name', 'Step protocol');set(fig6, 'Position', [200, 100, 400, 400]);set(gcf,'color','w');
subplot(2,1,1);
hold on;title('Contra','Color','b');
me=nanmean(reshape([red_c_norm{:}],11,length([red_c_norm{:}])/11),2);
st=nanstd(reshape([red_c_norm{:}],11,length([red_c_norm{:}])/11),1,2)/sqrt(length([red_c_norm{:}])/11);
itr=reshape([red_c_norm{:}],11,length([red_c_norm{:}])/11)
for i=1:length(itr)
exp=plot(itr(:,i),'-r')
exp.Color(4) = 0.05;
hold on;
end
hold on;
errorbar(me,st,'--or','MarkerSize',5,...
    'MarkerEdgeColor','k','MarkerFaceColor','r');hold on;
box off;ylabel('Peak normalized');xticks([1:1:11]);

% Contra green
hold on;
me=abs(nanmean(reshape([green_c_norm{:}],11,length([green_c_norm{:}])/11),2));
st=nanstd(reshape([green_c_norm{:}],11,length([green_c_norm{:}])/11),1,2)/sqrt(length([green_c_norm{:}])/11);
itg=abs(reshape([green_c_norm{:}],11,length([green_c_norm{:}])/11))
for i=1:length(itg)
exp=plot(itg(:,i),'-g')
exp.Color(4) = 0.2;
hold on;
end
errorbar(me,st,'--og','MarkerSize',5,...
    'MarkerEdgeColor','k','MarkerFaceColor','g');hold on;
box off;ylabel('Peak normalized');xticks([1:1:11]);
set(gca,'FontSize',10);

subplot(2,1,2);
hold on;title('Ipsi','Color','r');
%Ipsi red
me=abs(nanmean(reshape([red_i_norm{:}],11,length([red_i_norm{:}])/11),2));
st=nanstd(reshape([red_i_norm{:}],11,length([red_i_norm{:}])/11),1,2)/sqrt(length([red_i_norm{:}])/11);
hold on;
errorbar(me,st,'--or','MarkerSize',5,...
    'MarkerEdgeColor','k','MarkerFaceColor','r');hold on;
box off;ylabel('Peak normalized');xlabel('Step number');xticks([1:1:11]);

% Ipsi green
me=abs(nanmean(reshape([green_i_norm{:}],11,length([green_i_norm{:}])/11),2));
st=nanstd(reshape([green_i_norm{:}],11,length([green_i_norm{:}])/11),1,2)/sqrt(length([green_i_norm{:}])/11);
hold on;
errorbar(me,st,'--og','MarkerSize',5,...
    'MarkerEdgeColor','k','MarkerFaceColor','g');hold on;
box off;ylabel('Peak normalized');xlabel('Step number');xticks([1:1:11]);
legend('Chrimson','Chronos','Location','southeast');
legend boxoff;set(gca,'FontSize',10);

%Ipsi red individuals
itr_i=abs(reshape([red_i_norm{:}],11,length([red_i_norm{:}])/11))
for i=1:length(itr_i)
exp=plot(itr_i(:,i),'-r')
exp.Color(4) = 0.1;
hold on;
end
% Ipsi green individuals
itg_i=abs(reshape([green_i_norm{:}],11,length([green_i_norm{:}])/11))
for i=1:length(itg_i)
exp=plot(itg_i(:,i),'-g')
exp.Color(4) = 0.05;
hold on;
end

%% Dominat vs nondominat eye bino cells read out
for i=1:length(no_md3)
temp= data(no_md3(i)).step_red.neg_peak1(2,:);
temp2= data(no_md3(i)).step_blue.neg_peak2(2,:);
c_i(:,i)=data(no_md3(i)).brain_contra_ipsi;
av_peak1(:,i)=nanmean(temp(6:end));
av_peak2(:,i)=nanmean(temp2(6:end));
end
com_peak=[av_peak1;av_peak2];
[dom b]=find(com_peak==min(com_peak));
[ndom b]=find(com_peak==max(com_peak));
for i=1:length(no_md3);    
    DE(i)=com_peak(dom(i),i);
    NDE(i)=com_peak(ndom(i),i);
end
%% 
for i=1:length(no_md3)
    if c_i(i)==1
de_c_i(i)=dom(i)
    else
   de_c_i(i)=ndom(i)    
    end
    
end
%% Boxplot and plotspread plot for DE vs NDE
 %FIGURE for paper
figure; set(gcf,'color','w');
plotSpread(abs([DE' NDE'])/1000,'categoryIdx',[ones(1,length(DE))' ones(1,length(NDE))'*2], 'categoryMarkers',{'.','.'},'categoryColors',{'k','k'});
hold on;
con=[DE NDE];
grp=[ones(1,length(DE)) ones(1,length(NDE))*2];
boxplot(abs(con)/1000,grp,'Colors','k','OutlierSize',0.0001);
box off;xticklabels({'DE','NDE'});ylabel('Peak amplitude (nA)');
yticks([0:1:4]);
ylim([0 4]);
set(gca,'FontSize',10);
%% Scatter plot with refline
pointsize=25;
[cmap]=buildcmap('br');
figure;set(gcf,'color','w');scatter(abs(DE)/1000,abs(NDE)/1000,pointsize,de_c_i,'filled');
colormap(cmap)
xlabel('DE synaptic input (nA)');
ylabel('NDE synaptic input (nA)');
set(gca,'FontSize',12);
hold on;
r=refline(1,0);r.Color='k'
xlim([0 4]);
ylim([0 4]);
yticks([0:1:4]);
hold on;
scatter(0.5,3.8,'bo','filled');text(0.6,3.82,'contra');
hold on;
scatter(0.5,3.6,'ro','filled');text(0.6,3.62,'ipsi');
%% 
 edge=[0:0.25:3.5];
figure;h=histogram(abs(DE)/1000-abs(NDE)/1000,edge)
h.FaceColor=[0.5 0.5 0.5]
xlim([0 3.5])
xticks([0:0.5:3.5]);
set(gcf,'color','w');
box off;
ylabel('Counts');
ylim([0 20]);
yticks([0:5:20]);
xlabel('\Delta DE - NDE (nA)');
set(gca,'FontSize',12);
%% 

%% Chrimson only FIGURES: Data_Chrimson_only_MF.mat NEEDED

                         
 %% Sequential photostimulation prevents crosstalk: this uses 26 cells that MF patched with Chrimson only
%ChrimsonR only control
% 1 Blue ramp                                       establish blue threshold
% 2 red ramp no blue                           red should elicit a response but is then killed
% 3 red ramp with blue                         red should elicit a response but is not kiled (compare to 2)
% 4 double rep no blue                         This is the quantification of red rundown w/o blue
% 5 double rep with blue                      compare this to 4.
% 6 red ramp 10 no blue                       this ramp has 10ms red, no blue (response dies)
% 7 red ramp 10 with blue                     this ramp has 10ms red, 50ms blue (compare to 6)
% 8 red ramp 1 no blue                          this ramp has 1ms red, no blue (response dies)
% 9 red ramp 1 with blue                    this ramp has 1ms red, 50ms blue (compare to 8)
% 10 minimal red with blue                no explanation needed... should get minimal stim.
% 11 minimal red no blue                   no explanation needed, red should die, same power as 10
% 12 double rep with blue                 just to reset
% 13 minimal blue only                      using same irradiance for blue as used for red in 10
% 14 blue delay                                  250 red, 50 blue, window betwen red and blue for crosstalk depletion
% 15 red time                                      How long does the red need to be to get proper crosstalk depletion
str_LGN     = 'C:\Users\Simon-localadmin\Documents\MargrieLab\PhDprojects\LGN\AnalyzedData\';
folder_list = uipickfiles('FilterSpec',str_LGN);
load(char(folder_list));
%%  %FIGURE for paper: CROSSTALK PREVENTION
tr_l=1.2;
%linewidth of red blue indicator
i_l=5;
%which trace from the 11 steps 
trace_nr=7;
%Choose wether red should be contra or ipsi
contra_red=1;
base_start          =   1;
base_end            =   99;
redpeak_start       =   100;
redpeak_end         =   350;%looking for the peak shoud be the same as for blue even though the window is longer?
bluepeak_start      =   351;
bluepeak_end        =   400;
scale_x= 100;
scale_y= 400;
%Example trace for sequential photostimulation works in red only expressing
%cells
%good examples= 4, 7, 12, 16 maybe 1
cellnr=7;
if LGN(cellnr).experimentator=='SW'
    srF=1;
else LGN(cellnr).experimentator=='MF'  
    srF=2;
end
temp_size=size(LGN(cellnr).step_red.ephys_traces_70(:,:,:),3);
blue_o=LGN(cellnr).step_red.ephys_traces_70(:,11,1);
blue_red=LGN(cellnr).step_red.ephys_traces_70(:,3,3);
%Plotting
fig8= figure;set(fig8, 'Name', 'Step protocol');set(fig8, 'Position', [200, 100, 200, 400]);set(gcf,'color','w');
subplot(2,1,1);
ov_maxo=max(blue_o);
plot(blue_o,'Color','k','LineWidth',1);box off;axis off;
%hold on; x1= redpeak_start*srF;x2=redpeak_end*srF;p1=plot([x1 x2],[ov_maxo ov_maxo],'-','Color','r','LineWidth',i_l);
hold on; x1= bluepeak_start*srF;x2=bluepeak_end*srF; hold on;p1=plot([x1 x2],[ov_maxo*100 ov_maxo*100],'-','Color','b','LineWidth',i_l);
hold on;text(x1,ov_maxo*180,'472 nm','Color','b','FontSize',9);
ylim([min(blue_red) ov_maxo*180]);
    
subplot(2,1,2);
ov_maxo=max(blue_red);
plot(blue_red,'Color','k','LineWidth',1);box off;axis off;
hold on; x1= redpeak_start*srF;x2=redpeak_end*srF;p1=plot([x1 x2],[ov_maxo*100 ov_maxo*100],'-','Color','r','LineWidth',i_l);
hold on;text(x1,ov_maxo*180,'637','Color','r','FontSize',8);
hold on; x1= bluepeak_start*srF;x2=bluepeak_end*srF; hold on;p1=plot([x1 x2],[ov_maxo*100 ov_maxo*100],'-','Color','b','LineWidth',i_l);
hold on;text(x1,ov_maxo*180,'472 nm','Color','b','FontSize',8);
%ylim([min(blue_red) ov_maxo*180]);

ov_min=min(blue_red);
%scale barx
hold on;x1= 900*srF;x2=1000*srF;p1=plot([x1 x2],[ov_min*1.5 ov_min*1.5],'-','Color','k','LineWidth',1);
%scale bary
hold on;y2= ov_min+scale_y;y1=ov_min;p1=plot([x1 x1],[y1*1.5 y2*1.5],'-','Color','k','LineWidth',1); 
%
%% Compare red ramp with blue and no blue, red varies from 250, 10 and 1 ms, blue is always 50 ms
srF=2;
for m=1:length(LGN)
try
temp=LGN(m).step_red.ephys_traces_70(:,:,3);
red250b(:,m)=temp(:);
peak250b(:,m)=LGN(m).step_red.neg_peak1(3,:);
irr250b_r(:,m)=LGN(m).step_red.neg_irr_red(3,:);  
irr250b_b(:,m)=LGN(m).step_blue.neg_irr_blue(3,:); 
end
%10ms
try
temp2=LGN(m).step_red.ephys_traces_70(:,:,5);
red10b(:,m)=temp2(:);
peak10b(:,m)=LGN(m).step_red.neg_peak1(5,:);
irr10b_r(:,m)=(104.1*mean(LGN(m).step_red.bs_photodiodetraces_70(102*srF:105*srF,:,5))-3.467)/100;
irr10b_b(:,m)=(679.2*mean(LGN(m).step_red.bs_photodiodetraces_70(125*srF:150*srF,:,5))-26.82)/100;
end
try
temp3=LGN(m).step_red.ephys_traces_70(:,:,7);
red1b(:,m)=temp3(:);
peak1b(:,m)=LGN(m).step_red.neg_peak1(7,:);
irr1b_r(:,m)=(104.1*mean(LGN(m).step_red.bs_photodiodetraces_70(101*srF:101.5*srF,:,7))-3.467)/100;
irr1b_b(:,m)=(679.2*mean(LGN(m).step_red.bs_photodiodetraces_70(125*srF:150*srF,:,7))-26.82)/100;
 end
end
%% no blue
for m=1:length(LGN)
try
temp=LGN(m).step_red.ephys_traces_70(:,:,2);
red250bn(:,m)=temp(:);
peak250bn(:,m)=LGN(m).step_red.neg_peak1(2,:);
irr250bn_r(:,m)=LGN(m).step_red.neg_irr_red(2,:);  
irr250bn_b(:,m)=LGN(m).step_blue.neg_irr_blue(2,:);
end
try
temp2=LGN(m).step_red.ephys_traces_70(:,:,4);
red10bn(:,m)=temp2(:);
peak10bn(:,m)=LGN(m).step_red.neg_peak1(4,:);
irr10bn_r(:,m)=(104.1*mean(LGN(m).step_red.bs_photodiodetraces_70(102*srF:105*srF,:,4))-3.467)/100;
irr10bn_b(:,m)=(679.2*mean(LGN(m).step_red.bs_photodiodetraces_70(125*srF:150*srF,:,4))-26.82)/100;
end
try
temp3=LGN(m).step_red.ephys_traces_70(:,:,6);
red1bn(:,m)=temp3(:);
peak1bn(:,m)=LGN(m).step_red.neg_peak1(6,:);
irr1bn_r(:,m)=(104.1*mean(LGN(m).step_red.bs_photodiodetraces_70(101*srF:101.5*srF,:,6))-3.467)/100; 
irr1bn_b(:,m)=(679.2*mean(LGN(m).step_red.bs_photodiodetraces_70(125*srF:150*srF,:,6))-26.82)/100;
 end
end

%% READ OUT concatenated traces WITH BLUE
cell_ex=2;
trac_ex=[red250b(:,cell_ex) red10b(:,cell_ex) red1b(:,cell_ex)];
irr_exr=[irr250b_r(:,cell_ex) irr10b_r(:,cell_ex) irr1b_r(:,cell_ex)];
irr_exb=[irr250b_b(:,cell_ex) irr10b_b(:,cell_ex) irr1b_b(:,cell_ex)];
redpeak_start=100;
redpeak_end=[350 120 110];%SET TO LONGER VALUES SO THAT IT IS VISIBLE IN THE PLOT!!!!
bluepeak_end=[400 160 151];
tex={'red 250ms','red 10ms','red 1ms'};
tr1={'with blue',' ',' '};
%% Plotting
 %FIGURE for paper
for j=1:3
fig10= figure;set(fig10, 'Name', 'Control Chrimson');set(fig10, 'Position', [400+400*j, 400, 400, 300]);set(gcf,'color','w');
subplot(2,1,1);
plot(trac_ex(:,j),'Color','k','LineWidth',1);box off;axis off;%ylim([min(cat_bluesteps_a) abs(min(cat_bluesteps_a))/2]);
hold on; title(tex{j});
co=[0:1:10];
hold on;
for i=1:11
x1(i)=redpeak_end(j)*srF+length(trac_ex(:,j))/11*co(i);
x2(i)=bluepeak_end(j)*srF+length(trac_ex(:,j))/11*co(i);
p1=plot([x1(i) x2(i)],[abs(min(trac_ex(:,j)))/2-(abs(min(trac_ex(:,j)))/2)/3 abs(min(trac_ex(:,j)))/2-(abs(min(trac_ex(:,j)))/2)/3],'-','Color','b','LineWidth',4);
xr1(i)=redpeak_start*srF+length(trac_ex(:,j))/11*co(i);
xr2(i)=redpeak_end(j)*srF+length(trac_ex(:,j))/11*co(i);
p2=plot([xr1(i) xr2(i)],[abs(min(trac_ex(:,j)))/2-(abs(min(trac_ex(:,j)))/2)/3 abs(min(trac_ex(:,j)))/2-(abs(min(trac_ex(:,j)))/2)/3],'-','Color','r','LineWidth',4);
end
subplot(2,1,2);
scatter(xr1,irr_exr(:,j)','sr');axis off;
hold on;
scatter(x1,irr_exb(:,j)','sb');axis off;
hold on; title('with blue');
%'-or','MarkerSize',5,'MarkerEdgeColor','r','MarkerFaceColor','r');axis off
end
%% Choose an example cell WITHOUT BLUE
cell_ex=2;
trac_ex_n=[red250bn(:,cell_ex) red10bn(:,cell_ex) red1bn(:,cell_ex)];
irr_exr_n=[irr250bn_r(:,cell_ex) irr10bn_r(:,cell_ex) irr1bn_r(:,cell_ex)];
irr_exb_n=[irr250bn_b(:,cell_ex) irr10bn_b(:,cell_ex) irr1bn_b(:,cell_ex)];
redpeak_start=100;
redpeak_end=[350 120 110];%SET TO LONGER VALUES SO THAT IT IS VISIBLE IN THE PLOT!!!!
bluepeak_end=[400 160 151];
%% Plotting
 %FIGURE for paper
for j=1:3
fig10= figure;set(fig10, 'Name', 'Step protocol');set(fig10, 'Position', [400+400*j, 400, 400, 300]);set(gcf,'color','w');
hold on; title('without blue');
subplot(2,1,1);
plot(trac_ex_n(:,j),'Color','k','LineWidth',1);box off;axis off;%ylim([min(cat_bluesteps_a) abs(min(cat_bluesteps_a))/2]);
hold on;title(tex{j});
co=[0:1:10];
hold on;
for i=1:11
x1(i)=redpeak_end(j)*srF+length(trac_ex_n(:,j))/11*co(i);
x2(i)=bluepeak_end(j)*srF+length(trac_ex_n(:,j))/11*co(i);
%p1=plot([x1(i) x2(i)],[abs(min(trac_ex(:,j)))/2-(abs(min(trac_ex_n(:,j)))/2)/3 abs(min(trac_ex_n(:,j)))/2-(abs(min(trac_ex_n(:,j)))/2)/3],'-','Color','b','LineWidth',4);
xr1(i)=redpeak_start*srF+length(trac_ex_n(:,j))/11*co(i);
xr2(i)=redpeak_end(j)*srF+length(trac_ex_n(:,j))/11*co(i);
p2=plot([xr1(i) xr2(i)],[abs(min(trac_ex_n(:,j)))/2-(abs(min(trac_ex_n(:,j)))/2)/3 abs(min(trac_ex_n(:,j)))/2-(abs(min(trac_ex_n(:,j)))/2)/3],'-','Color','r','LineWidth',4);
end
subplot(2,1,2)
scatter(xr1,irr_exr_n(:,j)','sr');axis off;
hold on;
scatter(x1,irr_exb_n(:,j)','sb');axis off;
hold on; title('without blue');
%'-or','MarkerSize',5,'MarkerEdgeColor','r','MarkerFaceColor','r');axis off
end

%% One plot for an example cell: ALTERNATIVE TO THE ABOVE
%FIGURE for paper
ov_min=min(min(trac_ex));
ov_max=max(max(trac_ex)); 
fig11= figure;set(fig11, 'Name', 'Step protocol');set(fig11, 'Position', [200, 100, 800, 400]);set(gcf,'color','w');
pos=[1 3 5 7];
for j=1:3
hold on;
subplot(4,2,pos(j));
hold on;
%title(tr1{j});
plot(trac_ex(:,j),'Color','k','LineWidth',1);box off;axis off;%ylim([min(cat_bluesteps_a) abs(min(cat_bluesteps_a))/2]);
hold on;
co=[0:1:10];
for i=1:11
x1(i)=redpeak_end(j)*srF+length(trac_ex(:,j))/11*co(i);
x2(i)=bluepeak_end(j)*srF+length(trac_ex(:,j))/11*co(i);
p1=plot([x1(i) x2(i)],[abs(min(trac_ex(:,j)))/2-(abs(min(trac_ex(:,j)))/2)/3 abs(min(trac_ex(:,j)))/2-(abs(min(trac_ex(:,j)))/2)/3],'-','Color','b','LineWidth',4);
xr1(i)=redpeak_start*srF+length(trac_ex(:,j))/11*co(i);
xr2(i)=redpeak_end(j)*srF+length(trac_ex(:,j))/11*co(i);
p2=plot([xr1(i) xr2(i)],[abs(min(trac_ex(:,j)))/2-(abs(min(trac_ex(:,j)))/2)/3 abs(min(trac_ex(:,j)))/2-(abs(min(trac_ex(:,j)))/2)/3],'-','Color','r','LineWidth',4);
end
%title(tex{j},'Color','r');
end
%Irradiance
hold on;
subplot(4,2,pos(4));
for i=1:11
  plot(xr1(i),irr_exr(i,1)','sw','MarkerFaceColor','w','MarkerSize',5);box off
  hold on;
plot([xr1(i) xr1(i)],[0 irr_exr(i,1)'],'Color','r','LineWidth',1.5);axis off
end
hold on
plot(xr1(1:2),irr_exr(1:2,1)','sw','MarkerFaceColor','r','MarkerSize',5);box off
 hold on;
 for i=1:11
plot(x2(i),irr_exb(i,1)','-sb','MarkerFaceColor','w','MarkerEdgeColor','w','MarkerSize',5);box off
 hold on;
plot([x2(i) x2(i)],[0 irr_exb(i,1)'],'Color','b','LineWidth',1);axis off
end
% plot(xr1,irr_exr(:,1)','-sr','MarkerFaceColor','r','MarkerSize',5);
% %bar(xr1-1000,irr_exb(:,1)');
% %bar(xr1-1000,irr_exr(:,1)');
% hold on;
% plot(x1,irr_exb(:,1)','-sb','MarkerFaceColor','b','MarkerSize',5);box off
% ylabel('Irradiance (mW/mm^2)');
% xticks([x1(1):2000:x1(end)]);
% xticklabels({'1','2','3','4','5','6','7','8','9','10','11'});
% xlabel('Step');
% ylim([min(irr_exr(:,1)) 6]);
set(gca,'FontSize',10);  
hold on;
pos=[2 4 6 8];
for j=1:3
hold on;
subplot(4,2,pos(j));
plot(trac_ex_n(:,j),'Color','k','LineWidth',1);box off;axis off;%ylim([min(cat_bluesteps_a) abs(min(cat_bluesteps_a))/2]);
%title(tex{j},'Color','r');
hold on;
for i=1:11
x1(i)=redpeak_end(j)*srF+length(trac_ex_n(:,j))/11*co(i);
x2(i)=bluepeak_end(j)*srF+length(trac_ex_n(:,j))/11*co(i);
%p1=plot([x1(i) x2(i)],[abs(min(trac_ex_n(:,j)))/2-(abs(min(trac_ex_n(:,j)))/2)/3 abs(min(trac_ex_n(:,j)))/2-(abs(min(trac_ex_n(:,j)))/2)/3],'-','Color','b','LineWidth',4);
xr1(i)=redpeak_start*srF+length(trac_ex_n(:,j))/11*co(i);
xr2(i)=redpeak_end(j)*srF+length(trac_ex_n(:,j))/11*co(i);
p2=plot([xr1(i) xr2(i)],[abs(min(trac_ex_n(:,j)))/2-(abs(min(trac_ex_n(:,j)))/2)/3 abs(min(trac_ex_n(:,j)))/2-(abs(min(trac_ex_n(:,j)))/2)/3],'-','Color','r','LineWidth',4);
end
end
%scale barx
hold on;a1=x2(end)+500;a2=x2(end)+1000*srF;p1=plot([a2 a1],[ov_min ov_min],'-','Color','k','LineWidth',1);
%scale bary
scale_y=300;
hold on;y2=ov_min+scale_y;y1=ov_min;p1=plot([a1 a1],[y1 y2],'-','Color','k','LineWidth',1); 
hold on;
%Irradiance
subplot(4,2,pos(4));
for i=1:11
  plot(xr1(i),irr_exr(i,1)','sw','MarkerFaceColor','w','MarkerSize',5);box off
  hold on;
plot([xr1(i) xr1(i)],[0 irr_exr(i,1)'],'Color','r','LineWidth',1.5);axis off
end
hold on
plot(xr1(1:2),irr_exr(1:2,1)','sw','MarkerFaceColor','r','MarkerSize',5);box off
% plot(x1,irr_exr(:,1)','-sr','MarkerFaceColor','r','MarkerSize',5);box off
% ylabel('Irradiance (mW/mm^2)');
% xticks([x1(1):2000:x1(end)]);
% xticklabels({'1','2','3','4','5','6','7','8','9','10','11'});
% xlabel('Step');
% ylim([min(irr_exr(:,1)) 6]);
set(gca,'FontSize',10);  
%% Normalization to peak with blue FOR QUANTIFCATION
peakn_1=peak1b(:,[find(peak1b(1,:)~=0)])./min(peak1b(:,[find(peak1b(1,:)~=0)]));
peakn_10=peak10b(:,[find(peak10b(1,:)~=0)])./min(peak10b(:,[find(peak10b(1,:)~=0)]));
peakn_250=peak250b(:,[find(peak250b(1,:)~=0)])./min(peak250b(:,[find(peak250b(1,:)~=0)]));
peakn_1nb=peak1bn(:,[find(peak1bn(1,:)~=0)])./min(peak1b(:,[find(peak1b(1,:)~=0)]));
peakn_10nb=peak10bn(:,[find(peak10bn(1,:)~=0)])./min(peak10b(:,[find(peak10b(1,:)~=0)]));
peakn_250nb=peak250bn(:,[find(peak250bn(1,:)~=0)])./min(peak250b(:,[find(peak250b(1,:)~=0)]));
%% Plotting the ramps (250 , 10 , 1 ms red) with and without blue (50 ms)
%FIGURE for paper
fig9= figure;set(fig9, 'Name', 'Step protocol');set(fig9, 'Position', [200, 100, 600, 300]);set(gcf,'color','w');
subplot(1,2,2);
hold on;title('with 50 ms 472 nm','Color','b');
errorbar(nanmean(peakn_1,2),nanstd(peakn_1,1,2)/sqrt(size(peakn_1,2)),'-ob','MarkerSize',5,...
    'MarkerEdgeColor','k','MarkerFaceColor','r');hold on;
hold on;
errorbar(nanmean(peakn_10,2),nanstd(peakn_10,1,2)/sqrt(size(peakn_10,2)),'-^b','MarkerSize',5,...
    'MarkerEdgeColor','k','MarkerFaceColor','r');hold on;
errorbar(nanmean(peakn_250,2),nanstd(peakn_250,1,2)/sqrt(size(peakn_250,2)),'-sb','MarkerSize',6,...
    'MarkerEdgeColor','k','MarkerFaceColor','r');hold on;box off;xlabel('Step number');xticks([1:1:11]);
legend('1ms','10ms','250ms','Location','southeast');legend boxoff
set(gca,'FontSize',10);
subplot(1,2,1);
hold on;
title('637 nm only');
errorbar(nanmean(peakn_1nb,2),nanstd(peakn_1nb,1,2)/sqrt(size(peakn_1nb,2)),'-ok','MarkerSize',5,...
    'MarkerEdgeColor','k','MarkerFaceColor','r');hold on;
hold on;
errorbar(nanmean(peakn_10nb,2),nanstd(peakn_10nb,1,2)/sqrt(size(peakn_10nb,2)),'-^k','MarkerSize',5,...
    'MarkerEdgeColor','k','MarkerFaceColor','r');hold on;
errorbar(nanmean(peakn_250nb,2),nanstd(peakn_250nb,1,2)/sqrt(size(peakn_250nb,2)),'-sk','MarkerSize',6,...
    'MarkerEdgeColor','k','MarkerFaceColor','r');hold on;box off;xlabel('Step number');xticks([1:1:11]);ylabel('Peak normalized')
legend('1ms','10ms','250ms','Location','northwest');legend boxoff
set(gca,'FontSize',10);
%% Double rep with and without blue
for m=1:length(LGN)
try
temp=LGN(m).red_dr.ephys_traces_70(:,:,1);
red_nb(:,m)=temp(:);
peak_nb(:,m)=LGN(m).red_dr.neg_peak1(1,:);
irr_red_nb(:,m)=LGN(m).red_dr.neg_irr_red(1,:);  
irr_blue_nb(:,m)=LGN(m).blue_dr.neg_irr_blue(1,:); 
end
try
temp=LGN(m).red_dr.ephys_traces_70(:,:,2);
red_b(:,m)=temp(:);
peak_b(:,m)=LGN(m).red_dr.neg_peak1(2,:);
irr_red_b(:,m)=LGN(m).red_dr.neg_irr_red(2,:);  
irr_blue_b(:,m)=LGN(m).blue_dr.neg_irr_blue(2,:); 
end
end
%% Example cell for duble rep with and without blue
%FIGURE for paper
% Choose an example cell 
base_start          =   1;
base_end            =   99;
redpeak_start       =   100;
redpeak_end         =   350;%looking for the peak shoud be the same as for blue even though the window is longer?
bluepeak_start      =   351;
bluepeak_end        =   400;
scale_x= 100;
scale_y= 400;
xr1=[];
xr2=[];
x1=[];
x2=[];
cell_ex=4;
fig13= figure;set(fig13, 'Name', 'Step protocol');set(fig13, 'Position', [200, 100, 600, 400]);set(gcf,'color','w');
subplot(2,2,1);
title('637 nm only');hold on;
plot(red_nb(1:6001,cell_ex),'Color','k','LineWidth',1);box off;axis off;%ylim([min(cat_bluesteps_a) abs(min(cat_bluesteps_a))/2]);
hold on;
for i=1:3
x1(i)=redpeak_start*srF+length(red_nb(1:6000,cell_ex))/3*co(i);
x2(i)=redpeak_end*srF+length(red_nb(1:6000,cell_ex))/3*co(i);
p1=plot([x1(i) x2(i)],[abs(min(red_nb(1:6000,cell_ex)))/2-(abs(min(red_nb(1:6000,cell_ex)))/2)/3 abs(min(red_nb(1:6000,cell_ex)))/2-(abs(min(red_nb(1:6000,cell_ex)))/2)/3],'-','Color','r','LineWidth',4);
end
subplot(2,2,2);
hold on;
title('with 50 ms 472 nm');hold on;
plot(red_b(4000:end,cell_ex),'Color','k','LineWidth',1);box off;axis off;%ylim([min(cat_bluesteps_a) abs(min(cat_bluesteps_a))/2]);
co=[0:1:2];
hold on;
for i=1:3
x1(i)=redpeak_end*srF+length(red_nb(1:6000,cell_ex))/3*co(i);
x2(i)=bluepeak_end*srF+length(red_nb(1:6000,cell_ex))/3*co(i);
p1=plot([x1(i) x2(i)],[abs(min(red_b(4000:end,cell_ex)))/2-(abs(min(red_b(4000:end,cell_ex)))/2)/3 abs(min(red_b(4000:end,cell_ex)))/2-(abs(min(red_b(4000:end,cell_ex)))/2)/3],'-','Color','b','LineWidth',4);
xr1(i)=redpeak_start*srF+length(red_b(4000:end,cell_ex))/3*co(i);
xr2(i)=bluepeak_start*srF+length(red_b(4000:end,cell_ex))/3*co(i);
p2=plot([xr1(i) xr2(i)],[abs(min(red_b(4000:end,cell_ex)))/2-(abs(min(red_b(4000:end,cell_ex)))/2)/3 abs(min(red_b(4000:end,cell_ex)))/2-(abs(min(red_b(4000:end,cell_ex)))/2)/3],'-','Color','r','LineWidth',4);
end


ov_min=min(min(red_nb(1:6000,cell_ex)));
ov_max=max(max(red_nb(1:6000,cell_ex)));
%scale barx
hold on;a1=x2(end)+500;a2=x2(end)+600*srF;p1=plot([a2 a1],[ov_min ov_min],'-','Color','k','LineWidth',1);
%scale bary
scale_y=300;
hold on;y2=ov_min+scale_y;y1=ov_min;p1=plot([a1 a1],[y1 y2],'-','Color','k','LineWidth',1); 
hold on;

subplot(2,2,3);
irr_rb=irr_red_b(3:end,cell_ex)';
for i=1:3
  plot(xr1(i),irr_rb(i)','sw','MarkerFaceColor','w','MarkerSize',5);box off
  hold on;
plot([xr1(i) xr1(i)],[0 irr_rb(i)'],'Color','r','LineWidth',1.5);axis off
  hold on;
end
ylim([0 6]);
xticks([xr1(1):2000:xr1(end)]);
xlim([0 6000]);
xticklabels({'1','2','3'});
xlabel('Trial');
% plot(xr1,irr_red_nb(1:3,cell_ex)','-sr','MarkerFaceColor','r','MarkerSize',5);box off
% ylabel('Irradiance (mW/mm^2)');
% ylim([0 6]);
% xticks([xr1(1):2000:xr1(end)]);
% xlim([0 6000]);
% xticklabels({'1','2','3'});
% xlabel('Trial');

subplot(2,2,4);
%plot(xr1,irr_red_b(3:end,cell_ex)','-sr','MarkerFaceColor','r','MarkerSize',5);box off
% hold on;
% plot(x1,irr_blue_b(3:end,cell_ex)','-sb','MarkerFaceColor','b','MarkerSize',5);box off
% ylabel('Irradiance (mW/mm^2)');
% ylim([0 6]);
% xticks([xr1(1):2000:xr1(end)]);
% xlim([0 6000]);
% xticklabels({'1','2','3'});
% xlabel('Trial');
irr_rb=irr_red_b(3:end,cell_ex)';
irr_rnb=irr_blue_b(3:end,cell_ex)';
for i=1:3
  plot(xr1(i),irr_rb(i)','sw','MarkerFaceColor','w','MarkerSize',5);box off
  hold on;
plot([xr1(i) xr1(i)],[0 irr_rb(i)'],'Color','r','LineWidth',1.5);axis off
  hold on;
  plot(x1(i),irr_rnb(i)','sw','MarkerFaceColor','w','MarkerSize',5);box off
  hold on;
plot([x1(i) x1(i)],[0 irr_rnb(i)'],'Color','b','LineWidth',1.5);axis off
  hold on;
end
ylim([0 6]);
xticks([xr1(1):2000:xr1(end)]);
xlim([0 6000]);
xticklabels({'1','2','3'});
xlabel('Trial');

%% Normalization to peak with blue FOR QUANTIFCATION FOR DOUBLE REP
peak_b_norm=peak_b([3 4 5],[find(peak_b(1,:)~=0)])./peak_b(3,[find(peak_b(1,:)~=0)]);
peak_nb_norm=peak_nb([1 2 3],[find(peak_nb(1,:)~=0)])./peak_nb(1,[find(peak_nb(1,:)~=0)]);
%% Plotting
%FIGURE for paper
fig10= figure;set(fig10, 'Name', 'Step protocol');set(fig10, 'Position', [200, 100, 350, 200]);set(gcf,'color','w');
hold on;errorbar(nanmean(peak_nb_norm,2),nanstd(peak_nb_norm,1,2)/sqrt(size(peak_nb_norm,2)),'-ok','MarkerSize',5,...
    'MarkerEdgeColor','k','MarkerFaceColor','r');
hold on;errorbar(nanmean(peak_b_norm,2),nanstd(peak_b_norm,1,2)/sqrt(size(peak_b_norm,2)),'-ob','MarkerSize',5,...
    'MarkerEdgeColor','k','MarkerFaceColor','r');ylabel('Peak normalized');xlabel('Trial');xticks([1:1:3]);
legend('637 nm','637 + 472 nm','Location','southwest');legend boxoff                       
set(gca,'FontSize',10);
