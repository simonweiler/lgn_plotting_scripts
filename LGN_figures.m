%SW191127
%Script that needs the LGN structure with all extracted features
%Script displays one by one polished figures for the manuscript
%Additional functions needed
    %uipickfiles
 
%% Load LGN structure 
str_LGN     = 'C:\Users\Simon-localadmin\Documents\MargrieLab\PhDprojects\LGN';
folder_list = uipickfiles('FilterSpec',str_LGN);
load(char(folder_list));
data=data_mod;
%% Read out general features
for i=1:length(data)
    %Ocular category
    ocular_cat(i)=data(i).ocular_category;
    %Contra red or green
    contra_r(i)=data(i).brain_contra_ipsi;
    %Eye injection order
    eye(i)=data(i).eye_inj_ord;
    %Hemipshere
    hemi{i}=data(i).hemisphere;
    %MD
    md(i)=data(i).MD;
end
ocat1=find(ocular_cat==1);ocat2=find(ocular_cat==2);ocat3=find(ocular_cat==3);ocat4=find(ocular_cat==4);ocat5=find(ocular_cat==5);
c_red=find(contra_r==1);
i_red=find(contra_r==0);
md_n=find(md==0);
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
cat1=13; cat2=1; cat3=22; cat4=6 ; cat5=4;
%linewidth of traces
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
%Find the min/max for the traces for each category
sb_min=[min(data(r1(cat1)).step_red.ephys_traces_70(:,trace_nr,2)) min(data(r2(cat2)).step_red.ephys_traces_70(:,trace_nr,2))...
        min(data(r3(cat3)).step_red.ephys_traces_70(:,trace_nr,2)) min(data(r4(cat4)).step_red.ephys_traces_70(:,trace_nr,2))...
        min(data(r5(cat5)).step_red.ephys_traces_70(:,trace_nr,2))];
sb_max=[max(data(r1(cat1)).step_red.ephys_traces_40(:,trace_nr,2)) max(data(r2(cat2)).step_red.ephys_traces_40(:,trace_nr,2))...
        max(data(r3(cat3)).step_red.ephys_traces_40(:,trace_nr,2)) max(data(r4(cat4)).step_red.ephys_traces_40(:,trace_nr,2))...
        max(data(r5(cat5)).step_red.ephys_traces_40(:,trace_nr,2))];
ov_min=min(sb_min);
ov_max=max(sb_max);
figure;set(gcf,'color','w');plot(sb_min,'--or');hold on;plot(sb_max,'--ob');ylabel('Synaptic Input peak (pA)');xlabel('Category');box off
%The offset between the 70 and 40 mV
yoffset=ov_max/8;
ov_maxo=ov_max+ov_max/4;
%% Plot examples for the observed response types
%RED=Contra; GREEN=ipsi
fig1= figure;set(fig1, 'Name', 'Observed response types');set(fig1, 'Position', [400, 600, 1200, 200]);set(gcf,'color','w');
%Contra only
sb1=subplot(1,5,1);
if data(r1(cat1)).experimentator=='SW'
    srF=1;
else data(r1(cat1)).experimentator=='MF'  
    srF=2;
end
plot(data(r1(cat1)).step_red.ephys_traces_70(:,trace_nr,2),'Color','k','LineWidth',tr_l);box off;axis off
hold on;plot(data(r1(cat1)).step_red.ephys_traces_40(:,trace_nr,2)+yoffset,'Color','k','LineWidth',tr_l);box off;axis off
set(sb1,'ylim',[ov_min ov_maxo]);
hold on; x1= redpeak_start*srF;x2=redpeak_end*srF;p1=plot([x1 x2],[ov_maxo ov_maxo],'-','Color','r','LineWidth',i_l);
hold on; x1= bluepeak_start*srF;x2=bluepeak_end*srF; hold on;p1=plot([x1 x2],[ov_maxo ov_maxo],'-','Color','b','LineWidth',i_l);
%Ipsi only
sb2=subplot(1,5,2);
if data(r2(cat2)).experimentator=='SW'
    srF=1;
else data(r2(cat2)).experimentator=='MF'  
    srF=2;
end
plot(data(r2(cat2)).step_red.ephys_traces_70(:,trace_nr,2),'Color','k','LineWidth',tr_l);box off;axis off
hold on;plot(data(r2(cat2)).step_red.ephys_traces_40(:,trace_nr,2)+yoffset,'Color','k','LineWidth',tr_l);box off;axis off
set(sb2,'ylim',[ov_min ov_maxo]);
hold on; x1= redpeak_start*srF;x2=redpeak_end*srF;p1=plot([x1 x2],[ov_maxo ov_maxo],'-','Color','r','LineWidth',i_l);
hold on; x1= bluepeak_start*srF;x2=bluepeak_end*srF; hold on;p1=plot([x1 x2],[ov_maxo ov_maxo],'-','Color','b','LineWidth',i_l);


%Bino at 70 and 40 mV
sb3=subplot(1,5,3);  
if data(r3(cat3)).experimentator=='SW'
    srF=1;
else data(r3(cat3)).experimentator=='MF'  
    srF=2;
end
plot(data(r3(cat3)).step_red.ephys_traces_70(:,trace_nr,2),'Color','k','LineWidth',tr_l);box off;axis off
hold on;plot(data(r3(cat3)).step_red.ephys_traces_40(:,trace_nr,2)+yoffset,'Color','k','LineWidth',tr_l);box off;axis off
set(sb3,'ylim',[ov_min ov_maxo]);
hold on; x1= redpeak_start*srF;x2=redpeak_end*srF;p1=plot([x1 x2],[ov_maxo ov_maxo],'-','Color','r','LineWidth',i_l);
hold on; x1= bluepeak_start*srF;x2=bluepeak_end*srF; hold on;p1=plot([x1 x2],[ov_maxo ov_maxo],'-','Color','b','LineWidth',i_l);


%ipsi silent
sb4=subplot(1,5,4);  
if data(r4(cat4)).experimentator=='SW'
    srF=1;
else data(r4(cat4)).experimentator=='MF'  
    srF=2;
end
plot(data(r4(cat4)).step_red.ephys_traces_70(:,trace_nr,2),'Color','k','LineWidth',tr_l);box off;axis off
hold on;plot(data(r4(cat4)).step_red.ephys_traces_40(:,trace_nr,2)+yoffset,'Color','k','LineWidth',tr_l);box off;axis off
set(sb4,'ylim',[ov_min ov_maxo]);
hold on; x1= redpeak_start*srF;x2=redpeak_end*srF;p1=plot([x1 x2],[ov_maxo ov_maxo],'-','Color','r','LineWidth',i_l);
hold on; x1= bluepeak_start*srF;x2=bluepeak_end*srF; hold on;p1=plot([x1 x2],[ov_maxo ov_maxo],'-','Color','b','LineWidth',i_l);

%ipsi silent
sb5=subplot(1,5,5);  
if data(r5(cat5)).experimentator=='SW'
    srF=1;
else data(r5(cat5)).experimentator=='MF'  
    srF=2;
end
plot(data(r5(cat5)).step_red.ephys_traces_70(:,trace_nr,2),'Color','k','LineWidth',tr_l);box off;axis off
hold on;plot(data(r5(cat5)).step_red.ephys_traces_40(:,trace_nr,2)+yoffset,'Color','k','LineWidth',tr_l);box off;axis off
set(sb5,'ylim',[ov_min ov_maxo]);
hold on; x1= redpeak_start*srF;x2=redpeak_end*srF;p1=plot([x1 x2],[ov_maxo ov_maxo],'-','Color','r','LineWidth',i_l);
hold on; x1= bluepeak_start*srF;x2=bluepeak_end*srF; hold on;p1=plot([x1 x2],[ov_maxo ov_maxo],'-','Color','b','LineWidth',i_l);

%scale barx
hold on;x1= 900*srF;x2=1000*srF;p1=plot([x1 x2],[ov_min ov_min],'-','Color','k','LineWidth',2);
%scale bary
hold on;y2= ov_min+scale_y;y1=ov_min;p1=plot([x2 x2],[y1 y2],'-','Color','k','LineWidth',2); 
%% Example ramp for display purpose, atm binocular cells
%Blue AMPA/NMDA
%choose cell 
cellnr=7;
%Concatenate all 11 traces AMPA blue only
temp_cat=data(r_nm3(cellnr)).step_red.ephys_traces_70(:,:,1);
cat_bluesteps_a=temp_cat(:);
%Concatenate all 11 traces AMPA red with blue
temp_cat=data(r_nm3(cellnr)).step_red.ephys_traces_70(:,:,2);
cat_redsteps_a=temp_cat(:);
%Concatenate all 11 traces NMDA blue only
temp_cat=data(r_nm3(cellnr)).step_red.ephys_traces_40(:,:,1);
cat_bluesteps_n=temp_cat(:);
%Concatenate all 11 traces NMDA blue only
temp_cat=data(r_nm3(cellnr)).step_red.ephys_traces_40(:,:,2);
cat_redsteps_n=temp_cat(:);
%% Workflow for AMPA step
if data(r_nm3(cellnr)).experimentator=='SW'
    srF=1;
else data(r_nm3(cellnr)).experimentator=='MF'  
    srF=2;
end
fig7= figure;set(fig7, 'Name', 'Approach step display');set(fig7, 'Position', [200, 100, 600, 200]);set(gcf,'color','w');
%Plot blue only 
subplot(2,1,1);
plot(cat_bluesteps_a,'Color','k','LineWidth',1);box off;axis off;ylim([min(cat_bluesteps_a) abs(min(cat_bluesteps_a))/2]);
hold on;co=[0:1:10];
for i=1:11
x1=bluepeak_start*srF+length(cat_bluesteps_a)/11*co(i);
x2=bluepeak_end*srF+length(cat_bluesteps_a)/11*co(i);
p1=plot([x1 x2],[abs(min(cat_bluesteps_a))/2-(abs(min(cat_bluesteps_a))/2)/3 abs(min(cat_bluesteps_a))/2-(abs(min(cat_bluesteps_a))/2)/3],'-','Color','b','LineWidth',4);
end
text(0,abs(min(cat_bluesteps_a))/1.5,'472 nm','Color','b');

%Plot red/blue 
subplot(2,1,2);
plot(cat_redsteps_a,'Color','k','LineWidth',1);box off;axis off;ylim([min(cat_redsteps_a) abs(min(cat_redsteps_a))/2]);
hold on;co=[0:1:10];
for i=1:11
x1=bluepeak_start*srF+length(cat_bluesteps_a)/11*co(i);
x2=bluepeak_end*srF+length(cat_bluesteps_a)/11*co(i);
p1=plot([x1 x2],[abs(min(cat_redsteps_a))/2-(abs(min(cat_redsteps_a))/2)/3 abs(min(cat_redsteps_a))/2-(abs(min(cat_redsteps_a))/2)/3],'-','Color','b','LineWidth',4);
xr1=redpeak_start*srF+length(cat_redsteps_a)/11*co(i);
xr2=redpeak_end*srF+length(cat_redsteps_a)/11*co(i);
p2=plot([xr1 xr2],[abs(min(cat_redsteps_a))/2-(abs(min(cat_redsteps_a))/2)/3 abs(min(cat_redsteps_a))/2-(abs(min(cat_redsteps_a))/2)/3],'-','Color','r','LineWidth',4);
end
text(0,abs(min(cat_redsteps_a))/1.5,'637 nm','Color','r');text(1500*srF,abs(min(cat_redsteps_a))/1.5,'472 nm','Color','b');

%Scale bar
scale_x= 100;
scale_y= 400;
ov_min=min([min(cat_redsteps_a) min(cat_bluesteps_a)]);
%scale barx
hold on;x1= 10900*srF;x2=11000*srF;p1=plot([x1 x2],[ov_min ov_min],'-','Color','k','LineWidth',1.5);
%scale bary
hold on;y2= ov_min+scale_y;y1=ov_min;p1=plot([x2 x2],[y1 y2],'-','Color','k','LineWidth',1.5); 
%% ODI distributions and category plot
%Check the categories again: 2 bin cells have still ODI 1 and -1 %SW191118
%NO MD only atm
no_md=length([no_md1 no_md2 no_md3 no_md4 no_md5]);
fig2= figure;set(fig2, 'Name', 'ODI/category histogram');set(fig2, 'Position', [200, 600, 1000, 300]);set(gcf,'color','w');
%Step AMPA/NMDA
subplot(1,3,1);
b=bar([1 2 3 4 5],[length(no_md1)/no_md length(no_md2)/no_md length(no_md3)/no_md length(no_md4)/no_md length(no_md5)/no_md]*100);
box off;ylabel('Fraction of cells (%)');
b.FaceColor = 'flat';ylim([0 50]);
b.CData(2,:) = [1 0 0];b.CData(3,:) = [1 1 1]; b.CData(5,:) = [1 0 0];
yticks([0:10:50]); xticklabels({'contra','ipsi','bino','ipsi silent','contra silent'});
xtickangle(45);

subplot(1,3,2);
histogram([data([no_md1 no_md2 no_md3 no_md4 no_md5]).ODI_AMPA_step], 7, 'Normalization','probability','FaceColor',[0.7 0.7 0.7],'EdgeColor',[0 0 0]);
legend('AMPA');legend boxoff
ylim([0 0.8]);yticks([0:0.2:0.8]); 
yticklabels({'0','20','40','60','80'});ylabel('Fraction of cells (%)');xlabel('ODI');box off;

subplot(1,3,3);
histogram([data([no_md1 no_md2 no_md3 no_md4 no_md5]).ODI_NMDA_step], 7, 'Normalization','probability','FaceColor',[0.7 0.7 0.7],'EdgeColor',[0 0 0]);
legend('NMDA');legend boxoff
ylim([0 0.8]);yticks([0:0.2:0.8]); 
yticklabels({'0','20','40','60','80'});ylabel('Fraction of cells (%)');xlabel('ODI');box off;

% histogram([data(ocat3).ODI_AMPA_step],10,'FaceColor','k');hold on;histogram([data(ocat3).ODI_NMDA_step],10,'FaceColor','w');box off;xlabel('ODI');ylabel('Counts');
% ylim([0 50]);yticks([0:10:50]);legend('AMPA', 'NMDA');legend boxoff;

%% Sequential photostimulation prevents crosstalk
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
%% 
%Example trace for sequential photostimulation works in red only expressing
%cells
cellnr=4;
if LGN(cellnr).experimentator=='SW'
    srF=1;
else LGN(cellnr).experimentator=='MF'  
    srF=2;
end
temp_size=size(LGN(cellnr).step_red.ephys_traces_70(:,:,:),3);
blue_o=LGN(cellnr).step_red.ephys_traces_70(:,11,1);
blue_red=LGN(cellnr).step_red.ephys_traces_70(:,9,3);
%Plotting
fig8= figure;set(fig8, 'Name', 'Step protocol');set(fig8, 'Position', [200, 100, 200, 400]);set(gcf,'color','w');
subplot(2,1,1);
ov_maxo=max(blue_o);
plot(blue_o,'Color','k','LineWidth',1);box off;axis off;
%hold on; x1= redpeak_start*srF;x2=redpeak_end*srF;p1=plot([x1 x2],[ov_maxo ov_maxo],'-','Color','r','LineWidth',i_l);
hold on; x1= bluepeak_start*srF;x2=bluepeak_end*srF; hold on;p1=plot([x1 x2],[ov_maxo*5 ov_maxo*5],'-','Color','b','LineWidth',i_l);
subplot(2,1,2);
ov_maxo=max(blue_red);
plot(blue_red,'Color','k','LineWidth',1);box off;axis off;
hold on; x1= redpeak_start*srF;x2=redpeak_end*srF;p1=plot([x1 x2],[ov_maxo*30 ov_maxo*30],'-','Color','r','LineWidth',i_l);
hold on; x1= bluepeak_start*srF;x2=bluepeak_end*srF; hold on;p1=plot([x1 x2],[ov_maxo*30 ov_maxo*30],'-','Color','b','LineWidth',i_l);
ov_min=min(blue_red);
%scale barx
hold on;x1= 900*srF;x2=1000*srF;p1=plot([x1 x2],[ov_min*1.5 ov_min*1.5],'-','Color','k','LineWidth',2);
%scale bary
hold on;y2= ov_min+scale_y;y1=ov_min;p1=plot([x2 x2],[y1*1.5 y2*1.5],'-','Color','k','LineWidth',2); 
%
%% Compare red ramp with blue and no blue, red varies from 250, 10 and 1 ms, blue is always 50 ms
for m=1:length(LGN)
red250b=LGN(m).step_red.ephys_traces_70(:,:,3);
cat_bluesteps_a=temp_cat(:);
red10b=LGN(m).step_red.ephys_traces_70(:,:,7);
red1b=LGN(m).step_red.ephys_traces_70(:,:,9);
red250nb=LGN(m).step_red.ephys_traces_70(:,:,2);
red10nb=LGN(m).step_red.ephys_traces_70(:,:,6);
red1bn=LGN(m).step_red.ephys_traces_70(:,:,8);
end
%% Step protocol visualization prep
 rm={r_nm1 r_nm2 r_nm3 r_nm4 r_nm5};
 im={i_nm1 i_nm2 i_nm3 i_nm4 i_nm5};
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
%Contra red
fig6= figure;set(fig6, 'Name', 'Step protocol');set(fig6, 'Position', [200, 100, 400, 400]);set(gcf,'color','w');
subplot(2,1,1);
me=abs(nanmean(reshape([red_c_norm{:}],11,length([red_c_norm{:}])/11),2));
st=nanstd(reshape([red_c_norm{:}],11,length([red_c_norm{:}])/11),1,2)/sqrt(length([red_c_norm{:}])/11);
errorbar(me,st,'--or','MarkerSize',5,...
    'MarkerEdgeColor','k','MarkerFaceColor','r');hold on;
box off;ylabel('Peak normalized');xlabel('Step number');xticks([1:1:11]);
% Contra green
hold on;
me=abs(nanmean(reshape([green_c_norm{:}],11,length([green_c_norm{:}])/11),2));
st=nanstd(reshape([green_c_norm{:}],11,length([green_c_norm{:}])/11),1,2)/sqrt(length([green_c_norm{:}])/11);
errorbar(me,st,'--og','MarkerSize',5,...
    'MarkerEdgeColor','k','MarkerFaceColor','g');hold on;
box off;ylabel('Peak normalized');xlabel('Step number');xticks([1:1:11]);
legend('Contra Chrimson AMPA','Contra Chronos AMPA','Location','southeast');
legend boxoff;

subplot(2,1,2);
%Ipsi red
me=abs(nanmean(reshape([red_i_norm{:}],11,length([red_i_norm{:}])/11),2));
st=nanstd(reshape([red_i_norm{:}],11,length([red_i_norm{:}])/11),1,2)/sqrt(length([red_i_norm{:}])/11);
errorbar(me,st,'--or','MarkerSize',5,...
    'MarkerEdgeColor','k','MarkerFaceColor','r');hold on;
box off;ylabel('Peak normalized');xlabel('Step number');xticks([1:1:11]);
% Ipsi green
hold on;
me=abs(nanmean(reshape([green_i_norm{:}],11,length([green_i_norm{:}])/11),2));
st=nanstd(reshape([green_i_norm{:}],11,length([green_i_norm{:}])/11),1,2)/sqrt(length([green_i_norm{:}])/11);
errorbar(me,st,'--og','MarkerSize',5,...
    'MarkerEdgeColor','k','MarkerFaceColor','g');hold on;
box off;ylabel('Peak normalized');xlabel('Step number');xticks([1:1:11]);
legend('Ipsi Chrimson AMPA','Ipsi Chronos AMPA','Location','southeast');
legend boxoff;
%% %% Step protocol plots NMDA
tit={'contra only NMDA','ipsi only','bino','ipsi silent','contra silent'};
fig5= figure;set(fig5, 'Name', 'Step protocol');set(fig5, 'Position', [200, 100, 1400, 200]);set(gcf,'color','w');
 for j=1:5
 subplot(1,5,j);
 title(tit{j});
for i=1:length(rm{j})
    tem=rm{j};
    if j==1 | j==4
        try
    red_contra=data(tem(i)).step_red.pos_peak1(2,:);
    red_contra_norm=red_contra/max(abs(red_contra));
    r_c_n(:,i)=red_contra_norm;
    red_c_norm{:,j,i}=red_contra_norm;
    r_c_mean=nanmean(r_c_n,2);
    r_c_sem=nanstd(r_c_n,1,2)/sqrt(length(tem));
    all_rc(:,j)=r_c_mean;
    all_rc_sem(:,j)=r_c_sem;
        end
    %Plot individual curves that are contra and red
    plot(abs(red_contra_norm),'--r.');hold on;  box off;
    elseif j==2 | j==5
    green_ipsi=data(tem(i)).step_blue.pos_peak2(1,:); 
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
    red_contra=data(tem(i)).step_red.pos_peak1(2,:);
    red_contra_norm=red_contra/max(abs(red_contra));
    red_c_norm{:,j,i}=red_contra_norm;
    r_c_n(:,i)=red_contra_norm;
    r_c_mean=nanmean(r_c_n,2);
    r_c_sem=nanstd(r_c_n,1,2)/sqrt(length(tem));
    all_rc(:,j)=r_c_mean;
    all_rc_sem(:,j)=r_c_sem;
    plot(abs(red_contra_norm),'--r.');hold on;  box off;
    green_ipsi=data(tem(i)).step_blue.pos_peak2(1,:); 
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
    green_contra=data(tem_g(i)).step_blue.pos_peak2(1,:);
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
    red_ipsi=data(tem_g(i)).step_red.pos_peak1(2,:);
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
    green_contra=data(tem_g(i)).step_blue.pos_peak2(1,:);
    green_contra_norm=green_contra/max(abs(green_contra));
    g_c_n(:,i)=green_contra_norm;
     green_c_norm{:,j,i}=green_contra_norm;
    g_c_mean=nanmean(g_c_n,2);
    g_c_sem=nanstd(g_c_n,1,2)/sqrt(length(tem));
    all_gc(:,j)=g_c_mean;
    all_gc_sem(:,j)=g_c_sem;
    %Plot individual curves that are contra and green
    plot(abs(green_contra_norm),'--g.');hold on; box off;
    red_ipsi=data(tem_g(i)).step_red.pos_peak1(2,:);
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
%Contra red
fig6= figure;set(fig6, 'Name', 'Step protocol');set(fig6, 'Position', [200, 100, 400, 400]);set(gcf,'color','w');
subplot(2,1,1);
me=abs(nanmean(reshape([red_c_norm{:}],11,length([red_c_norm{:}])/11),2));
st=nanstd(reshape([red_c_norm{:}],11,length([red_c_norm{:}])/11),1,2)/sqrt(length([red_c_norm{:}])/11);
errorbar(me,st,'--or','MarkerSize',5,...
    'MarkerEdgeColor','k','MarkerFaceColor','r');hold on;
box off;ylabel('Peak normalized');xlabel('Step number');xticks([1:1:11]);
% Contra green
hold on;
me=abs(nanmean(reshape([green_c_norm{:}],11,length([green_c_norm{:}])/11),2));
st=nanstd(reshape([green_c_norm{:}],11,length([green_c_norm{:}])/11),1,2)/sqrt(length([green_c_norm{:}])/11);
errorbar(me,st,'--og','MarkerSize',5,...
    'MarkerEdgeColor','k','MarkerFaceColor','g');hold on;
box off;ylabel('Peak normalized');xlabel('Step number');xticks([1:1:11]);
legend('Contra Chrimson NMDA','Contra Chronos NMDA','Location','southeast');
legend boxoff;

subplot(2,1,2);
%Ipsi red
me=abs(nanmean(reshape([red_i_norm{:}],11,length([red_i_norm{:}])/11),2));
st=nanstd(reshape([red_i_norm{:}],11,length([red_i_norm{:}])/11),1,2)/sqrt(length([red_i_norm{:}])/11);
errorbar(me,st,'--or','MarkerSize',5,...
    'MarkerEdgeColor','k','MarkerFaceColor','r');hold on;
box off;ylabel('Peak normalized');xlabel('Step number');xticks([1:1:11]);
% Ipsi green
hold on;
me=abs(nanmean(reshape([green_i_norm{:}],11,length([green_i_norm{:}])/11),2));
st=nanstd(reshape([green_i_norm{:}],11,length([green_i_norm{:}])/11),1,2)/sqrt(length([green_i_norm{:}])/11);
errorbar(me,st,'--og','MarkerSize',5,...
    'MarkerEdgeColor','k','MarkerFaceColor','g');hold on;
box off;ylabel('Peak normalized');xlabel('Step number');xticks([1:1:11]);
legend('Ipsi Chrimson NMDA','Ipsi Chronos NMDA','Location','southeast');
legend boxoff;

%% Dominat vs nondominat eye bino cells
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
%% Boxplot and plotspread plot for DE vs NDE
figure; set(gcf,'color','w');
plotSpread(abs([DE' NDE'])/1000,'categoryIdx',[ones(1,length(DE))' ones(1,length(NDE))'*2], 'categoryMarkers',{'.','.'},'categoryColors',{'k','k'});
hold on;
con=[DE NDE];
grp=[ones(1,length(DE)) ones(1,length(NDE))*2];
boxplot(abs(con)/1000,grp,'Colors','k','OutlierSize',0.001);
box off;xticklabels({'DE','NDE'});ylabel('Peak amplitude (nA)');
yticks([0:1:4]);

%% 

                       
                        
                       
                        
                       
                        

