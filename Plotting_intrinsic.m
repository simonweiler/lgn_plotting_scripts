%% Load LGN structure with intrinsic and reps

%put folder for str_LGN
str_LGN     = 'C:\Users\Simon-localadmin\Documents\MargrieLab\PhDprojects\LGN';
folder_list = uipickfiles('FilterSpec',str_LGN);
load(char(folder_list));
data=LGN;
%% Plot the intrinsic traces: 30 ms
fig1= figure;set(fig1, 'Position', [400, 300, 800, 800]);set(gcf,'color','w');
for i=1:length(data)
hold on;
subplot(5,6,i);
plot(LGN(i).intr_1.traces1);
box off;
end

%% Plot the repetition traces
for i=1:length(data)
figure;set(gcf,'color','w');
hold on;
for k=1:size(data(i).red_dr.ephys_traces_40,3);
subplot(3,2,k);
plot(data(i).red_dr.ephys_traces_40(1:1000,:,k));
box off;
end
end