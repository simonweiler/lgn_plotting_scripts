function  [blue_ramp, red_ramp]=double_rep(list, idx, voltage, pathName, fc, show, ramp_rtrace, user, filterephys,adata_dir,injection_order);
%SW181229
%function to extract synaptic current peak, integral and photodiode signal
%for blue and red laser

%fct inputs
%   list= list of xsg files per cell
%   idx=  indices for ramp recodings
%   voltage= indices for voltage clamp (1 means -70 mV and 0 means 40 mV)
%   pathName=folder name of cell
%   fc= factor for threshold std
%   show= display plot or not (1 or 0)
%   ramp_rtrace= extract and save raw traces or not (1 or 0)

%define temporal windows
base_start          =   1;
base_end            =   99;

redpeak_start       =   100;
redpeak_end         =   350;%looking for the peak shoud be the same as for blue even though the window is longer?
bluepeak_start      =   351;
bluepeak_end        =   400;

base_start2          =   bluepeak_start-100;
base_end2            =   bluepeak_start-2;

BL_AMPA_P1 = [redpeak_start-28 redpeak_start-1];
BL_AMPA_P2 = [bluepeak_start-28 bluepeak_start-1];
BL_NMDA_P1 = [redpeak_start-99 redpeak_start-1];
BL_NMDA_P2 = [bluepeak_start-99 bluepeak_start-1];

RW_AMPA_P1 = [redpeak_start+3 redpeak_start+30];
RW_AMPA_P2 = [bluepeak_start+3 bluepeak_start+30];
RW_NMDA_P1 = [redpeak_start+3 redpeak_start+107];
RW_NMDA_P2 = [bluepeak_start+3 bluepeak_start+107];

% show = 1;

%% TR2019: filtering
% filterephys = 1;        % filtering yes/no?
cutoff      = 1000;     % Hz (use 500 Hz for mini event / amplitude detection and 1000Hz for max currents. Chen & Regehr 2000)
order       = 4;        % filter order ('pole'). (use 4 pole for minis and max current. Chen & Regehr 2000)
type        = 'Butter'; % filter type ('Bessel' or 'Butter' (for Butterworth -> ). Default: Bessel. Use Bessel at > 4 order to prevent ripples)

if filterephys
    disp('- - - - - - - -')
    disp(['Filtering: ' num2str(order) ' pole ' type '-Filter w/ ' num2str(cutoff) ' Hz cutoff']);
    disp('- - - - - - - -')
end

%% TR2019: plot specs
plotlength = 1; %seconds
savefig = 0; %save main figure

%create vector with start and end point for each ramp within the cell recording
%plot if wanted
if show==1
    try; clf(fig1); end
    try; clf(fitfig); end
    fig1 = figure;
    set(fig1, 'Name', char(pathName));
    set(fig1, 'Position', [200, 0, 1500, 1000]);
end

if ~user %SW
    runramp=1:11:length(idx);
    idm=voltage(idx);
    vclamp=idm(runramp);
    runramp=[runramp runramp(end)+11];
    numb_ramps = length(idx)/11; % how many ramps in total
else
    vclamp=voltage(idx);
    numb_ramps = length(idx); % how many ramps in total
end

%load each ramp per cell consecutively and extract relevant values such as
%snaptic current peak, integral and photodiode signal for blue and red
%temporal windows
for j=1:numb_ramps %loop across ramps per cell
    counter=1;

    
    if ~user %SW
        ramp_sequence = runramp(j):runramp(j+1)-1;
    else
        load([char(pathName) filesep list(idx(j)).name],'-mat');
        sr = header.ephys.ephys.sampleRate;%check sample rate
        srF = 1/(1000/sr);
        
        samples_per_sweep = header.ephys.ephys.traceLength*sr;
        timebase=1/sr:1/sr:samples_per_sweep/sr; %TR2019: timebase
        
        traces=data.ephys.trace_1;%raw ephys trace
        
        if filterephys % TR2019: filtering
            traces = lowpassfilt(traces, order, cutoff, sr, type);
        end
        
        photodiode=data.acquirer.trace_1;%photodiode (PD) signal
        ind_traces=reshape(traces,[length(traces)/(length(traces)/samples_per_sweep) length(traces)/samples_per_sweep]);
        photodiode=reshape(photodiode,[length(traces)/(length(traces)/samples_per_sweep) length(traces)/samples_per_sweep]);
        
        ramp_sequence = 1:size(ind_traces,2);
        
    end
    
    
    
    
    for i=ramp_sequence %within each ramp load xsg files (11 in total per ramp)
        
        if ~user % SW
            load([char(pathName) filesep list(idx(i)).name],'-mat');
            sr = header.ephys.ephys.sampleRate;%check sample rate
            srF = 1/(1000/sr);
            
            samples_per_sweep = header.ephys.ephys.traceLength*sr;
            timebase=1/sr:1/sr:samples_per_sweep/sr; %TR2019: timebase
            traces=data.ephys.trace_1;%raw ephys trace
            
            if filterephys % TR2019: filtering
                traces = lowpassfilt(traces, order, cutoff, sr, type);
            end
            
            photodiode=data.acquirer.trace_1;%photodiode (PD) signal
            
            try
                blue_amp(j,counter)=header.pulseJacker.pulseJacker.pulseDataMap{4,counter+1}.amplitude;%blue laser amplitude set in ephus
            catch
                blue_amp(j,counter)=0;
            end
            try
                red_amp(j,counter)=header.pulseJacker.pulseJacker.pulseDataMap{2,counter+1}.amplitude;%red laser amplitude set in ephus
            catch
                red_amp(j,counter)=0;
            end
            
            bs=traces(base_start*srF:base_end*srF,:);%first 100 ms baseline trace
            bs_std=std(bs);%std of baseline trace
            bs_traces=traces-mean(bs);%subtract baseline
            bs_photodiode=photodiode-mean(photodiode(base_start*srF:base_end*srF,:));
           
            
        else
            traces_clip=ind_traces(:,i);
            photodiode_clip=photodiode(:,i);
            try
            blue_amp(j,counter)=header.pulseJacker.pulseJacker.pulseDataMap{2,counter+1}.amplitude;%blue laser amplitude set in ephus
            catch
                blue_amp(j,counter)=0;
            end
            try
                red_amp(j,counter)=header.pulseJacker.pulseJacker.pulseDataMap{3,counter+1}.amplitude;%red laser amplitude set in ephus
            catch
                red_amp(j,counter)=0;
            end
            
            bs=traces_clip(base_start*srF:base_end*srF,:);%first 100 ms baseline trace
            bs_std=std(bs);%std of baseline trace
            bs_traces=traces_clip-mean(traces_clip(base_start*srF:base_end*srF,:));%subtract baseline
            bs_photodiode=photodiode_clip-mean(photodiode_clip(base_start*srF:base_end*srF,:)); 
        end        
        
        dp=mean(bs_photodiode((redpeak_start+50)*srF:(redpeak_end-50)*srF))>0.025;
        idx70=find(vclamp(j)==1 & dp==0);
        idx40=find(vclamp(j)==0 & dp==0);
        idx_sp70=find(vclamp(j)==1 & dp==1);
        idx_sp40=find(vclamp(j)==0 & dp==1);
        
        
        %photodiode
        PD1(j,counter)=mean(bs_photodiode(redpeak_start*srF:redpeak_end*srF,:));%max values of PD signal within the red stimulation window
        
        %%%%extract irradiance for red%%%%
        if ~user
            yirr_red(j,counter)=(12.19*PD1(j,counter)-0.4319)/100;
        else
            yirr_red(j,counter)=(104.1 *PD1(j,counter)-3.467)/100;
        end
        
        %for second window (same extraction as above for blue laser window
        integ2(j,counter)=trapz(bs_traces(bluepeak_start*srF:bluepeak_end*srF,:));
        
        %%%%extract irradiance for blue%%%%
        PD2(j,counter)=mean(bs_photodiode(bluepeak_start*srF:bluepeak_end*srF,:));
        if ~user
            yirr_blue(j,counter)=(7.232*PD2(j,counter)-0.9951)/100;%given in mW/mm2 compare to Klapoetke 2014
        else
            yirr_blue(j,counter)=(679.2*PD2(j,counter)-26.82)/100;
        end
        
        %% Test pulse for Rs,Rm and Cm calculation SW
        %peak1_test(j,counter)=min(bs_traces(testpulse_start*srF:testpulse_end*srF));
        
        if header.ephys.ephys.stimOnArray==0
           % header.acquirer.acquirer.triggerTime
            testpulse_start =header.ephys.ephys.pulseParameters{1, 1}.squarePulseTrainDelay*sr;
            
            testpulse_end       =   testpulse_start+srF*100;
            %             if header.ephys.ephys.pulseParameters{1, 1}.squarePulseTrainDelay==3
            %             testpulse_start     =   2980;
            %             testpulse_end       =   3080;
            %             else header.ephys.ephys.pulseParameters{1, 1}.squarePulseTrainDelay==2
            %             testpulse_start     =   1980;
            %             testpulse_end       =   2080;
            %             end
            test_window_ind=[(testpulse_start-1) : testpulse_end]';
            baseline_window_ind=[(base_start*srF+1) : base_end*srF]';
            
            numPointsTest=size(test_window_ind,1);
            numTraces=1;
            ydata=bs_traces;
            cellparam_ydata=zeros(numPointsTest, numTraces);
            cellparam_ydata=ydata(test_window_ind,:);
            magdY=abs(diff(cellparam_ydata(:,1)));
            divs3=find(magdY>0.5.*max(magdY));
            
            baseline_mean=median(ydata(baseline_window_ind,:));
            
            
            cellparam_ydata=cellparam_ydata(divs3(1)+round(srF./2):divs3(end)-round(srF./2),:);
            Slope1=cellparam_ydata(1,:)-cellparam_ydata(2,:);
            Slope2=cellparam_ydata(2,:)-cellparam_ydata(3,:);
            dSlope=Slope1-Slope2;
            SlopeEst1=Slope1+dSlope;
            SlopeEst2=Slope1+2.*dSlope;
            
            Rs=0;
            Rm=0;
            Cm=0;
            amp=-5;
            pretestbase=baseline_mean;
            
            temp_range = round(length(cellparam_ydata)*0.3):length(cellparam_ydata);
            midtestbase=mean(cellparam_ydata(temp_range)');
            
            peak=min(cellparam_ydata(:)');
            %improve estimate - first two points (sr=4 kHz) are missing!
            peak=peak+SlopeEst1.*1.3;
            %find tau
            r=find(cellparam_ydata(:) < (peak-(midtestbase)).*(exp(-1))+midtestbase);
            tau = length(r)./srF;   %time constant, in ms
            Rs = abs(amp./(peak-pretestbase)); %in Mohms
            %Rm = abs(amp./(midtestbase-pretestbase)); %YP's code
            Rm = abs((amp-(midtestbase-pretestbase).*Rs)./(midtestbase-pretestbase)); %also in Gohms
            Cm=(Rs+Rm).*tau./(Rs.*Rm); %in pF
            
            % read Rs header with time stamp for later plotting
            
            %%%%across holding potential and stepsper cell
            Rs_cell(j,counter)=Rs*1000;
            Rm_cell(j,counter)=Rm*1000;
            Cm_cell(j,counter)=Cm;
            try
            time_stamp(j,counter)=data.ephys.dataEventTimestamps_1(counter);
            catch
                time_stamp(j,counter)=data.ephys.dataEventTimestamps_1;
            end
        else
            Rs_cell(j,counter)=NaN;
            Rm_cell(j,counter)=NaN;
            Cm_cell(j,counter)=NaN;
            time_stamp(j,counter)=NaN;
        end
              
        %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%IMPLEMENTED AFTER MEETING FROM 190109%% neg_peak2(j,counter)=min(bs_traces(bluepeak_start*srF:bluepeak_end*srF,:));
        %neg peak2 is calculated using the current difference  between the last 10ms of
        %the first time window and the peak in the subsequent 2nd window to
        %correct for decay issues from the first pulse
        
        % first window
        integ1(j,counter)=trapz(bs_traces(redpeak_start*srF:redpeak_end*srF,:));%Integral within the red stimulation window
        
        base_neg_mean1(j,counter)=mean(bs_traces(BL_AMPA_P1(1)*srF:BL_AMPA_P1(2)*srF,:));
        neg_peak1(j,counter)=min(bs_traces(RW_AMPA_P1(1)*srF:RW_AMPA_P1(2)*srF,:));%negative peak within the red stimulation window
        neg_mean1(j,counter)=nanmean(bs_traces(RW_AMPA_P1(1)*srF:RW_AMPA_P1(2)*srF,:));%mean across redtrace
        neg_fail1(j,counter)=neg_peak1(j,counter)<fc*bs_std*(-1);%vector with binary values when neg peaks crossed definded std threshold
        
        % first window NMDA
        base_pos_mean1(j,counter)=mean(bs_traces(BL_NMDA_P1(1)*srF:BL_NMDA_P1(2)*srF,:));
        pos_peak1(j,counter)=max(bs_traces(RW_NMDA_P1(1)*srF:RW_NMDA_P1(2)*srF,:));%positive peak within the red stimulation window
        pos_mean1(j,counter)=nanmean(bs_traces(RW_NMDA_P1(1)*srF:RW_NMDA_P1(2)*srF,:));
        pos_fail1(j,counter)=pos_peak1(j,counter)>fc*bs_std;%vector with binary values when pos peaks crossed definded std threshold

        % second window AMPA
        base_neg_mean2(j,counter)=mean(bs_traces(BL_AMPA_P2(1)*srF:BL_AMPA_P2(2)*srF,:));
        neg_peak2(j,counter)=min(bs_traces(RW_AMPA_P2(1)*srF:RW_AMPA_P2(2)*srF,:));
        neg_mean2(j,counter)=nanmean(bs_traces(RW_AMPA_P2(1)*srF:RW_AMPA_P2(2)*srF,:));
        bs_std2 = std(bs_traces((bluepeak_start+7)*srF:(bluepeak_start+17)*srF,:))*-1;
        neg_fail2(j,counter)=neg_peak2(j,counter)<fc*bs_std2;
        
        
        %%%IMPLEMENTED AFTER MEETING FROM 190109%For NMDA: approach is to fit an expontial and then subtract this from
        %the actual curve to detect a second peak
        
        if ~isempty(idx70) %if j<=2
            base_pos_mean2(j,counter)=mean(bs_traces(BL_NMDA_P2(1)*srF:BL_NMDA_P2(2)*srF,:));
            pos_peak2(j,counter)=0;
            pos_mean2(j,counter)=NaN;
            pos_fail2(j,counter)=0;
                        
            yf=bs_traces;
            diff_bs_traces=bs_traces;
            
        elseif  ~isempty(idx_sp70)
            base_pos_mean2(j,counter)=mean(bs_traces(BL_NMDA_P2(1)*srF:BL_NMDA_P2(2)*srF,:));
            pos_peak2(j,counter)=0;
            pos_mean2(j,counter)=NaN;
            pos_fail2(j,counter)=0;
                        
            yf=bs_traces;
            diff_bs_traces=bs_traces;
            
        elseif ~isempty(idx40) % elseif j==3
            %INPORTANT SHIT HAPPENS HERE
            base_pos_mean2(j,counter)=mean(bs_traces(BL_NMDA_P2(1)*srF:BL_NMDA_P2(2)*srF,:));
            pos_peak2(j,counter)=max(bs_traces(RW_NMDA_P2(1)*srF:RW_NMDA_P2(2)*srF,:));
            pos_mean2(j,counter)=nanmean(bs_traces(RW_NMDA_P2(1)*srF:RW_NMDA_P2(2)*srF,:));
            pos_fail2(j,counter)=pos_peak2(j,counter)>fc*bs_std2;
                        
            yf=bs_traces;
            diff_bs_traces=bs_traces;
            
        elseif ~isempty(idx_sp40) %j==4;
            %currmaxpos(j,counter)=max(bs_traces(bluepeak_start*srF:(bluepeak_end+50)*srF,:));
            if ~pos_fail1(j,counter)% currmaxpos(j,counter)>pos_peak1(j,counter)
                base_pos_mean2(j,counter)=mean(bs_traces(BL_NMDA_P2(1)*srF:BL_NMDA_P2(2)*srF,:));
                pos_peak2(j,counter)=max(bs_traces(RW_NMDA_P2(1)*srF:RW_NMDA_P2(2)*srF,:));
                pos_mean2(j,counter)=nanmean(bs_traces(RW_NMDA_P2(1)*srF:RW_NMDA_P2(2)*srF,:));
                pos_fail2(j,counter)=pos_peak2(j,counter)>fc*bs_std2;
                
                yf=bs_traces;
                diff_bs_traces=bs_traces;
            else
                if ~user % SW
                    xt=1:50000;
                else
                    xt=1:200000;
                end
                
                A=pos_peak1(j,counter);
                t1=find(bs_traces==A)+(10*srF); % start 10 after peak
                t=t1:(redpeak_end-1)*srF;
                t=t';  % check dif SW MF
                curr_t=bs_traces(t);
                
                try
                    % fit with constraints
                    % single exponential
                    fo = fitoptions('Method','NonlinearLeastSquares', ...
                        'Lower',[0 -1], ...
                        'Upper', [pos_peak1(j,counter)*10 0],...
                        'StartPoint',[pos_peak1(j,counter)*3 -0.0005]);
                    ft = fittype('a*exp(b*(x))', 'options', fo);
                    [f1 gof1]=fit(t,curr_t,ft);
                    yf1=f1.a*exp(f1.b*xt);
                    %                     figure; plot(bs_traces); hold on; plot(yf)
                    
                    try % double exponential decay (if it fails or is worse exp1 is used)
                        fo = fitoptions('Method','NonlinearLeastSquares', ...
                            'MaxIter', 1000,...
                            'Lower',[0 -1 0 -1], ...
                            'Upper', [pos_peak1(j,counter)*10 0 pos_peak1(j,counter)*10 0],...
                            'StartPoint',[pos_peak1(j,counter)*1 -0.001 pos_peak1(j,counter)*1 -0.000005]);
                        ft = fittype('a*exp(b*(x)) + c*exp(d*(x))', 'options', fo);
                        [f2 gof2]=fit(t,curr_t,ft);
                        yf2=f2.a*exp(f2.b*xt) + f2.c*exp(f2.d*xt);
                        %                     figure; plot(bs_traces); hold on; plot(yf)
                    end
                    try
                        if gof2.adjrsquare>gof1.adjrsquare
                            gof = gof2;
                            f = f2;
                            yf = yf2;
                        else
                            gof = gof1;
                            f = f1;
                            yf = yf1;
                        end
                    catch
                        gof = gof1;
                        f = f1;
                        yf = yf1;
                    end
                    
                    
                    if ~user
                        for m=1:10000
                            diff_bs_traces(m,:)=bs_traces(m)-yf(m);
                        end
                    else
                        for m=1:40000
                            diff_bs_traces(m,:)=bs_traces(m)-yf(m);
                        end
                    end
                    
                    bs_diff_std2=std(diff_bs_traces(BL_NMDA_P2(1)*srF:BL_NMDA_P2(2)*srF,:));
                    
                    if show==1
                        fitfig = figure;
                        plot(bs_traces);hold on;plot(yf);plot(diff_bs_traces);
                        
                        set(fitfig, 'Name', ['FIT:' char(pathName) ]);
                        %%red vertical lines
                        hold on;
                        y1=get(gca,'ylim');
                        x1= redpeak_start*srF;
                        hold on;
                        p1=plot([x1 x1],y1,'--','Color','r');
                        p1.Color(4) = 0.3;
                        hold on;
                        y1=get(gca,'ylim');
                        x1=redpeak_end*srF;
                        hold on;
                        p2=plot([x1 x1],y1,'--','Color','r');
                        p2.Color(4) = 0.3;
                        hold on;
                        %%blue vertical lines
                        y1=get(gca,'ylim');
                        x1=bluepeak_start*srF;
                        hold on;
                        p3=plot([x1 x1],y1,'--','Color','b');
                        p3.Color(4) = 0.3;
                        hold on;
                        y1=get(gca,'ylim');
                        x1=bluepeak_end  *srF;
                        hold on;
                        p4=plot([x1 x1],y1,'--','Color','b');
                        p4.Color(4) = 0.3;
                        
                    end
                    
                    if gof.adjrsquare>0.9
                        base_pos_mean2(j,counter) = nanmean(diff_bs_traces(BL_NMDA_P2(1)*srF:BL_NMDA_P2(2)*srF,:));
                        pos_peak2(j,counter)=max(diff_bs_traces(RW_NMDA_P2(1)*srF:RW_NMDA_P2(2)*srF,:));
                        pos_mean2(j,counter)=nanmean(diff_bs_traces(RW_NMDA_P2(1)*srF:RW_NMDA_P2(2)*srF,:));
                        pos_fail2(j,counter)=pos_peak2(j,counter)>fc*1.5*bs_diff_std2;
                        
                    else
                        base_pos_mean2(j,counter) = nanmean(bs_traces(BL_NMDA_P2(1)*srF:BL_NMDA_P2(2)*srF,:));
                        pos_peak2(j,counter)=max(bs_traces(RW_NMDA_P2(1)*srF:RW_NMDA_P2(2)*srF,:));
                        pos_mean2(j,counter)=nanmean(bs_traces(RW_NMDA_P2(1)*srF:RW_NMDA_P2(2)*srF,:));
                        pos_fail2(j,counter)=pos_peak2(j,counter)>fc*1.5*bs_diff_std2;
                        
                    end
                catch
                    base_pos_mean2(j,counter) = nanmean(bs_traces(BL_NMDA_P2(1)*srF:BL_NMDA_P2(2)*srF,:));
                    pos_peak2(j,counter)=max(bs_traces(RW_NMDA_P2(1)*srF:RW_NMDA_P2(2)*srF,:));
                    pos_mean2(j,counter)=nanmean(bs_traces(RW_NMDA_P2(1)*srF:RW_NMDA_P2(2)*srF,:));
                    pos_fail2(j,counter)=pos_peak2(j,counter)>fc*bs_std2;
                    
                    yf=bs_traces;
                    diff_bs_traces=bs_traces;
                    
                end
            end
        end
        
        %ephys_traces
        ephys_traces(:,counter,j)=bs_traces;
%        fit_traces(:,counter,j)=yf;
        %diff_traces(:,counter,j)=diff_bs_traces;
        bs_photodiode_signal(:,counter,j)=bs_photodiode;
        yf=[];
        
        
        if exist('gof','var')==1
            gof_fit(counter,j)=gof.adjrsquare; clear gof.adjrsquare
            fit_param{counter,j} = f;
        else
            gof_fit(counter,j) = NaN; 
            fit_param{counter,j} = NaN;
        end
        
        counter=counter+1;
        traces=[];
        
        
        if show
            figure(fig1)
            subplot(2,numb_ramps-2,j);
            
            plot(bs_traces(1:plotlength*sr,:),'linewidth',1,'Color',[0 0 0]+0.05*counter);
            hold on;
            ylabel('Synaptic input (pA)');
            xlabel('Samples');

            %%red vertical lines
            hold on;
            y1=get(gca,'ylim');
            x1= redpeak_start*srF;
            hold on;
            p1=plot([x1 x1],y1,'--','Color','r');
            p1.Color(4) = 0.3;
            hold on;
            y1=get(gca,'ylim');
            x1=redpeak_end*srF;
            hold on;
            p2=plot([x1 x1],y1,'--','Color','r');
            p2.Color(4) = 0.3;
            hold on;
            %%blue vertical lines
            y1=get(gca,'ylim');
            x1=bluepeak_start*srF;
            hold on;
            p3=plot([x1 x1],y1,'--','Color','b');
            p3.Color(4) = 0.3;
            hold on;
            y1=get(gca,'ylim');
            x1=bluepeak_end  *srF;
            hold on;
            p4=plot([x1 x1],y1,'--','Color','b');
            p4.Color(4) = 0.3;
            drawnow
        end
    end
end

if savefig
    cd(adata_dir);
    saveas(fig1, [char(pathName) '.png'])
end

%%%%%%%%%%%%%%%%%%%%%  output %%%%%%%%%%%%%%%%%
red_ramp.response_windows = BL_AMPA_P1;
red_ramp.response_windows = BL_AMPA_P2;
red_ramp.response_windows = BL_NMDA_P1;
red_ramp.response_windows = BL_NMDA_P2;
red_ramp.response_windows = RW_AMPA_P1;
red_ramp.response_windows = RW_AMPA_P2;
red_ramp.response_windows = RW_NMDA_P1;
red_ramp.response_windows = RW_NMDA_P2;


%create structure with extracted parameters

%check whether VC at -70 or 40
idx70=find(vclamp==1);
idx40=find(vclamp==0);

red_ramp.neg_peak1=neg_peak1(idx70,:);
red_ramp.neg_mean1=neg_mean1(idx70,:);
red_ramp.neg_fail1=neg_fail1(idx70,:);
blue_ramp.neg_peak2=neg_peak2(idx70,:);
blue_ramp.neg_mean2=neg_mean2(idx70,:);
blue_ramp.neg_fail2=neg_fail2(idx70,:);
red_ramp.neg_integ1=integ1(idx70,:);
red_ramp.neg_PD=PD1(idx70,:);
blue_ramp.neg_PD=PD2(idx70,:);
red_ramp.neg_irr_red=yirr_red(idx70,:);
blue_ramp.neg_irr_blue=yirr_blue(idx70,:);
red_ramp.neg_laser_amp=red_amp(idx70,:);
blue_ramp.neg_laser_amp=blue_amp(idx70,:);
blue_ramp.neg_Rs=Rs_cell(idx70,:);
red_ramp.neg_Rs=Rs_cell(idx70,:);
red_ramp.neg_timestamp=time_stamp(idx70,:);
blue_ramp.neg_Cm=Cm_cell(idx70,:);
red_ramp.neg_Cm=Cm_cell(idx70,:);
red_ramp.neg_base_mean1=base_neg_mean1(idx70,:);
blue_ramp.neg_base_mean2=base_neg_mean2(idx70,:);


red_ramp.pos_peak1=pos_peak1(idx40,:);
red_ramp.pos_mean1=pos_mean1(idx40,:);
red_ramp.pos_fail1=pos_fail1(idx40,:);
blue_ramp.pos_peak2=pos_peak2(idx40,:);
blue_ramp.pos_mean2=pos_mean2(idx40,:);
blue_ramp.pos_fail2=pos_fail2(idx40,:);
blue_ramp.pos_integ2=integ2(idx40,:);
blue_ramp.pos_PD=PD2(idx40,:);
red_ramp.pos_PD=PD1(idx40,:);
red_ramp.pos_irr_red=yirr_red(idx40,:);
blue_ramp.pos_irr_blue=yirr_blue(idx40,:);
red_ramp.pos_laser_amp=red_amp(idx40,:);
blue_ramp.pos_laser_amp=blue_amp(idx40,:);
blue_ramp.pos_Rs=Rs_cell(idx40,:);
red_ramp.pos_Rs=Rs_cell(idx40,:);
blue_ramp.pos_Cm=Cm_cell(idx40,:);
red_ramp.pos_Cm=Cm_cell(idx40,:);
red_ramp.pos_timestamp=time_stamp(idx40,:);

red_ramp.pos_base_mean1=base_pos_mean1(idx40,:);
blue_ramp.pos_base_mean2=base_pos_mean2(idx40,:);

%save ephys traces if desired
sub_s=10;
if ramp_rtrace==1
    red_ramp.ephys_traces_70=ephys_traces(1:sub_s:plotlength*sr,:,idx70);
    %red_ramp.fit_traces_70=fit_traces(1:sub_s:plotlength*sr,:,idx70);
   % red_ramp.diff_traces_70=diff_traces(1:sub_s:plotlength*sr,:,idx70);
    red_ramp.ephys_traces_40=ephys_traces(1:sub_s:plotlength*sr,:,idx40);
   % red_ramp.fit_traces_40=fit_traces(1:sub_s:plotlength*sr,:,idx40);
   % red_ramp.diff_traces_40=diff_traces(1:sub_s:plotlength*sr,:,idx40);
    red_ramp.go_fit=gof_fit;
    red_ramp.fit_param=fit_param;
    red_ramp.bs_photodiodetraces_70=bs_photodiode_signal(1:sub_s:plotlength*sr,:,idx70);
    red_ramp.bs_photodiodetraces_40=bs_photodiode_signal(1:sub_s:plotlength*sr,:,idx40);
end

end