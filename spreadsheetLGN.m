function spreadsheetLGN(data,cellidx)

%Define figure window

fig1 = figure;
set(fig1, 'Name',  char([data(cellidx).patching_date data(cellidx).experimentator data(cellidx).cellname]));
set(fig1, 'Position', [0, 0, 800, 1000]);
field_number = (size(data,1)-1);
set(gcf,'color','w');
disp_ramp_anal = 1;
disp_min_anal = 1;
disp_morph_anal = 0;

if strcmp(data(cellidx).experimentator,'SW')
    srF = 1;
else
    srF = 2;
end

if disp_min_anal & disp_ramp_anal
    rowsn = 4;
    colsn = 5;
elseif disp_min_anal & ~disp_ramp_anal
    rowsn = 4;
    colsn = 4;
else ~disp_min_anal & disp_ramp_anal
    rowsn = 3;
    colsn = 4;
end

plot_n = 1;

%srF=0.1;

redpeak_start       =   100;
redpeak_end         =   350;%looking for the peak shoud be the same as for blue even though the window is longer?
bluepeak_start      =   351;
bluepeak_end        =   400;

% try; set(gcf,'Position',[-159 1094 1500 900]); end
% try; set(gcf,'Position',[-0 58 1065 900]); end
try; set(gcf,'Position',[81 0 1065 900]); end
% try; set(gcf,'Position',[1793 -33 1368 910]); end

% suptitle([data(cellidx).patching_date data(cellidx).experimentator data(cellidx).cellname ', ', ' Inj ord: ', num2str(data(cellidx).eye_inj_ord), ', Contra red: ', num2str(data(cellidx).brain_contra_ipsi) ', ', 'Hemi: ', data(cellidx).hemisphere]);
% set(fig1, 'WindowStyle', 'docked')
%% overview stack
subplot(rowsn,colsn,plot_n); plot_n = plot_n+1; hold off
try
    imshow(data(cellidx).Overview2p);
    title(data(cellidx).animal_name, 'Interpreter', 'none')
catch
    title(data(cellidx).animal_name, 'Interpreter', 'none')
end
%% Morpho;
if disp_morph_anal
    subplot(rowsn,colsn,plot_n); plot_n = plot_n+1; hold off
    try
        if ~isempty(data(cellidx).morphology)==1
            %extract traces for X Y Z
            traces=data(cellidx).morphology.traces;
            somaX=data(cellidx).morphology.somaX;
            somaY=data(cellidx).morphology.somaY;
            somaZ=data(cellidx).morphology.somaZ;
            plot(traces.X,traces.Y,'.','MarkerEdgeColor',[0 0 0],...
                'MarkerFaceColor',[0 0 0],...
                'LineWidth',0.1,'MarkerSize',0.01);
            hold on;
            plot(somaX,somaY,'.','MarkerEdgeColor',[1 0 0],...
                'MarkerFaceColor',[1 0 0],...
                'LineWidth',0.1,'MarkerSize',0.01);
            axis off;axis square
            xlim([-150 150]);ylim([-150 150]);
        else %just plot a point
            plot(1,1,'.','MarkerEdgeColor',[0 0 0],...
                'MarkerFaceColor',[0 0 0],...
                'LineWidth',0.1,'MarkerSize',0.01);
            axis off;axis square;xlim([-150 150]);ylim([-150 150]);
            
        end
        if ~isempty(data(cellidx).DOi)
            hold on;text(-120,-120, ['DOi=' num2str(round(data(cellidx).DOi,3))]);
        else
            %         hold on;text(-120,-120, ['DOi=' num2str(NaN)]);
        end
        title('Morph.');
    catch
        title('Morph.');
        box off
    end
end

%% Ephys traces AMPA
if disp_ramp_anal
    subplot(rowsn,colsn,plot_n); plot_n = plot_n+1; hold off
    try
        for i=1:11
            if ~isempty(data(cellidx).step_red.ephys_traces_70)==1
                plot(data(cellidx).step_red.ephys_traces_70(:,i,2),'linewidth',1,'Color',[0 0 0]+0.05*i);
                hold on;
            else
                plot(0,0);
            end
        end
        if length(data(cellidx).step_red.steps_use_AMPA)==1
            plot(data(cellidx).step_red.ephys_traces_70(:,data(cellidx).step_red.steps_use_AMPA,2),'r')
        elseif length(data(cellidx).step_red.steps_use_AMPA)>1
            plot(nanmean(data(cellidx).step_red.ephys_traces_70(:,data(cellidx).step_red.steps_use_AMPA,2)),'r')
        end
        plot(data(cellidx).step_red.ephys_traces_70(:,data(cellidx).step_red.steps_use_AMPA,2),'r')
        axis square;
        set(gca,'box','off'); set(gca,'TickDir','out');
        title('AMPA Step red, blue constant');
        ylabel('Syn. input (pA)') ;xlabel('Samples');
        
        %%red vertical lines
        hold on;y1=get(gca,'ylim');x1= redpeak_start*srF;hold on;p1=plot([x1 x1],y1,'--','Color','r');p1.Color(4) = 0.3;
        hold on;y1=get(gca,'ylim');x1=redpeak_end*srF;
        hold on;p2=plot([x1 x1],y1,'--','Color','r');p2.Color(4) = 0.3;hold on;
        
        %%blue vertical lines
        y1=get(gca,'ylim');x1=bluepeak_start*srF;hold on;p3=plot([x1 x1],y1,'--','Color','b');p3.Color(4) = 0.3;hold on;
        y1=get(gca,'ylim');x1=bluepeak_end  *srF;hold on;p4=plot([x1 x1],y1,'--','Color','b');
        p4.Color(4) = 0.3;
    catch
        title('AMPA Step red, blue constant'); box off
    end
end
%% Ephys traces NMDA
if disp_ramp_anal
    subplot(rowsn,colsn,plot_n); plot_n = plot_n+1; hold off
    try
        for i=1:11
            if ~isempty(data(cellidx).step_red.ephys_traces_40)==1
                plot(data(cellidx).step_red.ephys_traces_40(:,i,2),'linewidth',1,'Color',[0 0 0]+0.05*i);
                hold on;
            else
                plot(0,0);
            end
        end
        if length(data(cellidx).step_red.steps_use_NMDA)==1
            plot(data(cellidx).step_red.ephys_traces_40(:,data(cellidx).step_red.steps_use_NMDA,2),'r')
        elseif length(data(cellidx).step_red.steps_use_NMDA)>1
            plot(nanmean(data(cellidx).step_red.ephys_traces_40(:,data(cellidx).step_red.steps_use_NMDA,2)),'r')
        end
        axis square;
        set(gca,'box','off'); set(gca,'TickDir','out');
        title('NMDA Step red, blue constant');
        ylabel('Syn. input (pA)') ;xlabel('Samples');
        
        %%red vertical lines
        hold on;y1=get(gca,'ylim');x1= redpeak_start*srF;hold on;p1=plot([x1 x1],y1,'--','Color','r');p1.Color(4) = 0.3;
        hold on;y1=get(gca,'ylim');x1=redpeak_end*srF;
        hold on;p2=plot([x1 x1],y1,'--','Color','r');p2.Color(4) = 0.3;hold on;
        
        %%blue vertical lines
        y1=get(gca,'ylim');x1=bluepeak_start*srF;hold on;p3=plot([x1 x1],y1,'--','Color','b');p3.Color(4) = 0.3;hold on;
        y1=get(gca,'ylim');x1=bluepeak_end  *srF;hold on;p4=plot([x1 x1],y1,'--','Color','b');
        p4.Color(4) = 0.3;
    catch
        title('NMDA Step red, blue constant'); box off
    end
end
%% Ephys NMDA exp fits
if disp_ramp_anal
    subplot(rowsn,colsn,plot_n); plot_n = plot_n+1; hold off;
    try
        if length(data(cellidx).step_red.steps_use_NMDA)==1
            step_numb = data(cellidx).step_red.steps_use_NMDA;
            plot(data(cellidx).step_red.ephys_traces_40(:,step_numb,2),'k','linewidth',0.3); hold on
            plot(data(cellidx).step_red.fit_traces_40(:,step_numb,2),'r--','linewidth',0.3);
        elseif length(data(cellidx).step_red.steps_use_NMDA)>1
            c = colormap(lines);
            for step_numb = data(cellidx).step_red.steps_use_NMDA
                if ~isnan(nanmean(data(cellidx).step_red.go_fit(step_numb,:),2))
                    plot(data(cellidx).step_red.ephys_traces_40(:,step_numb,2),'','linewidth',0.3,'Color',c(step_numb,:)); hold on
                    plot(data(cellidx).step_red.fit_traces_40(:,step_numb,2),'--','linewidth',0.3,'Color',c(step_numb,:));
                end
            end
        end
        
        ylim([min(min(data(cellidx).step_red.ephys_traces_40(:,step_numb,2))) max(max(data(cellidx).step_red.ephys_traces_40(:,step_numb,2)))*1.3])
        axis square;
        set(gca,'box','off'); set(gca,'TickDir','out');
        title('NMDA Fit last 6 traces');
        ylabel('Syn. input (pA)') ;xlabel('Samples');
        
        %%red vertical lines
        y1=get(gca,'ylim'); x1= redpeak_start*srF; plot([x1 x1],y1,'--','Color','r');
        y1=get(gca,'ylim');  x1= redpeak_end*srF;   plot([x1 x1],y1,'--','Color','r');
        
        %%blue vertical lines
        y1=get(gca,'ylim'); x1=bluepeak_start*srF; plot([x1 x1],y1,'--','Color','b');
        y1=get(gca,'ylim'); x1=bluepeak_end*srF;   plot([x1 x1],y1,'--','Color','b');
        
        % disp fit param
        step_numb = data(cellidx).step_red.steps_use_NMDA;
        if length(data(cellidx).step_red.steps_use_NMDA)==1
            temp = data(cellidx).step_red.fit_param(step_numb,:);
            row_used = find(cell2mat(cellfun(@(x)any(isobject(x)),temp,'UniformOutput',false)));
            fit_param = data(cellidx).step_red.fit_param{step_numb,row_used};
            if strcmp(formula(fit_param),'a*exp(b*(x))')
                title({['NMDA red resp. step ' num2str(step_numb)];'single exp. decay fit'});
            else
                title({['NMDA red resp. ' num2str(step_numb)]; 'double exp. decay fit'});
            end
            
            temp_gof=nanmean(data(cellidx).step_red.go_fit(step_numb,row_used));
            text(mean(xlim),mean(ylim), ['gof=' num2str(round(temp_gof,3))]);
            
        elseif length(data(cellidx).step_red.steps_use_NMDA)>1 % needs to be tested in the case whre more than 1 fit was performed
            step_numb = find(~isnan(nanmean(data(cellidx).step_red.go_fit(step_numb,:),2)))+2;
            
            for i = 1:length(step_numb)
                temp = data(cellidx).step_red.fit_param(step_numb(i),:);
                row_used = find(cell2mat(cellfun(@(x)any(isobject(x)),temp,'UniformOutput',false)));
                fit_param = data(cellidx).step_red.fit_param{step_numb(i),row_used};
                temp_gof(i)=nanmean(data(cellidx).step_red.go_fit(step_numb,row_used));
                single_or_double_exp(i) = strcmp(formula(fit_param),'a*exp(b*(x))');
            end
            
            if length(step_numb)==1
                if all(single_or_double_exp)
                    title({['NMDA red resp. step ' num2str(step_numb)];'single exp. decay fit'});
                elseif all(~single_or_double_exp)
                    title({['NMDA red resp. ' num2str(step_numb)]; 'double exp. decay fit'});
                end
                text(500*srF,(max(max(data(cellidx).step_red.ephys_traces_40(:,step_numb,2)))-2)/2, ['gof=' num2str(round(nanmean(temp_gof),3))]);
            elseif length(step_numb)>1
                if all(single_or_double_exp)
                    title({['NMDA red resp. step ' num2str(step_numb)];'single exp. decay fits'});
                elseif all(~single_or_double_exp)
                    title({['NMDA red resp. ' num2str(step_numb)]; 'double exp. decay fits'});
                else
                    title({['NMDA red resp. ' num2str(step_numb)]; 'mixed exp. decay fits'});
                end
                text(mean(xlim),mean(ylim), ['mean gof=' num2str(round(nanmean(temp_gof),3))]);
            end
            
        end
        
        %     %Binary: Did last 6 fitted trace pass threshold?
        %     if ~isempty(data(cellidx).step_red.ephys_traces_40)==1
        %         ratio_suc=sum(data(cellidx).step_blue.pos_fail2(2,6:end))/6;
        %     else
        %         ratio_suc=NaN;
        %     end
        %     text(900*srF,(max(max(data(cellidx).step_red.fit_traces_40(:,:,2)))-2)/2, ['RT=' num2str(round(ratio_suc,3))]);
        
    catch
        title('NMDA Fit last 6 traces'); box off
    end
end
%% residuals during red pulse
if disp_ramp_anal
    subplot(rowsn,colsn,plot_n); plot_n = plot_n+1; hold off
    try
        step_numb = data(cellidx).step_red.steps_use_NMDA;
        if length(step_numb)>1
            for i=1:length(step_numb)
                if ~isempty(data(cellidx).step_red.fit_traces_40)==1
                    %po=plot(data(cellidx).step_red.fit_traces_40(:,i,2),'--','linewidth',0.3,'Color',[0 0 0]+0.05*i); po.Color(4) = 0.25;
                    
                    plot(data(cellidx).step_red.diff_traces_40(:,step_numb(i),2),'linewidth',1,'Color',[0 0 0]+0.05*step_numb(i)); hold on;
                    ymin(step_numb(i))=min(data(cellidx).step_red.diff_traces_40((redpeak_start+50)*srF:redpeak_end*srF,step_numb(i),2));
                    ymax(step_numb(i))=max(data(cellidx).step_red.diff_traces_40((redpeak_start+50)*srF:redpeak_end*srF,step_numb(i),2));
                else
                    plot(0,0);
                end
            end
            plot(nanmean(data(cellidx).step_red.diff_traces_40(:,step_numb,2),2),'r','linewidth',1)
            title({'Red resp selected traces: ' 'fit sub. resp (red: mean)'});
            
        else
            plot(data(cellidx).step_red.diff_traces_40(:,step_numb,2),'k','linewidth',1); hold on
            title({'Red resp selected trace: ' 'fit sub. resp'});
        end
        temp = data(cellidx).step_red.diff_traces_40((redpeak_start+50)*srF:(redpeak_end+50)*srF,step_numb,2);
        xlim([(redpeak_start)*srF , (redpeak_end+50)*srF])
        ylim([min(temp(:)) max(temp(:))]);
        
        axis square;
        set(gca,'box','off'); set(gca,'TickDir','out');
        hold on;
        y1=get(gca,'ylim');x1=(redpeak_start)*srF;hold on;p3=plot([x1 x1],y1,'--','Color','r');p3.Color(4) = 0.3;hold on;
        y1=get(gca,'ylim');x1=(redpeak_end)*srF;hold on;p4=plot([x1 x1],y1,'--','Color','r');
        p4.Color(4) = 0.3;
        xlabel('Samples')
        ylabel('Syn. input (pA)')
    catch
        title('Red resp fit sub. traces'); box off
    end
end
%% residuals during blue pulse
if disp_ramp_anal
    subplot(rowsn,colsn,plot_n); plot_n = plot_n+1; hold off
    try
        step_numb = data(cellidx).step_red.steps_use_NMDA;
        if length(step_numb)>1
            for i=1:length(step_numb)
                if ~isempty(data(cellidx).step_red.fit_traces_40)==1
                    %po=plot(data(cellidx).step_red.fit_traces_40(:,i,2),'--','linewidth',0.3,'Color',[0 0 0]+0.05*i); po.Color(4) = 0.25;
                    
                    plot(data(cellidx).step_red.diff_traces_40(:,step_numb(i),2),'linewidth',1,'Color',[0 0 0]+0.05*step_numb(i)); hold on;
                    ymin(step_numb(i))=min(data(cellidx).step_red.diff_traces_40((redpeak_start+50)*srF:redpeak_end*srF,step_numb(i),2));
                    ymax(step_numb(i))=max(data(cellidx).step_red.diff_traces_40((redpeak_start+50)*srF:redpeak_end*srF,step_numb(i),2));
                else
                    plot(0,0);
                end
            end
            plot(nanmean(data(cellidx).step_red.diff_traces_40(:,step_numb,2),2),'r','linewidth',1)
            title({'Blue resp selected traces: ' 'fit sub. resp (red: mean)'});
            
        else
            plot(data(cellidx).step_red.diff_traces_40(:,step_numb,2),'k','linewidth',1)
            title({'Blue resp selected trace: ' 'fit sub. resp'});
        end
        temp = data(cellidx).step_red.diff_traces_40((bluepeak_start-100)*srF:(bluepeak_start+300)*srF,step_numb,2);
        xlim([(bluepeak_start-100)*srF, (bluepeak_start+300)*srF])
        ylim([min(temp(:)) max(temp(:))]);
        
        axis square;
        set(gca,'box','off'); set(gca,'TickDir','out');
        %     title('Residuals blue pulse last 6 traces');
        hold on;
        y1=get(gca,'ylim');x1=(bluepeak_start)*srF;hold on;p3=plot([x1 x1],y1,'--','Color','b');p3.Color(4) = 0.3;hold on;
        y1=get(gca,'ylim');x1=(bluepeak_end)*srF;hold on;p4=plot([x1 x1],y1,'--','Color','b');
        p4.Color(4) = 0.3;
        xlabel('Samples')
        ylabel('Syn. input (pA)')
    catch
        title('Blue resp fit sub. traces'); box off
    end
end
%% peak amp vs step #
subplot(rowsn,colsn,plot_n); plot_n = plot_n+1; hold off;
try
    step_numb_A = data(cellidx).step_red.steps_use_AMPA; %% change to AMPA when error is fixed in str
    step_numb_N = data(cellidx).step_red.steps_use_NMDA;
    
    AMPA_1 = data(cellidx).step_red.neg_mean1(2,1:end)-data(cellidx).step_red.neg_base_mean1(2,1:end);
    AMPA_2 = data(cellidx).step_blue.neg_mean2(2,1:end)-data(cellidx).step_blue.neg_base_mean2(2,1:end);
    NMDA_1 = data(cellidx).step_red.pos_mean1(2,1:end)-data(cellidx).step_red.pos_base_mean1(2,1:end);
    NMDA_2 = data(cellidx).step_blue.pos_mean2(2,1:end)-data(cellidx).step_blue.pos_base_mean2(2,1:end);
    
    plot(1:11,AMPA_2./abs(min([AMPA_1, AMPA_2])).*100,'color', [0.2 1 0.2]); hold on;
    plot(1:11,AMPA_1./abs(min([AMPA_1, AMPA_2])).*100,'color', [1 0.2 0.2]);
    plot(1:11,NMDA_2./abs(max([NMDA_1, NMDA_2])).*100,'color', [0 0.5 0]);
    plot(1:11,NMDA_1./abs(max([NMDA_1, NMDA_2])).*100,'color', [0.5 0 0]);
    
    plot(step_numb_A,AMPA_2(step_numb_A)./abs(min([AMPA_1, AMPA_2])).*100,'o','color', [0.2 1 0.2])
    plot(step_numb_A,AMPA_1(step_numb_A)./abs(min([AMPA_1, AMPA_2])).*100,'o','color', [1 0.2 0.2])
    plot(step_numb_N,NMDA_2(step_numb_N)./abs(max([NMDA_1, NMDA_2])).*100,'o','color', [0 0.5 0]);
    plot(step_numb_N,NMDA_1(step_numb_N)./abs(max([NMDA_1, NMDA_2])).*100,'o','color', [0.5 0 0]);
    
    xlim([0 12]);
    axis square; set(gca,'TickDir','out');
    set(gca,'box','off');ylabel('Synaptic input (%)'); xlabel('Step Nr');
    %legend('AMPA b', 'AMPA r','NMDA b', 'NMDA r');
    title('Synaptic input vs step nr');
catch
    title('Synaptic input vs step nr');
end

%% Laser irradiance
subplot(rowsn,colsn,plot_n); plot_n = plot_n+1; hold off
try
    step_numb_A = data(cellidx).step_red.steps_use_AMPA; %% change to AMPA when error is fixed in str
    step_numb_N = data(cellidx).step_red.steps_use_NMDA;
    
    plot(1:11,data(cellidx).step_blue.neg_irr_blue(2,1:end),'color', [0.2 1 0.2]); hold on;
    plot(1:11,data(cellidx).step_red.neg_irr_red(2,1:end),'color', [1 0.2 0.2]);
    plot(1:11,data(cellidx).step_blue.pos_irr_blue(2,1:end),'color', [0 0.5 0]);
    plot(1:11,data(cellidx).step_red.pos_irr_red(2,1:end),'color', [0.5 0 0]);
    
    plot(step_numb_A,data(cellidx).step_blue.neg_irr_blue(2,step_numb_A),'o','color', [0.2 1 0.2]); hold on;
    plot(step_numb_A,data(cellidx).step_red.neg_irr_red(2,step_numb_A),'o','color', [1 0.2 0.2]);
    plot(step_numb_N,data(cellidx).step_blue.pos_irr_blue(2,step_numb_N),'o','color', [0 0.5 0]);
    plot(step_numb_N,data(cellidx).step_red.pos_irr_red(2,step_numb_N),'o','color', [0.5 0 0]);
    
    xlim([0 12]);
    axis square; box off; set(gca,'TickDir','out');
    set(gca,'box','off');ylabel('Irradiance (mW/mm2)');xlabel('Step Nr');
    %legend('AMPA b', 'AMPA r','NMDA b', 'NMDA r');
    title('Irradiance vs step nr');
catch
    title('Irradiance vs step nr');
end

%% automatic vs human classification missmatch
subplot(rowsn,colsn,plot_n); plot_n = plot_n+1; hold off;
try
    
    % region lines
    patch([-1 1 1 -1],[-1 -1 1 1]                       ,'b','FaceColor','none', 'EdgeColor', 'b', 'LineStyle', '--'); hold on %  true bino
    patch([1.4 1.6 1.6 1.4],[1.4 1.4 1.6 1.6]           ,'r','FaceColor','none', 'EdgeColor', 'r', 'LineStyle', '--') % contra only
    patch([-1.4 -1.6 -1.6 -1.4],[-1.4 -1.4 -1.6 -1.6]   ,'g','FaceColor','none', 'EdgeColor', 'g', 'LineStyle', '--') % ipsi only
    patch([1.4 1.6 1.6 1.4],[-1 -1 1 1]                 ,'m','FaceColor','none', 'EdgeColor', 'm', 'LineStyle', '--') % contra & silent ipsi
    patch([-1.4 -1.6 -1.6 -1.4],[-1 -1 1 1]             ,'c','FaceColor','none', 'EdgeColor', 'c', 'LineStyle', '--') % ipsi & contra ipsi
    patch([-1 1 1 -1],[-1.4 -1.4 -1.6 -1.6]             ,'k','FaceColor','none', 'LineStyle', '--') % ipsi & contra ipsi
    patch([-1 1 1 -1],[1.4 1.4 1.6 1.6]                 ,'k','FaceColor','none', 'LineStyle', '--') % ipsi & contra ipsi
    patch([-1.4 -1.6 -1.6 -1.4],[1.4 1.4 1.6 1.6]       ,'k','FaceColor','none', 'LineStyle', '--') % contra only
    patch([1.4 1.6 1.6 1.4],[-1.4 -1.4 -1.6 -1.6]       ,'k','FaceColor','none', 'LineStyle', '--') % contra only
    
    % make nice
    ylabel('Step prot. NMDA ODI'); xlabel('Step prot. AMPA ODI');
    xlim([-1.7 1.7]); ylim([-1.7 1.7]);
    yticks([-1.5 -1 -0.5 0 0.5 1 1.5])
    yticklabels({'-1' '-0.999' '-0.5' '0' '0.5' '0.999' '1'})
    xticks([-1.5 -1 -0.5 0 0.5 1 1.5])
    xticklabels({'-1' '-0.999' '-0.5' '0' '0.5' '0.999' '1'})
    set(gca,'TickDir','out'); box on
    set(gcf,'color','w');
    title('Human vs Automatic ODI class');
    
    % make temporary var
    clear data_temp
    if data(cellidx).ocular_category == 0 %%|| data(cellidx).ocular_category == 6
        data_temp.ocular_category = nan;
    else
        data_temp.ocular_category = data(cellidx).ocular_category;
    end
    
    % shift the monocular cells by 0.5 +/- a bit to more to clearly display them in the scatter plot
    if abs(data(cellidx).ODI_AMPA_step)<1
        data_temp.ODI_AMPA_step = data(cellidx).ODI_AMPA_step;
    elseif data(cellidx).ODI_AMPA_step == 1
        data_temp.ODI_AMPA_step = 1.5; %+(rand*0.2-0.1);
    elseif data(cellidx).ODI_AMPA_step == -1
        data_temp.ODI_AMPA_step = -1.5; %+(rand*0.2-0.1);
    elseif isnan(data(cellidx).ODI_AMPA_step)
        data_temp.ODI_AMPA_step = nan;
    end
    
    if abs(data(cellidx).ODI_NMDA_step)<1
        data_temp.ODI_NMDA_step = data(cellidx).ODI_NMDA_step;
    elseif data(cellidx).ODI_NMDA_step == 1
        data_temp.ODI_NMDA_step = 1.5; %+(rand*0.2-0.1);
    elseif data(cellidx).ODI_NMDA_step == -1
        data_temp.ODI_NMDA_step = -1.5; %+(rand*0.2-0.1);
    elseif isnan(data(cellidx).ODI_NMDA_step)
        data_temp.ODI_NMDA_step = nan;
    end
    
    % set color based on human annotation
    if data_temp.ocular_category == 1
        cell_cat_col = 'r'; % RED: contra_AMPA_contra_NMDA
    elseif data_temp.ocular_category == 2
        cell_cat_col = 'g'; % GREEN: ipsi_AMPA_ipsi_NMDA
    elseif data_temp.ocular_category == 3
        cell_cat_col = 'b'; % Blue: bino_AMPA_bino_NMDA
    elseif data_temp.ocular_category == 4
        cell_cat_col = 'm'; % Magenta: contra_AMPA_bino_NMDA_ipsi silent
    elseif data_temp.ocular_category == 5
        cell_cat_col = 'c'; % Cyan: ipsi_AMPA_bino_NMDA_contra silent
    elseif data_temp.ocular_category == 6
        cell_cat_col = 'k'; % Cyan: ipsi_AMPA_bino_NMDA_contra silent
    end
    
    % plot
    if ~isnan(data_temp.ocular_category) && ~isnan(data_temp.ODI_AMPA_step) && ~isnan(data_temp.ODI_NMDA_step)
        scatter([data_temp.ODI_AMPA_step],[data_temp.ODI_NMDA_step],20,cell_cat_col,'filled');
    elseif isnan(data_temp.ocular_category) && ~isnan(data_temp.ODI_AMPA_step) && ~isnan(data_temp.ODI_NMDA_step)
        scatter([data_temp.ODI_AMPA_step],[data_temp.ODI_NMDA_step],20,'k');
    elseif ~isnan(data_temp.ODI_AMPA_step)
        plot([data_temp.ODI_AMPA_step data_temp.ODI_AMPA_step], [-1.7 1.7],cell_cat_col)
    else
        text(-0.3,0,'NaN')
    end
    
    clear data_temp
    axis square;
    
catch
    title('Human vs Automatic ODI class');
    text(-0.3,0,'NaN')
    axis square;
end


%% AMPA to NMDA ratio
subplot(rowsn,colsn,plot_n); plot_n = plot_n+1; hold off
try
    clear cell_cat_col
    for i = 1:length(data)
        red_resp_AMPA_peak = [data(i).step_red.red_resp_AMPA_peak];
        blue_resp_AMPA_peak = [data(i).step_red.blue_resp_AMPA_peak];
        red_resp_NMDA_peak = [data(i).step_red.red_resp_NMDA_peak];
        blue_resp_NMDA_peak = [data(i).step_red.blue_resp_NMDA_peak];
        
        red_resp_AMPA_mean = [data(i).step_red.red_resp_AMPA];
        blue_resp_AMPA_mean = [data(i).step_red.blue_resp_AMPA];
        red_resp_NMDA_mean = [data(i).step_red.red_resp_NMDA];
        blue_resp_NMDA_mean = [data(i).step_red.blue_resp_NMDA];
        
        
        red_ratio = (red_resp_AMPA_peak - red_resp_NMDA_peak)./(red_resp_AMPA_peak + red_resp_NMDA_peak);
        blue_ratio = (blue_resp_AMPA_peak - blue_resp_NMDA_peak)./(blue_resp_AMPA_peak + blue_resp_NMDA_peak);
        
        if data(i).brain_contra_ipsi
            contra_ratio(i) = red_ratio;
            ipsi_ratio(i) = blue_ratio;
        else
            contra_ratio(i) = blue_ratio;
            ipsi_ratio(i) = red_ratio;
        end
        
        %         if data(i).ODI_AMPA_step == 1
        %             dominant_ratio(i) = contra_ratio(i);
        %             nondominant_ratio(i) = nan;
        %         elseif data(i).ODI_AMPA_step == -1
        %             dominant_ratio(i) = ipsi_ratio(i);
        %             nondominant_ratio(i) = nan;
        
        if data(i).ODI_AMPA_step >=0
            dominant_ratio(i) = contra_ratio(i);
            nondominant_ratio(i) = ipsi_ratio(i);
        elseif data(i).ODI_AMPA_step < 0
            dominant_ratio(i) = ipsi_ratio(i);
            nondominant_ratio(i) = contra_ratio(i);
        else
            dominant_ratio(i) = nan;
            nondominant_ratio(i) = nan;
        end
        
        % calculating Rs change here just to have it at hand
        Rs_step_AMPA = nanmean(data(i).step_blue.neg_Rs(2,:));
        if ~isempty(data(i).step_blue.pos_Rs)
            Rs_step_NMDA = nanmean(data(i).step_blue.pos_Rs(2,:));
        else
            Rs_step_NMDA = nan;
        end
        Rs_change(i) = (Rs_step_NMDA-Rs_step_AMPA)/(Rs_step_AMPA)*100;
        Rs_AMPA(i) = Rs_step_AMPA;
        Rs_NMDA(i) = Rs_step_NMDA;
        %     cell_cat_col(i) = log(Rs_change(i));
        cell_cat_col(i) = log(Rs_NMDA(i));
        %     cell_cat_col(i) = data(i).MD;
    end
    
    if ~isnan(dominant_ratio(cellidx)) & ~isnan(nondominant_ratio(cellidx))
        plotSpread([dominant_ratio(:) nondominant_ratio(:)]); hold on;
        scatter([1 2], [dominant_ratio(cellidx) nondominant_ratio(cellidx)],'r','filled')
        set(gca,'XTickLabel', {'dominant','non-dominant'});
        ylabel({'AMPA to NMDA response ratios' '(A-N)/(A+N)'})
        title({'AMPA to NMDA response ratio' ' for both eye inputs'})
        
    elseif isnan(dominant_ratio(cellidx)) & ~isnan(nondominant_ratio(cellidx))
        plotSpread([dominant_ratio(:) nondominant_ratio(:)]); hold on;
        scatter([1], [nondominant_ratio(cellidx)],'r','filled')
        set(gca,'XTickLabel', {'dominant','non-dominant'});
        ylabel({'AMPA to NMDA response ratios' '(A-N)/(A+N)'})
        title({'AMPA to NMDA response ratio' ' for both eye inputs'})
        
    elseif ~isnan(dominant_ratio(cellidx)) & isnan(nondominant_ratio(cellidx))
        plotSpread([dominant_ratio(:) nondominant_ratio(:)]); hold on;
        scatter(1,[dominant_ratio(cellidx)],'r','filled')
        set(gca,'XTickLabel', {'dominant','non-dominant'});
        ylabel({'AMPA to NMDA response ratios' '(A-N)/(A+N)'})
        title({'AMPA to NMDA response ratio' ' for both eye inputs'})
        
    elseif isnan(dominant_ratio(cellidx)) & isnan(nondominant_ratio(cellidx))
        %?
    end
    set(gca,'TickDir','out');
    axis square;
end

%% Failure Irradiance
if disp_min_anal
    subplot(rowsn,colsn,plot_n); plot_n = plot_n+1; hold off
    try
        if ~isempty(data(cellidx).red_failure_AMPA)
            
            % switch back to stim color instead of ips/contra
            % also plot std to see if their was an irradiance change
            % during the protocol.
            IR_red_AMPA_mean = nanmean(data(cellidx).red_failure_AMPA.IR_pulse1);
            IR_blue_AMPA_mean = nanmean(data(cellidx).blue_failure_AMPA.IR_pulse1);
            
            IR_red_NMDA_mean = nanmean(data(cellidx).red_failure_NMDA.IR_pulse1);
            IR_blue_NMDA_mean = nanmean(data(cellidx).blue_failure_NMDA.IR_pulse1);
            
            IR_red_AMPA_std = nanstd(data(cellidx).red_failure_AMPA.IR_pulse1);
            IR_blue_AMPA_std = nanstd(data(cellidx).blue_failure_AMPA.IR_pulse1);
            
            IR_red_NMDA_std = nanstd(data(cellidx).red_failure_NMDA.IR_pulse1);
            IR_blue_NMDA_std = nanstd(data(cellidx).blue_failure_NMDA.IR_pulse1);
            
            
            % plot
            hold off;
            b = bar([IR_red_AMPA_mean, IR_blue_AMPA_mean; IR_red_NMDA_mean, IR_blue_NMDA_mean]); hold on;
            barwidth = b(1).BarWidth/6;
            er = errorbar([1-barwidth 1+barwidth 2-barwidth 2+barwidth],...
                [IR_red_AMPA_mean, IR_blue_AMPA_mean, IR_red_NMDA_mean, IR_blue_NMDA_mean],...
                [IR_red_AMPA_std, IR_blue_AMPA_std, IR_red_NMDA_std, IR_blue_NMDA_std],'.');
            er.Color = [0 0 0];
            er.LineStyle = 'none';
            title('Failure Irradiance')
            
            % make nice
            b(1).FaceColor = 'r';
            b(2).FaceColor = 'g';
            set(gca,'XTickLabel',{'AMPA', 'NMDA'})
            box off; set(gca,'TickDir','out');
            ylabel(['Irradiance mean (mW/mm2)'])
            axis square;
        end
    catch
        title('Failure Irradiance'); box off
    end
end

%% Failure avg responses to AMPA & NMDA
if disp_min_anal
    
    for i = 1:4
        trace_cut = 500*srF;
        
        if  i == 1 &  ~isempty(data(cellidx).red_failure_AMPA) & ~isempty(data(cellidx).red_failure_AMPA.peaks)
            % create new subplot for average response and photodiode signal
            subplot(rowsn,colsn,plot_n); plot_n = plot_n+1; hold off
            % specify title
            title({'Minimal stim: red at -70mV'; 'and photodiode signal'}); hold on
            
            % define traces
            trace_mean = mean(data(cellidx).red_failure_AMPA.traces_all(1:trace_cut,:),2);
            trace_std = std(data(cellidx).red_failure_AMPA.traces_all(1:trace_cut,:)')'/sqrt(size(data(cellidx).red_failure_AMPA.traces_all,2));
            trace_bound_up = trace_mean+trace_std;
            trace_bound_bottom = trace_mean-trace_std;
            scalling_fac = max(abs(trace_mean))/max(data(cellidx).red_failure_AMPA.photodiode_mean_std(1,1:trace_cut)') ;
            photodiode_mean = data(cellidx).red_failure_AMPA.photodiode_mean_std(1,1:trace_cut)'.*scalling_fac;
            photodiode_std = data(cellidx).red_failure_AMPA.photodiode_mean_std(2,1:trace_cut)'.*scalling_fac;
            phototrace_bound_up = photodiode_mean+photodiode_std;
            phototrace_bound_bottom = photodiode_mean-photodiode_std;
            color1 = [1 0 0];
            color2 = [1 0.5 0.5];
            
        elseif i == 2 &  ~isempty(data(cellidx).blue_failure_AMPA) & ~isempty(data(cellidx).blue_failure_AMPA.peaks)
            subplot(rowsn,colsn,plot_n); plot_n = plot_n+1; hold off
            title({'Minimal stim: blue at -70mV'; 'and photodiode signal'}); hold on
            trace_mean = mean(data(cellidx).blue_failure_AMPA.traces_all(1:trace_cut,:),2);
            trace_std = std(data(cellidx).blue_failure_AMPA.traces_all(1:trace_cut,:)')'/sqrt(size(data(cellidx).blue_failure_AMPA.traces_all,2));
            trace_bound_up = trace_mean+trace_std;
            trace_bound_bottom = trace_mean-trace_std;
            scalling_fac = max(abs(trace_mean))/max(data(cellidx).blue_failure_AMPA.photodiode_mean_std(1,1:trace_cut)') ;
            photodiode_mean = data(cellidx).blue_failure_AMPA.photodiode_mean_std(1,1:trace_cut)'.*scalling_fac;
            photodiode_std = data(cellidx).blue_failure_AMPA.photodiode_mean_std(2,1:trace_cut)'.*scalling_fac;
            phototrace_bound_up = photodiode_mean+photodiode_std;
            phototrace_bound_bottom = photodiode_mean-photodiode_std;
            color1 = [0 1 0];
            color2 = [0.5 1 0.5];
            
        elseif  i == 3 &  ~isempty(data(cellidx).red_failure_NMDA) & ~isempty(data(cellidx).red_failure_NMDA.peaks)
            subplot(rowsn,colsn,plot_n); plot_n = plot_n+1; hold off
            title({'Minimal stim: red at +40mV'; 'and photodiode signal'}); hold on
            trace_mean = mean(data(cellidx).red_failure_NMDA.traces_all(1:trace_cut,:),2);
            trace_std = std(data(cellidx).red_failure_NMDA.traces_all(1:trace_cut,:)')'/sqrt(size(data(cellidx).red_failure_AMPA.traces_all,2));
            trace_bound_up = trace_mean+trace_std;
            trace_bound_bottom = trace_mean-trace_std;
            scalling_fac = max(abs(trace_mean))/max(data(cellidx).red_failure_NMDA.photodiode_mean_std(1,1:trace_cut)') ;
            photodiode_mean = data(cellidx).red_failure_NMDA.photodiode_mean_std(1,1:trace_cut)'.*scalling_fac;
            photodiode_std = data(cellidx).red_failure_NMDA.photodiode_mean_std(2,1:trace_cut)'.*scalling_fac;
            phototrace_bound_up = photodiode_mean+photodiode_std;
            phototrace_bound_bottom = photodiode_mean-photodiode_std;
            color1 = [1 0 0];
            color2 = [1 0.5 0.5];
            
        elseif i == 4 &  ~isempty(data(cellidx).blue_failure_NMDA) & ~isempty(data(cellidx).blue_failure_NMDA.peaks)
            subplot(rowsn,colsn,plot_n); plot_n = plot_n+1; hold off
            title({'Minimal stim: blue at +40mV'; 'and photodiode signal'}); hold on
            trace_mean = mean(data(cellidx).blue_failure_NMDA.traces_all(1:trace_cut,:),2);
            trace_std = std(data(cellidx).blue_failure_NMDA.traces_all(1:trace_cut,:)')'/sqrt(size(data(cellidx).blue_failure_NMDA.traces_all,2));
            trace_bound_up = trace_mean+trace_std;
            trace_bound_bottom = trace_mean-trace_std;
            scalling_fac = max(abs(trace_mean))/max(data(cellidx).blue_failure_NMDA.photodiode_mean_std(1,1:trace_cut)') ;
            photodiode_mean = data(cellidx).blue_failure_NMDA.photodiode_mean_std(1,1:trace_cut)'.*scalling_fac;
            photodiode_std = data(cellidx).blue_failure_NMDA.photodiode_mean_std(2,1:trace_cut)'.*scalling_fac;
            phototrace_bound_up = photodiode_mean+photodiode_std;
            phototrace_bound_bottom = photodiode_mean-photodiode_std;
            color1 = [0 1 0];
            color2 = [0.5 1 0.5];
            
        else
            continue
        end
        patch([[1:length(photodiode_mean)],[length(photodiode_mean):-1:1]]',...
            [phototrace_bound_up; flip(phototrace_bound_bottom)],[0.5 0.5 0.5],'EdgeColor','none'); hold on
        plot(photodiode_mean,'k')
        patch([[1:length(trace_mean)],[length(trace_mean):-1:1]]',...
            [trace_bound_up; flip(trace_bound_bottom)],color2,'EdgeColor','none'); hold on
        plot(trace_mean,'color',color1)
        
        ylabel({'Mean amplitude'; 'of response (pA)'});
        xlabel('Samples');
        set(gca,'box','off'); set(gca,'TickDir','out');
        axis square;
    end
end

%% Failure trial responses AMPA & NMDA
if disp_min_anal
    for i = 1:4
        if  i == 1 &  ~isempty(data(cellidx).red_failure_AMPA) & ~isempty(data(cellidx).red_failure_AMPA.peaks)
            try % Red AMPA
                % create new subplot for average response and photodiode signal
                subplot(rowsn,colsn,plot_n); plot_n = plot_n+1; hold off
                
                % extract values to be plotted
                resps = abs(data(cellidx).red_failure_AMPA.peaks);
                steady_state = data(cellidx).red_failure_AMPA.steady_state;
                if isnan(steady_state)
                    steady_state = 1;
                end
                threshold = data(cellidx).red_failure_AMPA.resp_thresh;
                data_color = [1 0 0];
                
                % calc response probability and average response amplitude
                resps_steady = resps(steady_state:end);
                resp_prob = sum(resps_steady>nanmean(threshold))/length(resps_steady)*100;
                avg_amp = mean(resps_steady(resps_steady>nanmean(threshold)));
                
                % specify title
                title({'Minimal stim: red at -70mV'; ['Rprob: ' num2str(resp_prob,3) '%, avg amp: ' num2str(avg_amp,3)]}); hold on
            end
        elseif i == 2 &  ~isempty(data(cellidx).blue_failure_AMPA) & ~isempty(data(cellidx).blue_failure_AMPA.peaks)
            try % Blue AMPA
                subplot(rowsn,colsn,plot_n); plot_n = plot_n+1; hold off
                resps = abs(data(cellidx).blue_failure_AMPA.peaks);
                steady_state = data(cellidx).blue_failure_AMPA.steady_state;
                if isnan(steady_state); steady_state = 1; end
                threshold = data(cellidx).blue_failure_AMPA.resp_thresh;
                data_color = [0 1 0];
                resps_steady = resps(steady_state:end);
                resp_prob = sum(resps_steady>nanmean(threshold))/length(resps_steady)*100;
                avg_amp = mean(resps_steady(resps_steady>nanmean(threshold)));
                title({'Minimal stim: blue at -70mV'; ['Rprob: ' num2str(resp_prob,3) '%, avg amp: ' num2str(avg_amp,3)]}); hold on
            end
        elseif  i == 3 &  ~isempty(data(cellidx).red_failure_NMDA) & ~isempty(data(cellidx).red_failure_NMDA.peaks)
            try % Red NMDA
                subplot(rowsn,colsn,plot_n); plot_n = plot_n+1; hold off
                resps = abs(data(cellidx).red_failure_NMDA.peaks);
                steady_state = data(cellidx).red_failure_NMDA.steady_state;
                if isnan(steady_state); steady_state = 1; end
                threshold = data(cellidx).red_failure_NMDA.resp_thresh;
                data_color = [1 0 0];
                resps_steady = resps(steady_state:end);
                resp_prob = sum(resps_steady>nanmean(threshold))/length(resps_steady)*100;
                avg_amp = mean(resps_steady(resps_steady>nanmean(threshold)));
                title({'Minimal stim: red at +40mV'; ['Rprob: ' num2str(resp_prob,3) '%, avg amp: ' num2str(avg_amp,3)]}); hold on
            end
        elseif i == 4 &  ~isempty(data(cellidx).blue_failure_NMDA) & ~isempty(data(cellidx).blue_failure_NMDA.peaks)
            try% Blue NMDA
                subplot(rowsn,colsn,plot_n); plot_n = plot_n+1; hold off
                resps = abs(data(cellidx).blue_failure_NMDA.peaks);
                steady_state = data(cellidx).blue_failure_NMDA.steady_state;
                if isnan(steady_state); steady_state = 1; end
                threshold = data(cellidx).blue_failure_NMDA.resp_thresh;
                data_color = [0 1 0];
                resps_steady = resps(steady_state:end);
                resp_prob = sum(resps_steady>nanmean(threshold))/length(resps_steady)*100;
                avg_amp = mean(resps_steady(resps_steady>nanmean(threshold)));
                title({'Minimal stim: blue at +40mV'; ['Rprob: ' num2str(resp_prob,3) '%, avg amp: ' num2str(avg_amp,3)]}); hold on
            end
        else
            continue
        end
        
        % plot
        try xlim([1 length(resps)]);
            plot(resps,'.','color',data_color); hold on
            plot([steady_state steady_state], ylim,'k--')
            plot(xlim, [mean(threshold) mean(threshold)],'k--')
            plot([threshold],'k--')
        end
        
        % make nice
        ylabel('Peak PSC amplitude (pA)');
        xlabel('Trial');
        set(gca,'box','off'); 
        set(gca,'TickDir','out');
        axis square;
    end
end

%% Series resistance change between protocols
subplot(rowsn,colsn,plot_n); plot_n = plot_n+1; hold off
try
    Rs_step_AMPA = nanmean(data(cellidx).step_blue.neg_Rs(2,:));
    if ~isempty(data(cellidx).step_blue.pos_Rs)
        Rs_step_NMDA = nanmean(data(cellidx).step_blue.pos_Rs(2,:));
    else
        Rs_step_NMDA = nan;
    end
    
    try % Red AMPA
        if ~isempty(data(cellidx).red_failure_AMPA) || ~isempty(data(cellidx).red_failure_AMPA.Rs_cell)
            Rs_fail_red_AMPA=nanmean(data(cellidx).red_failure_AMPA.Rs_cell);
        else
            Rs_fail_red_AMPA = nan;
        end
    catch
        Rs_fail_red_AMPA = nan;
    end
    
    try % Blue AMPA
        if ~isempty(data(cellidx).blue_failure_AMPA) || ~isempty(data(cellidx).blue_failure_AMPA.Rs_cell)
            Rs_fail_blue_AMPA=nanmean(data(cellidx).blue_failure_AMPA.Rs_cell);
        else
            Rs_fail_blue_AMPA = nan;
        end
    catch
        Rs_fail_blue_AMPA = nan;
    end
    
    try % Red NMDA
        if ~isempty(data(cellidx).red_failure_NMDA) || ~isempty(data(cellidx).red_failure_NMDA.Rs_cell)
            Rs_fail_red_NMDA=nanmean(data(cellidx).red_failure_NMDA.Rs_cell);
        else
            Rs_fail_red_NMDA = nan;
        end
    catch
        Rs_fail_red_NMDA = nan;
    end
    
    try % Blue NMDA
        if ~isempty(data(cellidx).blue_failure_NMDA) || ~isempty(data(cellidx).blue_failure_NMDA.Rs_cell)
            Rs_fail_blue_NMDA=nanmean(data(cellidx).blue_failure_NMDA.Rs_cell);
        else
            Rs_fail_blue_NMDA = nan;
        end
    catch
        Rs_fail_blue_NMDA = nan;
    end
    
    
    test1 = (Rs_step_NMDA-Rs_step_AMPA)/(Rs_step_AMPA)*100; % necessary for AMPA to NMDA ratio and ODI comparisons
    
    test2 = (Rs_fail_red_AMPA-Rs_step_AMPA)/(Rs_step_AMPA)*100; % necessary for fiber fraction calulcations
    test3 =  (Rs_fail_blue_AMPA-Rs_step_AMPA)/(Rs_step_AMPA)*100; % necessary for fiber fraction calulcations
    test4 =  (Rs_fail_blue_AMPA-Rs_fail_red_AMPA)/(Rs_fail_red_AMPA)*100; % necessary for fiber fraction calulcations
    
    test5 =  (Rs_fail_red_NMDA-Rs_step_NMDA)/(Rs_step_NMDA)*100;
    test6 =  (Rs_fail_blue_NMDA-Rs_step_NMDA)/(Rs_step_NMDA)*100;
    test7 =  (Rs_fail_blue_NMDA-Rs_fail_red_NMDA)/(Rs_fail_red_NMDA)*100;
    
    if all(isnan([test1 test2 test3 test4 test5 test6 test7]))
        error
    end
    
    bar([test1, test2, test3, test4, test5, test6, test7]); hold on
    title('Rs change');
    
    set(gca,'XTickLabel',{'StepA vs StepN', 'FailAcon vs StepA', 'FailAipsi vs StepA', 'FailAcon vs FailAipsi', 'FailNcon vs StepN' , 'FailNipsi vs StepN', 'FailNcon vs FailNipsi'})
    xtickangle(45)
    box off; set(gca,'TickDir','out');
    ylabel('% change in Rs')
    axis square;
catch
    title('Rs change');
    set(gca,'XTickLabel',{'StepA vs StepN', 'FailAcon vs StepA', 'FailAipsi vs StepA', 'FailAcon vs FailAipsi', 'FailNcon vs StepN' , 'FailNipsi vs StepN', 'FailNcon vs FailNipsi'})
    xtickangle(45)
    box off; set(gca,'TickDir','out');
    xlim([-0.2    8.2])
    ylabel('% change in Rs')
    text(3,0.5,'NaN')
    axis square;
end

%% Series resistance over time
% subplot(rowsn,colsn,plot_n); plot_n = plot_n+1; hold off
% try
% catch
% end

%%
%suptitle([data(cellidx).patching_date data(cellidx).experimentator data(cellidx).cellname ', ', ' Inj ord: ', num2str(data(cellidx).eye_inj_ord), ', Contra red: ', num2str(data(cellidx).brain_contra_ipsi) ', ', 'Hemi: ', data(cellidx).hemisphere]);


end