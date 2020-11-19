% Sample plot routines

Nres = length(Reservoir);
nModes = 2;
LW = 2;
cmap_perso = [51 51 255; 0 204 51; 204 0 0; 204 153 0; 153 0 102; 51 153 153; 204 102 204; 204 204 102]/255;

ModeName = {'Car','PT'};
%% Accumulation
figure
hold on
for i_m = 1:nModes
    for ires = 1:Nres
        modeindex = Reservoir(ires).ModeIndex{i_m};
        subplot(1,nModes,i_m)
        plot(CurrentTime,sum(Reservoir(ires).AccPerRoute(modeindex,:),1),'color',cmap_perso(ires,:),'linewidth',LW)
        title(ModeName{i_m})
        hold on
    end
end
if Nres == 2
    legend('R1','R2')
elseif Nres == 3
    legend('R1','R2','R3')
elseif Nres == 4
    legend('R1','R2','R3','R4')
end
set(gcf,'Position',[10 10 1400 2*300])

%% Inflow

figure
hold on
for i_m = 1:nModes
    demand = zeros(size(CurrentTime));
    for ires = 1:Nres
        modeindex = Reservoir(ires).ModeIndex{i_m};
        subplot(1,nModes,i_m)
        plot(CurrentTime,sum(Reservoir(ires).InflowPerRoute(modeindex,:),1),'color',cmap_perso(ires,:),'linewidth',LW)
        title(ModeName{i_m})
        hold on
        if ires == 1
            for im = modeindex
                demand = demand + Route(im).Demand;
            end
            plot(CurrentTime,demand,'k','linewidth',LW)
        end
    end
    
end
if Nres == 2
    legend('R1','Dem','R2')
elseif Nres == 3
    legend('R1','Dem','R2','R3')
elseif Nres == 4
    legend('R1','Dem','R2','R3','R4')
end
set(gcf,'Position',[10 10 1400 2*300])


%% Outflow

figure
hold on
for i_m = 1:nModes
    for ires = 1:Nres
        modeindex = Reservoir(ires).ModeIndex{i_m};
        subplot(1,nModes,i_m)
        plot(CurrentTime,sum(Reservoir(ires).OutflowPerRoute(modeindex,:),1),'color',cmap_perso(ires,:),'linewidth',LW)
        title(ModeName{i_m})
        hold on
    end  
end
if Nres == 2
    legend('R1','R2')
elseif Nres == 3
    legend('R1','R2','R3')
elseif Nres == 4
    legend('R1','R2','R3','R4')
end
set(gcf,'Position',[10 10 1400 2*300])

%% Inflow vs demand for each route in RES 1

% figure
% hold on
% for i = 1:length(Reservoir(1).InflowPerRoute(:,1))
%     plot(Reservoir(1).InflowPerRoute(i,:),'-','color',cmap_perso(i,:),'linewidth',LW)
%     plot(Route(i).Demand,'--','color',cmap_perso(i,:),'linewidth',LW)
% end
% legend('Inflow','Demand')
% set(gcf,'Position',[10 10 1400 2*300])

%% Compute Delay times
% figure
% hold on
% for i = 1:length(Reservoir(1).InflowPerRoute(:,1))
%     plot(Reservoir(1).InflowPerRoute(i,:),'-','color',cmap_perso(i,:))
%     plot(Route(i).Demand,'--','color',cmap_perso(i,:))
% end
% 
% figure
% hold on
% for i = 1:length(Reservoir(1).InflowPerRoute(:,1))
%     InflowCumulative(i,:) = cumsum(Reservoir(1).InflowPerRoute(i,:));
%     DemCumulative(i,:) = cumsum(Route(i).Demand);
%     nVeh = round(DemCumulative(i,end));
%     InflowTimesInterp = interp1(InflowCumulative(i,:),CurrentTime,1:nVeh);
%     DemTimesInterp = interp1(DemCumulative(i,:),CurrentTime,1:nVeh);
%     plot(DemTimesInterp,abs(DemTimesInterp-InflowTimesInterp))
% %     plot(InflowCumulative(i,:)/InflowCumulative(i,end),'color',cmap_perso(i,:))
% %     plot(DemCumulative(i,:)/DemCumulative(i,end),'--','color',cmap_perso(i,:))
% end
% aa = sum(InflowCumulative,1);bb = sum(DemCumulative,1);
% % plot(aa/aa(end),'color',cmap_perso(4,:))
% % plot(bb/bb(end),'--','color',cmap_perso(4,:))


% %% MFD surface
% 
% [kc,kb] = meshgrid(1:1000,1:50);
% Pc = zeros(size(kc));Pb = Pc;
% for iline = 1:size(kc,1)
%     k = [kc(iline,:);kb(iline,:)];
%     Temp_param = Reservoir(r).MFDfctParam;
%     [Temp_P] = MFDfct(k,Temp_param);
%     Pc(iline,:) = Temp_P(1,:);
%     Pb(iline,:) = Temp_P(2,:);
% end
% uc = Temp_param(1);ub  = Temp_param(2);bcc = Temp_param(3);bbc = Temp_param(4);
% bcb = Temp_param(5);bbb = Temp_param(6);
% b1  = 2*bcc; b2  = bbc + bcb; b3 = uc; % Critical acc. line for Global 3DMFD
% bc1 = 2*bcc; bc2 = bbc; bc3 = uc;      % Critical acc. line for Car 3DMFD
% bb1 = bcb;bb2 = 2*bbb;bb3 = ub;        % Critical acc. line for Bus 3DMFD
% ncline = [b1 b2 b3;...
%           bc1 bc2 bc3;...
%           bb1 bb2 bb3];
% for ires = 1:Nres
%     modeindexc = Reservoir(ires).ModeIndex{1};
%     nc = sum(Reservoir(ires).AccPerRoute(modeindexc,:),1);
%     modeindexb = Reservoir(ires).ModeIndex{2};
%     nb = sum(Reservoir(ires).AccPerRoute(modeindexb,:),1);
%     figure
%     subplot(1,2,1)
%     contour( kc, kb, Pc, 'ShowText','on' );
%     hold on
%     plot(-(ncline(2,3) + ncline(2,2)*[1:50])/ncline(2,1),[1:50],'k','LineWidth',LW)
%     plot(nc,nb,'linewidth',LW)
%     colormap('jet');colorbar
%     subplot(1,2,2)
%     contour( kc, kb, Pb, 'ShowText','on' );
%     hold on
%     plot(-(ncline(3,3) + ncline(3,2)*[1:50])/ncline(3,1),[1:50],'k','LineWidth',LW)
%     plot(nc,nb,'linewidth',LW)
%     colormap('jet');colorbar
%     set(gcf,'Position',[10 10 1400 2*300])
% end
