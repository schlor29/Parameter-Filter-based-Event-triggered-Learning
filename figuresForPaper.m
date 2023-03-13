%% ETL
clear all
close all
%%

load('.\figures\UpdateETL_004_tunedAndSmallEyeKFTuned\UpdateETL_004_tunedAndSmallEyeKFTuned.mat')
% load('.\figures\UpdateETL_004_tunedAndSmallEyeAlways\UpdateETL_004_tunedAndSmallEyeAlways.mat')
% load('.\figures\UpdateETL_004_tunedAndSmallNoUpdate\UpdateETL_004_tunedAndSmallNoUpdate.mat')


%% all triggering thresholds


% figure
figure('Renderer', 'painters', 'Position', [10 10 600 300])
hold on

%KFPartTuned
KFPartTunedchi2test = results.KFPartTunedchi2test.KFPartTunedchi2test(:);
p1 = plot(tMat,KFPartTunedchi2test,'LineWidth',2);
p2 = plot([tMat(1),tMat(end)],[1,1],'LineWidth',2);
for iSysChange = 1:numel(syst)
    try
        line([syst(iSysChange), syst(iSysChange)],get(gca,'YLim'),'Color',[0 0 0], 'Displayname', 'System change','LineWidth',2)
    catch
    end
end
for iLTrigger = 1:numel(lTriggerTimes)
    try
        line([lTriggerTimes(iLTrigger), lTriggerTimes(iLTrigger)],get(gca,'YLim'),'Color',[0 1 0],'LineWidth',2)
        f1 = fill([lTriggerTimes(iLTrigger), lTriggerTimes(iLTrigger), endExperimentTimes(iLTrigger), endExperimentTimes(iLTrigger)],[get(gca,'YLim'), fliplr(get(gca,'YLim'))] , 0.8*[0 1 0], 'EdgeColor','none','facealpha',.2);
    catch
    end
end
% title('KFPartTuned')

bT = KFPartTunedchi2test>1;
p3 = plot(tMat(bT), 1,'o','MarkerSize',10,...
    'MarkerEdgeColor',[1 0 0],'LineWidth',2);
%,...
%    'MarkerFaceColor',[0 1 0]

xlabel('t')
ylabel('test statistics')
legend([p1 p2 p3(1) f1],{'Test statistics','Threshold', 'Triger events', 'Experiment duration'})

ylim([0,1.2])
% xticks(0:10:50)

% linkaxes(sp,'x');
%%
cleanfigure('targetResolution',200,'scalePrecision',4)
matlab2tikz('standalone', true,'floatFormat','%.5g')


%% KFPartTuned
KFAllzCell = results.KFPartTunedz.KFPartTunedz(:)';
KFAllzMat = cell2mat(KFAllzCell);
KFAllPCell = results.KFPartTunedP.KFPartTunedP(:)';
KFAllt = results.KFPartTunedP.t(:)';

KFAllzModDiff = KFAllzMat - modzMatAllT;
zSysModDiff = syszMatAllT - modzMatAllT;
% %%%%% for always:
% KFAllzModDiff = KFAllzMat - KFAllzMat;
% zSysModDiff = syszMatAllT - KFAllzMat;
% %%%%% :for always

interval = sqrt( chi2inv(1-alpha,n*(n+m)) ./ cell2mat(cellfun(@(x) abs(diag(x^-1)),KFAllPCell,'UniformOutput',false)) );

triggerThresh = chi2inv(1-alpha,n*(n+m));

x = tMat;
ylims = [-0.6,0.6; -0.6,0.6; -1.5,1.5; -0.6,0.6];
% figure
figure('Renderer', 'painters', 'Position', [10 10 600 300])

iter = 0;
for iParam = [1,3,6,8]
    iter = iter+1;
%     subPlots(iParam) = subplot(n,n+m,iParam);
nexttile;
    hold on
%     %%%%% for always:
% plot([-1,0],[0,0])
% %%%%% :for always
    iParamModNM = iParam - (m+n).*floor(iParam./(m+n));
    if iParamModNM == 0
        iParamModNM = n+m;
    end
    fill([x,fliplr(x)], [KFAllzModDiff(iParam,:)+ interval(iParam,:),fliplr(KFAllzModDiff(iParam,:)- interval(iParam,:))], 0.8*[1 1 1], 'EdgeColor','none','facealpha',.5, 'Displayname', 'Projected erro covariance');
%     fill([x,fliplr(x)], [RLS.ThetaTildaVec(iParam,:)+ interval*cellfun(@(x) x(iParamModNM,iParamModNM),RLS.P),fliplr(RLS.ThetaTildaVec(iParam,:)- interval*cellfun(@(x) x(iParamModNM,iParamModNM),RLS.P))], 0.8*[1 1 1], 'EdgeColor','none','facealpha',.5, 'Displayname', ['variance param ' num2str(iParam)]);
    %     fill([x,fliplr(x)], [results.thetaTildaMean(iParam,:)+ interval*results.thetaTildaVar(iParam,:),fliplr(results.thetaTildaMean(iParam,:)- interval*results.thetaTildaVar(iParam,:))], 0.8*[1 1 1], 'facealpha',.25, 'Displayname', ['variance param ' num2str(iParam)]);
    
    h_ThetaTildaMean = plot(KFAllzModDiff(iParam,:)', 'DisplayName','Estimated Parameter error','LineWidth',1);
    h_trueError = plot(zSysModDiff(iParam,:)', 'DisplayName','Parameter error','LineWidth',1);
    ylim(ylims(iter,:));
    xlabel('t')
    grid on
    for iLTrigger = 1:numel(lTriggerTimes)
        try
            line([lTriggerTimes(iLTrigger), lTriggerTimes(iLTrigger)],get(gca,'YLim'),'Color',[0 1 0],'LineWidth',2)
            f1 = fill([lTriggerTimes(iLTrigger), lTriggerTimes(iLTrigger), endExperimentTimes(iLTrigger), endExperimentTimes(iLTrigger)],[get(gca,'YLim'), fliplr(get(gca,'YLim'))] , 0.8*[0 1 0], 'EdgeColor','none','facealpha',.2, 'displayname','Experiment duration');
        catch
        end
    end
    hold off
    legend off
    
cleanfigure('targetResolution',200,'scalePrecision',4)
end
legend([h_ThetaTildaMean, h_trueError, f1],'location', 'best')
% linkaxes(subPlots,'xy')
% sgtitle('KFTuned parameter estimate on partial data')
%%
% cleanfigure('targetResolution',10,'scalePrecision',4)
matlab2tikz('standalone', true,'floatFormat','%.5g')


%% Plot trace P
KFPartTunedPCell = results.KFPartTunedP.KFPartTunedP;
KFPartTunedPTrace = cellfun(@trace,KFPartTunedPCell)';
% KFAllPCell = results.KFAllP.KFAllP;
% KFAllPTrace = cellfun(@trace,KFAllPCell)';
figure('Renderer', 'painters', 'Position', [10 10 600 300])
hold on
plot(tMat,KFPartTunedPTrace)
% plot(tMat,KFAllPTrace)
ylim([0,3]);
for iSysChange = 1:numel(syst)
    try
        line([syst(iSysChange), syst(iSysChange)],get(gca,'YLim'),'Color',[0 0 0], 'LineWidth',2, 'Displayname', 'System change')
    catch
    end
end
% for iLTrigger = 1:numel(lTriggerTimes)
%     try
%         line([lTriggerTimes(iLTrigger), lTriggerTimes(iLTrigger)],get(gca,'YLim'),'Color',[0 1 0], 'LineWidth',2)
%         fill([lTriggerTimes(iLTrigger), lTriggerTimes(iLTrigger), endExperimentTimes(iLTrigger), endExperimentTimes(iLTrigger)],[get(gca,'YLim'), fliplr(get(gca,'YLim'))] , 0.8*[0 1 0], 'EdgeColor','none','facealpha',.2);
%     catch
%     end
% end
% title('Trace P')
xlabel('t')
ylabel('trace(P_{k|k})')

% ylim([0,3]);
%%
cleanfigure('targetResolution',200,'scalePrecision',4)
matlab2tikz('standalone', true,'floatFormat','%.5g')

%% plot motion of system
figure('Renderer', 'painters', 'Position', [10 10 600 600])
% nexttile
% plot(xVec(1,:),xVec(2,:))
T2D = [sys.kTheta  -sys.kTheta/sys.rho ];
TConstraint = 78.5398;
% TConstraint = 20;
XPoly2D = Polyhedron([T2D;-T2D], [TConstraint; TConstraint]);
XPoly2D = XPoly2D.intersect(Polyhedron.unitBox(2)*2);

nexttile
hold on
plot(XPoly2D)

xlim([-0.1 0.1])
ylim([-0.2 0.2])
plot(xVec(1,1:lTriggerTimes(1)),xVec(3,1:lTriggerTimes(1)))
title(['t \in [' num2str(1) ', ' num2str(lTriggerTimes(1)) ']'])
xlabel('x_1')
ylabel('x_3')
    
for iLTrigger = 1:numel(lTriggerTimes)
    try
        if iLTrigger == 2
            nexttile([2 1])
        else
            nexttile
        end
        
        hold on
        plot(XPoly2D)
        plot(xVec(1,lTriggerTimes(iLTrigger):endExperimentTimes(iLTrigger)),xVec(3,lTriggerTimes(iLTrigger):endExperimentTimes(iLTrigger)))
        if iLTrigger == 2
            xlim([-0.1 0.1])
            ylim([-0.55 0.85])
        else
            xlim([-0.1 0.1])
            ylim([-0.2 0.2])
        end
        title(['t \in [' num2str(lTriggerTimes(iLTrigger)) ', ' num2str(endExperimentTimes(iLTrigger)) ']'])
        
xlabel('x_1')
ylabel('x_3')
    catch
    end
    try
        nexttile
        hold on
        plot(XPoly2D)
xlim([-0.1 0.1])
ylim([-0.2 0.2])
        plot(xVec(1,endExperimentTimes(iLTrigger):lTriggerTimes(iLTrigger+1)),xVec(3,endExperimentTimes(iLTrigger):lTriggerTimes(iLTrigger+1)))
        title(['t \in [' num2str(endExperimentTimes(iLTrigger)) ', ' num2str(lTriggerTimes(iLTrigger+1)) ']'])
    
xlabel('x_1')
ylabel('x_3')
    catch
        delete(gca)
    end
end

nexttile
hold on
plot(XPoly2D)
xlim([-0.1 0.1])
ylim([-0.2 0.2])
plot(xVec(1,endExperimentTimes(end):end),xVec(3,endExperimentTimes(end):end))
title(['t \in [' num2str(endExperimentTimes(end)) ', ' num2str(tMat(end)) ']'])

xlabel('x_1')
ylabel('x_3')

% find all axes handle of type 'axes' and empty tag
% all_ha = findobj( gcf, 'type', 'axes', 'tag', '' );
% linkaxes( all_ha );

%%
% figure('Renderer', 'painters', 'Position', [10 10 600 600])
figure('Renderer', 'painters', 'Position', [10 10 300 300])
% nexttile
% plot(xVec(1,:),xVec(2,:))
T2D = [sys.kTheta  -sys.kTheta/sys.rho ];
TConstraint = 78.5398;
% TConstraint = 20;
XPoly2D = Polyhedron([T2D;-T2D], [TConstraint; TConstraint]);
XPoly2D = XPoly2D.intersect(Polyhedron.unitBox(2)*2);

nexttile
hold on
plot(XPoly2D)
plot(xVec(1,:),xVec(2,:))

xlim([-0.1 0.1])
ylim([-0.5 0.5])
title(['t \in [' num2str(1) ', ' num2str(tMat(end)) ']'])
xlabel('x_1')
ylabel('x_3')

%%
cleanfigure('targetResolution',200,'scalePrecision',4)
matlab2tikz('standalone', true,'floatFormat','%.5g')

%% Plot LQR stage cost
xCost = diag(xVec'*Q*xVec);
uCost = diag(uVec'*R*uVec);
cost = xCost + uCost;
figure
hold on
ylim([min(cost(2:end)),max(cost)])
plot(tMat,cost)
set(gca, 'YScale', 'log')
for iSysChange = 1:numel(syst)
    try
        line([syst(iSysChange), syst(iSysChange)],get(gca,'YLim'),'Color',[0 0 1], 'LineWidth',2, 'Displayname', 'System change')
    catch
    end
end
for iLTrigger = 1:numel(lTriggerTimes)
    try
        line([lTriggerTimes(iLTrigger), lTriggerTimes(iLTrigger)],get(gca,'YLim'),'Color',[1 0 0], 'LineWidth',2)
        fill([lTriggerTimes(iLTrigger), lTriggerTimes(iLTrigger), endExperimentTimes(iLTrigger), endExperimentTimes(iLTrigger)],[get(gca,'YLim'), fliplr(get(gca,'YLim'))] , 0.8*[0 1 0], 'EdgeColor','none','facealpha',.2);
    catch
    end
end
title('LQR stage cost')
averageCost = sum(cost)/numel(cost)



%% compare mse of model
clear all
load('.\figures\UpdateETL_004_tunedAndSmallEyeKFTuned\UpdateETL_004_tunedAndSmallEyeKFTuned.mat', 'results', 'modzMatAllT', 'syszMatAllT', 'lTriggerTimes', 'endExperimentTimes')
% load('.\figures\UpdateETL_004_tunedAndSmallEyeAlways\UpdateETL_004_tunedAndSmallEyeAlways.mat', 'results', 'modzMatAllT','syszMatAllT')
% load('.\figures\UpdateETL_004_tunedAndSmallNoUpdate\UpdateETL_004_tunedAndSmallNoUpdate.mat', 'results', 'modzMatAllT','syszMatAllT')

KFAllzCell = results.KFPartTunedz.KFPartTunedz(:)';
KFAllzMat = cell2mat(KFAllzCell);
KFAllPCell = results.KFPartTunedP.KFPartTunedP(:)';
KFAllt = results.KFPartTunedP.t(:)';

KFAllzModDiff = KFAllzMat - modzMatAllT;
zSysModDiff = syszMatAllT - modzMatAllT;
% %%%%% for always:
% KFAllzModDiff = KFAllzMat - KFAllzMat;
% zSysModDiff = syszMatAllT - KFAllzMat;
% %%%%% :for always

% figure
% plot(zSysModDiff')
mseSysMod_KFTuned = sum(sum(zSysModDiff'.^2))/3000/20;

experiments = sort([lTriggerTimes; endExperimentTimes]);
noExperiment = [1:experiments(1),experiments(2):experiments(3), experiments(4):3000];
zSysModDiff_exceptExperiment = zSysModDiff(:,noExperiment);
mseSysMod_KFTuned_exceptExperiment = sum(sum(zSysModDiff_exceptExperiment'.^2))/numel(zSysModDiff_exceptExperiment);

load('.\figures\UpdateETL_004_tunedAndSmallEyeAlways\UpdateETL_004_tunedAndSmallEyeAlways.mat', 'results', 'modzMatAllT','syszMatAllT')

KFAllzCell = results.KFPartTunedz.KFPartTunedz(:)';
KFAllzMat = cell2mat(KFAllzCell);
KFAllPCell = results.KFPartTunedP.KFPartTunedP(:)';
KFAllt = results.KFPartTunedP.t(:)';

% KFAllzModDiff = KFAllzMat - modzMatAllT;
% zSysModDiff = syszMatAllT - modzMatAllT;
%%%%% for always:
KFAllzModDiff = KFAllzMat - KFAllzMat;
zSysModDiff = syszMatAllT - KFAllzMat;
%%%%% :for always

% figure
% plot(zSysModDiff')
mseSysMod_Always = sum(sum(zSysModDiff'.^2))/3000/20;
zSysModDiff_exceptExperiment = zSysModDiff(:,noExperiment);
mseSysMod_Always_exceptExperiment = sum(sum(zSysModDiff_exceptExperiment'.^2))/numel(zSysModDiff_exceptExperiment);


load('.\figures\UpdateETL_004_tunedAndSmallNoUpdate\UpdateETL_004_tunedAndSmallNoUpdate.mat', 'results', 'modzMatAllT','syszMatAllT')

KFAllzCell = results.KFPartTunedz.KFPartTunedz(:)';
KFAllzMat = cell2mat(KFAllzCell);
KFAllPCell = results.KFPartTunedP.KFPartTunedP(:)';
KFAllt = results.KFPartTunedP.t(:)';

KFAllzModDiff = KFAllzMat - modzMatAllT;
zSysModDiff = syszMatAllT - modzMatAllT;
% %%%%% for always:
% KFAllzModDiff = KFAllzMat - KFAllzMat;
% zSysModDiff = syszMatAllT - KFAllzMat;
% %%%%% :for always

% figure
% plot(zSysModDiff')
mseSysMod_NoUpdate = sum(sum(zSysModDiff'.^2))/3000/20;
zSysModDiff_exceptExperiment = zSysModDiff(:,noExperiment);
mseSysMod_NoUpdate_exceptExperiment = sum(sum(zSysModDiff_exceptExperiment'.^2))/numel(zSysModDiff_exceptExperiment);


%% compare mse of model without Experiment
clear all
load('.\figures\UpdateETL_004_tunedAndSmallEyeKFTuned\UpdateETL_004_tunedAndSmallEyeKFTuned.mat', 'results', 'modzMatAllT', 'syszMatAllT', 'lTriggerTimes', 'endExperimentTimes')
% load('.\figures\UpdateETL_004_tunedAndSmallEyeAlways\UpdateETL_004_tunedAndSmallEyeAlways.mat', 'results', 'modzMatAllT','syszMatAllT')
% load('.\figures\UpdateETL_004_tunedAndSmallNoUpdate\UpdateETL_004_tunedAndSmallNoUpdate.mat', 'results', 'modzMatAllT','syszMatAllT')

KFAllzCell = results.KFPartTunedz.KFPartTunedz(:)';
KFAllzMat = cell2mat(KFAllzCell);
KFAllPCell = results.KFPartTunedP.KFPartTunedP(:)';
KFAllt = results.KFPartTunedP.t(:)';

KFAllzModDiff = KFAllzMat - modzMatAllT;
zSysModDiff = syszMatAllT - modzMatAllT;
% %%%%% for always:
% KFAllzModDiff = KFAllzMat - KFAllzMat;
% zSysModDiff = syszMatAllT - KFAllzMat;
% %%%%% :for always

% figure
% plot(zSysModDiff')
mseSysMod_KFTuned = sum(sum(zSysModDiff'.^2))/3000/20;

% lTriggerTimes(iLTrigger):endExperimentTimes(iLTrigger)
experiments = sort([lTriggerTimes; endExperimentTimes]);
noExperiment = [1:experiments(1),experiments(2):experiments(3), experiments(4):3000];
zSysModDiff_exceptExperiment = zSysModDiff(:,noExperiment);
mseSysMod_KFTuned_exceptExperiment = sum(sum(zSysModDiff_exceptExperiment'.^2))/numel(zSysModDiff_exceptExperiment);