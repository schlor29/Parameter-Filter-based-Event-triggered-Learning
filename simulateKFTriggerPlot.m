%%  implementation of event-triggered MPC
%   from Lehmann et al. (2013)
%
%   with assigned triggering sets and MPC control input
%
%   Sebastian Schlor
%   27.03.2020

%%
close all
clear all
clc

%% random generator
rng('default')

%% settings

% simulation length
T_sim = 3000;

TSysChange = 1000;


% Prediction horizon
N = 6;

%% system
% x_(t+1) = A_sys * x_t + B_sys * u_t + W * w_t

kindOfSystem = 'DCMotor';
% kindOfSystem = '2D1';
velocity = 2*pi*0.04;
switch kindOfSystem
    case 'DCMotor'
        
        % dimensions
        n = 4; % state
        m = 1; % input

        % % parameters
        sys.LS = 1.0;
        sys.dS = 0.02;
        sys.JM = 0.5;
        sys.betaM = 0.1;
        sys.R = 20;
        sys.KT = 10;
        sys.rho = 20;
        sys.kTheta = 1280.2;
        sys.JHatL = 20 * sys.JM;
        sys.betaL = 25;
        sys.Ts = 0.1;
        
        sys.factor = 20;
        sys.JL = sys.factor * sys.JM; % 20 * sys.JM: true system equals nominal model
        [Ad, Bd] = DCMotor(sys.JL);
        sys.A = Ad;
        sys.B = Bd;
        
%         sys.SigmaW = 0.01 * 0.005 * eye(n);
%         sys.SigmaW = 0.01 * eye(n);
        sys.SigmaW = 1* diag(1* [2*pi*1; velocity*1*sys.Ts; 2*pi*1*sys.rho; velocity*1*sys.rho*sys.Ts]/100) /16;
        
        %%%
        Q = diag([1e3, 1, 1e3, 1]);
        R = 10*eye(m);
        % Discrete-time LQR controller
        K = LQR(Q,R,sys);
        % LQR closed-loop dynamics
        A_K = sys.A+sys.B*K;
        %%%
        [V,D,W] = eig(sys.A);
% [V,D,W] = eig(A_K);
        wDirection = 0.01*real(W(1,:)); % actually not the left eigenvector, but produces nice disturance inside constraints 
%         wDirection = 0.05*real(W(:,4)); % actual right left eigenvector 
%         w = 0.01*wDirection' *mvnrnd(0,1,1); % noise in eigenvector direction 
% wDirection = ([3.9 0.016 78.5 0.31]* 1e-3).^.5; % as in our paper
        sys.SigmaW = diag((wDirection).^2);
        
        % disturbance
        wVec = zeros(n,T_sim);
        sysFactorVec = zeros(1,T_sim);
        for t = 1:T_sim
            w = mvnrnd(zeros(n,1),sys.SigmaW,1)'; % not noise in eigenvector direction 
            %     w = 0.0*mvnrnd(zeros(n,1),sys.SigmaW,1)' + 1*diag(sys.SigmaW);
%                 w = wDirection *mvnrnd(0,1,1); % eigenvector direction 
            wVec(:,t) = w;
            
            sysFactorVec(t) = 20+ 10 *2*(rand(1,1)-0.5);
        end
        sysFactorVec(1000) = 22;
        sysFactorVec(2000) = 19;
        
        %estimate system changes
%         [Al, Bl] = DCMotor(10* sys.JM);
%         [Am, Bm] = DCMotor(20* sys.JM);
%         [Au, Bu] = DCMotor(30* sys.JM);
%         ABDiffu = [Au, Bu] - [Am, Bm];
%         zDiffu = reshape(ABDiffu',[],1);
%         ABDiffl = [Al, Bl] - [Am, Bm];
%         zDiffl = reshape(ABDiffl',[],1);
%         zDiff = [zDiffu, zDiffl];
%         SigmaZDiff = zDiff*zDiff';
%         SigmaZTuned = 1/5000*SigmaZDiff;

        [Am, Bm] = DCMotor(20* sys.JM);
        zDiff = [];
        nPosibleChange = 2*n*(n+m);
        for iPosibleChange = 1:nPosibleChange
            factor = 20+ 10 *2*(rand(1,1)-0.5);
            [Au, Bu] = DCMotor(factor* sys.JM);
            ABDiffu = [Au, Bu] - [Am, Bm];
            zDiffu = reshape(ABDiffu',[],1);
            zDiff = [zDiff, zDiffu];
        end
        SigmaZDiff = zDiff*zDiff';
        SigmaZTuned = 1/(1000*nPosibleChange-1)*SigmaZDiff;
        SigmaZTuned = SigmaZTuned *10;
        SigmaZTuned = 0.000001*eye(n*(n+m)) + SigmaZTuned;
%         SigmaZTuned = 0.001*eye(n*(n+m)) + SigmaZTuned;
%         SigmaZTuned = 0.000001*ones(n*(n+m)) + SigmaZTuned;
        
%         SigmaZTuned = 0.001*eye(n*(n+m)) + 0.01*SigmaZDiff;
%         eig(SigmaZTuned)
        
        % constraints
        % State constraint set
        
        T = [sys.kTheta 0 -sys.kTheta/sys.rho 0];
        TConstraint = 78.5398;
        % TConstraint = 20;
        XPoly = Polyhedron([T;-T], [TConstraint; TConstraint]);
        XPoly = XPoly.intersect(Polyhedron.unitBox(n)*100);
        % XPoly = XPoly.intersect(Polyhedron.unitBox(n)*diag([10, 1, 150, 10]));
        
        % Input constraint set
        UPoly = 220 * Polyhedron.unitBox(m);
        % UPoly = 40000 * Polyhedron.unitBox(m);
        
    case '2D1'
        % dimensions
        n = 2; % state
        m = 1; % input
        
        sys.A = [1 1; 0 1];
        sys.B = [1;1];
        
        sys.SigmaW = 0.001*eye(n);
        
        SigmaZTuned = 0.01*eye(n*(n+m));
        
        % constraints
        % State constraint set
%         XPoly = 10* Polyhedron.unitBox(n);
        XPoly = diag([1,10])* Polyhedron.unitBox(n);
        
        % Input constraint set
        UPoly = 2 * Polyhedron.unitBox(m);
end

%% initial model
% x_m_(t+1) = A_mod * x_m_t + B_mod * u_t + w_m_t

mod = sys;





%% compute feasible x_ref

% % constant velocity
% % velocity = 2*pi*0.001; % sometimes infeasible
% % velocity = 2*pi*0.0005; %feasible
% % velocity = 2*pi*0.008;
% % velocity = 2*pi*0.01;
% % x_ref = @(t) 1* [velocity*t*sys.Ts; velocity*ones(1,numel(t)); velocity*t*sys.Ts*sys.rho; velocity*sys.rho*ones(1,numel(t))];
% 
% % sin
% velocity = 2*pi*0.04;
% x_ref = @(t) 1* [2*pi*sin(velocity/2/pi*t*sys.Ts); velocity*cos(velocity/2/pi*t*sys.Ts)*sys.Ts; 2*pi*sin(velocity/2/pi*t*sys.Ts)*sys.rho; velocity*cos(velocity/2/pi*t*sys.Ts)*sys.rho*sys.Ts];
% % x_ref = @(t) 1* [2*pi*sin(velocity/2/pi*t*sys.Ts); 0*velocity*cos(velocity/2/pi*t*sys.Ts)*sys.Ts; 0*2*pi*sin(velocity/2/pi*t*sys.Ts)*sys.rho; 0*velocity*cos(velocity/2/pi*t*sys.Ts)*sys.rho*sys.Ts];

% zeros reference
x_ref = @(t) zeros(n,numel(t));
% x_ref = @(t) 10*ones(n,numel(t));

% Initial condition
x1 = x_ref(1);



%% LQR
% Weighting matrices for LQR
% Q = eye(n);
% R = 10*eye(m);

Q = diag([1e3, 1, 1e3, 1]);
R = 10*eye(m);

% Discrete-time LQR controller
K = LQR(Q,R,mod);

% % LQR closed-loop dynamics
A_K = sys.A+sys.B*K;
eig(A_K)
eig(mod.A+mod.B*K)

%% Determine terminal cost
[QN,~,~] = dare(mod.A,mod.B,Q,R);

%% Initialize Controller
% TStart = 0.1 * Polyhedron.unitBox(n);
TStart = (chi2inv(0.9,4)* 16* diag(diag(mod.SigmaW)) + 0.1 * eye(n)) * Polyhedron.unitBox(n);
for t = 1:N
    TSets(t) = TStart;
end
TSets(N) = 0*TStart;


% triggerType = 'KFPart';
triggerType = 'KFPartTuned';
% triggerType = 'KFAll';
% triggerType = 'always';

%         estAlgo = 'LS';
%         estAlgo = 'KFAll';
%         estAlgo = 'KFPart';
        estAlgo = 'KFPartTuned';
%         estAlgo = 'true';

% setup.significanceLevel = 0.01;
% setup.doExperiment = true;
% setup.experimentDuration = 400;
setup.KFWeight = 1e6; % 100;
% setup.KFQ = 0.001*eye(n*(n+m));
setup.KFQ = 0.00001*eye(n*(n+m));
% alpha = 0.05;
alpha = 0.01;


%% simulation

% doExperiment = true;
doExperiment = false;
experimentDuration = 200;
% store results
results = resultsClass();
results.add(1,'sys',sys);
results.add(1,'mod',mod);

% initialize MPC
MPC = MPCClass(mod, N, Q, R, QN, XPoly, UPoly);

initZ = reshape([mod.A'; mod.B'],[],1);
P = 0.5*eye(numel(initZ));
KFAll = KalmanFilterSysEstClass(P, initZ, n, m, setup.KFQ, mod.SigmaW, 'KFAll', results,1);
% KFPart = KalmanFilterSysEstClass(P, initZ, n, m, setup.KFQ, mod.SigmaW, 'KFPart', results,1);
KFPartTuned = KalmanFilterSysEstClass(P, initZ, n, m, SigmaZTuned, mod.SigmaW, 'KFPartTuned', results,1);



tSinceLastTrigger = 0;
x = x1;
xData = [];
uData = [];
tauVec = [];

xPlan = zeros(n,N+1);
uPlan = zeros(m,N);
lTrigger = false;
modus = 'control';
identify = false;

d_tm1 = zeros(n+m,1);

% temporary different noise
% sys.SigmaW = 0.1*sys.SigmaW;

% traceAll= [];
% tracePartTuned= [];

for t = 1:T_sim
%     disp(t)
    
    if rem(t,TSysChange)==0
        switch kindOfSystem
            case 'DCMotor'
%                 sys.factor = 20+ 10 *2*(rand(1,1)-0.5);
                sys.factor = sysFactorVec(t);
                sys.JL = sys.factor * sys.JM; % 20 * sys.JM: true system equals nominal model
                [Ad, Bd] = DCMotor(sys.JL);
                sys.A = Ad;
                sys.B = Bd;
                
%                 SigmaWsqrt = sqrtm(sys.SigmaW) + 0.01* 2*(rand(n,n)-0.5);
%                 sys.SigmaW = SigmaWsqrt' * SigmaWsqrt;
                
                %             case '2D1'
            otherwise
                sys.A = sys.A + 0.1* 2*(rand(n,n)-0.5);
                sys.B = sys.B + 0.1* 2*(rand(n,m)-0.5);
                
                SigmaWsqrt = sqrtm(sys.SigmaW) + 0.01* 2*(rand(n,n)-0.5);
                sys.SigmaW = SigmaWsqrt' * SigmaWsqrt;
        end

        results.add(t,'sys',sys);
    end
    
    %% sensor side
    % measure x
    x_m = x;
    xData = [xData, x_m];
    if size(xData,2)>experimentDuration+1
        xData(:,1) = [];
    end
    
    % learning trigger
    t_KFAll = KFAll.update(x_m, d_tm1, mod, alpha);
    results.add(t,'t_KFAll',t_KFAll);
%     t_Res = Residual.update(xData, uData, mod, alpha);
%     results.add(t,'t_Res',t_Res);
    
    % evaluate event-trigger
    tSinceLastTrigger = tSinceLastTrigger + 1;
%     bTrigger = ~TSets(tSinceLastTrigger).contains(x_m - xPlan(:,tSinceLastTrigger+1)); % event-triggered
    bTrigger = true; % time-triggered
%     results.add(t,'bTrigger',bTrigger);
    
    if bTrigger
        disp(t)
        tau = tSinceLastTrigger;
        tauVec = [tauVec; tau];
        if numel(tauVec)>experimentDuration
            tauVec(1) = [];
        end
        results.add(t,'tau',tau);
        tSinceLastTrigger = 0;
        %% controller side
        
        % learning trigger
%         t_KFPart = KFPart.update(x_m, d_tm1, mod, alpha);
        t_KFPartTuned = KFPartTuned.update(x_m, d_tm1, mod, alpha);

        
    else
        % learning trigger without new data
%         t_KFPart = KFPart.update([], [], mod, alpha);
        t_KFPartTuned = KFPartTuned.update([], [], mod, alpha);

    end 
    
    
%     results.add(t,'t_KFPart',t_KFPart);
    results.add(t,'t_KFPartTuned',t_KFPartTuned);
%     traceAll(end+1) = trace(KFAll.P);
%     tracePartTuned(end+1) = trace(KFPartTuned.P);
    
    % learning trigger
    switch triggerType
%         case 'KFPart'
%             lTrigger = t_KFPart;
%             % wait until trigger 
%             lTrigger = lTrigger && size(xData,2)>=experimentDuration;
        case 'KFPartTuned'
            lTrigger = t_KFPartTuned;
            % wait until trigger 
            lTrigger = lTrigger && size(xData,2)>=experimentDuration;
        case 'KFAll'
            lTrigger = t_KFAll;
            % wait until trigger 
            lTrigger = lTrigger && size(xData,2)>=experimentDuration;
        case 'always'
            lTrigger = true;
            doExperiment = false;
    end
    
%     % wait until trigger 
%     lTrigger = lTrigger && size(xData,2)>=experimentDuration;
%     lTrigger = false;
%     lTrigger = (t>=200);
    
    switch modus
        case 'control'
            if lTrigger
                if ~isequal(triggerType,'always')
                    disp('trigger')
                    results.add(t,'lTrigger',1);
                end
                
                if doExperiment
                    modus = 'experiment';
                    experimentTimeCounter = 0;
                else
                    identify = true;
                end
                
            end
            
        case 'experiment'
            experimentTimeCounter = experimentTimeCounter + 1;
            
            if experimentTimeCounter == experimentDuration || any(x_m>100)
                identify = true;
            end
    end
    
    
    if false %identify
        % do identification based on data

        
        % learn new model
        switch estAlgo
            case 'LS'
                [A,B,w] = estSys(xData, [uData,zeros(m,1)]);
                mod.SigmaW = cov(w');
            case 'KFAll'
                [A, B] = KFAll.getAB;
                w = getw(xData, [uData,zeros(m,1)],A,B);
                mod.SigmaW = cov(w');
%             case 'KFPart'
%                 [A, B] = KFPart.getAB;
%                 w = getw(xData, [uData,zeros(m,1)],A,B);
%                 mod.SigmaW = cov(w');
            case 'KFPartTuned'
                [A, B] = KFPartTuned.getAB;
                w = getw(xData, [uData,zeros(m,1)],A,B);
                if ~any(any(isnan(cov(w'))))
                    mod.SigmaW = cov(w');
                end
            case 'true'
                A = sys.A;
                B = sys.B;
                mod.SigmaW = sys.SigmaW;
        end
        
        if true %~any(any(isnan([A,B])))
%             keyboard
            mod.A = A;
            mod.B = B;
            if ~isequal(triggerType,'always')
                results.add(t,'mod',mod);
            end
        
        
        % initialize controller new
        MPC = MPCClass(mod, N, Q, R, QN, XPoly, UPoly);
        KFAll = KalmanFilterSysEstClass(KFAll.P, KFAll.z, n, m, setup.KFQ, mod.SigmaW, 'KFAll', results,t);
%         KFPart = KalmanFilterSysEstClass(KFAll.P, KFAll.z, n, m, setup.KFQ, mod.SigmaW, 'KFPart', results,t);
%         KFPart = KalmanFilterSysEstClass(KFPart.P, KFPart.z, n, m, setup.KFQ, mod.SigmaW, 'KFPart', results,t);
        KFPartTuned = KalmanFilterSysEstClass(KFPartTuned.P, KFPartTuned.z, n, m, SigmaZTuned, mod.SigmaW, 'KFPartTuned', results,t);
        
        end
        % reset
        identify = false;
        xData = [];
        uData = [];
        tauVec = [];
        modus = 'control';
    end
    
    
    %% control
    if bTrigger
        %% controller side
        xRef = x_ref(t:t+N);
        % control modus
        switch modus
            case 'control'
                solutions = MPC.optimize(x_m,xRef);
%                 solutions = MPC.optimizeWithP(x_m, xRef, KFPartTuned.P, KFPartTuned.z, 0, SigmaZTuned, mod.SigmaW);
            case 'experiment'
                %                 solutions = MPC.optimize(x_m,xRef);
                switch estAlgo
                    case 'KFAll'
                        solutions = MPC.optimizeWithP(x_m, xRef, KFAll.P, KFAll.z, setup.KFWeight, setup.KFQ, mod.SigmaW);
%                     case 'KFPart'
%                         solutions = MPC.optimizeWithP(x_m, xRef, KFPart.P, KFPart.z, setup.KFWeight, setup.KFQ, mod.SigmaW);
                    case 'KFPartTuned'
                        solutions = MPC.optimizeWithP(x_m, xRef, KFPartTuned.P, KFPartTuned.z, setup.KFWeight, SigmaZTuned, mod.SigmaW);
                    otherwise
                        solutions = MPC.optimizeWithP(x_m, xRef, KFAll.P, KFAll.z, setup.KFWeight, setup.KFQ, mod.SigmaW);
                end
                
        end
        uPlan = solutions{1};
        xPlan = solutions{2};
    else
        % use predicted u
    end 
    u = uPlan(tSinceLastTrigger+1);
    
%     if t==100
%         u=230; % trace(P) drops instantly after high inputs. 
%     end
    
    uData = [uData, u];
    if size(uData,2)>experimentDuration
        uData(:,1) = [];
    end
        
    d_tm1 = [x;u];
    results.add(t,'d',d_tm1);
    
    
    %% system update
%     % disturbance
%     w = mvnrnd(zeros(n,1),sys.SigmaW,1)';
% %     w = 0.0*mvnrnd(zeros(n,1),sys.SigmaW,1)' + 1*diag(sys.SigmaW);
% %     [V,D,W] = eig(sys.A);
% %     wDirection = real(W(1,:));
% %     w = wDirection *mvnrnd(0,1,1);
    w = wVec(:,t);
    
    results.add(t,'w',w);
    x = sys.A * x + sys.B * u + w;
    
end


%% plots

%% time plot
tMat = results.d.t(:)';
T_sim = tMat(end);
dCell = results.d.d(:);
dMat = cell2mat(dCell');
xVec = dMat(1:n,:);
uVec = dMat(n+1:end,:);

% wCell = results.w.w(:);
% wMat = cell2mat(wCell');

syst = results.sys.t(:)';

try
    lTriggerTimes = results.lTrigger.t;
catch
    lTriggerTimes = [];
end
% lTriggerTimes = [];
endExperimentTimes = results.mod.t(2:end);
mAvlTriggerSections =  sort([lTriggerTimes; endExperimentTimes]);
triggerSections = [1; mAvlTriggerSections; T_sim+1];
% bTrigger = results.bTrigger.bTrigger';
% mAv = [];
% for iSection = 1:numel(triggerSections)-1
%     triggerInSection = bTrigger(triggerSections(iSection):triggerSections(iSection+1)-1);
%     mAv = [mAv, cumsum(triggerInSection)./[1:numel(triggerInSection)]];
% end
% mAv(triggerSections(2:end)-1)
% 1./mAv(triggerSections(2:end)-1)


% different learning triggers
t_KFAllVec = results.t_KFAll.t_KFAll;
% t_KFPartVec = results.t_KFPart.t_KFPart;
t_KFPartTunedVec = results.t_KFPartTuned.t_KFPartTuned;





% figure
% hold on
% plot(tMat,xVec, 'DisplayName','x')
% plot(tMat,uVec, 'DisplayName','u')
% % plot(tMat,wMat, 'DisplayName','w')
% legend

%% plot

tVec = tMat;

% close all
figure
% sfh1 = subplot(2,1,1);
% offset = 0.1;
% sfh1.Position = sfh1.Position + [0 -offset 0 offset];
hold on
% h_stateTrigger = plot(results.trigger + 10, 'DisplayName','state trigger');
% h_stateTriggerAverage = plot(cumsum(results.trigger)/[1:numel(results.trigger)]);
h_u = plot(tVec,uVec, 'DisplayName','u');
h_x = plot(tVec,xVec, 'DisplayName','x');
xrefVec = x_ref(tVec);
h_xref= plot(tVec, xrefVec, 'DisplayName','x reference');
% h_w = plot(tVec,wMat, 'DisplayName','w');

xlabel('t')

% learning trigger
% moving average
for iSysChange = 1:numel(syst)
    try
        line([syst(iSysChange), syst(iSysChange)],get(gca,'YLim'),'Color',[0 0 1], 'Displayname', 'System change')
    catch
    end
end

for iLTrigger = 1:numel(lTriggerTimes)
    try
        line([lTriggerTimes(iLTrigger), lTriggerTimes(iLTrigger)],get(gca,'YLim'),'Color',[1 0 0], 'Displayname', 'Learning trigger')
        fill([lTriggerTimes(iLTrigger), lTriggerTimes(iLTrigger), endExperimentTimes(iLTrigger), endExperimentTimes(iLTrigger)],[get(gca,'YLim'), fliplr(get(gca,'YLim'))] , 0.8*[0 1 0], 'EdgeColor','none','facealpha',.2);
    catch
    end
end
h = findobj(gca,'Type','line');

% for t = 1:length(results.triggerSets)
%     if ~isempty(results.triggerSets(t).V')
%         plot([t,t], results.triggerSets(t).V','r')
%     end
% end
% legend([h_stateTrigger, h_u, h_x, h_w, h_hoeffdingsTrigger, h_KSTrigger])
legend([h_u; h_x; h(1); h_xref], 'Location', 'best')



% sfh2 = subplot(2,1,2);
% sfh2.Position = sfh2.Position + [0 0 0 -offset];
% hold on
% h_stateTrigger = stem(tVec, bTrigger, 'DisplayName','state trigger', 'Marker', 'none');%plot(results.trigger, 'DisplayName','state trigger');
% h_stateTriggerAverage = plot(tVec, mAv, 'DisplayName','moving average of inter-event-times');
% 
% for iSysChange = 1:numel(syst)
%     try
%         line([syst(iSysChange), syst(iSysChange)],get(gca,'YLim'),'Color',[0 0 1], 'Displayname', 'System change')
%     catch
%     end
% end
% for iLTrigger = 1:numel(lTriggerTimes)
%     try
%         line([lTriggerTimes(iLTrigger), lTriggerTimes(iLTrigger)],get(gca,'YLim'),'Color',[1 0 0])
%         fill([lTriggerTimes(iLTrigger), lTriggerTimes(iLTrigger), endExperimentTimes(iLTrigger), endExperimentTimes(iLTrigger)],[get(gca,'YLim'), fliplr(get(gca,'YLim'))] , 0.8*[0 1 0], 'EdgeColor','none','facealpha',.2);
%     catch
%     end
% end
% legend([h_stateTrigger,h_stateTriggerAverage], 'Location', 'best')
% xlabel('t')
% ylabel('trigger')
% 
% linkaxes([sfh1,sfh2],'x');
% xlim(sfh2, [0 T_sim])



%%
% numberSubplots = 2;
% figure
% sp(1) = subplot(numberSubplots,1,1);
% hold on
% stem(tMat,t_KFAllVec)
% for iSysChange = 1:numel(syst)
%     try
%         line([syst(iSysChange), syst(iSysChange)],get(gca,'YLim'),'Color',[0 0 1], 'Displayname', 'System change')
%     catch
%     end
% end
% for iLTrigger = 1:numel(lTriggerTimes)
%     try
%         line([lTriggerTimes(iLTrigger), lTriggerTimes(iLTrigger)],get(gca,'YLim'),'Color',[1 0 0])
%         fill([lTriggerTimes(iLTrigger), lTriggerTimes(iLTrigger), endExperimentTimes(iLTrigger), endExperimentTimes(iLTrigger)],[get(gca,'YLim'), fliplr(get(gca,'YLim'))] , 0.8*[0 1 0], 'EdgeColor','none','facealpha',.2);
%     catch
%     end
% end
% title('KFAll')
% sp(2) =subplot(numberSubplots,1,2);
% hold on
% stem(tMat,t_KFPartVec)
% for iSysChange = 1:numel(syst)
%     try
%         line([syst(iSysChange), syst(iSysChange)],get(gca,'YLim'),'Color',[0 0 1], 'Displayname', 'System change')
%     catch
%     end
% end
% for iLTrigger = 1:numel(lTriggerTimes)
%     try
%         line([lTriggerTimes(iLTrigger), lTriggerTimes(iLTrigger)],get(gca,'YLim'),'Color',[1 0 0])
%         fill([lTriggerTimes(iLTrigger), lTriggerTimes(iLTrigger), endExperimentTimes(iLTrigger), endExperimentTimes(iLTrigger)],[get(gca,'YLim'), fliplr(get(gca,'YLim'))] , 0.8*[0 1 0], 'EdgeColor','none','facealpha',.2);
%     catch
%     end
% end
% 
% title('KFPart')
% 
% 
% linkaxes(sp,'x');

%% KF estimate

% system
sysCell = results.sys.sys(:)';
syst = results.sys.t(:)';
sysACell = {sysCell.A};
sysBCell = {sysCell.B};
syszCell = cellfun(@(A,B) reshape([A,B]',[],1),sysACell,sysBCell,'UniformOutput',false);
syszMat = cell2mat(syszCell);
syszMatAllT = [];
for iSys = 1:numel(syst)-1
    syszMatAllT = [syszMatAllT, repmat(syszMat(:,iSys), 1,syst(iSys+1)-syst(iSys))];
end
syszMatAllT = [syszMatAllT, repmat(syszMat(:,end), 1,tMat(end)+1 - syst(end))];

% model
modCell = results.mod.mod(:)';
modt = results.mod.t(:)';
modACell = {modCell.A};
modBCell = {modCell.B};
modzCell = cellfun(@(A,B) reshape([A,B]',[],1),modACell,modBCell,'UniformOutput',false);
modzMat = cell2mat(modzCell);
modzMatAllT = [];
for iSys = 1:numel(modt)-1
    modzMatAllT = [modzMatAllT, repmat(modzMat(:,iSys), 1,modt(iSys+1)-modt(iSys))];
end
modzMatAllT = [modzMatAllT, repmat(modzMat(:,end), 1,tMat(end)+1 - modt(end))];


%% KFAll
KFAllzCell = results.KFAllz.KFAllz(:)';
KFAllzMat = cell2mat(KFAllzCell);
KFAllPCell = results.KFAllP.KFAllP(:)';
KFAllt = results.KFAllP.t(:)';

KFAllzModDiff = KFAllzMat - modzMatAllT;
zSysModDiff = syszMatAllT - modzMatAllT;

interval = sqrt( chi2inv(1-alpha,n*(n+m)) ./ cell2mat(cellfun(@(x) abs(diag(x^-1)),KFAllPCell,'UniformOutput',false)) );

x = tMat;
figure
for iParam = 1:n*(n+m)
    subPlots(iParam) = subplot(n,n+m,iParam);
    hold on
    iParamModNM = iParam - (m+n).*floor(iParam./(m+n));
    if iParamModNM == 0
        iParamModNM = n+m;
    end
    fill([x,fliplr(x)], [KFAllzModDiff(iParam,:)+ interval(iParam,:),fliplr(KFAllzModDiff(iParam,:)- interval(iParam,:))], 0.8*[1 1 1], 'EdgeColor','none','facealpha',.5, 'Displayname', ['variance param ' num2str(iParam)]);
%     fill([x,fliplr(x)], [RLS.ThetaTildaVec(iParam,:)+ interval*cellfun(@(x) x(iParamModNM,iParamModNM),RLS.P),fliplr(RLS.ThetaTildaVec(iParam,:)- interval*cellfun(@(x) x(iParamModNM,iParamModNM),RLS.P))], 0.8*[1 1 1], 'EdgeColor','none','facealpha',.5, 'Displayname', ['variance param ' num2str(iParam)]);
    %     fill([x,fliplr(x)], [results.thetaTildaMean(iParam,:)+ interval*results.thetaTildaVar(iParam,:),fliplr(results.thetaTildaMean(iParam,:)- interval*results.thetaTildaVar(iParam,:))], 0.8*[1 1 1], 'facealpha',.25, 'Displayname', ['variance param ' num2str(iParam)]);
    
    h_ThetaTildaMean = plot(KFAllzModDiff(iParam,:)', 'DisplayName','ThetaTildaMean');
    h_trueError = plot(zSysModDiff(iParam,:)', 'DisplayName','trueTheta');
    ylim([-1,1]);
    xlabel('t')
    grid on
    for iLTrigger = 1:numel(lTriggerTimes)
        try
            line([lTriggerTimes(iLTrigger), lTriggerTimes(iLTrigger)],get(gca,'YLim'),'Color',[1 0 0])
            fill([lTriggerTimes(iLTrigger), lTriggerTimes(iLTrigger), endExperimentTimes(iLTrigger), endExperimentTimes(iLTrigger)],[get(gca,'YLim'), fliplr(get(gca,'YLim'))] , 0.8*[0 1 0], 'EdgeColor','none','facealpha',.2);
        catch
        end
    end
end
legend('location', 'best')
linkaxes(subPlots,'xy')
sgtitle('KF parameter estimate on all data')

%% KFPart
% KFAllzCell = results.KFPartz.KFPartz(:)';
% KFAllzMat = cell2mat(KFAllzCell);
% KFAllPCell = results.KFPartP.KFPartP(:)';
% KFAllt = results.KFPartP.t(:)';
% 
% KFAllzModDiff = KFAllzMat - modzMatAllT;
% % zSysModDiff = syszMatAllT - modzMatAllT;
% 
% interval = sqrt( chi2inv(1-alpha,n*(n+m)) ./ cell2mat(cellfun(@(x) abs(diag(x^-1)),KFAllPCell,'UniformOutput',false)) );
% 
% triggerThresh = chi2inv(1-alpha,n*(n+m));
% 
% x = tMat;
% figure
% for iParam = 1:n*(n+m)
%     subPlots(iParam) = subplot(n,n+m,iParam);
%     hold on
%     iParamModNM = iParam - (m+n).*floor(iParam./(m+n));
%     if iParamModNM == 0
%         iParamModNM = n+m;
%     end
%     fill([x,fliplr(x)], [KFAllzModDiff(iParam,:)+ interval(iParam,:),fliplr(KFAllzModDiff(iParam,:)- interval(iParam,:))], 0.8*[1 1 1], 'EdgeColor','none','facealpha',.5, 'Displayname', ['variance param ' num2str(iParam)]);
% %     fill([x,fliplr(x)], [RLS.ThetaTildaVec(iParam,:)+ interval*cellfun(@(x) x(iParamModNM,iParamModNM),RLS.P),fliplr(RLS.ThetaTildaVec(iParam,:)- interval*cellfun(@(x) x(iParamModNM,iParamModNM),RLS.P))], 0.8*[1 1 1], 'EdgeColor','none','facealpha',.5, 'Displayname', ['variance param ' num2str(iParam)]);
%     %     fill([x,fliplr(x)], [results.thetaTildaMean(iParam,:)+ interval*results.thetaTildaVar(iParam,:),fliplr(results.thetaTildaMean(iParam,:)- interval*results.thetaTildaVar(iParam,:))], 0.8*[1 1 1], 'facealpha',.25, 'Displayname', ['variance param ' num2str(iParam)]);
%     
%     h_ThetaTildaMean = plot(KFAllzModDiff(iParam,:)', 'DisplayName','ThetaTildaMean');
%     h_trueError = plot(zSysModDiff(iParam,:)', 'DisplayName','trueTheta');
%     ylim([-1,1]);
%     xlabel('t')
%     grid on
%     for iLTrigger = 1:numel(lTriggerTimes)
%         try
%             line([lTriggerTimes(iLTrigger), lTriggerTimes(iLTrigger)],get(gca,'YLim'),'Color',[1 0 0])
%             fill([lTriggerTimes(iLTrigger), lTriggerTimes(iLTrigger), endExperimentTimes(iLTrigger), endExperimentTimes(iLTrigger)],[get(gca,'YLim'), fliplr(get(gca,'YLim'))] , 0.8*[0 1 0], 'EdgeColor','none','facealpha',.2);
%         catch
%         end
%     end
% end
% legend('location', 'best')
% linkaxes(subPlots,'xy')
% sgtitle('KF parameter estimate on partial data')

%% KFPartTuned
KFAllzCell = results.KFPartTunedz.KFPartTunedz(:)';
KFAllzMat = cell2mat(KFAllzCell);
KFAllPCell = results.KFPartTunedP.KFPartTunedP(:)';
KFAllt = results.KFPartTunedP.t(:)';

KFAllzModDiff = KFAllzMat - modzMatAllT;
% zSysModDiff = syszMatAllT - modzMatAllT;

interval = sqrt( chi2inv(1-alpha,n*(n+m)) ./ cell2mat(cellfun(@(x) abs(diag(x^-1)),KFAllPCell,'UniformOutput',false)) );

triggerThresh = chi2inv(1-alpha,n*(n+m));

x = tMat;
figure
for iParam = 1:n*(n+m)
    subPlots(iParam) = subplot(n,n+m,iParam);
    hold on
    iParamModNM = iParam - (m+n).*floor(iParam./(m+n));
    if iParamModNM == 0
        iParamModNM = n+m;
    end
    fill([x,fliplr(x)], [KFAllzModDiff(iParam,:)+ interval(iParam,:),fliplr(KFAllzModDiff(iParam,:)- interval(iParam,:))], 0.8*[1 1 1], 'EdgeColor','none','facealpha',.5, 'Displayname', ['variance param ' num2str(iParam)]);
%     fill([x,fliplr(x)], [RLS.ThetaTildaVec(iParam,:)+ interval*cellfun(@(x) x(iParamModNM,iParamModNM),RLS.P),fliplr(RLS.ThetaTildaVec(iParam,:)- interval*cellfun(@(x) x(iParamModNM,iParamModNM),RLS.P))], 0.8*[1 1 1], 'EdgeColor','none','facealpha',.5, 'Displayname', ['variance param ' num2str(iParam)]);
    %     fill([x,fliplr(x)], [results.thetaTildaMean(iParam,:)+ interval*results.thetaTildaVar(iParam,:),fliplr(results.thetaTildaMean(iParam,:)- interval*results.thetaTildaVar(iParam,:))], 0.8*[1 1 1], 'facealpha',.25, 'Displayname', ['variance param ' num2str(iParam)]);
    
    h_ThetaTildaMean = plot(KFAllzModDiff(iParam,:)', 'DisplayName','ThetaTildaMean');
    h_trueError = plot(zSysModDiff(iParam,:)', 'DisplayName','trueTheta');
    ylim([-1,1]);
    xlabel('t')
    grid on
    for iLTrigger = 1:numel(lTriggerTimes)
        try
            line([lTriggerTimes(iLTrigger), lTriggerTimes(iLTrigger)],get(gca,'YLim'),'Color',[1 0 0])
            fill([lTriggerTimes(iLTrigger), lTriggerTimes(iLTrigger), endExperimentTimes(iLTrigger), endExperimentTimes(iLTrigger)],[get(gca,'YLim'), fliplr(get(gca,'YLim'))] , 0.8*[0 1 0], 'EdgeColor','none','facealpha',.2);
        catch
        end
    end
end
legend('location', 'best')
linkaxes(subPlots,'xy')
sgtitle('KFTuned parameter estimate on partial data')

%% factor in 4D
try
    factorvec = [results.sys.sys.factor]
catch
end

%% all triggering thresholds

numberSubplots = 2;
figure
sp(1) = subplot(numberSubplots,1,1);
hold on
% KFAll
KFAllchi2test = results.KFAllchi2test.KFAllchi2test(:);
plot(tMat,KFAllchi2test)
plot([tMat(1),tMat(end)],[1,1])
for iSysChange = 1:numel(syst)
    try
        line([syst(iSysChange), syst(iSysChange)],get(gca,'YLim'),'Color',[0 0 1], 'Displayname', 'System change')
    catch
    end
end
for iLTrigger = 1:numel(lTriggerTimes)
    try
        line([lTriggerTimes(iLTrigger), lTriggerTimes(iLTrigger)],get(gca,'YLim'),'Color',[1 0 0])
        fill([lTriggerTimes(iLTrigger), lTriggerTimes(iLTrigger), endExperimentTimes(iLTrigger), endExperimentTimes(iLTrigger)],[get(gca,'YLim'), fliplr(get(gca,'YLim'))] , 0.8*[0 1 0], 'EdgeColor','none','facealpha',.2);
    catch
    end
end
title('KFAll')
sp(2) =subplot(numberSubplots,1,2);
hold on
% %KFPart
% KFPartchi2test = results.KFPartchi2test.KFPartchi2test(:);
% plot(tMat,KFPartchi2test)
% plot([tMat(1),tMat(end)],[1,1])
% for iSysChange = 1:numel(syst)
%     try
%         line([syst(iSysChange), syst(iSysChange)],get(gca,'YLim'),'Color',[0 0 1], 'Displayname', 'System change')
%     catch
%     end
% end
% for iLTrigger = 1:numel(lTriggerTimes)
%     try
%         line([lTriggerTimes(iLTrigger), lTriggerTimes(iLTrigger)],get(gca,'YLim'),'Color',[1 0 0])
%         fill([lTriggerTimes(iLTrigger), lTriggerTimes(iLTrigger), endExperimentTimes(iLTrigger), endExperimentTimes(iLTrigger)],[get(gca,'YLim'), fliplr(get(gca,'YLim'))] , 0.8*[0 1 0], 'EdgeColor','none','facealpha',.2);
%     catch
%     end
% end
% title('KFPart')

% sp(3) =subplot(numberSubplots,1,3);
% hold on
%KFPartTuned
KFPartTunedchi2test = results.KFPartTunedchi2test.KFPartTunedchi2test(:);
plot(tMat,KFPartTunedchi2test)
plot([tMat(1),tMat(end)],[1,1])
for iSysChange = 1:numel(syst)
    try
        line([syst(iSysChange), syst(iSysChange)],get(gca,'YLim'),'Color',[0 0 1], 'Displayname', 'System change')
    catch
    end
end
for iLTrigger = 1:numel(lTriggerTimes)
    try
        line([lTriggerTimes(iLTrigger), lTriggerTimes(iLTrigger)],get(gca,'YLim'),'Color',[1 0 0])
        fill([lTriggerTimes(iLTrigger), lTriggerTimes(iLTrigger), endExperimentTimes(iLTrigger), endExperimentTimes(iLTrigger)],[get(gca,'YLim'), fliplr(get(gca,'YLim'))] , 0.8*[0 1 0], 'EdgeColor','none','facealpha',.2);
    catch
    end
end
title('KFPartTuned')

linkaxes(sp,'x');


%% Compare KFAll, KFPart, KFPartTuned
numberSubplots = 2;
figure
sp(1) = subplot(numberSubplots,1,1);
hold on
stem(tMat,t_KFAllVec)
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
title('KFAll trigger')
sp(2) =subplot(numberSubplots,1,2);
hold on
% stem(tMat,t_KFPartVec)
% for iSysChange = 1:numel(syst)
%     try
%         line([syst(iSysChange), syst(iSysChange)],get(gca,'YLim'),'Color',[0 0 1], 'LineWidth',2, 'Displayname', 'System change')
%     catch
%     end
% end
% for iLTrigger = 1:numel(lTriggerTimes)
%     try
%         line([lTriggerTimes(iLTrigger), lTriggerTimes(iLTrigger)],get(gca,'YLim'),'Color',[1 0 0], 'LineWidth',2)
%         fill([lTriggerTimes(iLTrigger), lTriggerTimes(iLTrigger), endExperimentTimes(iLTrigger), endExperimentTimes(iLTrigger)],[get(gca,'YLim'), fliplr(get(gca,'YLim'))] , 0.8*[0 1 0], 'EdgeColor','none','facealpha',.2);
%     catch
%     end
% end
% title('t_KFPart trigger')
% sp(3) =subplot(numberSubplots,1,3);
% hold on
stem(tMat,t_KFPartTunedVec)
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
title('t_KFPartTuned trigger')

linkaxes(sp,'x');

%% plot motion of system
figure
% nexttile
% plot(xVec(1,:),xVec(2,:))
T2D = [sys.kTheta  -sys.kTheta/sys.rho ];
TConstraint = 78.5398;
% TConstraint = 20;
XPoly2D = Polyhedron([T2D;-T2D], [TConstraint; TConstraint]);
XPoly2D = XPoly2D.intersect(Polyhedron.unitBox(2)*100);

nexttile
hold on
plot(XPoly2D)
xlim([-0.1 0.1])
ylim([-0.2 0.2])
plot(xVec(1,1:lTriggerTimes(1)),xVec(3,1:lTriggerTimes(1)))
for iLTrigger = 1:numel(lTriggerTimes)
    try
        nexttile
        hold on
        plot(XPoly2D)
xlim([-0.1 0.1])
ylim([-0.2 0.2])
        plot(xVec(1,lTriggerTimes(iLTrigger):endExperimentTimes(iLTrigger)),xVec(3,lTriggerTimes(iLTrigger):endExperimentTimes(iLTrigger)))
    catch
    end
    try
        nexttile
        hold on
        plot(XPoly2D)
xlim([-0.1 0.1])
ylim([-0.2 0.2])
        plot(xVec(1,endExperimentTimes(iLTrigger):lTriggerTimes(iLTrigger+1)),xVec(3,endExperimentTimes(iLTrigger):lTriggerTimes(iLTrigger+1)))
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

% find all axes handle of type 'axes' and empty tag
all_ha = findobj( gcf, 'type', 'axes', 'tag', '' );
linkaxes( all_ha );

%% Plot trace invP
KFPartTunedPCell = results.KFPartTunedP.KFPartTunedP;
KFPartTunedPInv = cellfun(@inv,KFPartTunedPCell,'UniformOutput',false);
KFPartTunedPInvTrace = cellfun(@trace,KFPartTunedPInv)';
KFAllPCell = results.KFAllP.KFAllP;
KFAllPInv = cellfun(@inv,KFAllPCell,'UniformOutput',false);
KFAllPInvTrace = cellfun(@trace,KFAllPInv)';
figure
hold on
plot(tMat,KFPartTunedPInvTrace)
plot(tMat,KFAllPInvTrace)
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

title('Trace inv(P)')

%% Plot trace P
KFPartTunedPCell = results.KFPartTunedP.KFPartTunedP;
KFPartTunedPTrace = cellfun(@trace,KFPartTunedPCell)';
KFAllPCell = results.KFAllP.KFAllP;
KFAllPTrace = cellfun(@trace,KFAllPCell)';
figure
hold on
plot(tMat,KFPartTunedPTrace)
plot(tMat,KFAllPTrace)
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
title('Trace P')


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

%%
keyboard
%% plot KFPartTunedP changing
% keyboard
% obj = VideoWriter('animation.avi');
% obj.Quality = 100;
% obj.FrameRate = 20;
% open(obj)

figure('Position', [10, 10, 1200, 1200])
for cnt=1:length(results.KFPartTunedP.KFPartTunedP)
heatmap(results.KFPartTunedP.KFPartTunedP{cnt,1});
caxis([0 0.5])
colormap jet;
colorbar;
title(num2str(cnt))
% drawnow;
pause(0.01)
% f = getframe(gcf);
% writeVideo(obj,f);
end
% obj.close();
%%
keyboard
%% plot KFAllP changing
% keyboard
figure
for cnt=1:length(results.KFAllP.KFAllP)
heatmap(results.KFAllP.KFAllP{cnt,1});
caxis([0 0.5])
colormap jet;
colorbar;
title(num2str(cnt))
% drawnow;
pause(0.01)
end
