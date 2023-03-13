classdef MPCClass < handle
    
    properties
        N
        Q
        R
        QN
        mod
        n
        m
        XPoly
        UPoly
        controller
        xTrajPrev
        uTrajPrev
        
        XPolyN
        
        identStruct
        identController
    end
    
    methods
        function obj = MPCClass(mod, N, Q, R, QN, XPoly, UPoly, identStruct)
            % Construct an instance of this class
            
            try
                obj.identStruct = identStruct;
            catch
                obj.identStruct = [];
            end
            
            obj.mod = mod;
            obj.N = N;
            obj.Q = Q;
            obj.R = R;
            obj.QN = QN;
            
            obj.n = size(Q,1);
            obj.m = size(R,1);
            
            if isempty(XPoly)
                obj.XPoly = Polyhedron.fullSpace(obj.n);
            else
                obj.XPoly = XPoly;
            end
            if isempty(UPoly)
                obj.UPoly = Polyhedron.fullSpace(obj.m);
            else
                obj.UPoly = UPoly;
            end
            
            %% controller
            % initialize YALMIP optimizer
            initializeController(obj)
            
        end
        
        function initializeController(obj)
            % terminal set
            % findmas (A First Course in Predictive Control, J. A. Rossiter)
            %%% Process is x(k+1)=A_K x(k)
            %%% Constraints at each sample are Cx <=f
            K = LQR(obj.Q,obj.R,obj.mod);
            C = [obj.XPoly.A; obj.UPoly.A* K];
            f = [obj.XPoly.b; obj.UPoly.b];
            A_K = obj.mod.A+obj.mod.B*K;
            [F,t]=findmas(A_K,C,f);
            obj.XPolyN = Polyhedron(F,t);
            obj.XPolyN.minHRep;
            
            
            % initialize YALMIP optimizer
            u = sdpvar(repmat(obj.m,1,obj.N),repmat(1,1,obj.N));
            x = sdpvar(repmat(obj.n,1,obj.N+1),repmat(1,1,obj.N+1));
            xRef = sdpvar(repmat(obj.n,1,obj.N+1),repmat(1,1,obj.N+1));
            
%             slackx = sdpvar(obj.n,obj.N); 10000*sum(sum(abs(slackx))) +
            slackXb = sdpvar(size(obj.XPoly.b,1),obj.N);
            slackXN = sdpvar(size(obj.XPolyN.b,1),1);
            
            slackObjective =  1000*sum(sum(abs(slackXb))) + sum(sum(slackXb.^2)) + 1000*sum(sum(abs(slackXN))) + sum(sum(slackXN.^2));
            constraints = [slackXb>=0, slackXN>=0];
            objective = 0;
            for k = 1:obj.N
                objective = objective + (x{k}-xRef{k})'*obj.Q*(x{k}-xRef{k}) + u{k}'*obj.R*u{k};
                constraints = [constraints, x{k+1} == obj.mod.A*x{k}+obj.mod.B*u{k}];% + slackx(:,k)];
                constraints = [constraints, obj.XPoly.A * x{k+1} <= obj.XPoly.b + slackXb(:,k)];
                constraints = [constraints, obj.UPoly.A * u{k} <= obj.UPoly.b];
                constraints = [constraints, obj.XPoly.Ae * x{k+1} == obj.XPoly.be];
                constraints = [constraints, obj.UPoly.Ae * u{k} == obj.UPoly.be];
                
            end
            objective = objective + (x{obj.N+1}-xRef{obj.N+1})'*obj.QN*(x{obj.N+1}-xRef{obj.N+1});
            constraints = [constraints, obj.XPolyN.A * x{obj.N+1} <= obj.XPolyN.b + slackXN];
            
            parameters_in = {x{1},[xRef{:}]};
            solutions_out = {[u{:}], [x{:}]};
            
            %             controller = optimizer(constraints, objective,sdpsettings('solver','gurobi'),parameters_in,solutions_out);
%             options = sdpsettings('verbose',2, 'debug',1, 'warning',1);
            options = sdpsettings('warning',1);
            obj.controller = optimizer(constraints, objective + slackObjective,options,parameters_in,solutions_out);
            
        end
        
        function initializeIdentController(obj)
            % initialize YALMIP optimizer
            u = sdpvar(repmat(obj.m,1,obj.N),repmat(1,1,obj.N));
            x = sdpvar(repmat(obj.n,1,obj.N+1),repmat(1,1,obj.N+1));
            xRef = sdpvar(repmat(obj.n,1,obj.N+1),repmat(1,1,obj.N+1));
            
            % variables for KF
            P = sdpvar(repmat(obj.n*(obj.n+obj.m),1,obj.N+1),repmat(obj.n*(obj.n+obj.m),1,obj.N+1));
            %             vecTheta =  sdpvar(repmat(obj.n*(obj.n+obj.m),1,obj.N+1),repmat(1,1,obj.N+1));
            vecTheta =  sdpvar(repmat(obj.n*(obj.n+obj.m),1,1),repmat(1,1,1));
            
            KF.SInv = sdpvar(repmat(obj.n,1,obj.N+1),repmat(obj.n,1,obj.N+1));
            
            KF.SigmaW = obj.mod.W*obj.mod.W';
            KF.Q = 0.0001*eye(obj.n*(obj.n+obj.m));
            
            constraints = [];
            objective = 0;
            for k = 1:obj.N
                objective = objective + (x{k}-xRef{k})'*obj.Q*(x{k}-xRef{k}) + u{k}'*obj.R*u{k};
                
                %                 constraints = [constraints, x{k+1} == obj.mod.A*x{k}+obj.mod.B*u{k}];
                KF.C = kron(eye(obj.n), [x{k}; u{k}]');
                constraints = [constraints, x{k+1} == KF.C * vecTheta];
                
                constraints = [constraints, obj.XPoly.A * x{k} <= obj.XPoly.b];
                constraints = [constraints, obj.UPoly.A * u{k} <= obj.UPoly.b];
                constraints = [constraints, obj.XPoly.Ae * x{k} == obj.XPoly.be];
                constraints = [constraints, obj.UPoly.Ae * u{k} == obj.UPoly.be];
                
                if ~isempty(obj.identStruct)
                    % add parameter covariance to cost function
                    
                    %%%
                    %Time update
                    %                   vecTheta = vecTheta;
                    P_kp1gk = P{k} + KF.Q;
                    
                    %Measurement update
                    KF.S = KF.C * P_kp1gk * KF.C' + KF.SigmaW;
                    constraints = [constraints, KF.SInv{k+1} * KF.S == eye(obj.n)];
                    KF.K = P_kp1gk * KF.C' * KF.SInv{k+1};
                    constraints = [constraints, P{k+1} == ( eye(obj.n*(obj.n+obj.m)) - KF.K * KF.C) * P_kp1gk * ( eye(obj.n*(obj.n+obj.m)) - KF.K * KF.C)' + KF.K * KF.SigmaW * KF.K'];
                    %%%
                    
                    objective = objective + trace(KF.P{k+1}); %minimize trace(P)
%                     objective = objective - trace(KF.P{k+1}^-1); %maximize trace(inv(P))
                end
            end
            objective = objective + (x{obj.N+1}-xRef{obj.N+1})'*obj.QN*(x{obj.N+1}-xRef{obj.N+1});
            
            parameters_in = {x{1},[xRef{:}]};
            solutions_out = {[u{:}], [x{:}]};
            
            %             controller = optimizer(constraints, objective,sdpsettings('solver','gurobi'),parameters_in,solutions_out);
            % options = sdpsettings('verbose',2, 'debug',1, 'warning',1);
            options = sdpsettings('warning',1);
            obj.identController = optimizer(constraints, objective,options,parameters_in,solutions_out);
            
        end
        
        function solutions = optimize(obj,x0,xRef)
            % Compute optimal open loop trajectory
            
            inputs = {x0, xRef};
            [solutions,diagnostics] = obj.controller{inputs};
            U = solutions{1};
            X = solutions{2};
            obj.uTrajPrev = U;
            obj.xTrajPrev = X;
            if diagnostics == 1
                error('The problem is infeasible');
            end
            
        end
        
        function solutions = optimizeWithP(obj, x0, xRef, P0, vecTheta, KFWeight, KFQ, KFSigmaW)
            algo = 'optimal';
%             algo = 'cheap';
            switch algo
                case 'cheap'
                    % cheap workaround with random noise
                    controlSolutions = optimize(obj,x0,xRef);
                    controlSolutions{1} = controlSolutions{1} + 10*randn(size(controlSolutions{1}));
                    for k = 1:obj.N
                        controlSolutions{2}(:,k+1) = obj.mod.A * controlSolutions{2}(:,k) + obj.mod.B * controlSolutions{1}(:,k);
                    end
                    solutions = controlSolutions;
                    
                case 'optimal'
                    % Compute optimal open loop trajectory
                    N = obj.N;
                    n = obj.n;
                    m = obj.m;
                    mod = obj.mod;
                    
%                     KFSigmaW = obj.mod.W*obj.mod.W';
%                     KFQ = 0.0001*eye(obj.n*(obj.n+obj.m));
                    
                    %%
                    % Define Matrcies for quadprog
                    % Set control and state constraints for time instance
                    
                    % Equality constraints
                    Aeq = full([kron(eye(N+1),eye(n)) + kron(spdiags(ones(N,1),-1,N+1,N+1),-mod.A), kron(spdiags(ones(N,1),-1,N+1,N),-mod.B)]);
                    beq = [x0; zeros(n*N,1)];
                    
                    
                    % Inequality constraints
                    Ax = [];
                    bx = [];
                    Au = [];
                    bu = [];
                    
                    for k = 0:N-1
                        
                        % new constraints
                        X_k = obj.XPoly;
                        U_k = obj.UPoly;
                        
                        Ax = blkdiag(Ax,X_k.A);
                        bx = [bx;X_k.b];
                        Au = blkdiag(Au,U_k.A);
                        bu = [bu;U_k.b];
                        
%                         disp(['k= ' num2str(k)])
                    end
                    
                    % Terminal constraint
%                     Ax = blkdiag(Ax,X_k.A);
%                     bx = [bx;X_k.b];
                    Ax = blkdiag(Ax,obj.XPolyN.A);
                    bx = [bx;obj.XPolyN.b];
                    
                    % concatenated
                    A_ineq = blkdiag(Ax,Au);
                    b_ineq = [bx;bu];
                    
                    % cost function
                    QBlkDiag = blkdiag(kron(eye(N),obj.Q),obj.QN);
                    RBlkDiag = kron(eye(N),obj.R);
                    QRBlkDiag = blkdiag(QBlkDiag, RBlkDiag);
                    
                    
                    
                    
                    % compute new MPC prediction
                    
                    % initial guess for optimization
                    K = LQR(obj.Q,obj.R,obj.mod);
                    x_MPC = zeros(n,N+1);
                    x_MPC(:,1) = x0;
                    u_MPC = zeros(m,N);
                    for tSim = 1:N
                        u_MPC(:,tSim) = K * (x_MPC(:,tSim) - xRef(tSim));
                        x_MPC(:,tSim+1) = mod.A * x_MPC(:,tSim) + mod.B * u_MPC(:,tSim);
                    end
                    xInit = reshape(x_MPC,[],1);
                    uInit = reshape(u_MPC,[],1);
                    
                    
                    %%%%%% fmincon
                    fminconoptions2 = optimoptions('fmincon','Display','notify', 'Algorithm', 'interior-point');
%                                     fminconoptions2 = optimoptions('fmincon','Display','notify', 'Algorithm', 'sqp');
%                                     fminconoptions2 = optimoptions('fmincon','Display','notify', 'Algorithm', 'active-set');
                    %                 solutionOL = fmincon(FUN,[xInit;uInit],A_ineq,b_ineq,Aeq,beq,[],[],[],fminconoptions2);
                    %
                    %                 %         f = zeros(1,n*(N+1)+m*N);
                    %                 %         solutionOL = quadprog(S,f,A_ineq,b_ineq,Aeq,beq,[],[],[xInit;uInit],options);
                    %
                    
%                     dref = [reshape(xRef,[],1); zeros(N*m,1)];
                    olutions = optimize(obj,x0,xRef);
%                     dref = [reshape(olutions{2},[],1); olutions{1}'];
                    dref = [reshape([olutions{2}(:,1:end-1), xRef(:,end)],[],1); olutions{1}']; % terminal point from reference.
                    
                    KFCostFunHan = @(d) KFCost(P0,  reshape(d(1:(N+1)*n),n,N+1), reshape(d((N+1)*n+1 : end),m,N), KFQ, KFSigmaW);
                    LQRFun = @(d) (d-dref)'*QRBlkDiag*(d-dref);
                    FUN = @(d) LQRFun(d) + KFWeight * KFCostFunHan(d);
                    
%                     [solutionOL,fval,exitflag,output] = fmincon(FUN,[xInit;uInit],A_ineq,b_ineq,Aeq,beq,[],[],[],fminconoptions2);
                    [solutionOL,fval,exitflag,output] = fmincon(FUN,dref,A_ineq,b_ineq,Aeq,beq,[],[],[],fminconoptions2);
                    
                    disp(['LQR cost: ' num2str(LQRFun(solutionOL)) ' KF cost: ' num2str(KFWeight * KFCostFunHan(solutionOL))])
                    
                    
                    if any(isnan(solutionOL)) || ~isequal(exitflag,1)
                        %                     keyboard
%                         fminconoptions2 = optimoptions('fmincon','Display','notify', 'Algorithm', 'sqp');
                        warning('nan')
%                         [solutionOL,fval,exitflag,output] = fmincon(FUN,[xInit;uInit],A_ineq,b_ineq,Aeq,beq,[],[],[],fminconoptions2);
%                         if any(isnan(solutionOL)) || ~isequal(exitflag,1)
%                             %                     keyboard
%                             fminconoptions2 = optimoptions('fmincon','Display','notify', 'Algorithm', 'active-set');
%                             warning('nan')
%                             [solutionOL,fval,exitflag,output] = fmincon(FUN,[xInit;uInit],A_ineq,b_ineq,Aeq,beq,[],[],[],fminconoptions2);
%                         end
%                         solutionOL = optimize(obj,x0,xRef);
                    end
                    
                    % extract solution
                    x_MPC = reshape(solutionOL(1:(N+1)*n),n,N+1);
                    u_MPC = reshape(solutionOL((N+1)*n+1 : end),m,N);
                    solutions{2} = x_MPC;
                    solutions{1} = u_MPC;
                    
                    
                    
                    %% fmincon
                    % beq = [x_m; zeros(n*N,1)];
                    % solutionOL = fmincon(FUN,[xInit;uInit],A_ineq,b_ineq,Aeq,beq,[],[],[],fminconoptions2);
                    %
                    
                otherwise
            end
            
        end
    end
end

