classdef KalmanFilterSysEstClass < handle
    
    
    properties
        name 
        results % to save data
        iter % numer of calls
        
        z % z is estimated vector
        P % P is estimation error covariance
        
        n % state dimension
        m % input dimension
        
        sigmaZ % covariance of system change process
        sigmaW % cocariance of measurement noise
        
    end
    
    methods
        function obj = KalmanFilterSysEstClass(initP, initZ, n, m, sigmaZ, sigmaW, name, results, iter)
            obj.P = initP;
            obj.z = initZ;
            
            obj.n = n;
            obj.m = m;
            
            obj.sigmaZ = sigmaZ;
            obj.sigmaW = sigmaW;
            
            try
                obj.name = name;
                obj.results = results;
                obj.iter = iter;
            catch
                obj.name = [];
                obj.results = [];
                obj.iter = 0;
            end
            
        end
        
        function estimate(obj, xMeas, xPrev, uPrev)
            
            obj.P = obj.P + obj.sigmaZ;
            
            if ~isempty(xMeas) % if measurement availabe
                C = kron(eye(obj.n), [xPrev', uPrev']);
                e = xMeas - C * obj.z;
                S = C * obj.P * C' + obj.sigmaW;
                K = obj.P * C' / S;
                obj.z = obj.z + K * e;
                obj.P = obj.P - K * C * obj.P;
            end
            
        end
        
        function lTrigger = update(obj, xMeas, xuPrev, mod, alpha)
            try
                xPrev = xuPrev(1:obj.n);
                uPrev = xuPrev(obj.n+1:end);
            catch
                xPrev = [];
                uPrev = [];
            end
            obj.estimate(xMeas, xPrev, uPrev);
            
            lTrigger = obj.trigger(mod, alpha);
            
            if ~isempty(obj.name) % save data
                obj.results.add(obj.iter,[obj.name,'obj'],obj)
                obj.results.add(obj.iter,[obj.name,'z'],obj.z)
                obj.results.add(obj.iter,[obj.name,'P'],obj.P)
                obj.results.add(obj.iter,[obj.name,'lTrigger'],lTrigger)
                obj.iter = obj.iter +1;
            end
            
        end
        
        function lTrigger = trigger(obj, mod, alpha)
            % KF Learning trigger
            zErr = obj.z - reshape([mod.A'; mod.B'],[],1);
            zErrSquare = zErr' * obj.P^-1 * zErr;
            PTrigger = 1-alpha;
            chi2test = zErrSquare / chi2inv(PTrigger,obj.n*(obj.n+obj.m));
            lTrigger = chi2test>1;
            
            if ~isempty(obj.name) % save data
                obj.results.add(obj.iter,[obj.name,'zErrSquare'],zErrSquare)
                obj.results.add(obj.iter,[obj.name,'chi2test'],chi2test)
            end
        end
        
        function [A, B] = getAB(obj)
            AB = reshape(obj.z,obj.n+obj.m,obj.n)';
            A = AB(:,1:obj.n);
            B = AB(:,obj.n+1:end);
        end
    end
end

