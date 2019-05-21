function [mk1k,Pk1k,vk1k,Vk1k] = factorizedGIWprediction(mkk,Pkk,vkk,Vkk,models)

% Models
Q = models.Q;
Ts = models.Ts;
d = models.d;
n = models.n;

% Motion model
[fm,Fm] = models.motionModel(mkk,Ts);

% Predicted mean
mk1k = fm;
% Predicted covariance
Pk1k = Fm*Pkk*Fm.' + Q;

% Transformation matrix
[Mm,M1,M2,I_x] = models.matrixTransformationFunction(mkk);

if ~isempty(I_x)
    
    if models.KLdiv_minimization_flag
        
        % Approximate the expected values
        [CmI,CmII] = approximateMatrixTransformationExpectedValues(Vkk,Pkk,Mm,M1,M2,I_x,d);
        
        % Compute the degrees of freedom by Newton's or Halley's method
        s_k = models.s_last;
        s_old = 0;
        sd1s_old = 0;
        iter = 1;
        while iter<1000
            % Increase iteration counter
            iter=iter+1;
            % Compute function h(s), and its first and second
            % order derivatives
            [h_k,hp_k,hb_k]= h_function(s_k,CmI,CmII,d);
            
            % Update the parameter
            % Halley's method
            s_k1 = s_k-(2*h_k*hp_k)/(2*hp_k^2-h_k*hb_k);
            if s_k1<d+2
                % Take Newton step instead
                s_k1 = s_k-h_k/hp_k;
            end
            
            % Make sure it is larger than minimum size
            s_k1=max(s_k1,d+2);
            
            % Check for convergence
            % 1) Difference between new parameter and last
            % parameter
            cc1 = abs(s_k1-s_k)<1e-3;
            % 2) Difference between new parameter and second to
            % last parameter. If they are the same, the
            % algorithm might be stuck going between two values
            cc2 = abs(s_k1-s_old)<1e-3;
            % 3) The parameter s occurs in the prediction of v
            % and V as 1-d/s-1/s. If this quantity has
            % converged, the parameter has converged.
            sd1s_new = 1-d/s_k-1/s_k;
            cc3 = abs(sd1s_new-sd1s_old)<1e-3;
            % Check if difference is smaller than threshold, or
            % if the difference itself is not changing
            if cc1 || cc2 || cc3
                s_k=s_k1;
                
                break
            else
                sd1s_old = sd1s_new;
                s_old = s_k;
                s_k=s_k1;
            end
        end
        
        s = s_k;
        % This is used to initiate s_k in the next prediction
        models.s_last = max(s_k,100);
        
        D = d+1;
        S = 1/(1-D/s);
        N = 1/(1-D/n);
        vpred = 2*D + (vkk-2*D)/(1 + (S*N-1)*(vkk/D-1));
        vpred = max(vpred,2*d+2+0.01);
        Vk1k = ((vpred-d-1)/(vkk-d-1))*(1-d/s-1/s)*(1-d/n-1/n)*CmII;
        vk1k = vpred;
        
    else
        
        % Compute approximations of expected values
        [~,E_X,E_invX] = approximateMatrixTransformationExpectedValues(Vkk,Pkk,Mm,M1,M2,I_x,d);
        
        % Wishart approximation by matching E[X] and E[X^-1]
        E_X_invX = E_X*E_invX;
        E_X_invX_I = E_X_invX - eye(d);
        % Least squares solution
        if norm(E_X_invX_I,'fro')<1e-6
            s = inf;
        else
            s = ((d+1)/d)*trace(E_X_invX_I\E_X_invX);
        end
        s = max(s,d+2);
        
        nu = 1 + (vkk-2*d-2)*(1/s + 1/n - (d+1)/n/s);
        
        vpred = d+1+(1/nu)*(vkk-d-1);
        vk1k = max(vpred,2*d+2+0.01);
        Vk1k = (1/nu)*(1-(d+1)/s)*(1-(d+1)/n)*E_X;
    end
    
else
    % If the matrix transformation function is not dependent on the
    % kinematic state, then the marginalisation does not incur any
    % uncertainty
    s_k = inf;
    CmII = Vkk;
    
    s = s_k;
    % This is used to initiate s_k in the next prediction
    models.s_last = max(s_k,100);
    
    D = d+1;
    S = 1/(1-D/s);
    N = 1/(1-D/n);
    vpred = 2*D + (vkk-2*D)/(1 + (S*N-1)*(vkk/D-1));
    vpred = max(vpred,2*d+2+0.01);
    Vk1k = ((vpred-d-1)/(vkk-d-1))*(1-d/s-1/s)*(1-d/n-1/n)*CmII;
    vk1k = vpred;
end



end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [h_k,hp_k,hb_k]= h_function(s_k,CmI,CmII,d)

% The parameter s is the solution to h(s)=0

% Function
h_k = d*log(s_k/2)-sum(psi(0,(s_k+1-(1:d))/2))+CmI-log(det(CmII));
% Differentiated once
hp_k = d/s_k-0.5*sum(psi(1,(s_k+1-(1:d))/2));
% Differentiated twice
hb_k = -d/s_k^2-0.25*sum(psi(2,(s_k+1-(1:d))/2));
end


