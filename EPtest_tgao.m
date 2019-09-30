 %function [alphaGD,pp]=GDtest_MIKE(alpha)

% This file should mimic the use of gradient descent, and other methods, to
% find the correct alpha value for the given test value
%
% I have written this using the outline given by Ting Gao
% The description of the problem can be found in a paper in this dropbox
% The numerical PDE solution is found in the mean_exit_time file

% Choose a true solution to the problem
alphaTrue = 1.5;
uTrue = escapenonsym(alphaTrue);

% Choose a norm in which to evaluate the error
normType = 2;

% Choose a boundary for alpha, beyond which we don't consider
alphaBOUNDS.min = .01;
alphaBOUNDS.max = 1.99;

% Create a function we are trying to minimize
f = @(u) norm(u-uTrue,normType).^2./norm(uTrue,2).^2;

% Choose tolerances for the gradient descent method
GDtol = 1e-5;
GDdiff = 1e-7;
maxIter = 30;
% Line search options:
% 0 - None, just take the full step
% 1 - Keep cutting in half until f is reduced
% 2 - Quadratic interpolant minimum
% 3 - Cubic interpolant minimum
% 4 - Scaled step size with iteration
lineSearch = 1;
lsMin = 1e-2; % minimum step size for option 1
stepPower = 2; % step size scaled by 1/iter^stepPower

% Initialize the gradient descent iteration
alpha = 1;
u = escapenonsym(alpha);
alphaHist = alpha;
fu = f(u);
fDiff = GDtol + 1; % Set to force entry into loop
iter = 1;

% Iterate until finding the correct alpha
while fDiff>GDtol && iter<=maxIter
    % Compute the finite difference gradient approximation
    alphaDiff = alpha + GDdiff;
    uPlusDelta = escapenonsym(alphaDiff);
    fuPlusDelta = f(uPlusDelta);
    alphaStep = (fuPlusDelta-fu)/GDdiff;
    
    % Perform the requested line search
    switch lineSearch
        % No line search, just take the full step
        case 0
            alphaTest = alpha - alphaStep;
            
            % Correct in case we are now outside the bounds
            alpha = alpha_check_tgao(alphaBOUNDS,alphaTest);
            
            unew = escapenonsym(alpha);
            funew = f(unew);
            
        % Cut in half line search, stops when norm decreases
        case 1
            ls = 1;
            while ls>lsMin
                alphaTest = alpha - ls*alphaStep;
                
                % Correct in case we are now outside the bounds
                alphaTest = alpha_check_tgao(alphaBOUNDS,alphaTest);
                
                unew = escapenonsym(alphaTest);
                funew = f(unew);
                if funew<fu   % Stop if there is a decrease
                    alpha = alphaTest;
                    ls = 0;
                else % Cut in half if no decrease occurs
                    ls = ls/2;
                end
            end
            
            if ls~=0
                funew = fu;
                alpha = alphaHist(iter);
                warning('Line search may have failed with option 1')
            end
            
        % Quadratic interpolation line search
        case 2
            % Evaluate the interpolant at a couple points
            alpha2 = alpha - alphaStep;
            alpha2 = alpha_check_tgao(alphaBOUNDS,alpha2);
            u2 = escapenonsym(alpha2);
            fu2 = f(u2);
            
            alpha1 = alpha - .5*alphaStep;
            if alpha2==alphaBOUNDS.min || alpha2==alphaBOUNDS.max
                alpha1 = alpha_check_tgao(alphaBOUNDS,alpha1,.2);
            end
            u1 = escapenonsym(alpha1);
            fu1 = f(u1);
            
            % Minimize the interpolant    
            alphaTest = fminbnd(@(aTest)polyval(polyfit([alpha,alpha1,alpha2],[fu,fu1,fu2],2),aTest),alphaBOUNDS.min,alphaBOUNDS.max);
            % Take the suggested minimum step that fits the boundary
            alphaBOUNDS.cushion = 0;
            alpha = alpha_check_tgao(alphaBOUNDS,alphaTest);
            unew = escapenonsym(alpha);
            funew = f(unew);
            
            % Check to see if we did better
            if funew>fu
                funew = fu;
                alpha = alphaHist(iter);
                warning('Line search may have failed with option 2')
            end
            
        % Cubic interpolation line search
        case 3
            % Evaluate the interpolant at a couple points
            alpha3 = alpha - alphaStep;
            alpha3 = alpha_check_tgao(alphaBOUNDS,alpha3);
            u3 = escapenonsym(alpha3);
            fu3 = f(u3);
            
            alpha1 = alpha - .25*alphaStep;
            alpha2 = alpha - .75*alphaStep;
            if alpha3==alphaBOUNDS.min || alpha3==alphaBOUNDS.max
                alpha1 = alpha_check_tgao(alphaBOUNDS,alpha1,.15);
                alpha2 = alpha_check_tgao(alphaBOUNDS,alpha2,.3);
            end
            u1 = escapenonsym(alpha1);
            fu1 = f(u1);
            u2 = escapenonsym(alpha2);
            fu2 = f(u2);
            
            % Minimize the interpolant    
            alphaTest = fminbnd(@(aTest)polyval(polyfit([alpha,alpha1,alpha2,alpha3],[fu,fu1,fu2,fu3],3),aTest),alphaBOUNDS.min,alphaBOUNDS.max);
            % Take the suggested minimum step that fits the boundary
            alpha = alpha_check_tgao(alphaBOUNDS,alphaTest);
            unew = escapenonsym(alpha);
            funew = f(unew);
            
            % Check to see if we did better
            if funew>fu
                funew = fu;
                alpha = alphaHist(iter);
                warning('Line search may have failed with option 3')
            end
            
        % Use the reduction in iteration per step
        case 4
            ls = iter^(-stepPower);
            alphaTest = alpha - ls*alphaStep;
            
            % Correct in case we are now outside the bounds
            alpha = alpha_check_tgao(alphaBOUNDS,alphaTest);
            
            unew = escapenonsym(alpha);
            funew = f(unew);
            
        otherwise
            error('Invalid line search option passed')
    end
    
    fDiff = abs(funew-fu);
    u = unew;
    fu = funew;
    iter = iter + 1;
    alphaHist(iter) = alpha;
end

% Save the gradient descent result
alphaGD = alpha;

% Find the result using built-in Matlab function
alphaFB = fminbnd(@(a)f(escapenonsym(a)),alphaBOUNDS.min,alphaBOUNDS.max);

pp.alpha = alphaTrue;
pp.normType = normType;
h = alpha_plots_tgao(3,pp);
hold on
%plot(alphaGD,f(escapenonsym(alphaGD)),'or')
plot(alphaFB,f(escapenonsym(alphaFB)),'xr')

%legend('Objective function f(\alpha)','Gradient descent','fminbnd','location','best')
legend('Objective function f(\alpha)','Estimation of \alpha')


