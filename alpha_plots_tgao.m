function h = alpha_plots_tgao(option,plotparams)
% function h = alpha_plots_MIKE(option)
% This function plots some useful info about the alpha parameter we are
% considering.  You can produce different plots depending on what option
% you pass and associated parameters.
% Inputs : option - Which plot you want produced
%            1: 2D plot of error of alpha guesses given true alphas
%                plotparams should be the norm (1,2,inf) for error eval
%                   The default value is 2, if plotparams not set
%            2: 2D plot of Mean Exit Times for all alphas
%                plotparams is unneeded
%            3: 1D plot of Mean Exit Time for a fixed alpha
%                plotparams.alpha should be the requested alpha value
%                plotparams.N should be the number of points in the plot
%                   The default value is 30, if plotparams.N not set
% Outputs : h - The figure handle of the requested plot

switch option
    % Plot the errors associated with all possible true alphas and alpha
    % guesses.  Hopefully, the best accuracy occurs when the right alpha is
    % selected
    case 1
        if nargin==2
            normType = plotparams;
        else
            normType = 2;
        end
        
        alphainvec = linspace(.05,1.95,20);
        alphaoutvec = linspace(.04,1.96,20);
        errmat = zeros(length(alphainvec),length(alphaoutvec));
        
        k = 1;
        for alphain=alphainvec
            uTrue = escapenonsym(alphain);
            f = @(u) norm(u-uTrue,normType); % Need the right uTrue
            m = 1;
            for alphaout=alphaoutvec
                uGuess = escapenonsym(alphaout);
                errmat(k,m) = f(uGuess);
                m = m + 1;
            end
            k = k + 1;
        end
        
        [AI,AO] = meshgrid(alphainvec,alphaoutvec);
        h = surf(AI,AO,errmat');
        shading interp
        xlabel('True \alpha')
        ylabel('Guess of \alpha')
        zlabel('Error')
    
    % Plot the Mean Exit Times for the range of alpha values and domain
    case 2
        x = -1:1/100:1;
        alphavec = linspace(.05,1.95,39);
        MET = zeros(length(alphavec),length(x));
        k = 1;
        for alpha=alphavec
            MET(k,:) = escapenonsym(alpha);
            k = k + 1;
        end
        [A,X] = meshgrid(alphavec,x);
        h = surf(A,X,MET');
        shading interp
        xlabel('\alpha')
        ylabel('x')
        zlabel('Mean Exit Time')
    case 3
    % Given a true alpha, plot the error of associated guesses
        alpha = plotparams.alpha;
        if alpha>=2 || alpha<=0
            error('passed alpha=%g; alpha must be in (0,2)',alpha)
        end
        if isfield(plotparams,'N')
            N = floor(plotparams.N);
        else
            N = 30;
        end
        if N<3
            error('passed N=%g unacceptable; N>3 required',N)
        end
        if isfield(plotparams,'normType')
            normType = plotparams.normType;
        else
            normType = 2;
        end
        if ~(normType==1 || normType==2 || normType==inf)
            error('normType must be either 1, 2, or inf')
        end
        
        uTrue = escapenonsym(alpha);
            f = @(u) norm(u-uTrue,normType).^2./norm(uTrue,2).^2; % Need the right uTrue
        alphavec = linspace(.05,1.95,39);
        for k=1:length(alphavec)
            uGuess = escapenonsym(alphavec(k));
            errvec(k) = f(uGuess);
        end
        h = plot(alphavec,errvec);
        xlabel('\alpha guess')
        ylabel('absolute error')
        title(sprintf('True \\alpha=%g, Norm Type: %g',alpha,normType))
    otherwise
        error('Invalid plot option passed')
end