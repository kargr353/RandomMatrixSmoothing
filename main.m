%% True track
plotFlag = 0;

% Number of measurements per time step (can be a vector)
NZ = 10; %[10 10];
% Probability of detection
PD = 0.49; %[0.25 0.5 0.75];

% Monte Carlo runs
MC = 100;

% Number of time steps
Nt = 100;

% Sample time
T = 1;
% Process noise, acceleration and turn-rate
siga = 10;
sigw = 1*pi/180;

% Target dimensions
target_width = 20; 
target_length = 50;

% Models
[models,models_conditional,Xbase] = ...
    initialize_tracking_parameters(T,siga,sigw,target_width,target_length);

% Dimension of extent
d=size(Xbase,1);

% Initial GIW parameters
m0 = [0 0 35 0 0].';
P0 = diag([10 10 50 50 1*pi/180].^2);
X0 = diag([25 25].^2);
v0 = 2*2+2+1;
V0 = (v0-2*2-2)*X0;


% True motion generation
TrueMotionModel = 1;

switch TrueMotionModel
    case 1
        % Coordinated turn, Cartesian velocity
        
        % Initial state
        x0 = [35; 0; 0*pi/180];
        
        % Maximum turn-rate
        maxTR = 10*pi/180;
        
        % Process noise, acceleration and turn-rate
        sigatrue = siga;
        sigwtrue = sigw;
        Gtrue = [0.5*T^2*eye(2) zeros(2,1); T*eye(2) zeros(2,1); zeros(1,2) T];
        Qatrue = blkdiag(eye(2)*0.5*sigatrue^2,sigwtrue^2);
        Qtrue = Gtrue*Qatrue*Gtrue';
        
        Qscale_cond = 1;
        Qscale_fact = 25;
        
        P0true = P0;
        
    case 2
        % Constant Cartesian velocity
        
        % Initial state
        x0 = [35; 0; 0*pi/180];
        
        % Maximum turn-rate
        maxTR = 10*pi/180;
        
        % Process noise, acceleration and turn-rate
        sigatrue = siga/100;
        sigwtrue = 0;
        Gtrue = [0.5*T^2*eye(2) zeros(2,1); T*eye(2) zeros(2,1); zeros(1,2) T];
        Qatrue = blkdiag(eye(2)*0.5*sigatrue^2,sigwtrue^2);
        Qtrue = Gtrue*Qatrue*Gtrue';
        
        Qscale_cond = 1;
        Qscale_fact = 25;
        
        P0true = P0;
end

for ipd = 1:length(PD)
    for inz = 1:length(NZ)
    
        
        % Number of measurement per time step
        Nz = NZ(inz);
        
        % Probability of detection
        p_D = PD(ipd);
        
        % Class instances
        clear GIWfilters
        GIWfilters = cell(1,0);
        
        % Conditional CV
        GIWfilters{end+1} = condGIW_forwardFilter_backwardSmoother;
        GIWfilters{end}.models = models_conditional;
        GIWfilters{end}.models.D = (Qscale_cond^2)*models_conditional.D;
        GIWfilters{end}.m0 = m0;
        GIWfilters{end}.P0 = P0+diag([10 10 0 0 0].^2);
        GIWfilters{end}.v0 = v0;
        GIWfilters{end}.V0 = V0;

        % Factorised CV
        GIWfilters{end+1} = factGIW_forwardFilter_backwardSmoother;
        GIWfilters{end}.models = models;
        GIWfilters{end}.models.motionModel = @ConstantVelocity;
        GIWfilters{end}.models.matrixTransformationFunction = @noMatrixTransformation;
        GIWfilters{end}.models.inverseMatrixTransformationFunction = @inverseNoMatrixTransformation;
        GIWfilters{end}.models.n_x = 4;
        GIWfilters{end}.models.H = [eye(GIWfilters{end}.models.d) zeros(GIWfilters{end}.models.d,GIWfilters{end}.models.n_x-GIWfilters{end}.models.d)];
        GIWfilters{end}.models.Q = (Qscale_fact^2)*models.Q(1:4,1:4);
        GIWfilters{end}.models.n = models_conditional.n;
        GIWfilters{end}.m0 = m0(1:4);
        GIWfilters{end}.P0 = P0(1:4,1:4)+diag([10 10 0 0].^2);
        GIWfilters{end}.v0 = v0;
        GIWfilters{end}.V0 = V0;
        
        % Factorised CT
        GIWfilters{end+1} = factGIW_forwardFilter_backwardSmoother;
        GIWfilters{end}.models = models;
        GIWfilters{end}.models.KLdiv_minimization_flag = 0;
        GIWfilters{end}.m0 = m0;
        GIWfilters{end}.P0 = P0+diag([10 10 0 0 0].^2);
        GIWfilters{end}.v0 = v0;
        GIWfilters{end}.V0 = V0;
        
        
        lgnd = {'CCV','FCV','FCT'};
        
        number_of_smoothers = numel(GIWfilters);
        
        % Run filters
        
        % Allocate memory
        GWD = zeros(Nt,MC,3,number_of_smoothers);
        KIN = zeros(Nt,MC,3,number_of_smoothers);
        EXT = zeros(Nt,MC,3,number_of_smoothers);
        
        runTime = 100;
        cycleTime = 0;
        ctr = 1;
        for imc = 1:MC
            runTime = runTime+cycleTime;
            if runTime>0.5
                fprintf([' ' num2str(imc)]);
                runTime = 0;
                ctr = ctr+1;
                if ctr>10
                    ctr = 1;
                    disp(' ')
                end
            end
            
            tic;
            
            %%%%%%%%%%%%%%%%%%%%%%%
            % Generate true track
            %%%%%%%%%%%%%%%%%%%%%%%
            [xtrue,Xtrue] = generateTrueTrack(Nt,x0,P0true,Xbase,Qtrue,T,maxTR);

            %%%%%%%%%%%%%%%%%%%%%%%%%
            % Generate measurements
            %%%%%%%%%%%%%%%%%%%%%%%%%
            Z = generateMeasurements(Nt,p_D,xtrue,Xtrue,Nz,models.R,d);

            
            if imc<=10
                figure(100),clf
                plot(xtrue(1,:),xtrue(2,:),'-','linewidth',2)
                hold on
                for tt = 1:Nt
                    plotCovariance(xtrue(1:2,tt),Xtrue(:,:,tt),'k',2);
                    hold on
                    plot(Z{tt}(1,:),Z{tt}(2,:),'r.')
                end
                axis equal
                xlabel('$p^{x}, \ [m]$','interpreter','latex','fontsize',15)
                ylabel('$p^{y}, \ [m]$','interpreter','latex','fontsize',15)
                drawnow
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%
            % Run filters
            %%%%%%%%%%%%%%%%%%%%%%%%%
            for ff = 1:number_of_smoothers
                GIWfilters{ff}.forwardFilter(Z);
                GIWfilters{ff}.backwardSmoother;
                GIWfilters{ff}.GWDmetric(xtrue,Xtrue);
                
                GWD(:,imc,:,ff) = GIWfilters{ff}.GWD;
                KIN(:,imc,:,ff) = GIWfilters{ff}.KIN;
                EXT(:,imc,:,ff) = GIWfilters{ff}.EXT;
            end
            
            cycleTime = toc;
        end
        disp(' ')

    end
end

%%
for ff = 1:number_of_smoothers
    figure(10+ff),clf
    for t=Nt:-1:1
        plotCovariance(xtrue(1:2,t),Xtrue(:,:,t),'k',2)
        hold on
        
        xpred = GIWfilters{ff}.mpred(1:2,t);
        Xpred = GIWfilters{ff}.Vpred(:,:,t)/(GIWfilters{ff}.vpred(t)-2*d-2);
        
        plotCovariance(xpred,Xpred,'g',2)
        
        xup = GIWfilters{ff}.mup(1:2,t);
        Xup = GIWfilters{ff}.Vup(:,:,t)/(GIWfilters{ff}.vup(t)-2*d-2);
        
        plotCovariance(xup,Xup,'r',2)
        
        xsm = GIWfilters{ff}.msm(1:2,t);
        Xsm = GIWfilters{ff}.Vsm(:,:,t)/(GIWfilters{ff}.vsm(t)-2*d-2);
        
        plotCovariance(xsm,Xsm,'b',2)
        
    end
    axis equal
    title(lgnd{ff})
end

%%

smcol = [clcol('b');clcol('p');clcol('o')];

ttl = {'$k|k-1$','$k|k$','$k|K$'};

meanKIN = mean(KIN,2);
meanEXT = mean(EXT,2);
meanGWD = mean(GWD,2);

figure(1),clf
for dd = 1:3
    subplot(3,3,dd)
    for ff = 1:number_of_smoothers
        semilogy(meanKIN(:,1,dd,ff),'linewidth',3,'color',smcol(ff,:))
        hold on
    end
    set(gca,'ylim',[min(meanKIN(:)) max(meanKIN(:))])
    grid on
    ylabel(['Kin ' ttl{dd}],'interpreter','latex')
    xlabel('Time','interpreter','latex')
    legend(lgnd,'location','best')
    
    subplot(3,3,3+dd)
    for ff = 1:number_of_smoothers
        semilogy(meanEXT(:,1,dd,ff),'linewidth',3,'color',smcol(ff,:))
        hold on
    end
    set(gca,'ylim',[min(meanEXT(:)) max(meanEXT(:))])
    grid on
    ylabel(['Ext ' ttl{dd}],'interpreter','latex')
    xlabel('Time','interpreter','latex')
    
    subplot(3,3,6+dd)
    for ff = 1:number_of_smoothers
        semilogy(meanGWD(:,1,dd,ff),'linewidth',3,'color',smcol(ff,:))
        hold on
    end
    set(gca,'ylim',[min(meanGWD(:)) max(meanGWD(:))])
    grid on
    ylabel(['GWD ' ttl{dd}],'interpreter','latex')
    xlabel('Time','interpreter','latex')
    
end

medianKIN = median(KIN,2);
medianEXT = median(EXT,2);
medianGWD = median(GWD,2);

figure(2),clf
for dd = 1:3
    subplot(3,3,dd)
    for ff = 1:number_of_smoothers
        semilogy(medianKIN(:,1,dd,ff),'linewidth',3,'color',smcol(ff,:))
        hold on
    end
    set(gca,'ylim',[min(medianKIN(:)) max(medianKIN(:))])
    grid on
    ylabel(['Kin ' ttl{dd}],'interpreter','latex')
    xlabel('Time','interpreter','latex')
    legend(lgnd,'location','best')
    
    subplot(3,3,3+dd)
    for ff = 1:number_of_smoothers
        semilogy(medianEXT(:,1,dd,ff),'linewidth',3,'color',smcol(ff,:))
        hold on
    end
    set(gca,'ylim',[min(medianEXT(:)) max(medianEXT(:))])
    grid on
    ylabel(['Ext ' ttl{dd}],'interpreter','latex')
    xlabel('Time','interpreter','latex')
    
    subplot(3,3,6+dd)
    for ff = 1:number_of_smoothers
        semilogy(medianGWD(:,1,dd,ff),'linewidth',3,'color',smcol(ff,:))
        hold on
    end
    set(gca,'ylim',[min(medianGWD(:)) max(medianGWD(:))])
    grid on
    ylabel(['GWD ' ttl{dd}],'interpreter','latex')
    xlabel('Time','interpreter','latex')
    
end

%%
smcol = [clcol('b');clcol('p');clcol('o')];

ttl = {'$k|k-1$','$k|k$','$k|K$'};

medianGWD = median(GWD,2);
titles = {'Prediction errors','Filtering errors','Smoothing errors'};

figure(3),clf
for dd = 1:3

    subplot(3,1,dd)
    for ff = 1:number_of_smoothers
        semilogy(medianGWD(:,1,dd,ff),'linewidth',3,'color',smcol(ff,:))
        hold on
    end
    %set(gca,'ylim',[min(medianGWD(:)) max(medianGWD(:))])
    %grid on
    ylabel(['$\Delta_{' ttl{dd}(2:end-1) '}$'],'interpreter','latex','fontsize',16)
    xlabel('Time','interpreter','latex','fontsize',16)
    title(titles{dd},'interpreter','latex','fontsize',16)
    
end
set(gcf,'color','white')


if MC>1
    number_of_bins = max(MC/25,10);
    
    figure(4),clf
    for dd = 1:3
        subplot(3,3,dd)
        for ff = 1:number_of_smoothers
            [N,X] = hist(mean(KIN(:,:,dd,ff),1),number_of_bins);
            plot(X,N/MC,'linewidth',3,'color',smcol(ff,:))
            hold on
        end
        ylabel(['Rel freq'],'interpreter','latex')
        xlabel(['Kin ' ttl{dd}],'interpreter','latex')
        
        subplot(3,3,3+dd)
        for ff = 1:number_of_smoothers
            [N,X] = hist(mean(EXT(:,:,dd,ff),1),number_of_bins);
            plot(X,N/MC,'linewidth',3,'color',smcol(ff,:))
            hold on
        end
        ylabel(['Rel freq'],'interpreter','latex')
        xlabel(['Ext ' ttl{dd}],'interpreter','latex')
        
        subplot(3,3,6+dd)
        for ff = 1:number_of_smoothers
            [N,X] = hist(mean(GWD(:,:,dd,ff),1),number_of_bins);
            plot(X,N/MC,'linewidth',3,'color',smcol(ff,:))
            hold on
        end
        ylabel(['Rel freq'],'interpreter','latex')
        xlabel(['GWD ' ttl{dd}],'interpreter','latex')
        
    end
    
    
    figure(5),clf
    for dd = 1:3
        subplot(3,3,dd)
        h = boxplot(squeeze(sum(KIN(:,:,dd,:),1)),'labels',lgnd,'color',smcol,'symbol','');
        set(h,{'linew'},{3});
        ylabel(['Kin ' ttl{dd}],'interpreter','latex')
        drawnow
        
        subplot(3,3,3+dd)
        h = boxplot(squeeze(sum(EXT(:,:,dd,:),1)),'labels',lgnd,'color',smcol,'symbol','');
        set(h,{'linew'},{3});
        ylabel(['Ext ' ttl{dd}],'interpreter','latex')
        drawnow
        
        subplot(3,3,6+dd)
        h = boxplot(squeeze(sum(GWD(:,:,dd,:),1)),'labels',lgnd,'color',smcol,'symbol','');
        set(h,{'linew'},{3});
        ylabel(['GWD ' ttl{dd}],'interpreter','latex')
        drawnow
    end
            
    figure(6),clf
    for dd = 1:3
        subplot(1,3,dd)
        boxplot(log(squeeze(sum(GWD(:,:,dd,:),1))),...
            'labels',lgnd,...
            'colors',repmat(smcol,[3 1]),'symbol','+');
        hnd = findobj(gca,'type','text');
        set(hnd,'interpreter','latex')
        a = get(get(gca,'children'),'children');
        %     t = get(a,'tag');
        for ia = 1:size(a,1)
            set(a(ia),'linewidth',2)
        end
        set(gca,'fontsize',12)
        ylabel(['$\log({\rm GWD})$'],'interpreter','latex','fontsize',16)
        title(ttl{dd},'interpreter','latex','fontsize',16)
        drawnow
    end
    
    
end
