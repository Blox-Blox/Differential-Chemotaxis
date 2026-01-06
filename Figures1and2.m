


%% Define some useful functions and values

% Concentration and derivatives thereof

conc = @(r,z,v,D) ( ...
    1./(4*pi*D.*sqrt(r.^2+z.^2)) .* exp( -v./(2*D).*(sqrt(r.^2+z.^2)-z) ) );

dc_dr = @(r,z,v,D) (...
    - r.*( 1./(r.^2+z.^2) + v./(2*D.*sqrt(r.^2+z.^2)) ) .* conc(r,z,v,D) );

dc_dz = @(r,z,v,D) (...
    ( - z./(r.^2+z.^2) + v./(2*D) - v.*z./(2*D.*sqrt(r.^2+z.^2)) ) ...
    .* conc(r,z,v,D) );


% Response function and derivatives thereof

X_ = @(r,z,vp,D1,D2,s1,s2,k1,k2,n1,n2) ( ...
      n1*log( 1 + s1.*conc(r,z,vp,D1)./k1 ) ...
    - n2*log( 1 + s2.*conc(r,z,vp,D2)./k2 ) );

dX_dpos = @(c1,c2,dc1_dpos,dc2_dpos,k1,k2,n1,n2) ( ...
      n1 .* dc1_dpos ./ (k1 + c1) ...
    - n2 .* dc2_dpos ./ (k2 + c2) );

dX_dr = @(r,z,vp,D1,D2,s1,s2,k1,k2,n1,n2) ( dX_dpos( ...
    s1.*conc(r,z,vp,D1),s2.*conc(r,z,vp,D2),...
    s1.*dc_dr(r,z,vp,D1),s2.*dc_dr(r,z,vp,D2),...
    k1,k2,n1,n2) );

dX_dz = @(r,z,vp,D1,D2,s1,s2,k1,k2,n1,n2) ( dX_dpos( ...
    s1.*conc(r,z,vp,D1),s2.*conc(r,z,vp,D2),...
    s1.*dc_dz(r,z,vp,D1),s2.*dc_dz(r,z,vp,D2),...
    k1,k2,n1,n2) );

dX_abs = @(r,z,vp,D1,D2,s1,s2,k1,k2,n1,n2) ( ...
    sqrt(  dX_dr(r,z,vp,D1,D2,s1,s2,k1,k2,n1,n2).^2 ...
         + dX_dz(r,z,vp,D1,D2,s1,s2,k1,k2,n1,n2).^2 ) );

     
% Microbe velocity (in the stationary frame)
     
v_abs = @(r,z,vp,D1,D2,s1,s2,k1,k2,n1,n2,vmax,kdX) ( ...
    vmax .* dX_abs(r,z,vp,D1,D2,s1,s2,k1,k2,n1,n2) ./ ...
    (kdX + dX_abs(r,z,vp,D1,D2,s1,s2,k1,k2,n1,n2)) );

vr_ = @(r,z,vp,D1,D2,s1,s2,k1,k2,n1,n2,vmax,kdX) ( ...
    v_abs(r,z,vp,D1,D2,s1,s2,k1,k2,n1,n2,vmax,kdX) .* ...
    dX_dr(r,z,vp,D1,D2,s1,s2,k1,k2,n1,n2) ./ ...
    dX_abs(r,z,vp,D1,D2,s1,s2,k1,k2,n1,n2) );

vz_ = @(r,z,vp,D1,D2,s1,s2,k1,k2,n1,n2,vmax,kdX) ( ...
    v_abs(r,z,vp,D1,D2,s1,s2,k1,k2,n1,n2,vmax,kdX) .* ...
    dX_dz(r,z,vp,D1,D2,s1,s2,k1,k2,n1,n2) ./ ...
    dX_abs(r,z,vp,D1,D2,s1,s2,k1,k2,n1,n2) );


% Endpoints for colormaps

paleblue = [0.95,0.95,1];
attractantgreen = [55,168,40]/255;
repellentred = [185,118,105]/255;


%% Define standard parameters

D1 = 1000;
D2 = 250;

s1 = 10^-12; % 1 fmol/sec
s2 = 10^-12;

k1 = 10^-7 * 10^-12; % (100 nM) * (10^-12 um^-3/M)
k2 = 10^-7 * 10^-12;

n1s = [6,4];
n2s = [0,2];

vmax = 6;
kdX = 0.004;

rp = 10;



%% %% %% %% %% %% %% %%          Figure 1B          %% %% %% %% %% %% %% %%

% Assume microbe swims at max speed

vb = vmax;

% Define example point and vp/vb ratio

r0 = -450;
z0 = -300;
vp0 = 0.75*vb;

% Define range of vp/vb ratios (for bottom panel)

vp = vb*linspace(0,1.25,126);

% Define range for plots

rplotlim = [-700,500];
zplotlim = [-800,400];
rs_plot = linspace(rplotlim(1),rplotlim(2),1201);
zs_plot = linspace(zplotlim(1),zplotlim(2),1201)';



% Calculate optimal and actual directions of swimming

i = 1j;

th_opt = -i.*log( ...
    (-i*vp.*r0-sqrt((vb-vp).*(vb+vp).*r0.^2+vb.^2.*z0.^2)) ./ ...
    (vb.*(r0-i*z0)));

th_opt0 = -i.*log( ...
    (-i*vp0.*r0-sqrt((vb-vp0).*(vb+vp0).*r0.^2+vb.^2.*z0.^2)) ./ ...
    (vb.*(r0-i*z0)));

th_actual = atan2(dc_dz(r0,z0,vp,D1),dc_dr(r0,z0,vp,D1));



% Plot the results

figure(1)
clf('reset')

subplot(2,1,1)

% Plot concentration as color map / contour plot

imCData = log10((s1/k1)*conc(rs_plot,zs_plot,vp0,D1));
contourf(rs_plot,zs_plot,imCData,[-inf,-1:(1/6):2],'LineColor','none');
caxis([-1,2])
cmap = makeUniformColormap([paleblue;attractantgreen],6*caxis*[-1;1]);
colormap(gca,cmap)

hold on
        
% Plot Particle (as a square with corners rounded in circle) 

rectangle('Position',[-rp,-rp,2*rp,2*rp],...
    'Curvature',[1,1],'FaceColor','w','EdgeColor','k')

%

plot(r0 + [0,300*cos(real(th_opt0))],...
    z0 + [0,300*sin(real(th_opt0))],...
    'LineWidth',3)
plot(r0 + [0,300*dc_dr(r0,z0,vp0,D1)./sqrt(dc_dr(r0,z0,vp0,D1).^2 + dc_dz(r0,z0,vp0,D1).^2)],...
    z0 + [0,300*dc_dr(r0,z0,vp0,D1)./sqrt(dc_dr(r0,z0,vp0,D1).^2 + dc_dz(r0,z0,vp0,D1).^2)],...
    'LineWidth',3)


set(gca,'Ydir','normal','FontSize',16)
axis equal tight

xlim([min(rs_plot(:)),max(rs_plot(:))])
ylim([min(zs_plot(:)),max(zs_plot(:))])


subplot(2,1,2)
plot(vp(imag(th_opt) < 1e-10)/vb,(180/pi)*(real(th_opt(imag(th_opt) < 1e-10)) - th_opt(1)))
hold on
plot(vp/vb,(180/pi)*(th_actual - th_opt(1)))
xlim([0,1.25])

drawnow



%% %% %% %% %% %% %% %%          Figure 1C          %% %% %% %% %% %% %% %%



rplotlim = [-600,100];
zplotlim = [-650,100];
rs_plot = linspace(rplotlim(1),rplotlim(2),401);
zs_plot = linspace(zplotlim(1),zplotlim(2),506)';


figure(2)
clf('reset')

for sp = 1:2
    
    if sp == 1
        
        [rs_calc,zs_calc] = meshgrid(...
            [-540:50:0,0],...
            -500:50:0);
        
        vps = repmat([0.75,1.25]*vb,numel(rs_calc),1);
        
    else
        
        [rs_calc,zs_calc] = meshgrid(...
            -399:5:0,...
            -599:5:0);
        
        vps = repmat(vp,numel(rs_calc),1);
        
    end
    
    calc_cut = log10((s1/k1)*conc(rs_calc,zs_calc,vp0,D1)) < -5/6;
    rs_calc(calc_cut) = nan;
    zs_calc(calc_cut) = nan;
    
    
    
    rs_calc = repmat(reshape(rs_calc,[],1),1,size(vps,2));
    zs_calc = repmat(reshape(zs_calc,[],1),1,size(vps,2));
    
    
    
    th_direct_calc = atan2(dc_dz(rs_calc,zs_calc,0,D1),dc_dr(rs_calc,zs_calc,0,D1));
    th_actual_calc = atan2(dc_dz(rs_calc,zs_calc,vps,D1),dc_dr(rs_calc,zs_calc,vps,D1));
    
    vr_opt_calc = (vps.*rs_calc.*zs_calc - rs_calc.*sqrt(-vps.^2.*rs_calc.^2+vb.^2.*(rs_calc.^2 + zs_calc.^2))) ./ ...
    (rs_calc.^2 + zs_calc.^2);
    vz_opt_calc = (-vps.*rs_calc.^2 - zs_calc.*sqrt(-vps.^2.*rs_calc.^2+vb.^2.*(rs_calc.^2 + zs_calc.^2))) ./ ...
        (rs_calc.^2 + zs_calc.^2);
    th_opt_calc = asin(real(vz_opt_calc./vb));
    
    th_opt_calc(real(th_opt_calc) == -pi/2 | (zs_calc > 0 & vps > vb)) = nan;
    th_opt_calc(any(isnan(th_opt_calc),2),:) = nan;
    
    
    if sp == 1
        
        calc_cut = ~any(isnan(th_opt_calc) | isnan(rs_calc),2);
        
        rs_calc = rs_calc(calc_cut,:);
        zs_calc = zs_calc(calc_cut,:);
        th_opt_calc = th_opt_calc(calc_cut,:);
        th_actual_calc = th_actual_calc(calc_cut,:);
        
    end
    
    
    subplot(2,1,sp)
    
    if sp == 1
        
        imCData = log10((s1/k1)*conc(rs_plot,zs_plot,vp0,D1));
        contourf(rs_plot,zs_plot,imCData,[-inf,-1:(1/6):2],'LineColor','none');
        caxis([-1,2])
        cmap = makeUniformColormap([paleblue;attractantgreen],6*caxis*[-1;1]);
        colormap(gca,cmap)
        colorbar


        hold on

        rectangle('Position',[-rp,-rp,2*rp,2*rp],...
            'Curvature',[1,1],'FaceColor','w','EdgeColor','k')

        
        quiver(rs_calc(:,1),zs_calc(:,1),cos(th_opt_calc(:,1)),sin(th_opt_calc(:,1)),...
            0.5,'Color',[0.6,0.6,0.6],'LineWidth',3)
        quiver(rs_calc(:,1),zs_calc(:,1),cos(th_actual_calc(:,1)),sin(th_actual_calc(:,1)),...
            0.5,'k','LineWidth',1)
        
        set(gca,'Ydir','normal','FontSize',16)
        axis equal tight

        xlim([min(rs_plot(:)),max(rs_plot(:))])
        ylim([min(zs_plot(:)),max(zs_plot(:))])
        
    else
        
        plot(vps(1,:)/vb,smooth((180/pi)*mean(th_opt_calc - th_direct_calc,1,'omitnan'),11),...
            'Color',[0.6,0.6,0.6],'LineWidth',2)
        hold on

        plot(vps(1,:)/vb,smooth((180/pi)*(mean(th_opt_calc - th_direct_calc,1,'omitnan') + ...
            std(th_opt_calc - th_direct_calc,[],1,'omitnan')),11),'Color',[0.6,0.6,0.6],'LineWidth',1)
        plot(vps(1,:)/vb,smooth((180/pi)*(mean(th_opt_calc - th_direct_calc,1,'omitnan') - ...
            std(th_opt_calc - th_direct_calc,[],1,'omitnan')),11),'Color',[0.6,0.6,0.6],'LineWidth',1)

        plot(vps(1,:)/vb,smooth((180/pi)*mean(th_actual_calc - th_direct_calc,1,'omitnan'),11),...
            'k','LineWidth',2)

        plot(vps(1,:)/vb,smooth((180/pi)*(mean(th_actual_calc - th_direct_calc,1,'omitnan') + ...
            std(th_actual_calc - th_direct_calc,[],1,'omitnan')),11),'k','LineWidth',1)
        plot(vps(1,:)/vb,smooth((180/pi)*(mean(th_actual_calc - th_direct_calc,1,'omitnan') - ...
            std(th_actual_calc - th_direct_calc,[],1,'omitnan')),11),'k','LineWidth',1)

        yticks(-105:15:45)
        xticks(0:0.25:1.25)
        
        
    end
    
    

end

drawnow



%% %% %% %% %% %% %% %%         Figure 2B,E         %% %% %% %% %% %% %% %%


rplotlim = [-650,650];
zplotlim = [-500,1300];
rs_plot = linspace(rplotlim(1),rplotlim(2),651);
zs_plot = linspace(zplotlim(1),zplotlim(2),901)';

rz_exgradpoint = [-200,100]; % For the microbe in Figure 2B

n1 = n1s(2);
n2 = n2s(2);

vp = 1.5*vmax;  

figure(3)
clf('reset')

for spi = 1:3
    
    subplot(1,3,spi)
    
    if spi == 1

        imCData = X_(rs_plot,zs_plot,vp,D1,D2,s1,s2,k1,k2,n1,0);
        contourf(rs_plot,zs_plot,imCData,[-inf,0:0.5:14],'LineColor','none');
        caxis([0,14])
        cmap = makeUniformColormap([paleblue;attractantgreen],2*caxis*[-1;1]);
        
    elseif spi == 2
        
        imCData = -X_(rs_plot,zs_plot,vp,D1,D2,s1,s2,k1,k2,0,n2);
        contourf(rs_plot,zs_plot,imCData,[-inf,0:0.5:9],'LineColor','none');
        caxis([0,9])
        cmap = makeUniformColormap([paleblue;repellentred],2*caxis*[-1;1]);
        
    else
        
        imCData = X_(rs_plot,zs_plot,vp,D1,D2,s1,s2,k1,k2,n1,n2);
        contourf(rs_plot,zs_plot,imCData,[-inf,-1:(1/6):6],'LineColor','none');
        caxis([-1,6])
        
        [cmap,inds] = makeUniformColormap([attractantgreen;paleblue;repellentred],10000);
        
        cdivs = 3;
        
        cmap = flipud([
            repelem(interp1(1:10000,cmap(:,1),...
                linspace(1,round(inds(2)*(max(caxis)+2)/max(caxis)),cdivs*(max(caxis)+2))),2);
            repelem(interp1(1:10000,cmap(:,2),...
                linspace(1,round(inds(2)*(max(caxis)+2)/max(caxis)),cdivs*(max(caxis)+2))),2);
            repelem(interp1(1:10000,cmap(:,3),...
                linspace(1,round(inds(2)*(max(caxis)+2)/max(caxis)),cdivs*(max(caxis)+2))),2)]');
            
        cmap = cmap([1:2:(4*cdivs),(4*cdivs+1):end],:);
        
        
    end
    
    colormap(gca,cmap)
    cb = colorbar;
    
    hold on
        
    rectangle('Position',[-rp,-rp,2*rp,2*rp],...
        'Curvature',[1,1],'FaceColor','w','EdgeColor','k')
    
    set(gca,'Ydir','normal','FontSize',16)
    axis equal tight
    
    if spi == 1
        plot(rz_exgradpoint(1) + (5e10)*[0,...
            dc_dr(rz_exgradpoint(1),rz_exgradpoint(2),vp,D1)],...
            rz_exgradpoint(2) + (5e10)*[0,...
            dc_dz(rz_exgradpoint(1),rz_exgradpoint(2),vp,D1)],'k-')
    elseif spi == 2
        plot(rz_exgradpoint(1) + (5e10)*[0,...
            dc_dr(rz_exgradpoint(1),rz_exgradpoint(2),vp,D2)],...
            rz_exgradpoint(2) + (5e10)*[0,...
            dc_dz(rz_exgradpoint(1),rz_exgradpoint(2),vp,D2)],'k-')
    end

    xlim([min(rs_plot(:)),max(rs_plot(:))])
    ylim([min(zs_plot(:)),max(zs_plot(:))])
    
end

drawnow



%% %% %% %% %% %% %% %%         Figure 1D and 2F         %% %% %% %% %% %% %% %%


% Define simulation parameters

rz0 = [-79.5; -400];

dt = 0.01;
tmax = 240;

rz_t = repmat([rz0,nan(2,round(tmax/dt))],[1,1,2]);
time2particle = nan(1,2);


% Run the simulation

for ni = 1:2
    n1 = n1s(ni);
    n2 = n2s(ni);
    
    for t = 1:round(tmax/dt)
        if ~isnan(time2particle(ni))
            break
        else
            rz_t(:,t+1,ni) = rz_t(:,t,ni) + ...
                [vr_(rz_t(1,t,ni),rz_t(2,t,ni)+vp*(t-1)*dt,...
                     vp,D1,D2,s1,s2,k1,k2,n1,n2,vmax,kdX);
                 vz_(rz_t(1,t,ni),rz_t(2,t,ni)+vp*(t-1)*dt,...
                     vp,D1,D2,s1,s2,k1,k2,n1,n2,vmax,kdX)] * dt;
        end
        if sum((rz_t(:,t+1,ni) - [0;-vp*(t-1)*dt]).^2) < rp^2
            time2particle(ni) = t*dt;
            break
        end
    end
end

% disp(time2particle)



% Plot Comparison of Differential and Purely Attractive Trajectories

ts_plot = 35.01 + 4.65*(-2:2);
rcrop = [-95,35];
zcrop = [-415,-185];
rs_plot_traj = linspace(rcrop(1),rcrop(2),151);
zs_plot_traj = linspace(zcrop(1),zcrop(2),526)';

figure(4)
clf('reset')

for tpi = length(ts_plot):-1:1
    t = round(ts_plot(tpi)/dt) + 1;
    for ni = 1:2
        if ni == 1
            subplot(2,length(ts_plot),tpi)
            
            im = image(rs_plot_traj,zs_plot_traj,...
                log10((s1/k1)*conc(rs_plot_traj,zs_plot_traj+vp*(t-1)*dt,vmax,D1)),...
                'CDataMapping','scaled');
            
            caxis([-1,2])
            cmap = makeUniformColormap([paleblue;attractantgreen],6*caxis*[-1;1]);
        else
            subplot(2,length(ts_plot),length(ts_plot)+tpi)
            
            
            imCData = ...
                X_(rs_plot_traj,zs_plot_traj+vp*(t-1)*dt,vp,D1,D2,s1,s2,k1,k2,n1s(ni),n2s(ni));
            contourf(rs_plot_traj,zs_plot_traj,imCData,[-inf,0:(1/6):6],'LineColor','none');
            
            caxis([0,6])
            cmap = makeUniformColormap([paleblue;attractantgreen],18);
        end
        
        colormap(gca,cmap)
        colorbar
        
        
        hold on
        
        rectangle('Position',[-rp,-rp-vp*(t-1)*dt,2*rp,2*rp],...
            'Curvature',[1,1],'FaceColor','w','EdgeColor','k')
        
        
        if ~isnan(rz_t(1,t,ni))
            plot(rz_t(1,1:t,ni),rz_t(2,1:t,ni),'k','LineWidth',1)
            plot(rz_t(1,t,ni),rz_t(2,t,ni),'k.')

        else
            plot(rz_t(1,1:t,ni),rz_t(2,1:t,ni),'y','LineWidth',1)
            plot(rz_t(1,find(~isnan(rz_t(1,:,ni)),1,'last'),ni),...
                rz_t(2,find(~isnan(rz_t(1,:,ni)),1,'last'),ni),'y.')

        end
        
        

        set(gca,'Ydir','normal','FontSize',16)
        axis equal tight
        
        xlim(rcrop)
        ylim(zcrop)
        
        drawnow
        
    end
end



% Plot Purely Attractive Trajectory

rcrop_attr = [-300,70];
zcrop_attr = [-485,70];
rs_plot_traj_attr = linspace(rcrop_attr(1),rcrop_attr(2),151);
zs_plot_traj_attr = linspace(zcrop_attr(1),zcrop_attr(2),526)';
ts_plot_attr = 15.5*(0:1:3);

figure(5)
clf('reset')

for tpi = length(ts_plot_attr):-1:1
    t = round(ts_plot_attr(tpi)/dt) + 1;
    for ni = 1
        subplot(1,length(ts_plot_attr),tpi)
        
        contourf(rs_plot_traj_attr,zs_plot_traj_attr,...
            log10((s1/k1)*conc(rs_plot_traj_attr,zs_plot_traj_attr+vp*(t-1)*dt,vmax,D1)),...
            [-inf,-1:(1/6):2],'LineColor','none')

        caxis([-1,2])
        cmap = makeUniformColormap([paleblue;attractantgreen],6*caxis*[-1;1]);
        
        colormap(gca,cmap)
        colorbar
        hold on
        
        rectangle('Position',[-rp,-rp-vp*(t-1)*dt,2*rp,2*rp],...
            'Curvature',[1,1],'FaceColor','w','EdgeColor','k')
        
        if ~isnan(rz_t(1,t,ni))
            plot(rz_t(1,1:t,ni),rz_t(2,1:t,ni),'k','LineWidth',1)
            plot(rz_t(1,t,ni),rz_t(2,t,ni),'k.')
            
        else
            plot(rz_t(1,1:t,ni),rz_t(2,1:t,ni),'y','LineWidth',1)
            plot(rz_t(1,find(~isnan(rz_t(1,:,ni)),1,'last'),ni),...
                rz_t(2,find(~isnan(rz_t(1,:,ni)),1,'last'),ni),'y.')
        end
        
        set(gca,'Ydir','normal','FontSize',16)
        axis equal tight
        
        xlim(rcrop_attr)
        ylim(zcrop_attr)
        
        drawnow
        
    end
end





