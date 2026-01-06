


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


% Define standard parameters

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
vp = 1.5*vmax;


%% %% %% %% %% %%        Figure 3A,B: Simulations         %% %% %% %% %% %%

z0 = -400;
r0 = -17.5*(0.5:16.5);

dt = 0.01;
tmax = 3600;

rz_t = cat(2,...
    repmat(cat(1,reshape(r0,1,1,[]),repmat(z0,1,1,length(r0))),[1,1,1,2]),...
    nan(2,round(tmax/dt),length(r0),2));

time2particle = nan(length(r0),2);

for ni = 1:2
    n1 = n1s(ni);
    n2 = n2s(ni);
    for r0i = 1:length(r0)
        for t = 1:round(tmax/dt)
            if ~isnan(time2particle(r0i,ni))
                break
            else
                rz_t(:,t+1,r0i,ni) = rz_t(:,t,r0i,ni) + ...
                    [vr_(rz_t(1,t,r0i,ni),rz_t(2,t,r0i,ni)+vp*(t-1)*dt,...
                         vp,D1,D2,s1,s2,k1,k2,n1,n2,vmax,kdX);
                     vz_(rz_t(1,t,r0i,ni),rz_t(2,t,r0i,ni)+vp*(t-1)*dt,...
                         vp,D1,D2,s1,s2,k1,k2,n1,n2,vmax,kdX)] * dt;
            end
            if sum((rz_t(:,t+1,r0i,ni) - [0;-vp*(t-1)*dt]).^2) < rp^2
                time2particle(r0i,ni) = t*dt;
                break
            end
        end
    end
end


%% %% %% %% %% %%      Figure 3A,B: Figure Generation     %% %% %% %% %% %%

rplotlim = [-120,50];
zplotlim = [-450,50];
rs_plot = linspace(rplotlim(1),rplotlim(2),151);
zs_plot = linspace(zplotlim(1),zplotlim(2),526)';

ts_plot = [0,22.5:6:42];

figure(5)
clf('reset')

for tpi = 1:length(ts_plot)
    t = round(ts_plot(tpi)/dt) + 1;
    for ni = 1:2
        if ni == 1
            subplot(2,length(ts_plot),tpi)
        else
            subplot(2,length(ts_plot),length(ts_plot)+tpi)
        end
        
        cmapsteps = [(5/4),(1/3)];
        cmapstep = cmapsteps(ni);
        
        if ni == 2
            imCData = X_(rs_plot,zs_plot+vp*(t-1)*dt,vp,D1,D2,s1,s2,k1,k2,n1s(ni),n2s(ni));
            contourf(rs_plot,zs_plot,imCData,[-inf,0:(1/3):(6)],'LineColor','none');
            
            caxis([0,6])
            
            cmap = makeUniformColormap([paleblue;attractantgreen],18);
        else
            imCData = X_(rs_plot,zs_plot+vp*(t-1)*dt,vp,D1,D2,s1,s2,k1,k2,n1s(ni),n2s(ni));
            contourf(rs_plot,zs_plot,imCData,[-inf,0:(1/2):(24)],'LineColor','none');
            
            caxis([0,24])
            
            cmap = makeUniformColormap([paleblue;attractantgreen],48);
            
            cmap = [cmap(1:2,:);...
                [repelem(cmap(4:2:end,1),2),repelem(cmap(4:2:end,2),2),repelem(cmap(4:2:end,3),2)]];
            
        end
        
        colorbar
        
        colormap(gca,cmap)
        
        hold on
        
        rectangle('Position',[-rp,-rp-vp*(t-1)*dt,2*rp,2*rp],...
            'Curvature',[1,1],'FaceColor','w','EdgeColor','k')
        
        for r0i = length(r0):-1:1
            
            if ~isnan(rz_t(1,t,r0i,ni))
                plot(rz_t(1,1:t,r0i,ni),rz_t(2,1:t,r0i,ni),'k','LineWidth',1)
                plot(rz_t(1,t,r0i,ni),rz_t(2,t,r0i,ni),'k.','MarkerSize',10)
                
                plot(-rz_t(1,1:t,r0i,ni),rz_t(2,1:t,r0i,ni),'k','LineWidth',1)
                plot(-rz_t(1,t,r0i,ni),rz_t(2,t,r0i,ni),'k.','MarkerSize',10)
            else
                col = [1,1,0];
                
                plot(rz_t(1,1:t,r0i,ni),rz_t(2,1:t,r0i,ni),'Color',col,'LineWidth',1)
                plot(rz_t(1,find(~isnan(rz_t(1,:,r0i,ni)),1,'last'),r0i,ni),...
                    rz_t(2,find(~isnan(rz_t(1,:,r0i,ni)),1,'last'),r0i,ni),'.',...
                    'Color',col,'MarkerSize',10)
                
                plot(-rz_t(1,1:t,r0i,ni),rz_t(2,1:t,r0i,ni),'Color',col,'LineWidth',1)
                plot(-rz_t(1,find(~isnan(rz_t(1,:,r0i,ni)),1,'last'),r0i,ni),...
                    rz_t(2,find(~isnan(rz_t(1,:,r0i,ni)),1,'last'),r0i,ni),'.',...
                    'Color',col,'MarkerSize',10)
            end
            
        end
        
        set(gca,'Ydir','normal','FontSize',16)
        axis equal tight
        
        xlim([min(rs_plot(:)),max(rs_plot(:))])
        ylim([min(zs_plot(:)),max(zs_plot(:))])
        
        drawnow
        
    end
end



%%



















