%% Figure 3 of Bloxham et al. "Repulsion improve chemotaxis" 2026
%
% This file contains code to generate Figure 3 in Bloxham, Lee, Gore
% "Repulsion from slow-diffusing nutrients improves microbial chemotaxis
% towards moving sources" Nature Communication 2026.
%
% NOTE: Results may vary if sections are run out of order as some
% parameters change between being fixed and being varied across a range of
% values for different subplots of the figure.
%
% This code was written using MATLAB R2020a.
%
% (c) Blox Bloxham 2026


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


% Define standard parameters. Units: um, s, M

D1 = 1000; % um^2/s
D2 = 250;

s1 = 10^-12; % 1 fmol/sec
s2 = 10^-12;

k1 = 10^-7 * 10^-12; % (100 nM) * (10^-12 um^-3/M)
k2 = 10^-7 * 10^-12;

n1s = [6,4]; % [attractive,differential], assuming receptor hexamers
n2s = [0,2];

vmax = 6;
kdX = 0.004; % Constant setting the initial steepness of the increase in
             % drift velocity as a function of response function gradient
             % according to: v(dX/dr) = vmax * (dX/dr) / (|dX/dr| + kdX)
rp = 10;

vp = 1.5*vmax; % Initially fixed. Later varied


% Turbo colormap (included for earlier versions of MATLAB)
turbo10 = [
    0.1900    0.0718    0.2322
    0.2733    0.3801    0.8404
    0.2138    0.6589    0.9796
    0.0996    0.8904    0.7239
    0.4278    0.9942    0.3857
    0.7661    0.9463    0.2031
    0.9732    0.7468    0.2254
    0.9723    0.4456    0.1080
    0.8161    0.1846    0.0181
    0.5199    0.0276    0.0078];


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









%% %% %% Figure 3C, D, and E - Colormapped Phase Spaces            %% %% %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

D2s = [250,250,500];

s2s = [1,5,1] * s1;

n1s = [6,4; 6,4; 6,3.6];
n2s = 6 - n1s;


dt = 0.01;
tmax = 600;

rmax = cell(1,length(length(D2s)));
improv_sp = cell(1,length(D2s));

redcreamblue_lab = [
    30+70*cos(2.7*linspace(-1,1,1001));
    10+45*cos(4*linspace(-1,1,1001)+2.2);
    -5+65*cos(3*linspace(-1,1,1001)+0.4) .* (1-0.4*exp(-100*linspace(-1,1,1001).^2))]';

redcreamblue_rgb = lab2rgb(redcreamblue_lab);

%
atanslopes = [2.6,1.6,1.2];


f = figure(6);
clf('reset')

tic

%resfactor = [10,5,5]; % Resolution used in paper
resfactor = [2,1,1];

for condi = 1:length(D2s)
    D2 = D2s(condi);
    s2 = s2s(condi);
    
    vps = vmax*logspace(0,log10(20),resfactor(condi)*13+1)';
    rps = logspace(0,log10(100),resfactor(condi)*20+1);
    
    rmax{condi} = nan(length(vps),length(rps),length(n1s));
    
    for ni = 1:2
     for rpi = 1:length(rps)
      for vpi = 1:length(vps)
           
        rp = rps(rpi);
        vp = vps(vpi);
        n1 = n1s(condi,ni);
        n2 = n2s(condi,ni);

        if vz_(0,rp,vp,D1,D2,s1,s2,k1,k2,n1,n2,vmax,kdX) + vp < 0

            az = linspace(1000,rp,1024);
            avz = vz_(0,az,vp,D1,D2,s1,s2,k1,k2,n1,n2,vmax,kdX);

            [~,zind] = max(avz + vp < 0);
            zend = az(zind);

            maxrz = [[0,0,0.01;rp,zend-0.01,zend-0.01],nan(2,round(tmax/dt))];

            i0 = 3;
            iend = (round(tmax/dt)+2);

        else

            th = linspace(0,pi,1024);

            pvr = vr_(rp*sin(th),rp*cos(th),...
                vp,D1,D2,s1,s2,k1,k2,n1,n2,vmax,kdX);
            pvz = vz_(rp*sin(th),rp*cos(th),...
                vp,D1,D2,s1,s2,k1,k2,n1,n2,vmax,kdX) + vp;

            [~,thind] = max(pvr.*sin(th) + pvz.*cos(th) < 0);

            thend = th(thind);
            rend = rp*sin(thend);
            zend = rp*cos(thend);

            maxrz = [[rend;zend],nan(2,round(tmax/dt))];
            i0 = 1;
            iend = round(tmax/dt);

        end

        if any(isnan(maxrz(:,i0)))
            warning(['maxrz(:,i0) is NaN at [condi,n1i,n2i,rpi,vpi] = [',...
                num2str(condi),',',num2str(n1i),',',num2str(n2i),...
                ',',num2str(rpi),',',num2str(vpi),']']);
        else
            for i = i0:iend
                maxrz(:,i+1) = maxrz(:,i) - [
                    vr_(maxrz(1,i),maxrz(2,i),vp,D1,D2,s1,s2,k1,k2,n1,n2,vmax,kdX);
                    vz_(maxrz(1,i),maxrz(2,i),vp,D1,D2,s1,s2,k1,k2,n1,n2,vmax,kdX) + vp] * dt;
                if sqrt(sum((maxrz(:,i+1) - maxrz(:,i)).^2)) < 0.01
                    maxrz(:,i+1) = maxrz(:,i) + 0.01 * ...
                        (maxrz(:,i+1) - maxrz(:,i)) / ...
                        sqrt(sum((maxrz(:,i+1) - maxrz(:,i)).^2));
                end
                
                if maxrz(2,i+1) < -1000 || (i == iend && maxrz(2,i+1) < -200)
                    rmax{condi}(vpi,rpi,ni) = ...
                        maxrz(1,find(~any(isnan(maxrz)),1,'last'));
                    break
                elseif maxrz(1,i+1) < rp/1000
                    rmax{condi}(vpi,rpi,ni) = 0;
                    break
                elseif isnan(any(maxrz(:,i+1)))
                    if maxrz(2,i) < -200
                        rmax{condi}(vpi,rpi,ni) = maxrz(1,i);
                        warning(['NaN produced at [condi,ni,rpi,vpi] = [',...
                            num2str(condi),',',num2str(ni),...
                            ',',num2str(rpi),',',num2str(vpi),']']);
                        break
                    else
                        warning(['NaN produced and no value assigned at [condi,ni,rpi,vpi] = [',...
                            num2str(condi),',',num2str(ni),...
                            ',',num2str(rpi),',',num2str(vpi),']']);
                        break
                    end
                elseif i == iend
                    warning(['No value at [condi,ni,rpi,vpi] = [',...
                        num2str(condi),',',num2str(ni),...
                        ',',num2str(rpi),',',num2str(vpi),']']);
                end
            end
        end
        
      end %vp loop
     end %rp loop
    end %n1 loop
    
    for ni = 0:2
    
        if ni == 0
            improv_sp{condi} = subplot(length(D2s),3,[1,2]+3*(condi-1));

            improv_im = image(100*(rmax{condi}(:,:,2).^2-rmax{condi}(:,:,1).^2)./...
                rmax{condi}(:,:,1).^2,'CDataMapping','scaled');
            improv_im.AlphaData = ~isnan(improv_im.CData);
            
            colormap(gca,min(max(redcreamblue_rgb,0),1))
            
            cb = colorbar;

            improv_im.Parent.CLim = [-1,1] * ...
                max(abs(improv_im.CData),[],'all','omitnan');

            cb.Limits = ...
                [min(min(improv_im.CData,[],'all','omitnan'),0),...
                 max(max(improv_im.CData,[],'all','omitnan'),0)];
             
        else
            
            subplot(2*length(D2s),3,3+6*(condi-1)+3*(ni-1))
            
            image(pi*rmax{condi}(:,:,ni).^2,'CDataMapping','scaled');
            colormap(gca,flipud(bone))
            colorbar
        
        end

        xticks(1 + (log([1,2,5,10,20,50,100])-log(rps(1))) ...
            *(length(rps) - 1) / (log(rps(end)) - log(rps(1))))
        xticklabels({'1','2','5','10','20','50','100'})
        if ni == 0 || (ni == 2 && condi == length(D2s))
            xlabel('Particle Radius (um)')
        end

        yticks(1 + (log([1,2,3,5,10,20]*vmax)-log(vps(1))) ...
            *(length(vps) - 1) / (log(vps(end)) - log(vps(1))))
        yticklabels({'1','2','3','5','10','20'})
        if ni == 0
            ylabel('Particle Velocity (Relative)')
        else
            ylabel('Particle Vel. (rel.)')
        end

        set(gca,'YDir','normal','FontSize',12,'PlotBoxAspectRatio',...
            [log(rps(end))-log(rps(1)),log(vps(end))-log(vps(1)),1])

        if ni == 0
            xlim(xlim + [0.5,-0.5])
            ylim(ylim + [0.5,-0.5])
            
            title(['Increase in Int. Kernel',...
                ' (DA/DR = ',num2str(D1/D2s(condi)),', sR/sA = ',num2str(s2s(condi)/s1),')'])
        elseif ni == 1
            title('Attractive Kernel')
        else
            title('Differential Kernel')
        end
        
    
    end
    
    drawnow
    
    toc
    
end



% Recolor

[redcreamblue_uniform,colptlocs] = ...
    makeUniformColormap(min(max(redcreamblue_rgb,0),1),1001);
midptloc = colptlocs((size(colptlocs,1)+1)/2);
if midptloc < (size(colptlocs,1)+1)/2
    redcreamblue = redcreamblue_uniform(1:(2*midptloc-1),:);
else
end

atanslope = 6;

redcreamblue_converging = min(max([
    interp1(linspace(-pi/2,pi/2,size(redcreamblue,1)),redcreamblue(:,1),...
    atan(linspace(-atanslope,atanslope,1001)),'pchip');
    interp1(linspace(-pi/2,pi/2,size(redcreamblue,1)),redcreamblue(:,2),...
    atan(linspace(-atanslope,atanslope,1001)),'pchip');
    interp1(linspace(-pi/2,pi/2,size(redcreamblue,1)),redcreamblue(:,3),...
    atan(linspace(-atanslope,atanslope,1001)),'pchip')]',0),1);

improve_minmax = [
    min(improv_sp{1}.Children(end).CData,[],'all'), ...
    max(improv_sp{1}.Children(end).CData,[],'all');
    min(improv_sp{2}.Children(end).CData,[],'all'), ...
    max(improv_sp{2}.Children(end).CData,[],'all');
    min(improv_sp{3}.Children(end).CData,[],'all'), ...
    max(improv_sp{3}.Children(end).CData,[],'all')];


colorbarlimits = nan([size(improv_sp{1}.Parent.Children,1),2]);
for i = 1:2:length(colorbarlimits)
    colorbarlimits(i,:) = improv_sp{1}.Parent.Children(i).Limits;
end

%

for i = 1:3
    improv_sp{i}.Colormap = redcreamblue_converging;
    improv_sp{i}.CLim = [-1,1]*max(abs(improve_minmax),[],'all');
end

%

for i = 1:2:length(colorbarlimits)
    improv_sp{1}.Parent.Children(i).Limits = colorbarlimits(i,:);
end



%%

%% %% %% Figure 3F                                                 %% %% %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

D2 = 250;

s2s = [0.1,0.2,0.5,1,2,5,10,20,50] * s1;

n1s = [6, 4];
n2s = [0, 2];

rps = logspace(log10(2.5),log10(20),10*3+1);
vps = 1.09 * rps.^2;

dt = 0.01;
tmax = 600;

rmax = nan(length(s2s),length(rps),length(n1s));

for s2i = 1:length(s2s)
  for rpi = 1:length(rps)
    for ni = 1:length(n1s)

        s2 = s2s(s2i);
        rp = rps(rpi);
        vp = vps(rpi);
        n1 = n1s(ni);
        n2 = n2s(ni);

        if vz_(0,rp,vp,D1,D2,s1,s2,k1,k2,n1,n2,vmax,kdX) + vp < 0

            az = linspace(1000,rp,1024);
            avz = vz_(0,az,vp,D1,D2,s1,s2,k1,k2,n1,n2,vmax,kdX);

            [~,zind] = max(avz + vp < 0);
            zend = az(zind);

            maxrz = [[0,0,0.01;rp,zend-0.01,zend-0.01],nan(2,round(tmax/dt))];

            i0 = 3;
            iend = (round(tmax/dt)+2);

        else

            th = linspace(0,pi,1024);

            pvr = vr_(rp*sin(th),rp*cos(th),...
                vp,D1,D2,s1,s2,k1,k2,n1,n2,vmax,kdX);
            pvz = vz_(rp*sin(th),rp*cos(th),...
                vp,D1,D2,s1,s2,k1,k2,n1,n2,vmax,kdX) + vp;

            [~,thind] = max(pvr.*sin(th) + pvz.*cos(th) < 0);

            thend = th(thind);
            rend = rp*sin(thend);
            zend = rp*cos(thend);

            maxrz = [[rend;zend],nan(2,round(tmax/dt))];
            i0 = 1;
            iend = round(tmax/dt);

        end

        if any(isnan(maxrz(:,i0)))
            %warning(['maxrz(:,i0) is NaN at [s2i,rpi,ni] = [',...
            %    num2str(s2i),',',num2str(rpi),',',num2str(ni),']']);
        else
            for i = i0:iend
                maxrz(:,i+1) = maxrz(:,i) - [
                    vr_(maxrz(1,i),maxrz(2,i),vp,D1,D2,s1,s2,k1,k2,n1,n2,vmax,kdX);
                    vz_(maxrz(1,i),maxrz(2,i),vp,D1,D2,s1,s2,k1,k2,n1,n2,vmax,kdX) + vp] * dt;
                if sqrt(sum((maxrz(:,i+1) - maxrz(:,i)).^2)) < 0.01
                    maxrz(:,i+1) = maxrz(:,i) + 0.01 * ...
                        (maxrz(:,i+1) - maxrz(:,i)) / ...
                        sqrt(sum((maxrz(:,i+1) - maxrz(:,i)).^2));
                end

                if maxrz(2,i+1) < -1000 || (i == iend && maxrz(2,i+1) < -200)
                    rmax(s2i,rpi,ni) = ...
                        maxrz(1,find(~any(isnan(maxrz)),1,'last'));
                    break
                elseif maxrz(1,i+1) < rp/1000
                    rmax(s2i,rpi,ni) = 0;
                    break
                elseif any(isnan(maxrz(:,i+1))) && maxrz(2,i) < -50
                    rmax(s2i,rpi,ni) = maxrz(1,i);
                    %warning(['NaN produced at [s2i,rpi,ni] = [',...
                    %    num2str(s2i),',',num2str(rpi),',',num2str(ni),']']);
                    break
                elseif any(isnan(maxrz(:,i+1)))
                    %warning(['NaN produced and no value assigned at [s2i,rpi,ni] = [',...
                    %    num2str(s2i),',',num2str(rpi),',',num2str(ni),']']);
                    break
                elseif i == iend
                    %error(['No value at [s2i,rpi,ni] = [',...
                    %    num2str(s2i),',',num2str(rpi),',',num2str(ni),']']);
                end
            end
        end
    end
  end
end

figure(7)



clf('reset')

subplot(3,1,1:2)

colororder(turbo10)
semilogx(rps,100*(rmax(:,:,2).^2-rmax(:,:,1).^2)./rmax(:,:,1).^2,...
    'LineWidth',2)
ylim([-100,500])
yticks(-100:100:500)

xlim([min(rps),max(rps)])
xticks(1:100)

legend('Location','northeast')

subplot(3,1,3)

loglog(rps,vps/vmax,'k','LineWidth',1)
xlim([min(rps),max(rps)])
ylim([min(vps),max(vps)]/vmax)
xticks(1:100)
yticks([1,2,5,10,20,50,100,200,500,1000,2000,5000,10000])

%% %% %% Figure 3G                                                 %% %% %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

D1 = 1000;
D2s = D1./[1.01,1.3,2,4,8,16];

s2 = s1;

n1s = 6*[0:0.1:0.4,0.49,linspace(0.5,1,26)];
n2s = 6 - n1s;

rps = logspace(log10(2.35),log10(20),5*3+1);
vps = 1.09 * rps.^2;

dt = 0.01;
tmax = 1200;

rmax = nan(length(D2s),length(n1s),length(rps));

for D2i = 1:length(D2s)
  for ni = 1:length(n1s)
    for rpi = 1:length(rps)

        D2 = D2s(D2i);
        n1 = n1s(ni);
        n2 = n2s(ni);
        rp = rps(rpi);
        vp = vps(rpi);

        if vz_(0,rp,vp,D1,D2,s1,s2,k1,k2,n1,n2,vmax,kdX) + vp < 0

            az = linspace(1000,rp,1024);
            avz = vz_(0,az,vp,D1,D2,s1,s2,k1,k2,n1,n2,vmax,kdX);

            [~,zind] = max(avz + vp < 0);
            zend = az(zind);

            maxrz = [[0,0,0.01;rp,zend-0.01,zend-0.01],nan(2,round(tmax/dt))];

            i0 = 3;
            iend = (round(tmax/dt)+2);

        else

            th = linspace(0,pi,1024);

            pvr = vr_(rp*sin(th),rp*cos(th),...
                vp,D1,D2,s1,s2,k1,k2,n1,n2,vmax,kdX);
            pvz = vz_(rp*sin(th),rp*cos(th),...
                vp,D1,D2,s1,s2,k1,k2,n1,n2,vmax,kdX) + vp;

            [~,thind] = max(pvr.*sin(th) + pvz.*cos(th) < 0);

            thend = th(thind);
            rend = rp*sin(thend);
            zend = rp*cos(thend);

            maxrz = [[rend;zend],nan(2,round(tmax/dt))];
            i0 = 1;
            iend = round(tmax/dt);

        end

        if any(isnan(maxrz(:,i0)))
            %warning(['maxrz(:,i0) is NaN at [s2i,rpi,ni] = [',...
            %    num2str(s2i),',',num2str(rpi),',',num2str(ni),']']);
        else
            for i = i0:iend
                maxrz(:,i+1) = maxrz(:,i) - [
                    vr_(maxrz(1,i),maxrz(2,i),vp,D1,D2,s1,s2,k1,k2,n1,n2,vmax,kdX);
                    vz_(maxrz(1,i),maxrz(2,i),vp,D1,D2,s1,s2,k1,k2,n1,n2,vmax,kdX) + vp] * dt;
                if sqrt(sum((maxrz(:,i+1) - maxrz(:,i)).^2)) < 0.01
                    maxrz(:,i+1) = maxrz(:,i) + 0.01 * ...
                        (maxrz(:,i+1) - maxrz(:,i)) / ...
                        sqrt(sum((maxrz(:,i+1) - maxrz(:,i)).^2));
                end

                if maxrz(2,i+1) < -1000 || (i == iend && maxrz(2,i+1) < -200)
                    rmax(D2i,ni,rpi) = ...
                        maxrz(1,find(~any(isnan(maxrz)),1,'last'));
                    break
                elseif maxrz(1,i+1) < rp/1000
                    rmax(D2i,ni,rpi) = 0;
                    break
                elseif any(isnan(maxrz(:,i+1))) && maxrz(2,i) < -50
                    rmax(D2i,ni,rpi) = maxrz(1,i);
                    %warning(['NaN produced at [s2i,rpi,ni] = [',...
                    %    num2str(s2i),',',num2str(rpi),',',num2str(ni),']']);
                    break
                elseif any(isnan(maxrz(:,i+1)))
                    %warning(['NaN produced and no value assigned at [s2i,rpi,ni] = [',...
                    %    num2str(s2i),',',num2str(rpi),',',num2str(ni),']']);
                    break
                elseif i == iend
                    %error(['No value at [s2i,rpi,ni] = [',...
                    %    num2str(s2i),',',num2str(rpi),',',num2str(ni),']']);
                end
            end
        end
    end
  end
end

%
figure(8)
clf('reset')
colororder(turbo10([2,3,4:2:end],:))
plot(n1s/6,geomean(pi*rmax.^2,3),'LineWidth',2)

