function ProcessDiskSims
% Load, plot and post analyse disk simulations

close all
addpath('Colormaps/');

% Video of cluster dynamics with rotation frequencies?
mk_vid = 0;

% Turn off irrelevant error message from loading the data
warning('off', 'MATLAB:dispatcher:UnresolvedFunctionHandle');

%%%%%%%%%%%%%%%%%%%%%%%% LOAD DATA %%%%%%%%%%%%%%%%%%%%%%%%
%%% Path
dat_path = 'Example/'; 

%%% Filenames
pmdata = 'Parameters_bare';
fndata = 'Simdata_bare';

%%% Load data
curr_data = load([dat_path,pmdata,'.mat']);
curr_Simdata = load([dat_path,fndata,'.mat']);

% How many disks?
N = curr_data.N;

% ODE solution
y = curr_data.y;

% Simulation times of that ODE model
t = curr_data.t;
       
% Disk rotation frequencies
omega_all = curr_Simdata.Omega_all;

% Domain size
L = 1.5*curr_data.L;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot all tracks
figure(1);
scatter(y(:,1),y(:,2),'filled');
hold on
for i = 2:N
    scatter(y(:,2*i-1),y(:,2*i),'filled');    
end
% Red squares to label starting points of trajectories
scatter(y(1,1:2:(2*N-1)),y(1,2:2:(2*N)),40,'r','square','filled');
axis equal
axis([-L,L,-L,L]);   
drawnow;

f2 = figure(2);
set(gcf,'color','w');
if mk_vid == 1
    v = VideoWriter([dat_path,fndata,'_PostVid'],'MPEG-4');
    v.Quality = 40;
    v.FrameRate = 20;
    open(v)        
end            

% Rotation frequencies of single embryos
omega_all_plot = [];        

    for k = 1:1:length(t) 
        % Load angular spinning frequency and positions from stored data (stored as ang. freq. w = 2*pi*f)
        omega_curr = omega_all(((k-1)*N+1):(k*N))/(2*pi); % (FREQUENCY IN HZ)       

        if k == 1 % First frame to plot
            set(0,'CurrentFigure',f2);
            sc = scatter(y(k,1:2:(2*N-1)),y(k,2:2:(2*N)),1.5e2,omega_curr','filled',...
                    'MarkerEdgeColor','k');
            hold on
            set(gca,'FontSize',15);

            axis equal
            axis([-L,L,-L,L]);           
            
            box on
            set(gca,'XTick',[])
            set(gca,'YTick',[])
            axis off

            caxis([0,0.7]); 
            cb = colorbar;               
            set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.1, 0.55, 0.5, 0.7]);                        
        else % Following loops just update axis content
            sc.XData = y(k,1:2:(2*N-1));
            sc.YData = y(k,2:2:(2*N));
            sc.CData = omega_curr';               
            caxis([0,0.7]); 
            cb = colorbar; 
         end                  
                         
        drawnow;              
        
        if (mk_vid == 1)
            set(0,'CurrentFigure',f2);
            frame = getframe(gcf);
            writeVideo(v,frame);
        end    
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot rotation frequencies of all embryos
figure(3)
title('Particle rotation frequencies')
plot(t(1:size(omega_all_plot,1)),omega_all_plot,'LineWidth',1.5);

if mk_vid == 1
    close(v);
end

end % Main

