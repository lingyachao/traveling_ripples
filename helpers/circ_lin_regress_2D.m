function [direction,spatial_frequency,sl,Rsquared] = circ_lin_regress_2D(circularV,linearV,clusteridx,fig)
% Code from Josh Jacobs lab https://github.com/jacobslab, download: 07092019
% Modified by Jon Kleen UCSF 2019
%%
% circ_lin_regress_2D is a function used to fit a 2D linear variables to a
% circular variable,handles wrapping-up cases. Can fit phase traveling wave
% from linear coordinates. ***Please note the direction of wave propagation will be 180 degree from this
% number as wave propagates to phase descending direction.
%
% For method details, see Zhang H, Watrous AJ, Patel A, Jacobs J. Theta and
% alpha oscillations are traveling waves in the human neocortex. Neuron. 2018.
%
% INPUTS
% circularV- (phase, in radians) is a N by 1 array, which we are trying to fit. 
% linearV- (coordinates) is a N by 2 matrix convey the linear x and y coordinate variables
    % Note: Handles wrapping-up cases on its own.
% clusteridx- index of which electrodes to use for this cluster (if empty, will use all available)
% fig- if 1 then generates 2 plots, the first showing the parameter space
% with residuals, and the other showing the fitting plan. If 0, no figures.
%
% OUTPUTS
% direction (?in radians?) of circular variable ascending with linear variables.
% spatial_frequency is the change rate of circular variable (phase) to
       % linear variables (distance) at the ascending direction
       % (units: radian / linear unit)
% sl is the slopes of phase change relative to each column in linearV.
% R-squared denotes how much variance of the circular variable is explained
    % by the regression model   ([Note - therefore, PGD in Zhang et al 2018?]). 
    % For statistical significance,  perform a permutation procedure: 
    %    shuffle circularV, & iterate x 1000 for 95% CI ([they take median PGD over all timepoints in a trial, and check vs shuffled])
%%

if nargin<2; disp('Incomplete input'); end
if ~exist('clusteridx','var'); clusteridx=ones(length(circularV)); end
if ~exist('fig','var'); fig=[0 0]; end
if size(circularV,1)==1; circularV=circularV'; end %expects phase vector as 1 column

circularV=circularV(clusteridx);
% linearV=linearV(clusteridx,:);

pos_x = linearV(:,1);
pos_y = linearV(:,2);
phase = mod(circularV,2*pi);
%helper function to calculate Residual resultant vector length between
%fit and actual phase
myfun = @(slope1,slope2)sqrt((sum(cos(phase-slope1*pos_x-slope2*pos_y)/length(phase)).^2 + (sum(sin(phase-(slope1*pos_x)-slope2*pos_y))/length(phase)).^2));


% arrange range and steplength for parameter space. angle ranges 2pi
% and spatial frequency range from 0 to 18 [90deg for LLNL microgrid -jk] degree per unit. the upper
% limit of spatial frequency is depend on the spatial nyquist
% frequency. You may change the steplength to anything but watch out
% for computation time.
steplength=1;
angle_range=pi*(0:steplength:360)/180; 
spatial_frequency_range=pi*(0:steplength:90)/180; %Zhang et al had 10mm spacing and was (0:.1:18)*pi/180;

[angleMatrix,spatial_frequency_Matrix]=meshgrid(angle_range,spatial_frequency_range); % make it to a matrix for arrayfun

% transfer angle and spatial_frequency into 2d linear slopes
slope1_Matrix=cos(angleMatrix).*spatial_frequency_Matrix; % slope for pos_x
slope2_Matrix=sin(angleMatrix).*spatial_frequency_Matrix; % slope for pos_y
Residual_resultant_vector_length=arrayfun(myfun,slope1_Matrix,slope2_Matrix); % calculate the resultant_vector_length for each possible slope combos
[row,column]=find(Residual_resultant_vector_length==max(max(Residual_resultant_vector_length))); % find the largest resultant_vector_length


%get the direction and spatial_frequency. If running traveling
%wave analysis, the direction should be flipped 180 degrees as waves
%propagates to the phase descending directions.
% jk 082019: added (1) since sometimes would arrive on two points 0 and 2pi simultaneously
direction=angleMatrix(row(1),column(1)); % find the best fit propagtion direction
spatial_frequency=spatial_frequency_Matrix(row(1),column(1)); % find the best fit spatial frequency

sl=spatial_frequency*[cos(direction) sin(direction)]; %transform to linear slopes

% calculate offset for fitted values
offs = atan2(sum(sin(phase-sl(1)*pos_x-sl(2)*pos_y)),sum(cos(phase-sl(1)*pos_x-sl(2)*pos_y))) ;

% circular-linear correlation:
pos_circ = mod(sl(1)*pos_x+sl(2)*pos_y+offs, 2*pi); % circular variable derived from the position
phase_mean = mod(angle(sum(exp(1i*phase))/length(phase)),2*pi); % circular mean of the theta phase
pos_circ_mean = mod(angle(sum(exp(1i*pos_circ))/length(phase)),2*pi); % circular mean of the circular position variable
cc = sum(sin(phase - phase_mean) .* sin(pos_circ - pos_circ_mean)) / sqrt( sum(sin(phase - phase_mean).^2) * sum(sin(pos_circ - pos_circ_mean).^2) );
% calculating the circ-correlation coefficient between measured phase and fitted
Rsquared=cc^2;
pos_circ=mod(pos_circ,2*pi);
if any(fig); set(gcf,'color','w','position',[100 100 1072 450]); subplot(1,1,1) %[770 50 500 650] for vertical plots
    if fig(1); subplot(1,2,1) % subplot 1 visualising resultant vector length on the parameter space.
    surf(slope1_Matrix,slope2_Matrix,Residual_resultant_vector_length) %make round parameter space.
    hold on
    scatter3(sl(1),sl(2),max(max(Residual_resultant_vector_length)),80,'y','filled')  % note the best fitting on parameter space
    plot3([sl(1) 0],[sl(2) 0],[max(max(Residual_resultant_vector_length)) max(max(Residual_resultant_vector_length))],'y:','LineWidth',2)
    shading('flat'); cm=jet(64); cm=cool*.95; colormap(gca,cm(50:end,:)); %colormap(gca,flipud(cm(11:end-15,[3 1 2]))); %cm=getjet; colormap(gca,flipud(cm(:,[3 1 2]))); 
    cb=colorbar;view(0,90); axis('off'); axis('equal'); caxis([0 1])
    xs=xlim; ys=ylim; zs=zlim; plot3([mean(ys) mean(ys);ys]',[xs;mean(xs) mean(xs)]',[zs(2) zs(2);zs(2) zs(2)],'k-'); axis([xs ys zs])
    axisshrinkfactor=1; axis(axis/axisshrinkfactor); 
    plot3(ang2rad(10+[10 30]),[ys(1) ys(1)]/axisshrinkfactor,[zs(2) zs(2)],'k-','LineWidth',3); text(ang2rad(20),min(ylim)-(diff(ylim)/20),zs(2),'20^o/mm')
    %for vertical plot: text([mean(xs)-.0125 xs(1)*.8],[ys(1)*1.15 mean(ys)],[zs(2) zs(2)],{'<-X->',''},'fontweight','bold','fontsize',12)
    %text([mean(xs)-.015 xs(1)*1.15],[ys(1)*.95 mean(ys)],[zs(2) zs(2)],{'<-X->',''},'fontweight','bold','fontsize',16)
    set(get(cb,'Title') ,'String','R','fontsize',14); %set(cb,'Position',[0.9    0.621    0.01    0.25])
    title({'Residual resultant vector length on parameter space',['R^2 = ' num2str(Rsquared) ', direction = ' num2str(rad2ang(direction)) ' degrees (' num2str(direction) ' radians)'],['spatial frequency (deg/mm): ' num2str(rad2ang(spatial_frequency)) ', wavelength: ' num2str(round(2*pi/spatial_frequency,1)) 'mm']}) %
    end
    if fig(2)
    subplot(1,2,2) %  subplot 2 visualising fitted plane and the residuals.
    [x_grid,y_grid]=meshgrid(linspace(min(pos_x),max(pos_x),50),linspace(min(pos_y),max(pos_y),50)); % make plane matrix
    mesh(x_grid,y_grid,mod(x_grid*sl(1)+y_grid*sl(2)+offs,2*pi),'facealpha',.5); % plot the fitted plane
    shading('flat');colormap(gca,hsv); cb=colorbar; set(cb,'ticks',[0 pi 2*pi],'ticklabels',{'0','pi','2*pi'},'fontsize',10)
    set(get(cb,'Title') ,'String','Phase','fontsize',15);
    hold on
    for i=1:length(phase)
        scatter3(pos_x(i),pos_y(i),phase(i),60,'r','filled') % plot the real phase values
        plot3([pos_x(i),pos_x(i)],[pos_y(i),pos_y(i)],[phase(i),pos_circ(i)],'k','linewidth',2) % lineup the real phase values and the plane
    end
    xlabel('X coordinates (mm)')
    ylabel('Y coordinates (mm)')
    zlabel('Circular values (Radians)'); caxis([0 2*pi])
    title('Circular-linear regression data fitting');
    [Y,X,Z]=meshgrid(0:2:14,0:2:6,0); Y=flipud(Y);
    plot3(X,Y,Z,'.','color',[.3 .3 .3],'markersize',20) %plot all
    surf(X,Y,Z,'facealpha',0,'edgealpha',.1); set(gca,'ztick',[0 pi 2*pi],'zticklabels',{'0','pi','2*pi'})
    axis equal; zlim([0 2*pi]); view(150,25); %view(255,35)
    end
end



