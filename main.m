clearvars; close all;
id = 'LA04_B4'; BLtime = [325,468];
% id = 'LA03_B2'; BLtime = [60,250];
Q = load(['./H-uG Data/' id '/MAT files/' id '_notched.mat']);


%% Define variables

sfx = Q.sfx;
d = Q.d';
BLtimeidx = BLtime*sfx;
d = d(BLtimeidx(1):BLtimeidx(2)-1,:);

[Nsamp,Nch] = size(d);
hwripbandfilt = ceil(1/50*(.5)*sfx); % half window in which to find min voltage in ripple-band-filtered trace
wd = ceil(sfx/2); % window where data around ripple is extracted - in reality it is wd+1
wd_center = wd/2+1; % idx for center of window
centered_t = ((-wd/2:wd/2)/sfx)';

[Y,X,~] = meshgrid(0:2:14,0:2:6,0);
x = X'; y = Y'; xy = [x(:) y(:)]; 

phzmap = colormap(flipud(cool)); phzmap(1,:) = [1 1 1];
close;

%% Filter and pre-process data
bc = isnan(d(1,:));
if any(bc); error(['channels ' find(bc) ' are bad']); end
[d_filt,d_hb,d_phz] = deal(nan(size(d)));
for i = find(~bc)
	d_filt(:,i) = jkfilt(d(:,i),sfx,50,150,'cheby2',8);
	d_hb(:,i) = abs(hilbert(d_filt(:,i)));
	d_hb(:,i) = smooth(d_hb(:,i),round(sfx*.025));
	d_phz(:,i) = angle(hilbert(d_filt(:,i)));
end

d_hb_z = nan(size(d_hb));
for i = 1:size(d_hb,2)
	[nm, ns, ~, ~] = robustMean(d_hb(:,i),[],4,[]);
	d_hb_z(:,i) = (d_hb(:,i)-nm)./ns; %z-score
	clear nm ns;
end

%% Ripples detection
% Method: Hilbert 50-150 Hz analytic amplitude envelope >1.8 robust-
% S.D. simultaneously on 9 or more electrodes

% H = d_hb_z'>1.8;
% HT = sum(H,1)>=9;
% eON = find(diff(HT)==1);
% eOFF = find(diff(HT)==-1);
% if isempty(eON)&&~isempty(eOFF)
% 	eON = 1;
% elseif isempty(eOFF)&&~isempty(eOFF)
% 	eOFF = length(HT);
% end
% if eOFF(1)<eON(1); eOFF(1)=[]; end
% if length(eOFF)<length(eON); eON(end)=[]; end
% if length(eOFF)~=length(eON); error('start and end of ripples is not matching up'); end
% ripples = [eON' eOFF'];
% 
% evmin = 0.02; evmax = 0.2;
% ripples = sortrows(ripples); % sorting chronologically is required for ANgetEventdata_... function
% evdur = diff(ripples,1,2)/sfx;
% evremove = evdur<=evmin|evdur>evmax|ripples(:,1)<wd|ripples(:,2)>Nsamp-wd;
% ripples = ripples(~evremove,:);

%% ALTERNATIVE: get ripples from file (post manual examination)

R = load(['./Ripples Timestamps/' id '_ripples.mat']);
ripples = round(R.Events(:,[1,2])*sfx/R.dsfx-BLtimeidx(1)+2);

%% Extract data around ripples

middlech = 12;
Nripples = size(ripples,1);

[R_filt,R_hb_z,R_phz] = deal(nan(Nripples,wd+1,Nch));
Rtimeofpeak = nan(Nripples, 1);

for i = 1:Nripples

	[~,r] = max(mean(d_hb_z(ripples(i,1):ripples(i,2),:),2)); r = r-1; % find peak of this ripple in time from onset to offset
	r = ripples(i,1)+r; % find moment of ripple peak in non-downsampled transform
	% r = ripples(ev,1)+r-1; % Jon's code

	[~,r2] = min(d(r-hwripbandfilt:r+hwripbandfilt,middlech)); r2 = r2-1; %find trough of ripple-band-filtered data around peak
	r = (r-hwripbandfilt)+r2;

	riprange = r-wd/2:r+wd/2;
	R_filt(i,:,:) = d_filt(riprange,:);
	R_hb_z(i,:,:) = d_hb_z(riprange,:);
	R_phz(i,:,:) = d_phz(riprange,:);
	Rtimeofpeak(i,1) = r;

	clear r r2 riprange;
end

%% Compute circular average across ripples
[R_phz_mean,R_phz_rtest] = deal(nan((wd+1),Nch));
for i = 1:(wd+1)
    for j = 1:Nch
        R_phz_mean(i,j) = circ_mean(R_phz(:,i,j));
        R_phz_rtest(i,j) = circ_rtest(R_phz(:,i,j));
    end
end; clear i j
R_phz_mean = mod(R_phz_mean,2*pi);
R_phz_mean(R_phz_rtest>.001) = nan;

%% Compute traveling direction for each Ripple AND each time point 

% centered by peak or trough
tw_filename = ['./tw_' id '_trough.mat'];

if isfile(tw_filename)
	load(tw_filename);
else
	fprintf('Computing traveling waves...');

	[tw_dirs,tw_dirs_Rsqr,tw_sf] = deal(nan(Nripples,wd+1));
	
	wd_true = wd+1;
	parfor i = 1:Nripples
		for j = 1:wd_true
			[dir,sf,~,Rsqr] = circ_lin_regress_2D(squeeze(R_phz(i,j,:)),xy,true(1,32),0);
			tw_dirs(i,j) = dir;
			tw_dirs_Rsqr(i,j) = Rsqr;
			tw_sf(i,j) = sf;
		end
		fprintf(['completed ripple ' num2str(i) '\n']);
	end

	save(tw_filename, 'tw_dirs', 'tw_dirs_Rsqr', 'tw_sf');
end

%% FIG 1: averaged across multiple ripples

fig1 = figure('color','w','position',[145 482 800 800], 'ToolBar','none'); 

s1 = subplot(3,1,1);
eegplotbytime2021(squeeze(mean(R_filt,1)),sfx,1500,-.25,s1);
xlim([-1 1]*.075); 
ylim([-33 0]); yticks(-32:8:-8); yticklabels({'32' '24' '16' '8'});

subplot(3,1,2);
imagesc(centered_t,1:32,R_phz_mean');
colormap(phzmap); caxis([0 2*pi]);
xlim([-1 1]*.075);
yticks(8:8:32); yticklabels({'8' '16' '24' '32'});

subplot(3,2,5);
v = R_phz_mean(wd_center,:); 
phzs = flip(flip(reshape(v,8,4)',1),2);
imagesc(phzs); hold on
contour(phzs,'color','k','LineWidth',1,'levelstep',ang2rad(10)); 
colormap(phzmap); caxis([0 2*pi]);
set(gca,'xtick',[],'ytick',[]);
c = colorbar; c.Ticks = [0 pi 2*pi]; c.TickLabels = {'0','\pi','2\pi'};


%% FIG 2: just looking at a single ripple

S = struct('sfx',sfx,'R_filt',R_filt,'R_phz',R_phz, ...
	'centered_t',centered_t,'phzmap',phzmap, ...
	'tw_dirs',tw_dirs,'tw_dirs_Rsqr',tw_dirs_Rsqr);

% fig2 = figure('color','w','position',[145 482 800 800], 'ToolBar','none'); 
% plot_ripple(10,fig2,S);

%% FIG 2 (INTERACTIVE): just looking at a single ripple

close(findall(0, 'type', 'figure'));

fig2_int = uifigure('color','w','Position',[145 482 800 800]);
fig2_p = uipanel(fig2_int,'Position',[15 100 770 650], ...
	'BackgroundColor',[0.97 0.97 0.97] ,'BorderType','none');
fig2_p.AutoResizeChildren = 'off';

warning('off','MATLAB:ui:Slider:fixedHeight');

lb = uilabel(fig2_int,'Position',[100 60 150 30]);
lb.FontSize = 12;
lb.FontWeight = 'bold';

sld = uislider(fig2_int,'Position',[15 50 770 50]);
sld.ValueChangedFcn = @(sld,event) plot_ripple(round(sld.Value),fig2_p,S,sld,lb);
sld.Limits = [1 Nripples]; sld.MajorTicks = 0:10:100; sld.MinorTicks = 0:1:Nripples;
sld.FontSize = 9;

plot_ripple(1,fig2_p,S,sld,lb);

btn_dec = uibutton(fig2_int,'Position',[15 60 30 30], 'Text','<');
btn_dec.FontSize = 16;
btn_dec.ButtonPushedFcn = @(btn_dec,event) plot_ripple(sld.Value-1,fig2_p,S,sld,lb);

btn_dec = uibutton(fig2_int,'Position',[55 60 30 30], 'Text','>');
btn_dec.FontSize = 16;
btn_dec.ButtonPushedFcn = @(btn_dec,event) plot_ripple(sld.Value+1,fig2_p,S,sld,lb);

%% FIG 3: traveling direction at center time BETWEEN RIPPLES 

fig3 = figure('color','w','position',[48 255 350 350],'ToolBar','none','MenuBar','none');
ax = polaraxes;

hist_range = wd_center;
tw_dirs_sig = tw_dirs;
% tw_dirs_sig(tw_dirs_Rsqr < 0.5) = nan;
all_dirs = tw_dirs_sig(:,hist_range);
polarhistogram(all_dirs(:),18,'facecolor',[0 0 0],'facealpha',0.2,'edgealpha',0.5);
ax.ThetaZeroLocation = 'bottom';



