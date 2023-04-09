function plot_ripple(ripple_num,fig,S,sld,lb)

	if ripple_num < 1; ripple_num = 1; end
	if ripple_num > size(S.R_filt,1); ripple_num = size(S.R_filt,1); end

	if exist('sld','var'); sld.Value = ripple_num; end
	if exist('lb','var'); lb.Text = ['Ripple #' num2str(ripple_num)]; end

	xlimit = [-1 1]*.075;

	s1 = subplot(3,1,1,'Parent',fig);

	% OPTION 1: plot all 32 channels
	eegplotbytime2021(squeeze(S.R_filt(ripple_num,:,:)),S.sfx,500,-.25,s1);
	xlim(s1,xlimit); ylim(s1,[-33 0]);
	
	% OPTION 2:average cross superior/inferior directions
% 	R_filt_avg = nan(size(S.R_filt,2),8);
% 	for i = 1:8
% 		R_filt_avg(:,i) = mean(S.R_filt(ripple_num,:,[i i+8 i+16 i+24]),3);
% 	end
% 	eegplotbytime2021(R_filt_avg,S.sfx,500,-.25,s1);
% 	xlim(s1,xlimit); ylim(s1,[-9 0]);

	yticks(s1,-32:8:-8); yticklabels(s1,{'32' '24' '16' '8'});

	s2 = subplot(3,1,2,'Parent',fig);
	imagesc(s2,S.centered_t,1:32,mod(squeeze(S.R_phz(ripple_num,:,:))',2*pi));
	colormap(s2, S.phzmap); caxis(s2,[0 2*pi]);
	xlim(s2,xlimit);
	yticks(s2,8:8:32); yticklabels(s2,{'8' '16' '24' '32'});

	s3 = subplot(3,1,3,'Parent',fig);
	
	t_sig = S.tw_dirs_Rsqr(ripple_num,:) > 0.50;
	scatter(s3, S.centered_t(t_sig),circ_rad2ang(S.tw_dirs(ripple_num,t_sig)), ...
		S.tw_dirs_Rsqr(ripple_num,t_sig)*250, ...
		'filled', ...
		'AlphaData',S.tw_dirs_Rsqr(ripple_num,t_sig), ...
		'MarkerFaceAlpha','flat');
	hold(s3,'on');
	scatter(s3, S.centered_t(~t_sig),circ_rad2ang(S.tw_dirs(ripple_num,~t_sig)), ...
		S.tw_dirs_Rsqr(ripple_num,~t_sig)*250, ...
		"red", ...
		'filled', ...
		'AlphaData',S.tw_dirs_Rsqr(ripple_num,~t_sig), ...
		'MarkerFaceAlpha','flat');
	hold(s3,'off');

	xlim(s3,xlimit); ylim(s3,[0 360]);
	yticks(s3,0:90:360); yticklabels(s3,{'0\circ' '90\circ' '180\circ' '270\circ' '360\circ'});
	xlabel(s3,'time (sec)'); ylabel(s3,'Directin of TW');
end