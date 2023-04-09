function [ha,hp,unwrhp]=jkhilbert(d,sfx,lo,hi,filttype,fig) % Jon Kleen 2016
% Hilbert transform of vector or matrix (of EEG/ECoG data), after bandpass-filter
% INPUTS: d is data as a vector or 2-d matrix
% sfx is sampling frequency
% lo and hi are the respective low and high frequency cutoffs for the band
% type is a string that specifies which of 2 filters to use: 
%         butterfield ('butter') or cheby type II ('cheby2' or 'cheby')
% OUTPUTS: hilbert amplitude (ha), phase (hp), and unwrapped phase (tohp)
% NOTE: You might want to downsample your data if high sampling rate.
if any(isnan(d)); disp('FYI: your data contains NaNs'); end
if length(size(d))>2; disp('too many dimensions. vector or 2-d matrix only'); return; end
flip=0; mat=0; 
if size(d,1)>1 && size(d,2)>1; mat=1; end %if d is a matrix
if size(d,1)<size(d,2); d=d'; flip=1; end %flip data to orient if needed
if nargin<5; filttype='butter'; end
[f] = jkfilt(d,sfx,lo,hi,filttype); % uses jkfilt on vector or matrix
if mat; for i=1:size(d,2); f(:,i)=f(:,i)-mean(f(:,i)); end
else f=f-mean(f);
end;
hilb=hilbert(f); %hilbert transform
ha=abs(hilb); %analytic amplitude
if nargout>=2; hp=angle(hilb); end %phase
if nargout>=3; unwrhp=unwrap(hp); end %unwrapped phase

if exist('fig','var') && fig; 
figure; x=0:1/sfx:size(d,1)/sfx-1/sfx; 
if mat
subplot(2,1,1); eegplot2022(f,sfx,125); xlabel('time (sec)'); ylabel('channels'); title(strcat('Filtered data (',num2str(lo),'-',num2str(hi),'hz)')); 
subplot(2,1,2); pcolor(x,1:size(d,2),ha'); xlabel('time (sec)'); ylabel('channels'); shading flat; set(gca,'YDir','reverse'); title('Analytic amplitude'); colormap(jet)
% if strcmp(questdlg('Graph phase data?','Graph phase data?','Yes','No','Yes'),'Yes'); 
figure; pcolor(x,1:size(d,2),hp'); xlabel('time (sec)'); ylabel('channels'); shading flat; set(gca,'YDir','reverse'); title('Phase'); colormap('hsv'); colorbar
cmg=colormap(gray); phzgray=cmg(1:2:end,:); phzgray(:,3)=phzgray(:,3)*-1+1; colormap([phzgray; (phzgray*-1)+1]);
% end
else hold on; in=ones(length(f),1); in([1:round(0.02*length(f)) end-round(0.02*length(f)):end])=0; in=logical(in);
% plot(x,(d/max(d))*max(ha),'color',[.4 .4 .4]); %scaled down amplitude
plot(x,d,'color',[.4 .4 .4]); 
plot(x(in),f(in),'color',[0 .7 0],'linewidth',1.5); 
% plot(x(in),ha(in),'b','linewidth',1.5); 
legend('Raw data (scaled down)', strcat('Filtered data (',num2str(lo),'-',num2str(hi),'hz)'),'Analytic amplitude');
end
end

if flip; %flip back if flipped
    ha=ha';
    if nargout>=2; hp=hp'; end
    if nargout>=3; unwrhp=unwrhp';  end
end 


% 
% function [amp,phz,totphz] = jkhilbert(d,sfx,lo,hi)
% % Gets hilbert transform amplitude (amp) and phase (phz) from EEG trace d,
% % first bandpass-filtered between frequencies lo and hi.
% % totphz is unwrapped phase.
% % YOU MIGHT WANT TO DOWNSAMPLE YOUR DATA if high sampling rate.
% 
% filterord=4; ripple=20;
% [b,a]=cheby2(filterord, ripple, [lo hi]/(sfx/2));
% h=filtfilt(b,a,d); 
% h=h-mean(h);
% hilb=hilbert(h)'; 
% amp=abs(hilb);
% if nargout>=2
%     phz=angle(hilb);
% end
% if nargout>=3
%     totphz=unwrap(phz);
% end
