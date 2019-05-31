% Att göra: 
% 1) Att läsa in och skapa ett powerspektrum för vart och ett av de
%    fyra fallen, samt olika kopplingar, precis som i figur 4a
% 2) Att läsa in och skapa DTF för vart och ett av de fyra fallen,
%    pröva några olika frekvenser


%  BADA: cue efter 2000 ms, 500 ms lång, 7000 ms långt delay
%  CUE PFC: cue efter 2000 ms, 7000 ms lång, inget delay
%  CUE PPC: cue efter 2000 ms, 7000 ms lång, inget delay
%  INGEN CUE: total tid 7000 ms.
%
%  Föreslår att vi skippar de första 500 ms på delayperioden, 
%  då får vi 4000 ms tid för varje typ




clear
filer{1}=extract_sim_dirs('MYSIMS\DTF\SIM_SERIER\SERIES_testDTF.txt');
filer{2}=extract_sim_dirs('MYSIMS\DTF\SIM_SERIER\SERIES_testDTF_CUE.txt');
titlename{1} = 'Delay';
titlename{2} = 'Cue till PFC';


% plots BUMP EEG. 
% dirname: The name of the simulation directory
% smooth =  degree of smoothing of histogram. Default = 2
%
% Version 2.0
% Author: Fredrik Edin, 2004
% Address: freedin@nada.kth.se

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% LOAD FILES %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
smooth = 2;

thisdir = pwd;

close all
LFP_ppc = cell(length(filer),1);
LFP_pfc = cell(length(filer),1);

for i = 1:length(filer)
    for j = 1:length(filer{i})
        
        cd( filer{i}{j} )
        
        % load data files
        load('Params.txt')
        version = Params(1);
        load('Q.txt')
        nmod = Params(2);
        if version == 3
            tStart = Params(4);
        	tmp = [ Params(1:6) ; 100 ; 100 ];
        	for k = 1:nmod
                tmp = [ tmp ; Params(3+4*k) ; 1000 ; Params(5+4*k) ];
    	        tmp = [ tmp ; Params(4+4*k) ; 1000 ; Params(6+4*k) ];
        	end
        	Params = tmp;
        elseif version == 4
            tStart = Params(4);
        	Params = [ Params(1:6) ; 100 ; 100 ; Params(7:end) ];
        elseif version == 5
            tStart = Params(4);
        end
        tStop = Params(5);
        
        
        % Name of the simulation
        filename = strcat( filer{i}{j}, '/', 'LFP_prox_dist.txt' );
        if exist(filename(1:end-4),'file')
            movefile(filename(1:end-4), filename);
        end
        tmp=load( filename );
        
        % plot the LFPgrams of the cell
        filename = strcat( filer{i}{j}, '/', 'Q.txt' );
        load( filename )
        if length(Q) > 0 %& Q(1,5) ~= 0
            % Normalize data to 0 mean and variance 1
            tmp2 = tmp(Q(1,1)-tStart+1000:Q(1,1)-tStart+5000,2);
            tmp2 = (tmp2 - mean(tmp2)) / std(tmp2);
            tmp3 = tmp(Q(1,1)-tStart+1000:Q(1,1)-tStart+5000,3);
            tmp3 = (tmp3 - mean(tmp3)) / std(tmp3);
            % Store as LFPs
            LFP_ppc{i} = [LFP_ppc{i} tmp2]; % OBS, har ar ordningen ratt igen eftersom jag anvander nya simulatorn
            LFP_pfc{i} = [LFP_pfc{i} tmp3]; % OBS, har ar ordningen ratt igen eftersom jag anvander nya simulatorn
        end
        cd( thisdir )
    end
    
    LFP_ppc{i} = LFP_ppc{i} - ones(size(LFP_ppc{i},1),1)*mean(LFP_ppc{i});
    LFP_pfc{i} = LFP_pfc{i} - ones(size(LFP_ppc{i},1),1)*mean(LFP_pfc{i});
    
end
dt = tmp(2,1)-tmp(1,1);
fS = 1000/dt;
fN = fS/2;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% BANDPASS %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fl = 0;  % Lower cut-off in Hz
fu = 55; % Higher cut-off in Hz

Wp = [fl fu]/fN;
fwinL = 0.005; % The width of the difference between lower pass- and stopbands
fwinU = 0.05;  % The width of the difference between upper pass- and stopbands
Ws = [];
Rp = 1;
Rs = 60;

% Create a bandpass filter
if fu == 100
    Wp(2) = [];
    Ws = min(Wp-fwinL,0.5*Wp);
    [n,Wn] = ellipord(Ws, Wp, Rp, Rs);
    [b,a] = ellip(n,1,60,Wn,'high');
else
    Ws = min(Wp(end)+fwinU,Wp(end)+0.25*(1-Wp(end)));
    if fl > 0
        Ws = [max(Wp(1)-fwinL,0.5*Wp(1)) Ws];
    else
        Wp(1) = [];
    end
    [n,Wn] = ellipord(Wp, Ws, Rp, Rs);
    [b,a] = ellip(n,1,60,Wn);
end 



       
leg = [];
wa = [];
Pp = [];
Pf = [];

for i = 1:length(LFP_pfc)
    Pp{i} = []; % Power spectrum of the 
    Pf{i} = [];
    for j = 1:size(LFP_pfc{i},2)
        filt_Lp{i}(:,j) = filter(b,a,reshape(LFP_ppc{i}(:,j),[],1));
        filt_Lf{i}(:,j) = filter(b,a,reshape(LFP_pfc{i}(:,j),[],1));
        [Pp{i}(:,j), Ppc{i}(:,2*j-1:2*j), f] = psd(filt_Lp{i}(:,j),1000,fS);
        [Pf{i}(:,j), Pfc{i}(:,2*j-1:2*j), f] = psd(filt_Lf{i}(:,j),1000,fS);
        Pp_o{i}(:,j) = psd(reshape(LFP_ppc{i}(:,j),[],1),1000,fS);
        Pf_o{i}(:,j) = psd(reshape(LFP_pfc{i}(:,j),[],1),1000,fS);
    end
end
cd(thisdir)



dt = tmp(2,1)-tmp(1,1);
% t = 0:dt:dt*(size(filt_Lp{1},1)-1);
f = 0:1:500;
color = 'rgb';
for i = 1:length(LFP_pfc)
        meanPf_o(:,i) = [mean(Pf_o{i}')]';
        meanPp_o(:,i) = [mean(Pp_o{i}')]';
end
figure(1)
clf
subplot(1,2,1)
color = 'bgrkm';
for i = 1:length(Pf_o)
    plot(f,log10(Pf_o{i})*10,'color',color(i))
    hold on
end
hold on
xlim([0 60])
ylim([-20 20])
legend({'0','100'},2)
title('PFC')
subplot(1,2,2)
for i = 1:length(Pp_o)
    plot(f,log10(Pp_o{i})*10,'color',color(i))
    hold on
end
hold on
xlim([0 60])
ylim([-20 20])
legend({'0','100'},3)
title('PPC')


cd(thisdir)




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% DTF %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
thisdir = pwd;



pmin = 1;
pmax = 100;
orsel = 'fpe';


for i = 1:length(filer)
    
    A0 = [];
    A25 = [];
    
    for j = 1:1
        
        % Do ARfit for prox_dist
%         A0 = [A0 ; filt_Lp{i}(:,j) filt_Lf{i}(:,j)];
%         A25 = [A25 ; filt_Lp{i}(:,j+1) filt_Lf{i}(:,j+1)];
%         A50 = [A50 ; filt_Lp{i}(:,j+2) filt_Lf{i}(:,j+2)];
        A0 = [A0 ; LFP_ppc{i}(:,j) LFP_pfc{i}(:,j)];
        A25 = [A25 ; LFP_ppc{i}(:,j+1) LFP_pfc{i}(:,j+1)];
        
    end
    
    [w0{i},AR0{i},C0{i},SBC0,FPE0,th0]=arfit(A0,pmin,pmax,orsel,'zero');
    order0(i) = size(AR0{i},2)/2;
    [w25{i},AR25{i},C25{i},SBC25,FPE25,th25]=arfit(A25,pmin,pmax,orsel,'zero');
    order25(i) = size(AR25{i},2)/2;
%     [w50{i},AR50{i},C50{i},SBC50,FPE50,th50]=arfit(A50,pmin,pmax,orsel,'zero');
%     order50(i) = size(AR50{i},2)/2;
    
    % results of ARFIT
    disp(['Order ' titlename{i} '  0%:' int2str(order0(i))])
    disp(['Order ' titlename{i} ' 100%:' int2str(order25(i))])
%     disp(['Order ' titlename{i} ' 50%:' int2str(order50(i))])
    
end 

% Do directed transfer function
f = 0:1:50;
dt = 0.001;

for i = 1:length(filer)
    
    for j = 1:length(f)
        
        [DTF0{i}{j},H0]=DirTransFunc(AR0{i},f(j),dt);
        wp0{i}(j) = DTF0{i}{j}(2,1) - DTF0{i}{j}(1,2);
        [DTF25{i}{j},H25]=DirTransFunc(AR25{i},f(j),dt);
        wp25{i}(j) = DTF25{i}{j}(2,1) - DTF25{i}{j}(1,2);
%         [DTF50{i}{j},H50]=DirTransFunc(AR50{i},f(j),dt);
%         wp50{i}(j) = DTF50{i}{j}(2,1) - DTF50{i}{j}(1,2);
        
    end
end


style = '--';
tit = [];
for i = 1:length(filer)
    
    figure(2)
    plot(f,wp0{i},style(1:i),f,wp25{i},style(1:i))%,f,wp50{i})
    hold on
    tit = [ tit titlename{i} ' (' style(1:i) ')  ' ];
    title(tit)
    legend('0','100')%,'50')
    ylim([-1 1])
    line(xlim, [0 0], 'Color', 'k', 'Linestyle', '--')
    yl = ylim;
    xl = xlim;
    text(xl(1) + 0.9*diff(xl), yl(1) + 0.25*diff(yl), 'P->F < F->P','HorizontalAlignment', 'right')
    text(xl(1) + 0.9*diff(xl), yl(1) + 0.20*diff(yl), 'DTF > 0 means netflow goes in direction P->F','HorizontalAlignment', 'right')
    
    disp(titlename{i})
    disp(['0: ' num2str(sum(wp0{i})) ' 100: ' num2str(sum(wp25{i}))])
    
end

for i = 1:length(filer)
        
    disp([titlename{i} ' 40Hz'])
    disp(['0: ' num2str(wp0{i}(41)) ' 100: ' num2str(wp25{i}(41))])
    
end

for i = 1:length(filer)
        
    disp([titlename{i} ' Gamma'])
    disp(['0: ' num2str(sum(wp0{i}(31:end))) ' 100: ' num2str(sum(wp25{i}(31:end)))])
    
end


