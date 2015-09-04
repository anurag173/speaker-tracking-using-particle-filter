function [res,data] = IMRevTimeAnalysis(SetupStr,varargin)
%IMRevTimeAnalysis  reverberation time analysis of image method setup
%
% [res] = IMRevTimeAnalysis(SetupStr,'ArgName',ArgVal,...)
%
% This function performs an analysis of the reverberation time resulting
% from image method computations for a specific user-defined configuration.
% The input string 'SetupStr' corresponds to the name (without '.m'
% extension) of the m-function containing the desired setup definitions,
% including parameters such as sampling frequency, room dimensions and
% reflection coefficients (see 'IMSetup.m' for an example). This function
% then simulates a source of white noise stopping abruptly in the defined
% environment. Measuring the slope of the decay of sound pressure level
% in the room after the source is stopped then provides a direct measurement
% of the T60 value. This function provides another two T60 estimates by 
% Schroeder's integration method on the audio data, and by determining the
% energy decay rate from the computed transfer function itself (provides the
% most accurate T60 estimates). This process is then repeated for a series
% of source-receiver configuration to obtain a statistically representative
% set of measurements.
%
% This function also accepts a series of 'ArgName' - 'ArgVal' parameters
% corresponding to the user-definable parameters below:
%
%         PlotRes: set to 1 if plots of intermediate results desired 
%                  (execution will be paused). Disabled by default.
%        LimAtten: maximum attenuation in dB in IM transfer function
%                  computations, see function 'MakeIMResp.m'. Default 45.
%   WriteRes2File: set to 0 if analysis results are not to be appended to
%                  the end of setup m-file 'SetupStr'. Enabled by default.
%       NumConfig: number of source-receiver configurations. Default 50.
%
% Returns the following parameters:
%
%    res.medMeasT60: median T60 value (over source-receiver configurations) 
%                    obtained from slope measurements in the SPL decay curve
% res.medSchMeasT60: median T60 value (over source-receiver configurations) 
%                    obtained from slope measurements in the audio data using 
%                    Schroeder's integration method
%      res.medIMT60: median T60 value (over source-receiver configurations) 
%                    obtained directly from the computed transfer functions 
%                    using Schroeder's integration method

% Release date: September 2006
% Author: Eric A. Lehmann, WATRI, Perth, Australia (www.watri.org.au)
%         Eric.Lehmann@watri.org.au -- Tel. +61 (0)8 6488 4642
%
% Copyright 2006 Eric A. Lehmann
%
% This program is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation; either version 2 of the License, or (at your
% option) any later version.
% 
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
% Public License for more details.
% 
% You should have received a copy of the GNU General Public License along
% with this program; if not, write to the Free Software Foundation, Inc.,
% 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

% User input variables:
VarList = {'PlotRes'        0;        % set to 1 if plots of intermediate results desired
           'LimAtten'       45;       % maximum attenuation in IM transfer function computation
           'NumConfig'      50;       % total number of src-rcv configurations considered
           'WriteRes2File'  1};       % set to 1 if analysis results are to be appended to setup file
eval(SetUserVars(VarList,varargin));   % set user-definable variables

nNoiseT60period = 7;    % number of periods (T60 length) of simulated noise, must be bigger than 2
nSilT60period = 7;      % number of periods (T60 length) of simulated silence
LineTolNumSTD = 8;      % number of STD used as tolerance to check for SPL decrease slope
NumSlopeVals = 14;      % approximate number of points during the SPL decrease period

eval(['SetupStruc = ' SetupStr ';']);    % fetch parameters for desired setup
Fs = SetupStruc.Fs;
cc = SetupStruc.c;
room = SetupStruc.room;
beta = SetupStruc.beta;
if isequal(beta(:),zeros(6,1)),
    error('Reverberation time analysis cannot be carried out with all ''beta'' parameters equal to zero!'); 
end

% Determine rough estimate of reverberation time:
V = prod(room);
S = room(1)*room(2)*2 + room(1)*room(3)*2 + room(2)*room(3)*2;
aa_sab = ((2-beta(1)^2-beta(2)^2)*room(2)*room(3) + (2-beta(3)^2-beta(4)^2)*room(1)*room(3) + (2-beta(5)^2-beta(6)^2)*room(1)*room(2))/S;
T60sab = 0.161 * V/(S*aa_sab);        % Sabine's reverberation time (in s).
T60sab_samp = round(T60sab*Fs);       % Sabine's reverberation time (in samples).
aa_eyr = -log(1-aa_sab);
T60eyr = 0.161 * V/(S*aa_eyr);      % Eyring's reverberation time (in samples).

% Length of chunks considered when determining envelope of received SPL
maxlen = ceil(T60sab_samp*LimAtten/60/NumSlopeVals);    % allows approximately NumSlopeVals values in the SPL decrease slope

% source signal: noise followed by silence:
src_data = [randn(1,nNoiseT60period*T60sab_samp) randn(1,nSilT60period*T60sab_samp)*0.000001]; 

% simulate a few audio samples with different source-mic configurations
MeasT60vec = zeros(1,NumConfig);
IMT60vec = zeros(1,NumConfig);
SchMeasT60vec = zeros(1,NumConfig);

PrintLoopPCw('   [IMRevTimeAnalysis] Computing sample TFs. ');
for ii=1:NumConfig,
    PrintLoopPCw(ii,NumConfig);
    
    X_src = rand(1,3).*room*.8 + .1*room;   % avoid positions close to walls
    X_rcv = rand(1,3).*room*.8 + .1*room;
    while norm(X_src-X_rcv)<1,    % choose new points if they're too close to each other...
        X_src = rand(1,3).*room*.8 + .1*room;
        X_rcv = rand(1,3).*room*.8 + .1*room;
    end
    
    TFcoeffs = MakeIMResp(Fs,beta,X_src,X_rcv,room,cc,LimAtten,[],1);  % Compute TF for given src/rcv positions
    TFlen = length(TFcoeffs);
    DPdel = norm(X_rcv-X_src)/cc;   % Direct path delay in s
    DPdel_samp = ceil(DPdel*Fs);    % Direct path delay in samples
    
    %-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    % T60 measurement from resulting TF using Schroeder integration method
    %---- this method checks for a slope fit starting from when the overall TF 
    %---- energy becomes below 1/3 of the total energy, and ending at the index 
    %---- which generates the lowest median squared fitting error.
    EnDecayVec = zeros(1,TFlen);
    for kk=1:TFlen,
        EnDecayVec(kk) = sum(TFcoeffs(kk:end).^2);  % Energy decay using Schroeder's integration method
    end
    foo = find(EnDecayVec>=EnDecayVec(1)/3);    % Discard direct path and early reflections
    intstart = foo(end);                        % Index from which to start fitting decay line
    EnDecayVec(EnDecayVec==0) = eps;
    EnDecayVec = 10*log10(EnDecayVec);          % Decay curve in dB.

    slopevec = NaN*ones(1,TFlen);
    errorvec = NaN*ones(1,TFlen);
    for kk=intstart+1:TFlen,
        polyp = polyfit([intstart:kk],EnDecayVec(intstart:kk),1);   % Fit decay line on considered part of the energy decay curve
        slopevec(kk) = 60/abs(polyp(1))/Fs;                         % Slope of the decay line, i.e. T60 value
        errorvec(kk) = median((EnDecayVec - polyval(polyp,[1:TFlen])).^2);    % Fitting error for current slope based on entire decay curve
    end                                                             % using the error median avoids influence of end part of the curve
    foo = find(errorvec==nanmin(errorvec));     % T60 estimate as decay line minimising squared error
    IMT60vec(ii) = slopevec(foo(1));     % T60 estimate from TF computations
    %-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    
    rcv_data = conv(TFcoeffs,src_data);     % signal at receiver
    rcv_data = rcv_data(1:end-TFlen+1);     % crop to useful length

    %-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    % T60 measurement from sound field simulation using Schroeder integration method
    DecayStartInd = floor((nNoiseT60period*T60sab + DPdel)*Fs);
    DecayStopInd = DecayStartInd + TFlen - DPdel_samp;%ceil(T60sab*Fs*1.3);          % integrate for one T60 max
    DecayIndLen = DecayStopInd-DecayStartInd+1;
    EnDecayVec = zeros(1,DecayIndLen);
    for kk=DecayStartInd:DecayStopInd,
        EnDecayVec(kk-DecayStartInd+1) = sum(rcv_data(kk:DecayStopInd).^2);
    end
    foo = find(EnDecayVec>=EnDecayVec(1)/2);    % Discard direct path and early reflections
    FitStartInd = foo(end)+DecayStartInd-1;                     % Index from which to start fitting decay line
    EnDecayVec(EnDecayVec==0) = eps;
    EnDecayVec = 10*log10(EnDecayVec);          % Decay curve in dB.

    slopevec = zeros(1,DecayIndLen);
    if PlotRes, polyvec = zeros(2,DecayIndLen); end;
    errorvec = NaN*ones(1,DecayIndLen);
    warning('off','MATLAB:polyfit:RepeatedPointsOrRescale');
    for kk=FitStartInd+1:DecayStopInd,
        polyp = polyfit([FitStartInd:kk],EnDecayVec(FitStartInd-DecayStartInd+1:kk-DecayStartInd+1),1);   % Fit decay line on considered part of the energy decay curve
        slopevec(kk-DecayStartInd+1) = 60/abs(polyp(1))/Fs;                         % Slope of the decay line, i.e. T60 value
        if PlotRes, polyvec(:,kk-DecayStartInd+1) = [polyp(1); polyp(2)]; end;
        errorvec(kk-DecayStartInd+1) = median((EnDecayVec - polyval(polyp,[DecayStartInd:DecayStopInd])).^2);    % Fitting error for current slope based on entire decay curve
    end                                                             % using the error median avoids influence of end part of the curve
    warning('on','MATLAB:polyfit:RepeatedPointsOrRescale');
    fooind = find(errorvec==nanmin(errorvec));     % T60 estimate as decay line minimising squared error
    FitStopInd = fooind(1)+DecayStartInd-1;
    if PlotRes, polyvec = polyvec(:,fooind(1)); end;
    SchMeasT60vec(ii) = slopevec(fooind(1));
    %-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    
    % determine envelope of received SPL:
    nrch = floor(length(rcv_data)/maxlen);
    maxindvec = zeros(1,nrch); 
    maxvalvec = zeros(1,nrch);
    for i=0:nrch-1,
        foo = (rcv_data(i*maxlen+1:(i+1)*maxlen)).^2;       % determine energy within chunk
        foo = 10*log10( sum(foo)/maxlen );
        maxindvec(i+1) = round((i+.5)*maxlen);
        maxvalvec(i+1) = foo;
    end
    
    rcv_data(rcv_data==0) = eps;            % avoid log of 0
    rcv_spl = 10*log10((rcv_data).^2);      % transfer to sound pressure level

    % determine average SPL during noise and silence periods:
    startind = find(maxindvec>=(2*T60sab_samp));
    startind = startind(1);
    endind = find(maxindvec>(nNoiseT60period*T60sab_samp));
    fooendind = endind(1) - 1;      % index where source stops emitting noise
    NoiseAvSPL = median(maxvalvec(startind:fooendind)); % median discard potential settling values at start
    LineTol = LineTolNumSTD*std(abs(maxvalvec(startind:fooendind)-NoiseAvSPL)); % tolerance for SPL decrease slope points below
    startind = find(maxindvec>=((nNoiseT60period+2)*T60sab_samp + DPdel_samp));
    startind = startind(1);
    SilAvSPL = mean(maxvalvec(startind:end));

    % find slope best fitting the data (i.e. containing the most data points):
    foo = nNoiseT60period*T60sab*Fs + TFlen;   % where the useful data ends!... (TFlen samples after noise stops)
    endind = find(maxindvec<=foo);
    endind = endind(end);
    
    SlopeStopTime = maxindvec(endind)/Fs;
    SlopeStartTime = nNoiseT60period*T60sab + DPdel;    % slope start: when last samples of noise reach receiver (direct path)
    startind = find(maxindvec>(SlopeStartTime*Fs));
    startind = startind(1);
    SlopeValVec = maxvalvec(startind:endind);
    SlopeIndVec = maxindvec(startind:endind);

    SlopeStartPt = [SlopeStartTime NoiseAvSPL];     % condition: slope must go through SlopeStartPt
    SlopeValVec_len = length(SlopeValVec);
    slopevec = zeros(1,SlopeValVec_len);
    PtsWithinTolVec = zeros(1,length(SlopeValVec));
    for ll=1:SlopeValVec_len,       % Determine candidate slope, then see how many line points are included. Best slope is that containing most points.
        slope = [SlopeIndVec(ll)/Fs SlopeValVec(ll)] - SlopeStartPt;
        slopevec(ll) = slope(2)/slope(1);   % current slope, deltaY/deltaX
        avslope = median(slopevec);         % compute average slope of all the points so far (median discards extremes!!)
        bb = SlopeStartPt(2) - avslope*SlopeStartPt(1); % line through SlopeStartPt with desired slope

        nPtsWithinTol = 0;          % for current average slope, check how many line points are included (i.e. within tolerance)
        for mm=1:length(SlopeValVec),
            foo = avslope * SlopeIndVec(mm)/Fs + bb;    % line point at current time
            if (SlopeValVec(mm)<=(foo+LineTol)) && (SlopeValVec(mm)>=(foo-LineTol)),
                nPtsWithinTol = nPtsWithinTol + 1;      % current point is within tolerance from line point
            end
        end

        PtsWithinTolVec(ll) = nPtsWithinTol;    % number of points within tolerance for this slope
    end

    % refine slope estimate for line not necessarily passing through SlopeStartPt
    maxind = find(PtsWithinTolVec==max(PtsWithinTolVec));
    maxind = max([maxind(1) max(PtsWithinTolVec)]);  % pick index yielding the maximum number of points
    if maxind>1,
        RealSlopeStopTime = SlopeIndVec(maxind)/Fs;  % for plots below
        polyxvec = SlopeIndVec(1:maxind)/Fs;
        polyyvec = SlopeValVec(1:maxind);
        polyp = polyfit(polyxvec,polyyvec,1);       % 1st order polynomial approximation of line points
        ResReverbTime = 60/abs(polyp(1));           % measured T60 derived from slope of line approximation
        MeasT60vec(ii) = ResReverbTime;
    else
        RealSlopeStopTime = NaN;    % no line points available for slope estimate
        ResReverbTime = NaN;
    end

        plot([1:length(rcv_data)]/Fs,rcv_data);
        axis tight; xlabel('time (s)'); title('Received data');
    % results plot if desired
    if PlotRes,
        reusefig('T60 Analysis of Image Method computations'); clf;
        subplot(2,2,1);
        ymin = SilAvSPL-20; ymax = NoiseAvSPL+20;
        plot([1:length(rcv_spl)]/Fs,rcv_spl,'color',[.75 .75 .75]); hold on;       % SPL data
        plot(maxindvec/Fs,maxvalvec,'color','r','marker','o','linestyle','-','markersize',3);   % envelope
        plot(nNoiseT60period*T60sab*ones(1,2),[ymin ymax],'b--');       % end of noise line
        plot(SlopeStartTime*ones(1,2),[ymin ymax],'b--');               % start of slope line
        plot(SlopeStopTime*ones(1,2),[ymin ymax],'b--');                % max end of slope line
        plot(RealSlopeStopTime*ones(1,2),[ymin ymax],'b:');             % considered end of SPL decrease
        plot([1 length(rcv_spl)]/Fs,NoiseAvSPL*ones(1,2),'g-');         % average noise SPL
        plot([1 length(rcv_spl)]/Fs,SilAvSPL*ones(1,2),'g-');           % average silence SPL
        if maxind>1,
            plot(([ymin ymax]-polyp(2))/polyp(1),[ymin ymax],'k');      % estimated line
            plot(([ymin ymax]-polyp(2))/polyp(1),[ymin ymax]+LineTol,'k:');
            plot(([ymin ymax]-polyp(2))/polyp(1),[ymin ymax]-LineTol,'k:');
        end
        text(nNoiseT60period*T60sab,ymin,'  Noise end','verticalalignment','bottom','color','k','rotation',90);
        text(SlopeStartTime,ymin,'  Decay start','verticalalignment','top','color','k','rotation',90);
        text(RealSlopeStopTime,ymin,'  Decay end','verticalalignment','bottom','color','k','rotation',90);
        text(SlopeStopTime,ymin,'  Max decay end','verticalalignment','top','color','k','rotation',90);
        text((nNoiseT60period+2)*T60sab,NoiseAvSPL,'Noise SPL  ','verticalalignment','bottom','horizontalalignment','right','color','k');
        text((nNoiseT60period-1)*T60sab,SilAvSPL,'  Silence SPL','verticalalignment','bottom','color','k');
        title(sprintf('Measured T60: %.4f (s)',ResReverbTime));
        axis tight;
        xlim([(nNoiseT60period-1)*T60sab (nNoiseT60period+2)*T60sab]);
        ylim([ymin ymax]);
        xlabel('time (s)'); ylabel('SPL (dB)');

        subplot(2,2,2);
        indvec = [DecayStartInd:DecayStopInd];
        plot(indvec/Fs,EnDecayVec,'r'); hold on;
        plot(indvec/Fs,polyval(polyvec,indvec),'k');
        axis tight; yl = ylim; yl = [yl(1)/2+yl(2) yl(2)]; ylim(yl);
        plot([FitStartInd FitStartInd]/Fs,[yl(1) yl(2)],'b--');
        plot([FitStopInd FitStopInd]/Fs,[yl(1) yl(2)],'b--');
        text(FitStartInd/Fs,yl(1),'  Slope fitting start','verticalalignment','bottom','color','k','rotation',90);
        text(FitStopInd/Fs,yl(1),'  Slope fitting stop','verticalalignment','top','color','k','rotation',90);
        xlabel('time (s)'); ylabel('Energy (dB)');
        title(sprintf('T60 via Schroeder integration method: %.4f (s)',SchMeasT60vec(end)));

        subplot(2,2,3);
        plot([1:length(rcv_data)]/Fs,rcv_data);
        axis tight; xlabel('time (s)'); title('Received data');

        subplot(2,2,4);
        plot([1:length(TFcoeffs)]/Fs,TFcoeffs,'r');
        axis tight; xlabel('time (s)'); ylabel('TF coeffs');
        title(sprintf('TF''s T60 (Schroeder integration): %.4f (s)',IMT60vec(ii)));

        if NumConfig-ii>=1, pause; end;
    end

end
 
medMeasT60 = median(MeasT60vec);
medIMT60 = median(IMT60vec);
iqrMeasT60 = iqr(MeasT60vec);
iqrIMT60 = iqr(IMT60vec);
medSchMeasT60 = median(SchMeasT60vec);
iqrSchMeasT60 = iqr(SchMeasT60vec);

res.medMeasT60 = medMeasT60;
res.medSchMeasT60 = medSchMeasT60;
res.medIMT60 = medIMT60;
res.T60sab = T60sab;
res.T60eyr = T60eyr;
data.MeasT60vec = MeasT60vec;
data.SchMeasT60vec = SchMeasT60vec;
data.IMT60vec = IMT60vec;

if WriteRes2File,
    fid = fopen(which(SetupStr),'a');
elseif nargout==0
    fid = 1;
else
    return
end

fprintf(fid,'\nOutStruc.measT60 = %.5f;\n\n',medIMT60);
fprintf(fid,'\n%%-=:=- -=:=- -=:=- -=:=- -=:=- -=:=- -=:=- -=:=- -=:=- -=:=- -=:=- -=:=- -=:=-');
fprintf(fid,'\n%%-=:=- Results from reverberation time analysis (''IMRevTimeAnalysis.m'')  -=:=-');
fprintf(fid,'\n%%\n%% Simulation parameters:  Fs =   %.5g',Fs);
fprintf(fid,'\n%%                         c =    %.5g',cc);
fprintf(fid,'\n%%                         room = [%.5g  %.5g  %.5g]',room(1),room(2),room(3));
fprintf(fid,'\n%%                         beta = [%.5g  %.5g  %.5g  %.5g  %.5g  %.5g]',beta(1),beta(2),beta(3),beta(4),beta(5),beta(6));
fprintf(fid,'\n%%\n%% Median T60 measured in room, over %d configurations:',length(MeasT60vec));
fprintf(fid,'\n%%  > T60 = %.5f s (inter-quartile range: %.5g s)',medMeasT60,iqrMeasT60);
fprintf(fid,'\n%%\n%% Median T60 measured in room (Schroeder integration method), over %d configurations:',length(SchMeasT60vec));
fprintf(fid,'\n%%  > T60 = %.5f s (inter-quartile range: %.5g s)',medSchMeasT60,iqrSchMeasT60);
fprintf(fid,'\n%%\n%% Median T60 estimated from transfer function, over %d configurations: ',length(IMT60vec));
fprintf(fid,'\n%%  > T60 = %.5f s (inter-quartile range: %.5g s)',medIMT60,iqrIMT60);
fprintf(fid,'\n%%\n%% Sabine''s T60 estimate: T60 = %.5g s',T60sab);
fprintf(fid,'\n%% Eyring''s T60 estimate: T60 = %.5g s',T60eyr);
fprintf(fid,'\n%%-=:=- -=:=- -=:=- -=:=- -=:=- -=:=- -=:=- -=:=- -=:=- -=:=- -=:=- -=:=- -=:=-\n\n');

if WriteRes2File,
    fclose(fid);
    fprintf('   [IMRevTimeAnalysis] Analysis results written to file ''%s''\n',which(SetupStr));
end

