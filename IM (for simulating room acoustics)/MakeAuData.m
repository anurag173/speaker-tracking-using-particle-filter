function [AuData] = MakeAuData(TFFileName,AuFileName,SrcSignal,varargin)
%MakeAuData  creates audio samples from TF bank
%
% [AuData] = MakeAuData(TFFileName,AuFileName,SrcSignal)
% [AuData] = MakeAuData( ... ,'SNRval',noiseSNR)
% [AuData] = MakeAuData( ... ,'Dir',DirString)
% [AuData] = MakeAuData( ... ,'SaveRes',SaveFlag)
% [AuData] = MakeAuData( ... ,'AddNoise',NoiseFlag)
% 
% This function generates samples of audio data based on a bank of
% pre-computed transfer functions. The input variable 'TFFileName'
% determines the name of the .mat file where the TF bank is stored
% (typically as a result of using the function 'MakeTrajIMTFs.m'). The
% suffix '.mat' is also automatically appended to the given file name,
% which can be a full access path to the desired file. If no access path is
% given, the function looks for the desired file in the current working
% directory.
%
% The specific simulation parameters, such as room dimensions, microphone
% positions, number of microphones, source trajectory, etc., are
% implicitely defined by the set of transfer functions in 'TFFileName'. It
% is hence up to the user to ensure consistency between the resulting audio
% data and the environmental setup parameters that were used to compute the
% impulse responses.
%
% The resulting audio data, together with various simulation parameters, is
% stored on file under the name defined by the input variable 'AuFileName'.
% The suffix '.mat' is also automatically appended to the given file name,
% which can be a full access path to the desired file. If no access path is
% given, the function saves the resulting audio data in the current working
% directory. The resulting audio samples are stored in a single matrix of
% multi-channel data 'AuData', each column containing the data generated
% for the corresponding microphone.
%
% The input variable 'SrcSignal' is a 1-dimensional vector of audio data
% corresponding to the signal emitted by the source. Note that the sampling
% frequency of this signal must be in agreement with the sampling frequency
% used when computing the transfer functions. The length of the source
% audio sample, together with the source trajectory defined when computing
% the impulse responses, define the velocity of the speaker across the
% environment.
%
% This function also accepts an optional parameter to set the level of
% white noise added to the microphone signals. Set this parameter by
% defining a value (in dB) for 'noiseSNR', otherwise it will default to 20
% dB.
%
% This function creates the audio data by splitting the source data
% 'SrcSignal' into as many (non-overlapping) frames as the number of source
% trajectory points. Each frame of signal is then convolved with the
% corresponding transfer function, and the convolution results are then
% combined additively to generate the microphone signal. This process is
% repeated for each trajectory point and each microphone. The assumption is
% that the source remains stationary during each frame.
%
% You can influence the source path by setting 'DirString' to one of the 
% following strings:
%   'SE':  (default) the source simply follows the path defined in 
%          'TFFileName' from start to end
%   'ES':  the source follows the path backwards, from end to start
%   'SMS': the source follows half the path from the start to middle point, 
%          then returns back to the start position
%   'EME': the source follows half the path from the end to middle point, 
%          then returns back to the end position
%
% Set 'SaveFlag' value to 0 to disable saving the results on file (enabled
% by default). If set to 0, the parameter 'AuFileName' is simply ignored.
%
% Set 'NoiseFlag' value to 0 to avoid additive noise in the resulting 
% signals (enabled by default). If set to 0, the 'SNRval' parameter is 
% ignored.

% Release date: 28 October 2005
% Author: Eric A. Lehmann, WATRI, Perth, Australia (www.watri.org.au)
%         Eric.Lehmann@watri.org.au -- Tel. +61 (0)8 6488 4642
%
% Copyright 2005 Eric A. Lehmann
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

VarList = {'SNRval'     20;            % Desired SNR in the resulting signals (in dB)
           'AddNoise'   1;             % Set to 0 if you don't want additive noise
           'Dir'        'SE';          % direction of the source
           'SilentFlag' 0;             % set to 1 for silent behaviour.
           'SaveRes'    1};            % set to 0 to disable saving to file
eval(SetUserVars(VarList,varargin));   % set user-definable variables

if min(size(SrcSignal))~=1,
    error('Source signal must be one-dimensional (single channel).');
end
SrcSignal = SrcSignal(:);   % make sure vector is column vector

if isempty(find(TFFileName=='\')) & isempty(find(TFFileName=='/')),     % TFFileName doesn't have access path (full or partial)
%     foo = which(TFFileName);
%     if strcmp(foo,''),
%         foo = which([TFFileName '.mat']);
%         if strcmp(foo,''),
%             error('Cannot locate TFFileName!');
%         end
%     end
%     TFFileName = foo;      % complete TFFileName with full access path 
    foo = pwd;
    if ~isempty(find(foo=='\')),      % complete TFFileName with full access path to current directory
        TFFileName = [foo '\' TFFileName];
    else
        TFFileName = [foo '/' TFFileName];
    end
end

foo = load(TFFileName);       % load TF bank and other parameters from file
TFcell = foo.TFcell;
SetupStruc = foo.SetupStruc;
LimAtten = foo.LimAtten;
Fs = SetupStruc.Fs;
clear foo;

nSamp = length(SrcSignal);  % total number of samples in the audio data
nMics = size(TFcell,1);     % number of mics
nFrames = size(TFcell,2);   % number of trajectory points

if strcmp(Dir,'ES'),
    TFcell = TFcell(:,end:-1:1);
elseif strcmp(Dir,'SMS'),
    TFcell = TFcell(:,[1:floor(nFrames/2) ceil(nFrames/2):-1:1]);
elseif strcmp(Dir,'EME'),
    TFcell = TFcell(:,[end:-1:ceil(nFrames/2+1) floor(nFrames/2+1):end]);
end

nSampPerFrame = floor(nSamp/nFrames);	% Discard the last sample values if not enough for one frame.
nSamp = nFrames*nSampPerFrame; 		% Redefine nSamp to be the exact number of sample values used.
SrcSignal = SrcSignal(1:nSamp);

AuData = zeros(nSamp,nMics);
if ~SilentFlag, PrintLoopPCw('   [MakeAuData] Computing audio data. '); end;
for tt=1:nFrames,
   FrameStartInd = (tt-1)*nSampPerFrame+1;	% Start/end indices of the current frame in the overall audio sample.
   FrameStopInd = tt*nSampPerFrame;
   FrameData = SrcSignal(FrameStartInd:FrameStopInd);   % get one frame of data
   for mm=1:nMics,    % Compute the received signal by convolving the source signal with the IR
      if ~SilentFlag, PrintLoopPCw((tt-1)*nMics+mm,nFrames*nMics); end;
      hh = TFcell{mm,tt};
      TFlen = length(hh);       % length of current TF (can be variable!)
      EndIndex = FrameStopInd+TFlen-1;  % max length for the current convolution.
      if EndIndex<=nSamp,       % current end signal index is within AuData limits
          AuData(FrameStartInd:EndIndex,mm) = AuData(FrameStartInd:EndIndex,mm) + conv(hh,FrameData);
      else                      % otherwise, only used convolution results within limits
          foo = conv(hh,FrameData);
          AuData(FrameStartInd:end,mm) = AuData(FrameStartInd:end,mm) + foo(1:end-EndIndex+nSamp);
      end
   end
end

if AddNoise,
    av_pow = mean( sum(AuData.^2,1)/nSamp );	% Average mic power across all received signals.
    sigma_noise = sqrt( av_pow/(10^(SNRval/10)) );		% STD of white noise component to achieve desired SNR.
    AuData = AuData + sigma_noise*randn(size(AuData));	% Add some random noise
end

%-=:=- Save results into .mat file -=:=-
if SaveRes,
    if isempty(find(AuFileName=='\')) & isempty(find(AuFileName=='/')),     % AuFileName doesn't have access path (full or partial)
        foo = pwd;
        if ~isempty(find(foo=='\')),      % complete AuFileName with full access path to current directory
            AuFileName = [foo '\' AuFileName];
        else
            AuFileName = [foo '/' AuFileName];
        end
    end

    ii = 1;         % checks for existing file and alter name if necessary
    tmp = AuFileName;
    while exist([tmp '.mat'],'file')==2,
        tmp = [AuFileName '_' num2str(ii)];   % if existing file detected, increment and add suffix
        ii = ii+1;
    end
    if ii~=1,
        fprintf('   [MakeAuData] *Warning* File name altered to avoid overwriting existing results\n');
    end
    AuFileName = [tmp '.mat'];

    save(AuFileName,'AuData','TFFileName','SetupStruc','LimAtten','SNRval','Dir','AddNoise');
    if ~SilentFlag, fprintf('   [MakeAuData] Audio data saved in file ''%s''\n',AuFileName); end;
end
