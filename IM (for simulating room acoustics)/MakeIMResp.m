function [FilterCoeffs,SpectrumVec] = MakeIMResp(Fs,beta,X_src,X_rcv,room,cc,LimAtten_dB,measT60,silentflag)
%MakeIMResp  impulse response based on image method calculations
%
% [h,H] = MakeIMRes(Fs,beta,source,sink,room,c,LimAtten,measT60)
% 
% This function generates the impulse response between a sound source and
% an acoustic receiver, based on various environmental parameters such as
% source and sensor positions, enclosure's dimension and reflection
% coefficients, etc. Input parameters are defined as follows:
%
%       Fs: scalar, sampling frequency in hertz. Eg: 8000.
%     beta: vector of dimension 6, corresponding to each wall's reflection
%           coefficient: [x1 x2 y1 y2 z1 z2]. Index 1 indicates wall closest
%           to origin. Set to [0 0 0 0 0 0] to obtain anechoic response (only 
%           direct path). Eg: [0.75 0.75 0.85 0.25 0.3 0.9].
%   source: vector dimension 3, indicating the location of the source in
%           space: [x y z]. Eg: [1 1 1.5].
%     sink: vector of dimension 3, indicating the location of the microphone
%           in space: [x y z]. Eg: [2 2 1.5].
%     room: vector of dimension 3, indicating the rectangular room dimensions:
%           [x_length y_length z_length]. Eg: [4 4 3].
%        c: scalar, speed of acoustic waves. Eg: 342.
% LimAtten: scalar, parameter in dB determining how much of the resulting
%           transfer function is cropped: the impulse response is computed
%           until the time index where its overall energy content has
%           decreased by 'LimAtten' dB, after which the computations stop.
%           Not relevant if beta=zeros(1,6). Eg: 50.
%  measT60: scalar, measured reverberation time in seconds. The above beta
%           parameters typically generate a reverberation time in the TF
%           that is slightly different from the Sabine or Eyring formula
%           using the given beta values. If the T60 value has been
%           practically measured (e.g. using IMRevTimeAnalysis.m) for the
%           current setup, set 'measT60' to this value for a more accurate
%           computation of the TF. Otherwise leave empty, i.e. define as [].
% 
% This function returns the time coefficients of the filter (transfer
% function) in 'h' and its corresponding frequency response 'H'. The filter
% coefficients are real and non-normalised. The first value in vector 'h',
% h(1), corresponds to time t=0. The number of coefficients returned is
% variable and results from the value of 'LimAtten' defined by the user:
% the filter length will be as large as necessary to capture all the
% relevant highest-order reflections. 
% 
% Implementation based on Allen and Berkley's Image Method for Efficiently
% Simulating Small-room Acoustics. See: J.B. Allen and D.A. Berkley, "Image
% method for efficiently simulating small-room acoustics", J. Acoust. Soc.
% Am., issue 65, number 4, April 1979. The computations are carried out in
% the frequency domain, which allows fractional delays for all reflections.

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


% Explanations for the following code -------------------------------------
% The following implementation of the image method principle has been
% speficically optimised for execution speed. The following code is based
% on the observation that for a particular dimension, the delays from the
% image sources to the receiver increases monotonically as the absolute
% value of the image index (m, n, or l) increases. Hence, all image sources
% for whose indices are above or equal to a specific limit index (for which
% the received delay is above the relevant cut-off value) can be discarded.
% The following code checks, for each dimension, the delay of each received
% path and automatically determines when to stop, thus avoiding unnecessary
% computations.
% This approach is more effective than a similar approach based on the
% amplitude level of the received image source signal. The latter approach
% will typically cut off image source that are below a given amplitude
% threshold, but it is possible that two such sources would still produce a
% combined amplitude larger than the threshold in the resulting TF (this
% might typically leads to significant errors in the tail of the TF). The
% technique used in the code below effectively crops the end of the TF,
% leaving the uncropped part of the TF unaltered (the amount of TF cropped
% depends on the 'LimAtten' parameter).
% The resulting number of considered image sources hence automatically
% results from environmental factors, such as the room dimensions, the
% source and sensor positions, and the walls' reflection coefficients. As a
% result, the length of the computed transfer function has an optimally
% minimum length (no extra padding with negligibly small values).
%--------------------------------------------------------------------------

global SpectrumVec FreqPoints      % not too pretty, but this avoids passing potentially large vectors to frequently called functions...

%---- Check user input:
if X_rcv(1)>=room(1) | X_rcv(2)>=room(2) | X_rcv(3)>=room(3) | X_rcv(1)<=0 | X_rcv(2)<=0 | X_rcv(3)<=0,
    error('Receiver must be within the room boundaries.');
elseif X_src(1)>=room(1) | X_src(2)>=room(2) | X_src(3)>=room(3) | X_src(1)<=0 | X_src(2)<=0 | X_src(3)<=0,
    error('Source must be within the room boundaries.');
elseif ~isempty(find(beta>=1)) | ~isempty(find(beta<0)),
    error('Parameter ''beta'' must be in range [0...1[.');
end

if nargin<9,
    silentflag = 0;     % set to 1 to disable on-screen messages
end

X_src = X_src(:);       % Source location
X_rcv = X_rcv(:);       % Receiver location
beta = beta(:);
Rr = 2*room(:);         % Room dimensions
DPdel = norm(X_rcv - X_src)/cc;         % direct path delay in [s]

%---- Define enough frequency points for resulting time impulse response:
if ~isequal(beta,zeros(6,1)),
    if isempty(measT60),    % if no practical T60 measurement available, use Sabine estimate
        V = prod(room);
        aa_sab = (2-beta(1)^2-beta(2)^2)*room(2)*room(3) + (2-beta(3)^2-beta(4)^2)*room(1)*room(3) + (2-beta(5)^2-beta(6)^2)*room(1)*room(2);
        T60val = 0.161*V/aa_sab;      % Sabine's reverberation time in [s]
    else
        T60val = measT60;   % practical T60 measurement determines real energy decay in TF!
    end
    foo = LimAtten_dB * T60val / 60;    % desired length of TF (TF decays by 60dB for T60 seconds after direct path delay)
    MaxDelay = DPdel + foo;             % maximum delay in TF: direct path plus TForder
else
    MaxDelay = 2*DPdel;         % anechoic case: allow for 2 times direct path in TF 
end
TForder = ceil(MaxDelay*Fs);    % total TF length

FreqPoints = linspace(0,Fs/2,TForder).';
SpectrumVec = zeros(TForder,1);

%---- summation over room dimensions:
if ~silentflag, fprintf('   [MakeIMResp] Computing transfer function '); end;
for a = 0:1
    for b = 0:1
        for d = 0:1
            if ~silentflag, fprintf('.'); end;
            
            m = 1;              % Check delay values for m=1 and above
            FoundLValBelowLim = Check_lDim(a,b,d,m,X_rcv,X_src,Rr,cc,MaxDelay,beta);
            while FoundLValBelowLim==1,
                m = m+1;
                FoundLValBelowLim = Check_lDim(a,b,d,m,X_rcv,X_src,Rr,cc,MaxDelay,beta);
            end
            
            m = 0;              % Check delay values for m=0 and below
            FoundLValBelowLim = Check_lDim(a,b,d,m,X_rcv,X_src,Rr,cc,MaxDelay,beta);
            while FoundLValBelowLim==1,
                m = m-1;
                FoundLValBelowLim = Check_lDim(a,b,d,m,X_rcv,X_src,Rr,cc,MaxDelay,beta);
            end

        end
    end
end
if ~silentflag, fprintf('\n'); end;

%---- Inverse Fourier transform:
SpectrumVec(1) = SpectrumVec(1)/2;      % remove DC component in resulting time coefficients.
freqvec = -i*2*pi*linspace(0,1/2,TForder);
FilterCoeffs = zeros(TForder,1);
for ii=1:TForder,
    freq = exp((ii-1)*freqvec);
    FilterCoeffs(ii) = real(freq*SpectrumVec);
end
FilterCoeffs = FilterCoeffs/TForder;


%============
function [FoundLValBelowLim] = Check_lDim(a,b,d,m,X_rcv,X_src,Rr,cc,MaxDelay,beta)

FoundLValBelowLim = 0;

l = 1;              % Check delay values for l=1 and above
FoundNValBelowLim = Check_nDim(a,b,d,l,m,X_rcv,X_src,Rr,cc,MaxDelay,beta);
while FoundNValBelowLim==1,
    l = l+1;
    FoundNValBelowLim = Check_nDim(a,b,d,l,m,X_rcv,X_src,Rr,cc,MaxDelay,beta);
end
if l~=1, FoundLValBelowLim = 1; end;

l = 0;              % Check delay values for l=0 and below
FoundNValBelowLim = Check_nDim(a,b,d,l,m,X_rcv,X_src,Rr,cc,MaxDelay,beta);
while FoundNValBelowLim==1,
    l = l-1;
    FoundNValBelowLim = Check_nDim(a,b,d,l,m,X_rcv,X_src,Rr,cc,MaxDelay,beta);
end
if l~=0, FoundLValBelowLim = 1; end;


%============
function [FoundNValBelowLim] = Check_nDim(a,b,d,l,m,X_rcv,X_src,Rr,cc,MaxDelay,beta)

global SpectrumVec FreqPoints

FoundNValBelowLim = 0;

n = 1;          % Check delay values for n=1 and above
dist = norm( [2*a-1; 2*b-1; 2*d-1].*X_src + X_rcv - Rr.*[n;l;m] );
foo_time = dist/cc;
while foo_time<=MaxDelay,    % if delay is below TF length limit for n=1, check n=2,3,4...
    foo_amplitude = prod(beta.^abs([n-a; n; l-b; l; m-d; m])) / (4*pi*dist);
    SpectrumVec = SpectrumVec + foo_amplitude * exp(i*2*pi*foo_time*FreqPoints);
    n = n+1;
    dist = norm( [2*a-1; 2*b-1; 2*d-1].*X_src + X_rcv - Rr.*[n;l;m] );
    foo_time = dist/cc;
end
if n~=1, FoundNValBelowLim = 1; end;

n = 0;          % Check delay values for n=0 and below
dist = norm( [2*a-1; 2*b-1; 2*d-1].*X_src + X_rcv - Rr.*[n;l;m] );
foo_time = dist/cc;
while foo_time<=MaxDelay,    % if delay is below TF length for n=0, check n=-1,-2,-3...
    foo_amplitude = prod(beta.^abs([n-a; n; l-b; l; m-d; m])) / (4*pi*dist);
    SpectrumVec = SpectrumVec + foo_amplitude * exp(i*2*pi*foo_time*FreqPoints);
    n = n-1;
    dist = norm( [2*a-1; 2*b-1; 2*d-1].*X_src + X_rcv - Rr.*[n;l;m] );
    foo_time = dist/cc;
end
if n~=0, FoundNValBelowLim = 1; end;

