function [TFcell] = MakeTrajIMTFs(SetupStr,TFFileName,varargin)
%MakeTrajIMTFs  generate transfer functions for a specific setup using the Image Method
%
% [TFcell] = MakeTrajIMTFs(SetupStr,TFFileName)
%
% This function generates a bank of transfer functions for a particular
% user-defined room setup using the image method for simulating small room
% acoustics. The input string 'SetupStr' corresponds to the name (without
% '.m' extension) of the m-function containing the desired definitions for
% the environmental setup and source trajectory, including parameters such
% as sampling frequency, room dimensions and reflection coefficients (see
% 'IMSetup.m' for an example).
%
% This function returns a 2-dimensional cell array 'TFcell' containing the
% transfer functions for each source trajectory points and each
% microphones, organised as follows: TFcell{mic_index,traj_index}. 
% The filter length might differ in each transfer function.
%
% This function also saves the computation results on file. The argument
% 'TFFileName' determines the name of the .mat file where the variable
% 'TFcell' is to be saved. If a file already exists with the same name
% given as argument, a suffix index will be added in order to avoid
% overwriting existing files. The suffix '.mat' is also automatically
% appended to the given file name, which can be a full access path to the
% desired file. If no access path is given, the file is saved in the
% current working directory. 
%
% [TFcell] = MakeTrajIMTFs(SetupStr,TFFileName,'LimAtten',LimAttenVal)
%
% Passing the value 'LimAttenVal' to the function allows to set the maximum
% attenuation in dB (compared to the direct path amplitude) determining
% which image sources to include in the TF computations (see function
% 'MakeIMResp.m'). Defaults to 45 dB if no value is defined.
%
% For reference, the parameters 'SetupStr' and 'LimAtten' are also saved in
% the resulting data file 'TFFileName' (together with the TF bank variable 
% 'TFcell').

% Release date: 17 August 2005
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

VarList = {'LimAtten'  45};       % maximum attenuation in IM transfer function computation
eval(SetUserVars(VarList,varargin));   % set user-definable variables

eval(['SetupStruc = ' SetupStr ';']);    % fetch parameters for desired setup
Fs = SetupStruc.Fs;
cc = SetupStruc.c;
room = SetupStruc.room;
micpos = SetupStruc.micpos;
beta = SetupStruc.beta;
straj = SetupStruc.straj;
if isfield(SetupStruc,'measT60'),
    measT60 = SetupStruc.measT60;
else
    measT60 = [];
end

nMics = size(micpos,1);     % number of microphones
nSPts = size(straj,1);      % number of source trajectory points

TFcell = cell(nMics,nSPts); % pre-allocate cell array
PrintLoopPCw('   [MakeTrajIMTFs] Computing transfer functions. ');
for mm=1:nMics,
    X_rcv = micpos(mm,:);
    for tt=1:nSPts,         % compute IM transfer function for each source-receiver combinations
        PrintLoopPCw((mm-1)*nSPts+tt,nMics*nSPts);
        X_src = straj(tt,:);
        TFcell{mm,tt} = MakeIMResp(Fs,beta,X_src,X_rcv,room,cc,LimAtten,measT60,1);
    end
end

%-=:=- Save results into .mat file -=:=-
if isempty(find(TFFileName=='\')) & isempty(find(TFFileName=='/')),     % TFFileName doesn't have access path (full or partial)
    foo = pwd;
    if ~isempty(find(foo=='\')),      % complete TFFileName with full access path to current directory
        TFFileName = [foo '\' TFFileName];
    else
        TFFileName = [foo '/' TFFileName];
    end
end

ii = 1;         % checks for existing file and alter name if necessary
tmp = TFFileName;
while exist([tmp '.mat'],'file')==2,
    tmp = [TFFileName '_' num2str(ii)];   % if existing file detected, increment and add suffix
    ii = ii+1;
end
if ii~=1,
    fprintf('   [MakeTrajIMTFs] *Warning* File name altered to avoid overwriting existing results *\n');
end
TFFileName = [tmp '.mat'];
    
save(TFFileName,'TFcell','SetupStruc','LimAtten');
fprintf('   [MakeTrajIMTFs] TF bank saved in file ''%s''\n',TFFileName);
