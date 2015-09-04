function [OutStruc] = IMSetup()
%[OutStruc] = IMSetup()  environmental parameters for image method simulation
%
% Edit this function example to define the different parameters of your
% image method simulation. Returns the structure 'OutStruc' with the 
% following fields:
%
%     Fs: sampling frequency
%      c: sound velocity
%   room: 3 element vector of enclosure dimensions [x_length y_length z_length]
% micpos: N by 3 matrix, [x y z] positions of N microphones
%   beta: 6 element vector, wall reflection coefficients
%  straj: M by 3 matrix, [x y z] positions of M source trajectory points

OutStruc.Fs = 16000;                % sampling frequency in Hz
OutStruc.c = 343;                   % propagation speed of acoustic waves in m/s

OutStruc.room = [5  5  2.7];        % room dimensions in m
OutStruc.micpos = [1  2.3  0.1;     % [x y z] position of array microphones in m
                   1  2.7  0.1;
                   2.3 4 0.1;
                   2.7 4 0.1;
                   4 2.7 0.1;
                   4 2.3 0.1;
                   2.7 1 0.1;
                   2.3 1 0.1];
      
OutStruc.beta = ones(1,6)*0.1;    % reflection coefficients in range [0 ... 1[

OutStruc.straj = [linspace(1.85,3.05,50).'  linspace(1.85,3.05,50).'  0.1*ones(50,1)];     
                                    % [x y z] positions in source trajectory in m
                                    % straight line in front of mic pair, 50 source points
OutStruc.measT60 = 0.10987;
%-=:=- -=:=- -=:=- -=:=- -=:=- -=:=- -=:=- -=:=- -=:=- -=:=- -=:=- -=:=- -=:=-
%-=:=- Results from reverberation time analysis ('IMRevTimeAnalysis.m')  -=:=-
%
% Simulation parameters:  Fs =   16000
%                         c =    343
%                         room = [5  5  2.7]
%                         beta = [0.5  0.5  0.5  0.5  0.5  0.5]
%
% Median T60 measured in room, over 50 configurations:
%  > T60 = 0.09956 s (inter-quartile range: 0.0041759 s)
%
% Median T60 measured in room (Schroeder integration method), over 50 configurations:
%  > T60 = 0.10757 s (inter-quartile range: 0.0051386 s)
%
% Median T60 estimated from transfer function, over 50 configurations: 
%  > T60 = 0.10987 s (inter-quartile range: 0.0058195 s)
%
% Sabine's T60 estimate: T60 = 0.13933 s
% Eyring's T60 estimate: T60 = 0.075377 s
%-=:=- -=:=- -=:=- -=:=- -=:=- -=:=- -=:=- -=:=- -=:=- -=:=- -=:=- -=:=- -=:=-


OutStruc.measT60 = 0.10972;


%-=:=- -=:=- -=:=- -=:=- -=:=- -=:=- -=:=- -=:=- -=:=- -=:=- -=:=- -=:=- -=:=-
%-=:=- Results from reverberation time analysis ('IMRevTimeAnalysis.m')  -=:=-
%
% Simulation parameters:  Fs =   16000
%                         c =    343
%                         room = [5  5  2.7]
%                         beta = [0.5  0.5  0.5  0.5  0.5  0.5]
%
% Median T60 measured in room, over 50 configurations:
%  > T60 = 0.10087 s (inter-quartile range: 0.0052682 s)
%
% Median T60 measured in room (Schroeder integration method), over 50 configurations:
%  > T60 = 0.11051 s (inter-quartile range: 0.007997 s)
%
% Median T60 estimated from transfer function, over 50 configurations: 
%  > T60 = 0.10972 s (inter-quartile range: 0.0053356 s)
%
% Sabine's T60 estimate: T60 = 0.13933 s
% Eyring's T60 estimate: T60 = 0.075377 s
%-=:=- -=:=- -=:=- -=:=- -=:=- -=:=- -=:=- -=:=- -=:=- -=:=- -=:=- -=:=- -=:=-
OutStruc.measT60 = 0.10906;
%-=:=- -=:=- -=:=- -=:=- -=:=- -=:=- -=:=- -=:=- -=:=- -=:=- -=:=- -=:=- -=:=-
%-=:=- Results from reverberation time analysis ('IMRevTimeAnalysis.m')  -=:=-
%
% Simulation parameters:  Fs =   16000
%                         c =    343
%                         room = [5  5  2.7]
%                         beta = [0.5  0.5  0.5  0.5  0.5  0.5]
%
% Median T60 measured in room, over 50 configurations:
%  > T60 = 0.09854 s (inter-quartile range: 0.0039346 s)
%
% Median T60 measured in room (Schroeder integration method), over 50 configurations:
%  > T60 = 0.10687 s (inter-quartile range: 0.0074981 s)
%
% Median T60 estimated from transfer function, over 50 configurations: 
%  > T60 = 0.10906 s (inter-quartile range: 0.0063426 s)
%
% Sabine's T60 estimate: T60 = 0.13933 s
% Eyring's T60 estimate: T60 = 0.075377 s
%-=:=- -=:=- -=:=- -=:=- -=:=- -=:=- -=:=- -=:=- -=:=- -=:=- -=:=- -=:=- -=:=-
OutStruc.measT60 = 0.10837;
%-=:=- -=:=- -=:=- -=:=- -=:=- -=:=- -=:=- -=:=- -=:=- -=:=- -=:=- -=:=- -=:=-
%-=:=- Results from reverberation time analysis ('IMRevTimeAnalysis.m')  -=:=-
%
% Simulation parameters:  Fs =   16000
%                         c =    343
%                         room = [5  5  2.7]
%                         beta = [0.5  0.5  0.5  0.5  0.5  0.5]
%
% Median T60 measured in room, over 50 configurations:
%  > T60 = 0.10178 s (inter-quartile range: 0.0050022 s)
%
% Median T60 measured in room (Schroeder integration method), over 50 configurations:
%  > T60 = 0.10900 s (inter-quartile range: 0.0061421 s)
%
% Median T60 estimated from transfer function, over 50 configurations: 
%  > T60 = 0.10837 s (inter-quartile range: 0.0067357 s)
%
% Sabine's T60 estimate: T60 = 0.13933 s
% Eyring's T60 estimate: T60 = 0.075377 s
%-=:=- -=:=- -=:=- -=:=- -=:=- -=:=- -=:=- -=:=- -=:=- -=:=- -=:=- -=:=- -=:=-


OutStruc.measT60 = 0.10868;


%-=:=- -=:=- -=:=- -=:=- -=:=- -=:=- -=:=- -=:=- -=:=- -=:=- -=:=- -=:=- -=:=-
%-=:=- Results from reverberation time analysis ('IMRevTimeAnalysis.m')  -=:=-
%
% Simulation parameters:  Fs =   16000
%                         c =    343
%                         room = [5  5  2.7]
%                         beta = [0.5  0.5  0.5  0.5  0.5  0.5]
%
% Median T60 measured in room, over 50 configurations:
%  > T60 = 0.09792 s (inter-quartile range: 0.0035388 s)
%
% Median T60 measured in room (Schroeder integration method), over 50 configurations:
%  > T60 = 0.10433 s (inter-quartile range: 0.0040081 s)
%
% Median T60 estimated from transfer function, over 50 configurations: 
%  > T60 = 0.10868 s (inter-quartile range: 0.0052061 s)
%
% Sabine's T60 estimate: T60 = 0.13933 s
% Eyring's T60 estimate: T60 = 0.075377 s
%-=:=- -=:=- -=:=- -=:=- -=:=- -=:=- -=:=- -=:=- -=:=- -=:=- -=:=- -=:=- -=:=-


OutStruc.measT60 = 0.11002;


%-=:=- -=:=- -=:=- -=:=- -=:=- -=:=- -=:=- -=:=- -=:=- -=:=- -=:=- -=:=- -=:=-
%-=:=- Results from reverberation time analysis ('IMRevTimeAnalysis.m')  -=:=-
%
% Simulation parameters:  Fs =   16000
%                         c =    343
%                         room = [5  5  2.7]
%                         beta = [0.5  0.5  0.5  0.5  0.5  0.5]
%
% Median T60 measured in room, over 50 configurations:
%  > T60 = 0.10376 s (inter-quartile range: 0.0048619 s)
%
% Median T60 measured in room (Schroeder integration method), over 50 configurations:
%  > T60 = 0.11520 s (inter-quartile range: 0.0083013 s)
%
% Median T60 estimated from transfer function, over 50 configurations: 
%  > T60 = 0.11002 s (inter-quartile range: 0.0089278 s)
%
% Sabine's T60 estimate: T60 = 0.13933 s
% Eyring's T60 estimate: T60 = 0.075377 s
%-=:=- -=:=- -=:=- -=:=- -=:=- -=:=- -=:=- -=:=- -=:=- -=:=- -=:=- -=:=- -=:=-


OutStruc.measT60 = 0.11039;


%-=:=- -=:=- -=:=- -=:=- -=:=- -=:=- -=:=- -=:=- -=:=- -=:=- -=:=- -=:=- -=:=-
%-=:=- Results from reverberation time analysis ('IMRevTimeAnalysis.m')  -=:=-
%
% Simulation parameters:  Fs =   16000
%                         c =    343
%                         room = [5  5  2.7]
%                         beta = [0.5  0.5  0.5  0.5  0.5  0.5]
%
% Median T60 measured in room, over 50 configurations:
%  > T60 = 0.10774 s (inter-quartile range: 0.0065984 s)
%
% Median T60 measured in room (Schroeder integration method), over 50 configurations:
%  > T60 = 0.12663 s (inter-quartile range: 0.032545 s)
%
% Median T60 estimated from transfer function, over 50 configurations: 
%  > T60 = 0.11039 s (inter-quartile range: 0.0045311 s)
%
% Sabine's T60 estimate: T60 = 0.13933 s
% Eyring's T60 estimate: T60 = 0.075377 s
%-=:=- -=:=- -=:=- -=:=- -=:=- -=:=- -=:=- -=:=- -=:=- -=:=- -=:=- -=:=- -=:=-


OutStruc.measT60 = 0.11078;


%-=:=- -=:=- -=:=- -=:=- -=:=- -=:=- -=:=- -=:=- -=:=- -=:=- -=:=- -=:=- -=:=-
%-=:=- Results from reverberation time analysis ('IMRevTimeAnalysis.m')  -=:=-
%
% Simulation parameters:  Fs =   16000
%                         c =    343
%                         room = [5  5  2.7]
%                         beta = [0.5  0.5  0.5  0.5  0.5  0.5]
%
% Median T60 measured in room, over 50 configurations:
%  > T60 = 0.11041 s (inter-quartile range: 0.018473 s)
%
% Median T60 measured in room (Schroeder integration method), over 50 configurations:
%  > T60 = 0.14719 s (inter-quartile range: 0.050365 s)
%
% Median T60 estimated from transfer function, over 50 configurations: 
%  > T60 = 0.11078 s (inter-quartile range: 0.0068627 s)
%
% Sabine's T60 estimate: T60 = 0.13933 s
% Eyring's T60 estimate: T60 = 0.075377 s
%-=:=- -=:=- -=:=- -=:=- -=:=- -=:=- -=:=- -=:=- -=:=- -=:=- -=:=- -=:=- -=:=-


OutStruc.measT60 = 0.11061;


%-=:=- -=:=- -=:=- -=:=- -=:=- -=:=- -=:=- -=:=- -=:=- -=:=- -=:=- -=:=- -=:=-
%-=:=- Results from reverberation time analysis ('IMRevTimeAnalysis.m')  -=:=-
%
% Simulation parameters:  Fs =   16000
%                         c =    343
%                         room = [5  5  2.7]
%                         beta = [0.5  0.5  0.5  0.5  0.5  0.5]
%
% Median T60 measured in room, over 50 configurations:
%  > T60 = 0.11195 s (inter-quartile range: 0.0082793 s)
%
% Median T60 measured in room (Schroeder integration method), over 50 configurations:
%  > T60 = 0.11826 s (inter-quartile range: 0.010749 s)
%
% Median T60 estimated from transfer function, over 50 configurations: 
%  > T60 = 0.11061 s (inter-quartile range: 0.0084476 s)
%
% Sabine's T60 estimate: T60 = 0.13933 s
% Eyring's T60 estimate: T60 = 0.075377 s
%-=:=- -=:=- -=:=- -=:=- -=:=- -=:=- -=:=- -=:=- -=:=- -=:=- -=:=- -=:=- -=:=-


OutStruc.measT60 = 0.15436;


%-=:=- -=:=- -=:=- -=:=- -=:=- -=:=- -=:=- -=:=- -=:=- -=:=- -=:=- -=:=- -=:=-
%-=:=- Results from reverberation time analysis ('IMRevTimeAnalysis.m')  -=:=-
%
% Simulation parameters:  Fs =   16000
%                         c =    343
%                         room = [5  5  2.7]
%                         beta = [0.1  0.1  0.1  0.1  0.1  0.1]
%
% Median T60 measured in room, over 50 configurations:
%  > T60 = 0.06201 s (inter-quartile range: 0.014684 s)
%
% Median T60 measured in room (Schroeder integration method), over 50 configurations:
%  > T60 = 0.12591 s (inter-quartile range: 0.039555 s)
%
% Median T60 estimated from transfer function, over 50 configurations: 
%  > T60 = 0.15436 s (inter-quartile range: 0.067565 s)
%
% Sabine's T60 estimate: T60 = 0.10555 s
% Eyring's T60 estimate: T60 = 0.022691 s
%-=:=- -=:=- -=:=- -=:=- -=:=- -=:=- -=:=- -=:=- -=:=- -=:=- -=:=- -=:=- -=:=-

