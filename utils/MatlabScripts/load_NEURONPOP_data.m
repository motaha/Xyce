function [data] = load_NEURONPOP_data(filename)
%
% *************************************************************************
% load_NEURONPOP_data.m  --  Programmer:  Stephen J. Verzi (03-15-11).
% *************************************************************************
%
% [data] = load_NEURONPOP_data(filename)
%
% Load output data from a file resulting from a NEURONPOP netlist (Xyce).
%
% Input Arguments:
%   filename - Name of the NEURONPOP output file.
%
% Output Arguments:
%   data     - MATLAB data structure containing data in the output file of
%              a NEURONPOP netlist.
%

data = [];

msg = nargchk(1, 1, nargin);
if (~isempty(msg) || isempty(filename))
   warning(['load_NEURONPOP_data: ', msg]);
   help load_NEURONPOP_data;
   return;
end;

% NOTE: both of the following need to be coordinated with NEURONPOP (Xyce)
MAX_NIC = 2; % maximum number of internal connections allowed per neuron
MAX_NEC = 2; % maximum number of external connections allowed per neuron

% NOTE: the following loop has not been optimized for speed
data.time = 0;
data.number_neurons = 0;
data.x = [];
data.y = [];
data.voltage = [];
data.ics = [];
data.ecs = [];
kk = 0;
fid = fopen(filename, 'r');
   while (~feof(fid))
      fline = fgetl(fid);
      metadata = sscanf(fline, '%g,', 2);
      t = metadata(1);
      N = metadata(2);
      ldata = sscanf(fline, '%g,', ((3 + MAX_NIC + MAX_NEC)*N + 2));
      ldata = reshape(ldata(3:end), N, (3 + MAX_NIC + MAX_NEC));
      kk = kk + 1;
      data(kk).time = t;
      data(kk).number_neurons = N;
      data(kk).x = ldata(:,1);
      data(kk).y = ldata(:,2);
      data(kk).voltage = ldata(:,3);
      data(kk).ics = ldata(:,4:(4+MAX_NIC-1));
      data(kk).ecs = ldata(:,(4+MAX_NIC):end);
   end;
fclose(fid);

end
