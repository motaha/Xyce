function [H] = plot_NEURONPOP_data(data, dt)
%
% *************************************************************************
% plot_NEURONPOP_data.m  --  Programmer:  Stephen J. Verzi (03-15-11).
% *************************************************************************
%
% [H] = plot_NEURONPOP_data(data, {dt})
%
% Plot output data from NEURONPOP netlist (Xyce) sequentially (i.e., in a
% movie-like fashion).
%
% Input Arguments:
%   data     - MATLAB data structure containing data in the output file of
%              a NEURONPOP netlist (see load_NEURONPOP_data).
%   dt       - delay time between each plot in sequence.  Optional, if not
%              present or empty, 1/4 second is used.
%
% Output Arguments:
%   H        - handle to MATLAB plot window (optional).
%
% Notes:
%   Currently time is assumed to be in seconds (for title display).
%

H = [];

msg = nargchk(1, 2, nargin);
if (~isempty(msg) || isempty(data))
   warning(['plot_NEURONPOP_data: ', msg]);
   help plot_NEURONPOP_data;
   return;
end;

if ((nargin < 2) || isempty(dt))
   dt = 0.25;
end;

for ii=1:length(data)
   plot3(data(ii).x, data(ii).y, data(ii).voltage, '*');
   tstr = ['voltages at time ', num2str(data(ii).time), ' (sec) - ', ...
       num2str(data(ii).number_neurons), ' neurons'];
   title(tstr);
   xlabel('x position');
   ylabel('y position');
   zlabel('voltage');
   pause(dt);
end;

end

