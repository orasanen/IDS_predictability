function drawstds(handle,xvals,means,devs,barwidth,linewidth,color,horshift)
% function drawstds(handle,xvals,means,devs,barwidth,linewidth,color)
%
% Draws standard deviation bars around the given mean curve
%
% Inputs:
%
%   handle  : image handle
%   xvals   : x-axis values
%   means   : y-axis mean values
%   devs    : STDs corresponding each mean value
%   barwidth: width of the horizontal STD bars (default 2)
%   linewidth: thickness of the STD markings (default 2)
%   color   : color of the STD markings (default black)

if nargin < 5
    barwidth = 2;
end
if nargin <6
    linewidth = 2;
end

if nargin <7
    color = 'black';
end

if nargin <8
    horshift = 0;
end
  

figure(handle);
hold on;

for k = 1:length(means)    
    
    line([xvals(k)+horshift xvals(k)+horshift],[means(k)-devs(k) means(k)+devs(k)],'LineWidth',linewidth,'Color',color);    
    line([xvals(k)-barwidth+horshift xvals(k)+barwidth+horshift],[means(k)-devs(k) means(k)-devs(k)],'LineWidth',linewidth,'Color',color);
    line([xvals(k)-barwidth+horshift xvals(k)+barwidth+horshift],[means(k)+devs(k) means(k)+devs(k)],'LineWidth',linewidth,'Color',color);
end

