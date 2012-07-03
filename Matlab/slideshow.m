
function slideshow(directory,firstFrame,lastFrame)
% SLIDESHOW(DIRECTORY,FIRSTFRAME,LASTFRAME) shows the 3D-PTV images in the directory in 4
% subplots between first and last frames, or all.
%
% Example:
%       directory = '/Users/alex/Desktop/GUI/pyptv/experiments/exp1/img';
%       firstFrame = 10160;
%       lastFrame = 10003;
%       slideshow(directory, firstFrame, lastFrame);
%
%

nCam = 4;


if ~nargin
   directory = 'c:\PTV\origo\working_folder_Dumbbell_b_10_08\img_35\';
   firstFrame =351106;
   lastFrame = 351374;
end

d = dir(fullfile(directory,'cam*.*'));
disp d
% clean out _targets files
ind_to_remove = zeros(1,length(d));
for i = 1:length(d)
   if findstr(d(i).name,'_targets')
       ind_to_remove(i) = 1;
   end
end
d = d(~ind_to_remove);




nFrames = floor(length(d)/nCam);

figure;
set(gcf, 'Position', get(0,'Screensize'));
ha = tight_subplot(2,2,[.035 .035],[.035 .035],[.035 .035])
for i = 1:nFrames
   for j = 0:nCam-1
       fileName = fullfile(directory,d(i+j*nFrames).name);
       frameNum = str2double(fileName(findstr(fileName,'.')+1:end));
       if frameNum >= firstFrame && frameNum <= lastFrame
           im = imread(fileName);
           axes(ha(j+1)); subimage(im); axis ij %#ok<LAXES>
           if j < 2, title(sprintf('Frame %d',frameNum)); end
       end
   end
   drawnow
end


function ha = tight_subplot(Nh, Nw, gap, marg_h, marg_w)

% tight_subplot creates "subplot" axes with adjustable gaps and margins
%
% ha = tight_subplot(Nh, Nw, gap, marg_h, marg_w)
%
%   in:  Nh      number of axes in hight (vertical direction)
%        Nw      number of axes in width (horizontaldirection)
%        gap     gaps between the axes in normalized units (0...1)
%                   or [gap_h gap_w] for different gaps in height and width
%        marg_h  margins in height in normalized units (0...1)
%                   or [lower upper] for different lower and upper margins
%        marg_w  margins in width in normalized units (0...1)
%                   or [left right] for different left and right margins
%
%  out:  ha     array of handles of the axes objects
%                   starting from upper left corner, going row-wise as in
%                   going row-wise as in
%
%  Example: ha = tight_subplot(3,2,[.01 .03],[.1 .01],[.01 .01])
%           for ii = 1:6; axes(ha(ii)); plot(randn(10,ii)); end
%           set(ha(1:4),'XTickLabel',''); set(ha,'YTickLabel','')

% Pekka Kumpulainen 20.6.2010   @tut.fi
% Tampere University of Technology / Automation Science and Engineering


if nargin<3; gap = .02; end
if nargin<4 || isempty(marg_h); marg_h = .05; end
if nargin<5; marg_w = .05; end

if numel(gap)==1;
   gap = [gap gap];
end
if numel(marg_w)==1;
   marg_w = [marg_w marg_w];
end
if numel(marg_h)==1;
   marg_h = [marg_h marg_h];
end

axh = (1-sum(marg_h)-(Nh-1)*gap(1))/Nh;
axw = (1-sum(marg_w)-(Nw-1)*gap(2))/Nw;

py = 1-marg_h(2)-axh;

ha = zeros(Nh*Nw,1);
ii = 0;
for ih = 1:Nh
   px = marg_w(1);

   for ix = 1:Nw
       ii = ii+1;
       ha(ii) = axes('Units','normalized', ...
           'Position',[px py axw axh], ...
           'XTickLabel','', ...
           'YTickLabel','');
       px = px+axw+gap(2);
   end
   py = py-axh-gap(1);
end

function [ax,h]=suplabel(text,whichLabel,supAxes)
% PLaces text as a title, xlabel, or ylabel on a group of subplots.
% Returns a handle to the label and a handle to the axis.
%  [ax,h]=suplabel(text,whichLabel,supAxes)
% returns handles to both the axis and the label.
%  ax=suplabel(text,whichLabel,supAxes)
% returns a handle to the axis only.
%  suplabel(text) with one input argument assumes whichLabel='x'
%
% whichLabel is any of 'x', 'y', 'yy', or 't', specifying whether the
% text is to be the xlable, ylabel, right side y-label,
% or title respectively.
%
% supAxes is an optional argument specifying the Position of the
%  "super" axes surrounding the subplots.
%  supAxes defaults to [.08 .08 .84 .84]
%  specify supAxes if labels get chopped or overlay subplots
%
% EXAMPLE:
%  subplot(2,2,1);ylabel('ylabel1');title('title1')
%  subplot(2,2,2);ylabel('ylabel2');title('title2')
%  subplot(2,2,3);ylabel('ylabel3');xlabel('xlabel3')
%  subplot(2,2,4);ylabel('ylabel4');xlabel('xlabel4')
%  [ax1,h1]=suplabel('super X label');
%  [ax2,h2]=suplabel('super Y label','y');
%  [ax3,h2]=suplabel('super Y label (right)','yy');
%  [ax4,h3]=suplabel('super Title'  ,'t');
%  set(h3,'FontSize',30)
%
% SEE ALSO: text, title, xlabel, ylabel, zlabel, subplot,
%           suptitle (Matlab Central)

% Author: Ben Barrowes <barrowes@alum.mit.edu>

%modified 3/16/2010 by IJW to make axis behavior re "zoom" on exit same as
%at beginning. Requires adding tag to the invisible axes


currax=findobj(gcf,'type','axes','-not','tag','suplabel');

if nargin < 3
 supAxes=[.08 .08 .84 .84];
 ah=findall(gcf,'type','axes');
 if ~isempty(ah)
 supAxes=[inf,inf,0,0];
 leftMin=inf;  bottomMin=inf;  leftMax=0;  bottomMax=0;
 axBuf=.04;
 set(ah,'units','normalized')
 ah=findall(gcf,'type','axes');
 for ii=1:length(ah)
  if strcmp(get(ah(ii),'Visible'),'on')
   thisPos=get(ah(ii),'Position');
   leftMin=min(leftMin,thisPos(1));
   bottomMin=min(bottomMin,thisPos(2));
   leftMax=max(leftMax,thisPos(1)+thisPos(3));
   bottomMax=max(bottomMax,thisPos(2)+thisPos(4));
  end
 end
 supAxes=[leftMin-axBuf,bottomMin-axBuf,leftMax-leftMin+axBuf*2,bottomMax-bottomMin+axBuf*2];
 end
end
if nargin < 2, whichLabel = 'x';  end
if nargin < 1, help(mfilename); return; end

if ~isstr(text) | ~isstr(whichLabel)
 error('text and whichLabel must be strings')
end
whichLabel=lower(whichLabel);

ax=axes('Units','Normal','Position',supAxes,'Visible','off','tag','suplabel');
if strcmp('t',whichLabel)
 set(get(ax,'Title'),'Visible','on')
 title(text);
elseif strcmp('x',whichLabel)
 set(get(ax,'XLabel'),'Visible','on')
 xlabel(text);
elseif strcmp('y',whichLabel)
 set(get(ax,'YLabel'),'Visible','on')
 ylabel(text);
elseif strcmp('yy',whichLabel)
 set(get(ax,'YLabel'),'Visible','on')
 ylabel(text);
 set(ax,'YAxisLocation','right')
end

for k=1:length(currax), axes(currax(k));end % restore all other axes

if (nargout < 2)
 return
end
if strcmp('t',whichLabel)
 h=get(ax,'Title');
 set(h,'VerticalAlignment','middle')
elseif strcmp('x',whichLabel)
 h=get(ax,'XLabel');
elseif strcmp('y',whichLabel) | strcmp('yy',whichLabel)
 h=get(ax,'YLabel');
end

%%%ah=findall(gcf,'type','axes');
%%%'sssssssss',kb