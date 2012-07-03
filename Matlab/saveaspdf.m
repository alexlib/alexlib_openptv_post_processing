function saveaspdf(h,filename)
%SAVEASPDF Save a figure as a clean pdf file ready for publication.
%   saveaspdf saves the current figure to the path defined in the global
%   SAVEASPDF_PATH, or to your Desktop if no path is specified.
%
%   saveaspdf(filename) saves the current figure as filename to the path
%   defined in the global SAVEASPDF_PATH, or to your Desktop if no path is
%   specified.
%
%   saveaspdf(h,filename) saves the figure with handle h as filename to the
%   path defined in the global SAVEASPDF_PATH, or to your Desktop if no
%   path is specified.
%
%   Before printing the pdf, lines are scaled using the global scale factor
%   SAVEASPDF_SCALELINEWIDTH (default 1.4), and font sizes are scaled using
%   the global scale factor SAVEASPDF_SCALEFONTSIZE (default 1.3).
%
%   See also saveas, print.

%   Version: 23/09/10
%   Authors: Laurent Sorber (Laurent.Sorber@cs.kuleuven.be)

% Define globals.
global SAVEASPDF_PATH;
global SAVEASPDF_SCALELINEWIDTH;
global SAVEASPDF_SCALEFONTSIZE;

% Set up default values.
if nargin == 1 && ischar(h)
    filename = h;
    h = gcf;
elseif nargin == 0
    h = gcf;
    filename = get(h,'Name');
    if isempty(filename)
        filename = ['Figure ' int2str(gcf)];
    end
end
if ~isempty(SAVEASPDF_PATH)
    path = SAVEASPDF_PATH;
else
    if ispc, path = 'USERPROFILE'; else path = 'HOME'; end
end
path = fullfile(getenv(path),'Desktop');
filename = fullfile(path,regexprep(filename,'\.[a-zA-Z]+','.pdf'));
if ~isempty(SAVEASPDF_SCALELINEWIDTH)
    scalelw = SAVEASPDF_SCALELINEWIDTH;
else
    scalelw = 1.4;
end
if ~isempty(SAVEASPDF_SCALEFONTSIZE)
    scalefs = SAVEASPDF_SCALEFONTSIZE;
else
    scalefs = 1.3;
end

% Backup current settings.
oldUnits = get(h,'Units');
oldPaperPosition = get(h,'PaperPosition');
oldPaperSize = get(h,'PaperSize');
lines = findall(h,'type','line');
text = [findall(h,'type','text'); findall(h,'type','axes')];
oldLineWidth = get(lines,'LineWidth');
oldFontSize = get(text,'FontSize');
if ~iscell(oldLineWidth), oldLineWidth = num2cell(oldLineWidth); end
if ~iscell(oldFontSize), oldFontSize = num2cell(oldFontSize); end

% Scale the line width and font size, reset units and update the positions.
set(lines,{'LineWidth'},num2cell(cellfun(@(x)x*scalelw,oldLineWidth)));
set(text,{'FontSize'},num2cell(cellfun(@(x)round(x*scalefs),oldFontSize)));
set(h,'Units','inches');
curPosition = get(h,'Position');
set(h,'PaperPosition',[0,0,curPosition(3:4)]);
set(h,'PaperSize',curPosition(3:4));

% Print the pdf.
print(h,'-dpdf',filename);

% Restore the previous settings.
set(h,'Units',oldUnits);
set(h,'PaperPosition',oldPaperPosition);
set(h,'PaperSize',oldPaperSize);
set(text,{'FontSize'},oldFontSize);
set(lines,{'LineWidth'},oldLineWidth);
