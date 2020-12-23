function photo_helper(image_path)
% FERRAMENTA PARA MEDIÇÂO DE DISTANCIA EM IMAGENS
% REF: https://www.mathworks.com/help/images/measure-distances-in-images.html
% Dezembro 2020 - @UFSM

%%% TUTORIAL ------------------------------------------------------------->
% ~ Ao abrir a imagem, clique arraste e solte para gerar um segmento de reta,
%   o segmento mostrará a distancia em pixels.
%
% ~ Crie um segmento sobre um objeto na imagem com tamanho real conhecido.
%
% ~ Clique duplo sobre o segmento, na caixa de dialogo atualize o valor em
%   pixels para o correspondente comprimento em milimetros.
%
% ~ Os demais segmentos serão medidos com o valor corrigido!
% 
% Enjoy

%% Carregar imagem 
im = imread(image_path);

%% Propriedades da imagem 
sz = size(im);
myData.Units = 'pixels';
myData.MaxValue = hypot(sz(1),sz(2));
myData.Colormap = hot;
myData.ScaleFactor = 1;

%% Visualizar e extrair medidas (manualmente)
figure()
hIm = imshow(im);
    
%%% Inicializar função de desenho sobre a imagem 
% warning('off')
hIm.ButtonDownFcn = @(~,~) startDrawing(hIm.Parent,myData);


end






%% %%%%%%%%%%%%%%%%%%%%%%% FUNÇÕES INTERNAS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function startDrawing(hAx,myData)
% Create a line ROI object. Specify the initial color of the line and
% store the |myData| structure in the |UserData| property of the ROI.
h = images.roi.Line('Color',[0, 0, 0.5625],'UserData',myData);

% Set up a listener for movement of the line ROI. When the line ROI moves,
% the |updateLabel| callback updates the text in the line ROI label and
% changes the color of the line, based on its length.
addlistener(h,'MovingROI',@updateLabel);

% Set up a listener for clicks on the line ROI. When you click on the line
% ROI, the |updateUnits| callback opens a GUI that lets you specify the
% known distance in real-world units, such as, meters or feet.
addlistener(h,'ROIClicked',@updateUnits);

% Get the current mouse location from the |CurrentPoint| property of the
% axes and extract the _x_ and _y_ coordinates.
cp = hAx.CurrentPoint;
cp = [cp(1,1) cp(1,2)];

% Begin drawing the ROI from the current mouse location. Using the
% |beginDrawingFromPoint| method, you can draw multiple ROIs.
h.beginDrawingFromPoint(cp);

% Add a custom option to the line ROI context menu to delete all existing
% line ROIs.
c = h.UIContextMenu;
uimenu(c,'Label','Delete All','Callback',@deleteAll);
end


%%%% Update Labels %%%%
function updateLabel(src,evt)

% Get the current line position.
pos = evt.Source.Position;

% Determine the length of the line.
diffPos = diff(pos);
mag = hypot(diffPos(1),diffPos(2));

% Choose a color from the color map based on the length of the line. The
% line changes color as it gets longer or shorter.
color = src.UserData.Colormap(ceil(64*(mag/src.UserData.MaxValue)),:);

% Apply the scale factor to line length to calibrate the measurements.
mag = mag*src.UserData.ScaleFactor;

% Update the label.
set(src,'Label',[num2str(mag,'%30.1f') ' ' src.UserData.Units],'Color',color);
end


%%%% Atualizar unidades %%%%
function updateUnits(src,evt)
% When you double-click the ROI label, the example opens a popup dialog box
% to get information about the actual distance. Use this information to
% scale all line ROI measurements.
if strcmp(evt.SelectionType,'double') && strcmp(evt.SelectedPart,'label')

    % Display the popup dialog box.
    answer = inputdlg({'Known distance','Units'},...
        'Specify known distane',[1 20],{'10','mm'});

    % Determine the scale factor based on the inputs.
    num = str2double(answer{1});

    % Get the length of the current line ROI.
    pos = src.Position;
    diffPos = diff(pos);
    mag = hypot(diffPos(1),diffPos(2));

    % Calculate the scale factor by dividing the known length value by the
    % current length, measured in pixels.
    scale = num/mag;

    % Store the scale factor and the units information in the |myData|
    % structure.
    myData.Units = answer{2};
    myData.MaxValue = src.UserData.MaxValue;
    myData.Colormap = src.UserData.Colormap;
    myData.ScaleFactor = scale;

    % Reset the data stored in the |UserData| property of all existing line
    % ROI objects. Use |findobj| to find all line ROI objects in the axes.
    hAx = src.Parent;
    hROIs = findobj(hAx,'Type','images.roi.Line');
    set(hROIs,'UserData',myData);

    % Update the label in each line ROI object, based on the information
    % collected in the input dialog.
    for i = 1:numel(hROIs)
        pos = hROIs(i).Position;
        diffPos = diff(pos);
        mag = hypot(diffPos(1),diffPos(2));
        set(hROIs(i),'Label',[num2str(mag*scale,'%30.1f') ' ' answer{2}]);
    end

    % Reset the |ButtonDownFcn| callback function with the current |myData|
    % value.
    hIm = findobj(hAx,'Type','image');
    hIm.ButtonDownFcn = @(~,~) startDrawing(hAx,myData);
end
end


%%%% Apagar todos os segmentos %%%%
function deleteAll(src,~)
hFig = ancestor(src,'figure');
hROIs = findobj(hFig,'Type','images.roi.Line');
delete(hROIs)
end