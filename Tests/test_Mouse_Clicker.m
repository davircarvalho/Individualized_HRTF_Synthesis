clear all; clc;
Obj = SOFAload('individuo_140.sofa');
posi = Obj.SourcePosition;

%% create test figure
f = figure(1);
scatter(posi(:,1), posi(:,2), 26,'square', 'filled')
axis tight

%% set function to call on mouse click
set(f,'WindowButtonMotionFcn',@mouse_mov);
set(f,'WindowButtonDownFcn',@mouse_down);
set(f,'WindowButtonUpFcn',@mouse_up);

% function called on mouse click in the figure
function mouse_down(~,~)
    global state 
    state = true;
    get(gca, 'CurrentPoint')   
end

function mouse_up(~,~)
    global state 
    state = false;
end

function mouse_mov(~,~)
global state
    while state            
        get(gca, 'CurrentPoint')  
        pause(0)
     end
end