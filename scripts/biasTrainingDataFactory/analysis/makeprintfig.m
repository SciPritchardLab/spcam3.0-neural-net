function [hf] = makeprintfig (w,h)
% Creates a new figure with visual size in the same proportions as its
% print will be, and with the 

    if (h/w > 800/1340) 
        fheight = 900;
        fwidth = w/h*fheight;
    else
        fwidth = 1400;
        fheight = h/w*fwidth;
    end
%    set (gcf, 'PaperSize',[w h],'PaperPosition',[0 0 w h],'PaperUnits','inches');
    hf = figure,clf;  
    aux = get (hf,'Position'); screenaspect = aux(3)/aux(4); 
set (gcf, 'OuterPosition',[aux(1) aux(2) fwidth*screenaspect fheight]);     
    set (hf, 'Papersize',[w h],'PaperPosition',[0 0 w h]);
