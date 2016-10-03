x= [1:4];
y = [1:4];
[X,Y] = meshgrid(x,y)
%Note that size(Z) is the same as size(x) and size(y)
%Z = [1 0 -2 1 ;1 0 -2 1;1 0 -2 1 ;1 0 -2 1];
Z=1*ones(4);
 Z=[ [0.1 0.2 0.3; 0.4 0.5 0.6;0.7 0.8 0.9] zeros(3,1);
        zeros(1,3) 0];
Z_gt=[eye(3) zeros(3,1);
       zeros(1,3) 0]; % create a colormap having RGB values of dark green,
%light green, white, dark red and light red.
%map2 = [0 1 0;1 0 0 ];
% %use the user defined colormap for figure.

%plot the figure
figure(1)
colormap('parula(10)');
set(gca,'clim',[0 1]);
pcolor(X,Y,Z);
axis ij
figure(2)
colormap('parula(10)');
set(gca,'clim',[0 1]);
pcolor(X,Y,Z_gt);
axis ij
%set the x and y labels
%set(gca,'XTick',[1 2 3 4 5 6],'YTick',[1 2 3 4],'XTicklabel',[' ';'a';'b'; 'c'; 'd';'e'],'YTicklabel',[' ';'f';'g';'h']);
%set the color limits
%caxis([-2 2])