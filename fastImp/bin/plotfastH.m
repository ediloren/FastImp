function plotfastH(name);
% Load and plot a 3D fastimp structure produced with "fasthenry/bin/zbuf -m name" 
% where name is without the ".mat"

xt = []; yt = []; zt =  []; 
xq = []; yq = []; zq = []; 

eval(['load ' name]);
fprintf(1, 'loaded %d panels\n', length(xt) + length(xq)); 
hold off; 
if length(xt) > 0, 
  ht = fill3(xt, yt, zt, 'r'); 
  hold on; 
end
X = max([max(xt) max(yt) max(zt)]); 
Y = min([min(xt) min(yt) min(zt)]); 

if length(xq) > 0, 
  hq = fill3(xq, yq, zq, 'y'); 
end; 

X = max([X max(xq) max(yq) max(zq)]); 
Y = min([Y min(xq) min(yq) min(zq)]); 

axis([Y X Y X Y X]); 
fprintf(1, 'finished filling polygons\n'); 

axis('square');
set(ht, 'FaceColor', 'w')
set(ht, 'EdgeColor', 'k')
set(hq, 'FaceColor', 'w')
set(hq, 'EdgeColor', 'k')
axis('square'); 
axis('off'); 

f = gcf; 
set(f, 'Color', [1 1 1]); 
set(f, 'PaperOrientation', 'portrait'); 
print -deps panels.eps;
