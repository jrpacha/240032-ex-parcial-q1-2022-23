% Problema 2 #5
clearvars
close all
format short e

eval('meshHole');

%
% Part A
%
vertexsT1 = nodes(elem(218,:),:);
x_barycenter = sum(vertexsT1(:,1))/3;

fprintf(['Part A\n',...
    'The x-coordinate of the barycenter of the triangle T1 is ',...
    'x = %.6f\n'],x_barycenter)
%
% Part B
%
vertexsT2 = nodes(elem(214,:),:);
alphas = [0.2,0.5,0.3];
y = alphas*vertexsT2(:,2);
fprintf(['Part B\n',...
    'The y-coordinate of the point inside the triangle T2 with ',...
    'barycentric coordinates\n\t alpha = [0.2,0.5,0.3] is y = %.6f\n'],y)
%
% Part C
%
fprintf('Part C\n')
nods = [elem(218,:);elem(214,:)];% Check which are the global nodes of
                                 % T1 and T2 (see result at command window)
fprintf('\tNodes T1:%5d%5d%5d\n',nods(1,:));
fprintf('\tNodes T2:%5d%5d%5d\n',nods(2,:));
vertexsQuad = [nodes(58,:); ...  % 1st node -> global node 58
    nodes(67,:); ....            % 2nd node -> global node 67
    nodes(84,:); ....            % 3rd node -> global node 84
    nodes(75,:)];                % 4th node -> global node 75
vertexsQuadPlot = [vertexsQuad; vertexsQuad(1,:)];
plot(vertexsQuad(:,1),vertexsQuad(:,2),...
    'Marker','o',...
    'MarkerFaceColor','green',...
    'MarkerEdgeColor','black',...
    'MarkerSize',10)
hold on
plot(vertexsQuadPlot(:,1),vertexsQuadPlot(:,2),...
    'Color','blue',...
    'LineWidth',2)
hold off

q = [0.625,0.555];
[alphas,isInside] = baryCoordQuad(vertexsQuad,q);
fprintf(['The third barycentric coordinate of Q wrt the quadrilateral ',...
    'is alpha(3) = %.6f\n'],alphas(3))
%
% Part D 
%
f = @(x,y) x + cos(y);
xi = vertexsQuad(:,1); yi = vertexsQuad(:,2);
fi = f(xi,yi);
f_interp = alphas*fi;
fprintf(['Part D\n',...
    'The interpolated value of f(Q) is f_interp(Q) = %.6f\n'],f_interp)

