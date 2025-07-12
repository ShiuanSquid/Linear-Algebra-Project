%% Cubic Spline Interpolation for CO2 Data
% Project: Linear Algebra Project
% Author: SOH SHIUAN MING MAT2409043
% Date: 2025-07-05
% Description: Implements natural cubic and parabolic runout spline interpolation for atmospheric CO2 concentration data from 1850-2000

clc; clear; close all;

%Method Name.
method1 = 'Natural Cubic Spline';
method2 = 'Parabolic Runout Spline';

%Matrix for Natural Cubic Spline.
A1 =[4 1 0 0 0;
    1 4 1 0 0;
    0 1 4 1 0;
    0 0 1 4 1;
    0 0 0 1 4];

%Matrix for Parabolic Runout Spline.
A2 =[5 1 0 0 0;
    1 4 1 0 0;
    0 1 4 1 0;
    0 0 1 4 1;
    0 0 0 1 5];

%Graph title name.
titlename1 = 'Estimated $CO_2$ Concentration (ppm) using Natural Cubic Spline Interpolation';
titlename2 = 'Estimated $CO_2$ Concentration (ppm) using Parabolic Runout Spline Interpolation';
titlename3 = 'Comparison of Natural and Parabolic Runout Spline Estimates';

%Color for Graph.
GraphColor1 = turbo(6);
GraphColor2 = turbo(6);
GraphColor3 = [ 0, 1, 0;
    0, 1, 0;
    0, 1, 0;
    0, 1, 0;
    0, 1, 0;
    0, 1, 0];
GraphColor4 = [ 1, 0, 0;
    1, 0, 0;
    1, 0, 0;
    1, 0, 0;
    1, 0, 0;
    1, 0, 0];

%Draw Graph 1, 2 and 3.
DrawSpline(A1, 1, GraphColor1, titlename1, method1)
DrawSpline(A2, 2, GraphColor2, titlename2, method2)
DrawSpline(A1, 3, GraphColor3, titlename3, method1)
DrawSpline(A2, 3, GraphColor4, titlename3, method2)

function DrawSpline(A, figureNumber, GraphColor, titlename, method)
%Given data.
Actual1990 = 354.29;
Actual2010 = 389.21;

Year = 1850:25:2000;
CO2 = [285.2, 288.6, 295.7, 305.3,311.3, 331.36, 369.64];

x = Year;
y = CO2;
%Compute h.
h = x(2) - x(1);

%Create data base for coefficients and M.
M = [0 0 0 0 0 0 0];
a = [0 0 0 0 0 0];
b = [0 0 0 0 0 0];
c = [0 0 0 0 0 0];
d = [0 0 0 0 0 0];

%Data base for point of curve.
xx = zeros(6, 1000);
g = zeros(6, 1000);

B = (6./(h.^2)).*[y(1) - 2.*y(2) + y(3);
    y(2) - 2.*y(3) + y(4);
    y(3) - 2.*y(4) + y(5);
    y(4) - 2.*y(5) + y(6);
    y(5) - 2.*y(6) + y(7)];

%Find M2 to M6.
X = A\B;

%Saves data to M.
for i = 2:6
    M(i) = X(i-1);
end

%Endpoint Constraints.
if A(1,1) == 5 
    M(1) = M(2);
    M(7) = M(6);
else 
    M(1) = 0;
    M(7) = 0;
end

%Compute and save coefficient.
for i = 1:6
    a(i) = (M(i+1) - M(i))./(6.*h);
    b(i) = M(i)./2;
    c(i) = (y(i+1)-y(i))./h - ((M(i+1) + 2.*M(i)).*h )./6;
    d(i) = y(i);
end

for i = 1:6
    xx(i, :) = linspace(x(i),x(i+1),1000);
    g(i, : ) = a(i).* (xx(i, :) - x(i)).^3 +b(i).*(xx(i, :) - x(i)).^2 +c(i).*(xx(i, :) - x(i)) + d(i);
end

figure(figureNumber)
hold on
%Plot the curve.
for i = 1:6
    plot(xx(i, :),g(i, : ),'-','LineWidth',0.2, 'Color', GraphColor(i,:))
end

%Plot the given data point.
plot(x,y,'o','MarkerFaceColor','white', ...
    'MarkerEdgeColor','black', ...
    'MarkerSize',5, ...
    'LineWidth',1)

%Compute the cubic spline estimate for the CO2 concentration in 1990 and 2010.
Spline1990 = a(6).* (1990 - x(6)).^3 +b(6).*(1990 - x(6)).^2 +c(6).*(1990 - x(6)) + d(6);
Spline2010 = a(6).* (2010 - x(6)).^3 +b(6).*(2010 - x(6)).^2 +c(6).*(2010 - x(6)) + d(6);

%Absolute error and relative error.
absErr1990 = abs(Spline1990 - Actual1990);
absErr2010 = abs(Spline2010 - Actual2010);

relErr1990 = absErr1990 / Actual1990 * 100;
relErr2010 = absErr2010 / Actual2010 * 100;

%Print all the result.
Year = [1990; 2010];
Estimate = [Spline1990; Spline2010];
Actual = [Actual1990; Actual2010];
AbsErr = [absErr1990; absErr2010];
RelErr = [relErr1990; relErr2010];

ResultTable = table(Year, Actual, Estimate, AbsErr, RelErr);
ResultTable.Properties.VariableNames = {'Year','Actual_CO2', 'Estimated_CO2',  'Absolute_Error', 'Relative_Error_Percent'};

if figureNumber == 1 || figureNumber == 2
fprintf('Using %s:\n', method);
disp(ResultTable);
fprintf('In 1990, the spline estimate is %.4f ppm, while the actual CO2 concentration is %.2f ppm.\n', ...
    Spline1990, Actual1990);
fprintf('This results in an absolute error of %.4f ppm and a relative error of %.4f%%.\n\n', ...
    absErr1990, relErr1990);
fprintf('In 2010, the spline estimate is %.4f ppm, while the actual CO2 concentration is %.2f ppm.\n', ...
    Spline2010, Actual2010);
fprintf('This results in an absolute error of %.4f ppm and a relative error of %.4f%%.\n\n', ...
    absErr2010, relErr2010);
end

if figureNumber == 1 || figureNumber == 2
% Display B
Bi_Index = (1:5)';    
Bi_Values = B(:);      
B_Table = table(Bi_Index, Bi_Values);  
B_Table.Properties.VariableNames = {'i', 'B_i'}; 
disp('Values B_i:');
disp(B_Table);

% Display M1 to M7.
Mi_Index = (1:7)';
Mi_Values = M(:);
M_Table = table(Mi_Index, Mi_Values);
M_Table.Properties.VariableNames = {'i', 'M_i'};
disp(['Values M_i for ', method, ':']);
disp(M_Table);

% Display spline coefficients a_i, b_i, c_i, d_i for i = 1 to 6.
Interval = (1:6)';
a_col = a(:);
b_col = b(:);
c_col = c(:);
d_col = d(:);
Coeff_Table = table(Interval, a_col, b_col, c_col, d_col);
Coeff_Table.Properties.VariableNames = {'i', 'a_i', 'b_i', 'c_i', 'd_i'};
disp(['Spline coefficients for ', method, ':']);
disp(Coeff_Table);
end

%Legend for compare.
if figureNumber == 3 || figureNumber == 4
    LegendColors = [0, 1, 0; 1, 0, 0];  
    LegendNames = {'Natural Spline', 'Parabolic Runout Spline'};
hold on;
for i = 1:2
    k(i) = plot(nan, nan, '-', 'Color', LegendColors(i,:), ...
        'LineWidth', 0.2, 'DisplayName', LegendNames{i});
end
legend([k(1), k(2)], LegendNames,'Location','best','FontName','Times New Roman', 'FontSize',10);
end

h = gca;
h.XTick = 1850:25:2000;
xtickangle(0)

title(titlename, ...
      'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', 10);
xlabel('Year', 'Interpreter','latex')
ylabel('$CO_2$ (ppm)', 'Interpreter','latex')
axis tight 
box off
grid on
hold off

end

% Save figure 1, 2, and 3 as PNG.
% f1 = figure(1);
% exportgraphics(f1, 'NaturalCubicSpline.png');
% 
% f2 = figure(2);
% exportgraphics(f2, 'ParabolicRunoutSpline.png');
% 
% f3 = figure(3);
% exportgraphics(f3, 'Comparison.png');
