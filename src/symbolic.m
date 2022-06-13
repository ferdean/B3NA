%% % Simple case

syms E I K L w(x) k L

Dw      = diff(w, x);
DDw     = diff(w, x, 2);
DDDw    = diff(w, x, 3);

eqn     = E*I * diff(w, x, 4) == k * x;

cond_1  = w(0) == 0;
cond_2  = Dw(0) == 0;
cond_3  = DDw(L) == 0;
cond_4  = DDDw(L) == 0;

solution = dsolve(eqn, cond_1, cond_2, cond_3, cond_4);

disp(solution)

%% % Frame

syms E I M P w(x)


Dw      = diff(w, x);
DDw     = diff(w, x, 2);

eqn     = E*I*DDw == (P);

cond_1  = w(0) == 1/3;
cond_2  = Dw(0) == 0;

solution = dsolve(eqn, cond_1, cond_2);

disp(solution)

%% Basic frame
%%% Unitary characteristics
P = 1; 
E = 1; 
I = 1;

%%% Line definition
x     = 0:0.05:2;
x_2   = 0:0.05:1;

def_1 = (x.*(- P.*x.^2 + 3.*P.*x))./(6*E*I);
def_2 = -(P.*x_2.^2)./(2*E*I);
def_3 = (x_2.*(- P.*x_2.^2 + 3.*P.*x_2))./(6*E*I);

%%% Plot body
fig = figure(1);

set(fig, 'Units', 'centimeters')
pos = get(fig,'Position');
set(fig,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])

bar_1 = plot(def_1, x, 'LineWidth', 1.5, 'Color', 'k');
hold on
bar_2 = plot(x_2 + 2/3, def_2 + 2, 'LineWidth', 1.5, 'Color', 'k');
hold on
bar_3 = plot(def_3 + 4/3, x_2 + 0.5, 'LineWidth', 1.5, 'Color', 'k');

hold on
plot([0 0 1 1], [0 2 2 1], 'LineWidth', 1.25, 'LineStyle', '--', 'Color', '#828282')

p1 = [1 1];                         
p2 = [1.2 1];                         
dp = p2-p1;                         
quiver(p1(1),p1(2),dp(1),dp(2),0, 'Color', 'r', 'LineWidth', 1.2, 'MaxHeadSize', 1.5 )

p1 = [4/3 0.5];                         
p2 = [4/3 + 0.2 0.5];                         
dp = p2-p1;                         
quiver(p1(1),p1(2),dp(1),dp(2),0, 'Color', 'r', 'LineWidth', 1.2, 'MaxHeadSize', 1.5 )

yline(0, 'LineWidth', 5)

xlim([-1.5 2.5])
ylim([-0.5 3.5])

bar_2.Annotation.LegendInformation.IconDisplayStyle = 'off';
bar_3.Annotation.LegendInformation.IconDisplayStyle = 'off';

legend('deformed', 'original', 'interpreter', 'latex')

xlabel('x (-)', 'Interpreter', 'latex')
ylabel('y (-)', 'Interpreter', 'latex')

