x0 = [3;4;6;3;2;1;0;2;1;2;6;5;4;3;2;1;6;5;4;2;0;3;2;4;3;3;2;0]; % Punto factible 

% Función objetivo
f0 = @(x11,x21,x31,x12,x22,x32,x13,x23,x33,x14,x24,x34,x15,x25,x35,x16,x26,x36,x17,x27,x37,L1,L2,L3,L4,L5,L6,L7) -((200*x11 ...
      + 300*x21 + 500*x31 + 200*x12 + 300*x22 + 500*x32 ... 
      + 200*x13 + 300*x23 + 500*x33 + 200*x14 + 300*x24 + 500*x34 ...
      + 200*x15 + 300*x25 + 500*x35 + 200*x16 + 300*x26 + 500*x36 ...
      + 200*x17 + 300*x27 + 500*x37) - ((3*(x11+x21+x31)+6)*L1 ...
      + (3*(x12+x22+x32)+6)*L2 + (3*(x13+x23+x33)+6)*L3 ...
      + (3*(x14+x24+x34)+6)*L4 + (3*(x15+x25+x35)+6)*L5 ...
      + (3*(x16+x26+x36)+6)*L6 + (3*(x17+x27+x37)+6)*L7) ...
      - (20*(x21*x31) + 20*(x22*x32) + 20*(x23*x33) ...
      + 20*(x24*x34) + 20*(x25*x35) + 20*(x26*x36) ...
      + 20*(x27*x37)) - 5000);


%    x11,x21,x31,L1,x12,x22,x32,L2,x13,x23,x33,L3,x14,x24,x34,L4,x15,x25,x35,L5,x16,x26,x36,L6,x17,x27,x37,L7
A = [-1, -1, -1, 1, 0,  0,  0,  0, 0,  0,  0,  0, 0,  0,  0,  0, 0,  0,  0,  0, 0,  0,  0,  0, 0,  0,  0,  0;
     0,  0,  0,  0, -1, -1, -1, 1, 0,  0,  0,  0, 0,  0,  0,  0, 0,  0,  0,  0, 0,  0,  0,  0, 0,  0,  0,  0;
     0,  0,  0,  0, 0,  0,  0,  0, -1, -1, -1, 1, 0,  0,  0,  0, 0,  0,  0,  0, 0,  0,  0,  0, 0,  0,  0,  0;
     0,  0,  0,  0, 0,  0,  0,  0, 0,  0,  0,  0, -1, -1, -1, 1, 0,  0,  0,  0, 0,  0,  0,  0, 0,  0,  0,  0;
     0,  0,  0,  0, 0,  0,  0,  0, 0,  0,  0,  0, 0,  0,  0,  0, -1, -1, -1, 1, 0,  0,  0,  0, 0,  0,  0,  0;
     0,  0,  0,  0, 0,  0,  0,  0, 0,  0,  0,  0, 0,  0,  0,  0, 0,  0,  0,  0, -1, -1, -1, 1, 0,  0,  0,  0;
     0,  0,  0,  0, 0,  0,  0,  0, 0,  0,  0,  0, 0,  0,  0,  0, 0,  0,  0,  0, 0,  0,  0,  0, -1, -1, -1, 1;
     1,  0,  0,  0, 0,  0,  0,  0, 0,  0,  0,  0, 0,  0,  0,  0, 0,  0,  0,  0, 0,  0,  0,  0,  0,  0,  0, 0;
     0,  1,  0,  0, 0,  0,  0,  0, 0,  0,  0,  0, 0,  0,  0,  0, 0,  0,  0,  0, 0,  0,  0,  0,  0,  0,  0, 0;
     0,  0,  1,  0, 0,  0,  0,  0, 0,  0,  0,  0, 0,  0,  0,  0, 0,  0,  0,  0, 0,  0,  0,  0,  0,  0,  0, 0;
     0,  0,  0,  0, 1,  0,  0,  0, 0,  0,  0,  0, 0,  0,  0,  0, 0,  0,  0,  0, 0,  0,  0,  0,  0,  0,  0, 0;
     0,  0,  0,  0, 0,  1,  0,  0, 0,  0,  0,  0, 0,  0,  0,  0, 0,  0,  0,  0, 0,  0,  0,  0,  0,  0,  0, 0;
     0,  0,  0,  0, 0,  0,  1,  0, 0,  0,  0,  0, 0,  0,  0,  0, 0,  0,  0,  0, 0,  0,  0,  0,  0,  0,  0, 0;
     0,  0,  0,  0, 0,  0,  0,  0, 1,  0,  0,  0, 0,  0,  0,  0, 0,  0,  0,  0, 0,  0,  0,  0,  0,  0,  0, 0;
     0,  0,  0,  0, 0,  0,  0,  0, 0,  1,  0,  0, 0,  0,  0,  0, 0,  0,  0,  0, 0,  0,  0,  0,  0,  0,  0, 0;
     0,  0,  0,  0, 0,  0,  0,  0, 0,  0,  1,  0, 0,  0,  0,  0, 0,  0,  0,  0, 0,  0,  0,  0,  0,  0,  0, 0;
     0,  0,  0,  0, 0,  0,  0,  0, 0,  0,  0,  0, 1,  0,  0,  0, 0,  0,  0,  0, 0,  0,  0,  0,  0,  0,  0, 0;
     0,  0,  0,  0, 0,  0,  0,  0, 0,  0,  0,  0, 0,  1,  0,  0, 0,  0,  0,  0, 0,  0,  0,  0,  0,  0,  0, 0;
     0,  0,  0,  0, 0,  0,  0,  0, 0,  0,  0,  0, 0,  0,  1,  0, 0,  0,  0,  0, 0,  0,  0,  0,  0,  0,  0, 0;
     0,  0,  0,  0, 0,  0,  0,  0, 0,  0,  0,  0, 0,  0,  0,  0, 1,  0,  0,  0, 0,  0,  0,  0,  0,  0,  0, 0;
     0,  0,  0,  0, 0,  0,  0,  0, 0,  0,  0,  0, 0,  0,  0,  0, 0,  1,  0,  0, 0,  0,  0,  0,  0,  0,  0, 0;
     0,  0,  0,  0, 0,  0,  0,  0, 0,  0,  0,  0, 0,  0,  0,  0, 0,  0,  1,  0, 0,  0,  0,  0,  0,  0,  0, 0;
     0,  0,  0,  0, 0,  0,  0,  0, 0,  0,  0,  0, 0,  0,  0,  0, 0,  0,  0,  0, 1,  0,  0,  0,  0,  0,  0, 0;
     0,  0,  0,  0, 0,  0,  0,  0, 0,  0,  0,  0, 0,  0,  0,  0, 0,  0,  0,  0, 0,  1,  0,  0,  0,  0,  0, 0;
     0,  0,  0,  0, 0,  0,  0,  0, 0,  0,  0,  0, 0,  0,  0,  0, 0,  0,  0,  0, 0,  0,  1,  0,  0,  0,  0, 0;
     0,  0,  0,  0, 0,  0,  0,  0, 0,  0,  0,  0, 0,  0,  0,  0, 0,  0,  0,  0, 0,  0,  0,  0,  1,  0,  0, 0;
     0,  0,  0,  0, 0,  0,  0,  0, 0,  0,  0,  0, 0,  0,  0,  0, 0,  0,  0,  0, 0,  0,  0,  0,  0,  1,  0, 0;
     0,  0,  0,  0, 0,  0,  0,  0, 0,  0,  0,  0, 0,  0,  0,  0, 0,  0,  0,  0, 0,  0,  0,  0,  0,  0,  1, 0;
     -1, 0,  0,  0, 0,  0,  0,  0, 0,  0,  0,  0, 0,  0,  0,  0, 0,  0,  0,  0, 0,  0,  0,  0,  0,  0,  0, 0;
     0, -1,  0,  0, 0,  0,  0,  0, 0,  0,  0,  0, 0,  0,  0,  0, 0,  0,  0,  0, 0,  0,  0,  0,  0,  0,  0, 0;
     0,  0, -1,  0, 0,  0,  0,  0, 0,  0,  0,  0, 0,  0,  0,  0, 0,  0,  0,  0, 0,  0,  0,  0,  0,  0,  0, 0;
     0,  0,  0, -1, 0,  0,  0,  0, 0,  0,  0,  0, 0,  0,  0,  0, 0,  0,  0,  0, 0,  0,  0,  0,  0,  0,  0, 0;
     0,  0,  0,  0,-1,  0,  0,  0, 0,  0,  0,  0, 0,  0,  0,  0, 0,  0,  0,  0, 0,  0,  0,  0,  0,  0,  0, 0;
     0,  0,  0,  0, 0, -1,  0,  0, 0,  0,  0,  0, 0,  0,  0,  0, 0,  0,  0,  0, 0,  0,  0,  0,  0,  0,  0, 0;
     0,  0,  0,  0, 0,  0, -1,  0, 0,  0,  0,  0, 0,  0,  0,  0, 0,  0,  0,  0, 0,  0,  0,  0,  0,  0,  0, 0;
     0,  0,  0,  0, 0,  0,  0, -1, 0,  0,  0,  0, 0,  0,  0,  0, 0,  0,  0,  0, 0,  0,  0,  0,  0,  0,  0, 0;
     0,  0,  0,  0, 0,  0,  0,  0,-1,  0,  0,  0, 0,  0,  0,  0, 0,  0,  0,  0, 0,  0,  0,  0,  0,  0,  0, 0;
     0,  0,  0,  0, 0,  0,  0,  0, 0, -1,  0,  0, 0,  0,  0,  0, 0,  0,  0,  0, 0,  0,  0,  0,  0,  0,  0, 0;
     0,  0,  0,  0, 0,  0,  0,  0, 0,  0, -1,  0, 0,  0,  0,  0, 0,  0,  0,  0, 0,  0,  0,  0,  0,  0,  0, 0;
     0,  0,  0,  0, 0,  0,  0,  0, 0,  0,  0, -1, 0,  0,  0,  0, 0,  0,  0,  0, 0,  0,  0,  0,  0,  0,  0, 0;
     0,  0,  0,  0, 0,  0,  0,  0, 0,  0,  0,  0,-1,  0,  0,  0, 0,  0,  0,  0, 0,  0,  0,  0,  0,  0,  0, 0;
     0,  0,  0,  0, 0,  0,  0,  0, 0,  0,  0,  0, 0, -1,  0,  0, 0,  0,  0,  0, 0,  0,  0,  0,  0,  0,  0, 0;
     0,  0,  0,  0, 0,  0,  0,  0, 0,  0,  0,  0, 0,  0, -1,  0, 0,  0,  0,  0, 0,  0,  0,  0,  0,  0,  0, 0;
     0,  0,  0,  0, 0,  0,  0,  0, 0,  0,  0,  0, 0,  0,  0, -1, 0,  0,  0,  0, 0,  0,  0,  0,  0,  0,  0, 0;
     0,  0,  0,  0, 0,  0,  0,  0, 0,  0,  0,  0, 0,  0,  0,  0,-1,  0,  0,  0, 0,  0,  0,  0,  0,  0,  0, 0;
     0,  0,  0,  0, 0,  0,  0,  0, 0,  0,  0,  0, 0,  0,  0,  0, 0, -1,  0,  0, 0,  0,  0,  0,  0,  0,  0, 0;
     0,  0,  0,  0, 0,  0,  0,  0, 0,  0,  0,  0, 0,  0,  0,  0, 0,  0, -1,  0, 0,  0,  0,  0,  0,  0,  0, 0;
     0,  0,  0,  0, 0,  0,  0,  0, 0,  0,  0,  0, 0,  0,  0,  0, 0,  0,  0, -1, 0,  0,  0,  0,  0,  0,  0, 0;
     0,  0,  0,  0, 0,  0,  0,  0, 0,  0,  0,  0, 0,  0,  0,  0, 0,  0,  0,  0,-1,  0,  0,  0,  0,  0,  0, 0;
     0,  0,  0,  0, 0,  0,  0,  0, 0,  0,  0,  0, 0,  0,  0,  0, 0,  0,  0,  0, 0, -1,  0,  0,  0,  0,  0, 0;
     0,  0,  0,  0, 0,  0,  0,  0, 0,  0,  0,  0, 0,  0,  0,  0, 0,  0,  0,  0, 0,  0, -1,  0,  0,  0,  0, 0;
     0,  0,  0,  0, 0,  0,  0,  0, 0,  0,  0,  0, 0,  0,  0,  0, 0,  0,  0,  0, 0,  0,  0, -1,  0,  0,  0, 0;
     0,  0,  0,  0, 0,  0,  0,  0, 0,  0,  0,  0, 0,  0,  0,  0, 0,  0,  0,  0, 0,  0,  0,  0, -1,  0,  0, 0;
     0,  0,  0,  0, 0,  0,  0,  0, 0,  0,  0,  0, 0,  0,  0,  0, 0,  0,  0,  0, 0,  0,  0,  0,  0, -1,  0, 0;
     0,  0,  0,  0, 0,  0,  0,  0, 0,  0,  0,  0, 0,  0,  0,  0, 0,  0,  0,  0, 0,  0,  0,  0,  0,  0, -1, 0;
     0,  0,  0,  0, 0,  0,  0,  0, 0,  0,  0,  0, 0,  0,  0,  0, 0,  0,  0,  0, 0,  0,  0,  0,  0,  0,  0,-1];
%    x11,x21,x31,L1,x12,x22,x32,L2,x13,x23,x33,L3,x14,x24,x34,L4,x15,x25,x35,L5,x16,x26,x36,L6,x17,x27,x37,L7

b = [-5;-5;-5;-5;-5;-5;-5;6;6;6;6;6;6;6;6;6;6;6;6;6;6;6;6;6;6;6;6;6;
      0; 0; 0; 0; 0; 0; 0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0];

syms x11 x21 x31 L1 x12 x22 x32 L2 x13 x23 x33 L3 x14 x24 x34 L4 x15 x25 x35 L5 x16 x26 x36 L6 x17 x27 x37 L7
% Calculamos el gradiente
grad = gradient(f0, [x11,x21,x31,L1,x12,x22,x32,L2,x13,x23,x33,L3,x14,x24,x34,L4,x15,x25,x35,L5,x16,x26,x36,L6,x17,x27,x37,L7]);
fgrad = inline(grad);


% Utilizamos gradiente proyectado y redondeamos
xmin = GradienteProyectado(f0,fgrad,A, 0.0001,x0, b);
xminr = round(xmin)

% Graficamos
x = [xminr(1) + xminr(2) + xminr(3), xminr(5) + xminr(6) + xminr(7), xminr(9) + xminr(10) + xminr(11), xminr(13) + xminr(14) + xminr(15), xminr(17) + xminr(18) + xminr(19), xminr(21) + xminr(22) + xminr(23), xminr(25) + xminr(26) + xminr(27)];
y = x+5;
yl = 2*x;
plot(x,y,'m',x,yl,'g')
legend('Guardias','Ladrones')
