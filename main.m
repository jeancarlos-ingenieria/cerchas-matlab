clc
clf
clear

% Introduccion de datos

% Seccion y modulo elastico

A = 1;
E = 1;

% Escalas

grafDeformada = 1;
escalaDef = 0.005;

% Nudos

nudos = [1 0 0;
    2 10 0;
    3 10 10;
    4 0 10];

% Elementos

elementos = [1 1 2 A E;
        2 1 3 A E;
        3 1 4 A E;
        4 3 4 A E;
        5 2 3 A E];

% Apoyos

apoyos = [2 1 0;
    3 1 1];

% Cargas

cargas = [4 2 -4];

% ===================================
% Inicio del programa
% ===================================

temp = size(elementos);
nroElem = temp(1);

% Longitudes y angulos

for i = 1:nroElem

    nudo1 = elementos(i, 2);
    nudo2 = elementos(i, 3);

    x1 = nudos(nudo1, 2);
    x2 = nudos(nudo2, 2);
    y1 = nudos(nudo1, 3);
    y2 = nudos(nudo2, 3);

    longitud = ((x2 - x1)^2 + (y2 - y1)^2)^0.5;
    longitudes(i) = longitud;

    angulo = acosd((x2 - x1) / longitud);

    if (y2 - y1) < 0
        angulo = -angulo;
    endif

    angulos(i) = angulo;
endfor

% Matriz de rigidez de elementos

for i = 1:nroElem

    alfa = angulos(i);
    cs = cosd(alfa);
    sn = sind(alfa);
    longitud = longitudes(i);

    kLocal(:, :, i) = E * A / longitud * [1 -1; -1 1];
    T = [cs sn 0 0; 0 0 cs sn];
    kGlobal(:, :, i) = T' * kLocal(:, :, i) * T;

endfor

% Grados de libertad

for i = 1:nroElem
    gdl(i, 1) = i;
    nudo1 = elementos(i, 2);
    nudo2 = elementos(i, 3);
    gdl(i, 2) = nudo1 * 2 - 1;
    gdl(i, 3) = nudo1 * 2;
    gdl(i, 4) = nudo2 * 2 - 1;
    gdl(i, 5) = nudo2 * 2;

endfor

% Matriz de rigidez de la estructura

temp = size(nudos);
nroNudos = temp(1);

kEstructura = zeros(nroNudos * 2, nroNudos * 2);

for i = 1:nroElem
    gdlBarra = gdl(i, 2:5);
    kEstructura(gdlBarra, gdlBarra) = kGlobal(:, :, i) + kEstructura(gdlBarra, gdlBarra);
endfor

% Apoyos de estructura

temp = size(apoyos);
nroApoyos = temp(1);

contador = 0;

for i = 1:nroApoyos

    if (apoyos(i, 2) == 1)
        contador = contador +1;
        nudoApoyado = apoyos(i, 1);
        gdlRestringidos(contador) = nudoApoyado * 2 - 1;
    endif

    if (apoyos(i, 3) == 1)
        contador = contador +1;
        nudoApoyado = apoyos(i, 1);
        gdlRestringidos(contador) = nudoApoyado * 2;
    endif

endfor

% Matriz de rigidez reducida

kReducida = kEstructura;

kReducida(gdlRestringidos, :) = [];
kReducida(:, gdlRestringidos) = [];

% Vector de cargas

F = zeros(nroNudos * 2, 1);

temp = size(cargas);
nroCargas = temp(1);

for i = 1:nroCargas
    nudoCargado = cargas(i, 1);
    Fx = cargas(i, 2);
    Fy = cargas(i, 3);

    F(nudoCargado * 2 - 1, 1) = Fx;
    F(nudoCargado * 2, 1) = Fy;
endfor

% Vector de cargas reducido

FReducido = F;
FReducido(gdlRestringidos, :) = [];

%% Resolucion del sistema de ecuaciones

UReducido = linsolve(kReducida, FReducido);

% GDL libres
gdlLibres = 1:1:nroNudos * 2;
gdlLibres(gdlRestringidos) = [];

% Completando el vector de desplazamiento

U = zeros(nroNudos * 2, 1);
U(gdlLibres) = UReducido;

% Calculando reacciones de la estructura

% [k][U]=[F]+[R]
% [R]=[k][U]-[F]

R = kEstructura * U - F;

% Ploteando

hold on

% Estructura sin deformar

for i = 1:nroElem
    nudo1 = elementos(i, 2);
    nudo2 = elementos(i, 3);

    x1 = nudos(nudo1, 2);
    x2 = nudos(nudo2, 2);
    y1 = nudos(nudo1, 3);
    y2 = nudos(nudo2, 3);

    plot([x1 x2], [y1 y2], "linewidth", 1, "color", "k")
endfor



for i = 1:nroNudos
    nudosDesplazados(i, 1) = nudos(i, 1);
    nudosDesplazados(i, 2) = nudos(i, 2) + U(i * 2 - 1,1) * escalaDef;
    nudosDesplazados(i, 3) = nudos(i, 3) + U(i * 2,1) * escalaDef;
endfor

% Estructura deformada


for i = 1:nroElem
    nudo1 = elementos(i, 2);
    nudo2 = elementos(i, 3);

    x1 = nudosDesplazados(nudo1, 2);
    x2 = nudosDesplazados(nudo2, 2);
    y1 = nudosDesplazados(nudo1, 3);
    y2 = nudosDesplazados(nudo2, 3);

    plot([x1 x2], [y1 y2], "linewidth", 1, "color", "r")
endfor

axis equal

