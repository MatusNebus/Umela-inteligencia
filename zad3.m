clc; clear; close all;

%Nastavenia GA:
runs    = 10;               % pocet behov
popsize = 60;               % velkost populacie
ngen    = 1000;             % pocet generacii
eliteCount = 3;             % elitizmus (koľko najlepších ide priamo ďalej)
selectionName = "seltourn"; % metody vyberu: selsus, selbest, selrand, seltourn
pmut  = 0.5;                % pravdepodobnosť mutácie

%Súradnice bodov (matica B) zo zadania:
B = [ ...
    0, 0;     17, 100;  51, 15;   70, 62;   42, 25;   32, 17;   51, 64; ...
    39, 45;   68, 89;   20, 19;   12, 87;   80, 37;   35, 82;   2, 15; ...
    38, 95;   33, 50;   85, 52;   97, 27;   99, 10;   37, 67;   20, 82; ...
    49, 0;    62, 14;   7, 60;    0, 0  ...
];

nPoints = size(B, 1);           % cize 25 bodov
startIdx = 1;                   % prvý index
endIdx   = startIdx;            % koncime zas v bode 1
innerIdx = 2:(nPoints-1);       % posledny ignorujeme
nInner   = numel(innerIdx);     % pocet vnutornych

%Úložiská výsledkov pre graf a vyhodnotenie:
bestHistAll = nan(runs, ngen);    % (run, generácia) = najlepšia fitness v populácii
bestXAll    = nan(runs, nInner);  % uložíme len vnútornú permutáciu (2..24)
bestFAll    = nan(runs, 1);       % finálna najlepšia dĺžka pre každý beh

%HLAVNY CYKLUS:
for r = 1:runs
    pop = gen_perm_pop(popsize, innerIdx);      %Inicializacia populacie
    f = evaluate_pop(B, pop, startIdx, endIdx); %Vyhodnotenie fitness pre celu populaciu

    for g = 1:ngen
        elites = selbest(pop, f, eliteCount);   %eliteCount najlepsich ide priamo dalej
        nKids = popsize - eliteCount;           %pocet zmenenych potomkov

        %vyber:
        parents = do_selection(selectionName, pop, f, nKids);

        %krizenie:
        kids = crosord(parents, 1);

        %mutacia:
        kids = invord(kids, pmut);

        pop = [elites; kids]; % Nova populacia
        f = evaluate_pop(B, pop, startIdx, endIdx); % Vyhodnotenie novej populacie
 
        bestHistAll(r, g) = min(f); % Ulozime najlepsiu hodnotu v generacii
    end

    % Najlepsi jedinec po poslednej generacii v tomto behu
    [bestF, bestIdx] = min(f);
    bestFAll(r) = bestF;
    bestXAll(r,:) = pop(bestIdx,:);

    if bestF <= 480
        if round(bestF) == 468
            fprintf('Run %02d: bestF = %.6f <<OK>> (BEST)\n', r, bestF);
        else
            fprintf('Run %02d: bestF = %.6f <<OK>>\n', r, bestF);
        end
    else
        fprintf('Run %02d: bestF = %.6f\n', r, bestF);
    end
end

%  VYHODNOTENIE + GRAFY
[globalBestF, ridx] = min(bestFAll);
globalBestInner = bestXAll(ridx,:);
globalBestRoute = [startIdx, globalBestInner, endIdx];

nSuccess = sum(bestFAll <= 480);
successPct = 100 * nSuccess / runs;
fprintf("\nuspesnost (<=480): (%.1f %%)\n", successPct);
fprintf("\n====================\n");
fprintf("GLOBAL BEST (z %d behov):\n", runs);
fprintf("bestF = %.6f\n", globalBestF);
fprintf("bestRoute = ["); fprintf(" %d", globalBestRoute); fprintf(" ]\n");
fprintf("====================\n");

% Graf konvergencie: kazdy beh + priemer
figure; hold on; grid on;
plot(1:ngen, bestHistAll', 'LineWidth', 1);
plot(1:ngen, mean(bestHistAll, 1, 'omitnan'), 'k', 'LineWidth', 2); % priemer ciernou
xlabel("Generácia");
ylabel("Najlepšia dĺžka trasy v populácii");
title("Genetický algoritmus");
legendStrings = [compose("Run %d", 1:runs), "Priemer"];
legend(legendStrings, 'Location', 'northeastoutside');
hold off;

% Vizualizacia najlepsiej trasy v rovine (body + spojnice)
figure; hold on; grid on; axis equal;
xlim([0 100]);
ylim([0 100]);
scatter(B(:,1), B(:,2), 40, 'filled');
bestPts = B(globalBestRoute, :);
plot(bestPts(:,1), bestPts(:,2), 'r-', 'LineWidth', 1);
xlabel("x");
ylabel("y");
title(sprintf("Najlepsia trasa (dĺžka = %.3f)", globalBestF));
hold off;


%pomocne funkcie:

%vyber podla zvoleneho typu:
function parents = do_selection(name, pop, f, nSel)
    switch lower(name)
        case "selsus"
            parents = selsus(pop, f, nSel);
        case "selbest"
            parents = selbest(pop, f, nSel);
        case "selrand"
            parents = selrand(pop, f, nSel);
        case "seltourn"
            parents = seltourn(pop, f, nSel);
        otherwise
            error("Neznáma selekcia: %s", name);
    end
end

%vygenerovanie novej populacie:
function pop = gen_perm_pop(popsize, innerIdx)
    nInner = numel(innerIdx);
    pop = zeros(popsize, nInner);
    for i = 1:popsize
        pop(i,:) = innerIdx(randperm(nInner));
    end
end

%dlzka trasy medzi dvomi bodmi:
function f = route_length(B, route)
    pts = B(route, :);                % vyber body v poradí trasy
    dxy = diff(pts, 1, 1);            % rozdiely po krokoch
    f = sum( sqrt(sum(dxy.^2, 2)) );  % euklidovská dĺžka lomennej čiary
end

%fitnes:
function f = evaluate_pop(B, pop, startIdx, endIdx)
    n = size(pop, 1);
    f = zeros(n, 1);
    for i = 1:n
        route = [startIdx, pop(i,:), endIdx];
        f(i) = route_length(B, route);
    end
end