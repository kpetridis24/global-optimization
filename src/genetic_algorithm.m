
%% Necessary tools initialization

% Unknown function to approximate, only used to generate input-output values.
f = @(u1, u2) sin(u1 + u2) .* sin(u2.^2);

% Gaussian function used to approximate the unknown function.
gauss = @(u1, u2, c1, c2, sigma1, sigma2) exp(-((u1 - c1).^2 ./...
    (2 * sigma1^2) + (u2 - c2).^2 ./(2 * sigma2^2)));

%% Input ranges

a1 = -1;
b1 = 2;
a2 = -2;
b2 = 1;

%% Plot the objective, unknown function

plotObjective();
figure;

%% Generate random population

n = 70;
bestToChoose = 28;
chromosomes = randomPopulation(n);
numOfGenerations = 1000;
threshold = 0.001;

%% Approximation with the Genetic Algorithm

[u1, u2] = meshgrid(-1:0.01:2, -2:0.01:1);
mutationProbability = 0.2;

for generation = 1:numOfGenerations
    chromosomes = mutate(n, chromosomes, mutationProbability);
    
    [mse, isBest, index] = totalFitness(chromosomes, n, u1, u2, gauss, threshold, f);
    if isBest == true
        plotApproximated(index, chromosomes, gauss);
        break;
    end
    
    % Sort (DES) the MSEs so we can pic the first k, thus the best chromos.
    [mse, indices] = sort(abs(mse));
    chromosomes = chooseBest(mse, indices, chromosomes, bestToChoose, n);
   
    % All the posible combinations between the best chromos.
    combinations = nchoosek(1:bestToChoose, 2);
    chromosomes = intersection(chromosomes, bestToChoose, combinations, n);
    
end

%% Plot the approximation function.

plotApproximated(1, chromosomes, gauss);

%% Generates random number in the specified range

function n = random(a, b, numOfValues)
    n = (b - a) .* rand(numOfValues, 1) + a;
end

%% Each chromosome has 15 gaussians with parameters c1, c2, s1, s2

function chromosome = randomChromosome()
    chromosome = zeros(15, 4);
    for i = 1:15
        chromosome(i, :) = [random(-4, 4, 1), random(-4, 4, 1), rand(), rand()];
    end
end

%% Generates random population of chromosomes

function chromosomes = randomPopulation(size)
    chromosomes = zeros(size, 15, 4);
    for i = 1:size
        chromosomes(i, :, :) = randomChromosome();
    end
end

%% MSE of the approximation output, relative to the real output

function mse = fitness(y, y_real)
    len = size(y);
    y = y - y_real;
    y = y.*2;
    temp = sum(y);
    mse = sum(temp) / len(1);
end

%% Combines characteristics of two parents and produces children

function [child1, child2] = intersect(parent1, parent2)
    child1 = zeros(15, 4);
    child2 = zeros(15, 4);
    for i = 1:15
        if rand() < 0.5
            child1(i, :) = parent1(1, i, :);
            child2(i, :) = parent2(1, i, :);
        else
            child1(i, :) = parent2(1, i, :);
            child2(i, :) = parent1(1, i, :);
        end
    end
end

%% Mutates the population with a low probability

function population = mutate(n, population, probability)
    for i = 1:n
        if rand() < probability
            for j = 1:15
                if rand() <= 0.9
                    population(i, j, :) = [random(-4, 4, 1),...
                        random(-4, 4, 1), rand(), rand()];
                end
            end
        end
    end
end

%% Performs intersection amongst the best chromosomes of the population

function chromosomes = intersection(chromosomes, bestToChoose, combinations, n)
    index = bestToChoose + 1;
    numOfCombinations = size(combinations);
    for i = 1:numOfCombinations
%         if rand() <= 0.7
            if rand() <= 0.5
                parent1 = chromosomes(combinations(i, 1), :, :);
                parent2 = chromosomes(combinations(i, 2), :, :);
                [child1, child2] = intersect(parent1, parent2);

                chromosomes(index, :, :) = child1;
                chromosomes(index + 1, :, :) = child2;
                index = index + 2;
                if index >= n
                    index = index - 2;
                    break;
                end
            end
%         end
    end
    
    % Fill the new population if needed.
    if index < n
        while index <= n
            chromosomes(index, :, :) = randomChromosome();
            index = index + 1;
        end
    end
end

%% MSE calculation for every chromosome based on input-output

function [mse, isBest, index] = totalFitness(chromosomes, n, u1, u2, gauss, threshold, f)
    mse = zeros(1, n);
    index = 1;
    isBest = false;
    for chromo = 1:n
        res = 0;
        for i = 1:15
            res = res - 0.08 + gauss(u1, u2, chromosomes(chromo, i, 1),...
                chromosomes(chromo, i, 2), chromosomes(chromo, i, 3),...
                chromosomes(chromo, i, 4));
        end
        mse(1, chromo) = fitness(res/15, f(u1, u2));
        if abs(mse(1, chromo)) <= threshold
            isBest = true;
            break;
        end
        index = index + 1;
    end
end

%% Plots the unknown function

function plotObjective()
    [X, Y] = meshgrid(-1:0.05:2, -2:0.05:1);
	Z = @(u1, u2) sin(u1 + u2) .* sin(u2.^2);
    surf(X, Y, Z(X, Y));
end

%% Plots the approximated function

function plotApproximated(index, chromosomes, gauss)
    [X, Y] = meshgrid(-1:0.05:2, -2:0.05:1);
    res = 0;
    for i = 1:15
        res = res - 0.08 + gauss(X, Y, chromosomes(index, i, 1),...
            chromosomes(index, i, 2), chromosomes(index, i, 3),...
            chromosomes(index, i, 4));
    end
    surf(X, Y, res);
end

%%

function chromosomes = chooseBest(mse, indices, chromosomes, bestToChoose, n)
    added = zeros(1, n);
    k = 1;
    for i = 1:bestToChoose
        if ismember(mse(1, i), added) == true
            continue;
        end
        added(k) = mse(1, i);
        chromosomes(i, :, :) = chromosomes(indices(i), :, :);
        k = k + 1;
    end
end



