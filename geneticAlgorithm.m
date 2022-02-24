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

% Input values used to compute the real vs approximate output.
numOfInputs = 300;

%% Plot the objective, unknown function

plotObjective();
figure;

%% Generate random population

n = 12;
bestToChoose = 4;
chromosomes = randomPopulation(n);
numOfGenerations = 2000;

%% Approximation with the Genetic Algorithm

for generation = 1:numOfGenerations
    u1 = random(a1, b1, numOfInputs);
    u2 = random(a2, b2, numOfInputs);
    
    mutationProbability = 0.2;
    chromosomes = mutate(n, chromosomes, mutationProbability);
    
    mse = totalFitness(chromosomes, n, u1, u2, gauss, f);
 
    % Sort (DES) the MSEs so we can pic the first k, thus the best chromos.
    [mse, indices] = sort(mse);
    for i = 1:bestToChoose
        chromosomes(i, :, :) = chromosomes(indices(i), :, :);
    end
      
    % All the posible combinations between the best chromos.
    combinations = nchoosek(1:bestToChoose, 2);
    chromosomes = intersection(chromosomes, bestToChoose, combinations, n);
    
end

%% Plot the approximation function.

plotApproximated(chromosomes, gauss);

%% Generates random number in the specified range

function n = random(a, b, numOfValues)
    n = (b - a) .* rand(numOfValues, 1) + a;
end

%% Each chromosome has 15 gaussians with parameters c1, c2, s1, s2

function chromosome = randomChromosome()
    chromosome = zeros(15, 4);
    for i = 1:15
        chromosome(i, :) = [(rand()-0.5)*2, (rand()-0.5)*2, rand(), rand()];
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
    mse = 0;
    for i = 1:size(y)
        mse = mse + (y(i) - y_real(i)).^2;
    end
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
                if rand() <= 0.8
                    population(i, j, :) = [(rand()-0.5)*2, (rand()-0.5)*2,...
                        rand(), rand()];
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
    
    % Fill the new population if needed.
    if index < n
        while index <= n
            chromosomes(index, :, :) = randomChromosome();
            index = index + 1;
        end
    end
end

%% MSE calculation for every chromosome based on input-output

function mse = totalFitness(chromosomes, n, u1, u2, gauss, f)
    mse = zeros(1, n);
    for chromo = 1:n
        res = 0;
        for i = 1:15
            res = res + (1/15) * gauss(u1, u2, chromosomes(chromo, i, 1),...
                chromosomes(chromo, i, 2), chromosomes(chromo, i, 3),...
                chromosomes(chromo, i, 4));
        end
        mse(1, chromo) = fitness(res, f(u1, u2));
    end
end

%% Plots the unknown function

function plotObjective()
    [X, Y] = meshgrid(-1:0.1:2, -2:0.1:1);
	Z = @(u1, u2) sin(u1 + u2) .* sin(u2.^2);
    surf(X, Y, Z(X, Y));
end

%% Plots the approximated function

function plotApproximated(chromosomes, gauss)
    [X, Y] = meshgrid(-1:0.1:2, -2:0.1:1);
    res = 0;
    for i = 1:15
        res = res + gauss(X, Y, chromosomes(1, i, 1), chromosomes(1, i, 2),...
            chromosomes(1, i, 3), chromosomes(1, i, 4));
    end
    surf(X, Y, res);
end


