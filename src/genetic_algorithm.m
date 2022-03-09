
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

n = 100;
bias = - 0.1;
bestToChoose = 30;
numOfGaussians = 35;
chromosomes = randomPopulation(n, numOfGaussians);
numOfGenerations = 5000;
threshold = 0.000000008;
% threshold = 0.0001;

%% Approximation with the Genetic Algorithm

[u1, u2] = meshgrid(-1:0.05:2, -2:0.05:1);
% [u1, u2] = meshgrid(-1:0.1:2, -2:0.1:1);
mutationProbability = 0.1;

for generation = 1:numOfGenerations
    chromosomes = mutate(n, chromosomes, mutationProbability, numOfGaussians);
    
    [mse, isBest, index] = totalFitness(chromosomes, n, u1, u2, gauss,...
                            threshold, numOfGaussians, bias, f);
    if isBest == true
        plotApproximated(index, numOfGaussians, bias, chromosomes, gauss);
        break;
    end
    
    % Sort (DES) the MSEs so we can pic the first k, thus the best chromos.
    [mse, indices] = sort(abs(mse)); mse(1:bestToChoose)
    chromosomes = chooseBest(mse, indices, chromosomes, bestToChoose, n);
   
    % All the posible combinations between the best chromos.
    combinations = nchoosek(1:bestToChoose, 2);
    chromosomes = intersection(chromosomes, bestToChoose, combinations, n, numOfGaussians);
    
end

%% Plot the approximation function.

plotApproximated(1, numOfGaussians, bias, chromosomes, gauss);

%% Generates random number in the specified range

function n = random(a, b, numOfValues)
    n = (b - a) .* rand(numOfValues, 1) + a;
end

%% Each chromosome has 'numOfGaussians' gaussians with parameters c1, c2, s1, s2

function chromosome = randomChromosome(numOfGaussians)
    chromosome = zeros(numOfGaussians, 4);
    for i = 1:numOfGaussians
        chromosome(i, :) = [random(-4, 4, 1), random(-4, 4, 1),...
            random(0.1, 5, 1), random(0.1, 5, 1)];
    end
end

%% Generates random population of chromosomes

function chromosomes = randomPopulation(size, numOfGaussians)
    chromosomes = zeros(size, numOfGaussians, 4);
    for i = 1:size
        chromosomes(i, :, :) = randomChromosome(numOfGaussians);
    end
end

%% MSE of the approximation output, relative to the real output

function mse = fitness(y, y_real)
% mse = immse(y, y_real);
    len = size(y);
    y = y - y_real;
    y = y.*2;
    temp = sum(y);
    mse = sum(temp) / len(1);
end

%% Combines characteristics of two parents and produces children

function [child1, child2] = intersect(parent1, parent2, numOfGaussians)
    child1 = zeros(numOfGaussians, 4);
    child2 = zeros(numOfGaussians, 4);
    for i = 1:numOfGaussians
        if rand() < 0.2
            child1(i, :) = parent1(1, i, :);
            child2(i, :) = parent2(1, i, :);
        else
            for j = 1:4
                parameter1 = parent1(1, i, j);
                parameter2 = parent2(1, i, j);
                parameter1_binary = floatTobinary(parameter1, 3, 30);
                parameter2_binary = floatTobinary(parameter2, 3, 30);
                for k = 1:length(parameter1_binary)
                    if rand() <= 0.8 
                        temp = parameter1_binary(k);
                        parameter1_binary(k) = parameter2_binary(k);
                        parameter2_binary(k) = temp;
                    end
                end
                child1(i, j) = binaryToFloat(parameter2_binary, 3, 30);
                child2(i, j) = binaryToFloat(parameter1_binary, 3, 30);
            end
%             child1(i, :) = parent2(1, i, :);
%             child2(i, :) = parent1(1, i, :);
        end
    end
end

%% Mutates the population with a low probability

function population = mutate(n, population, probability, numOfGaussians)
    for i = 1:n
        if rand() < probability
            for j = 1:numOfGaussians
                    for z = 1:4
                        parameter = population(i, j, z);
                        parameter_binary = floatTobinary(parameter, 3, 30);
                        for k = 1:length(parameter_binary)
                            if rand() <= 0.05
                                if parameter_binary(k) == 0
                                    parameter_binary(k) = 1;
                                else
                                    parameter_binary(k) = 0;
                                end
                            end
                        end
                        population(i, j, z) = binaryToFloat(parameter_binary, 3, 30);
                    end
%                     population(i, j, :) = [random(-4, 4, 1),...
%                         random(-4, 4, 1), random(0.1, 5, 1), random(0.1, 5, 1)];
            end
        end
    end
end

%% Performs intersection amongst the best chromosomes of the population

function chromosomes = intersection(chromosomes, bestToChoose, combinations,...
    n, numOfGaussians)
    index = bestToChoose + 1;
    numOfCombinations = size(combinations);
    for i = 1:numOfCombinations
            if rand() <= 0.5
                parent1 = chromosomes(combinations(i, 1), :, :);
                parent2 = chromosomes(combinations(i, 2), :, :);
                [child1, child2] = intersect(parent1, parent2, numOfGaussians);

                chromosomes(index, :, :) = child1;
                chromosomes(index + 1, :, :) = child2;
                index = index + 2;
                if index >= n
                    index = index - 2;
                    break;
                end
            end
    end
    
    % Fill the new population if needed.
    if index < n
        while index <= n
            chromosomes(index, :, :) = randomChromosome(numOfGaussians);
            index = index + 1;
        end
    end
end

%% MSE calculation for every chromosome based on input-output

function [mse, isBest, index] = totalFitness(chromosomes, n, u1, u2, gauss,...
            threshold, numOfGaussians, bias, f)
    mse = zeros(1, n);
    index = 1;
    isBest = false;
    for chromo = 1:n
        res = 0;
        for i = 1:numOfGaussians
            res = res + bias + gauss(u1, u2, chromosomes(chromo, i, 1),...
                chromosomes(chromo, i, 2), chromosomes(chromo, i, 3),...
                chromosomes(chromo, i, 4));
        end
        mse(1, chromo) = fitness(res/numOfGaussians, f(u1, u2));
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

function plotApproximated(index, numOfGaussians, bias, chromosomes, gauss)
    [X, Y] = meshgrid(-1:0.05:2, -2:0.05:1);
    res = 0;
    for i = 1:numOfGaussians
        res = res + bias + gauss(X, Y, chromosomes(index, i, 1),...
            chromosomes(index, i, 2), chromosomes(index, i, 3),...
            chromosomes(index, i, 4));
    end
    surf(X, Y, res);
end

%% Picks the k-best chromosomes, based on their MSE

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


%% Converts float number to binary array

function bin = floatTobinary(a, integerBits, fractionBits)
    bin = [ fix(rem(fix(a) * pow2(-(integerBits - 1):0),2)),...
            fix(rem( rem(a, 1) * pow2(1:fractionBits), 2))];
end

%% Converts binary array to the float number that produced it

function float = binaryToFloat(bin, integerBits, fractionBits)
    float = bin * pow2([integerBits-1:-1:0 -(1:fractionBits)].');
end



