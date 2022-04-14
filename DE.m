function [speak, sallpeak] = DE(fn, run)
    if nargin == 0
        fn = 1;
        run = 1;
    end

    algRand = RandStream.create('mt19937ar','seed', run);
    RandStream.setGlobalStream(algRand);
    
    % initialize the problem
    pro = DMMOP(fn);
    
    N = 100;
    F = 0.5;
    CR = 0.5;
    D = pro.D;
    pop = rand(N, D) .* (pro.upper - pro.lower) + pro.lower;
    fits = pro.GetFits(pop);
    
    while ~pro.Terminate()
        disp(pro.evaluated);
        % mutation, DE/rand/1
        index = zeros(N, 3);
        for i = 1:N
            index(i, :) = randperm(algRand, N-1, 3);
            index(i, index(i, :) >= i) = index(i, index(i, :) >= i) + 1;
        end
        
        mutant = pop(index(:, 1), :) - F .* (pop(index(:, 2), :) - pop(index(:, 3), :));
        
        % crossover 
        cross = rand(algRand, N, D) < CR;
        parent_only_idx = find(sum(cross, 2) == 0);
        for i = parent_only_idx
            cross(i, randi(algRand, D)) = true;
        end
        off = cross .* mutant + (1-cross) .* pop;
        
        % bound check
        off = boundary_check(off, pro.lower, pro.upper);
        
        % evaluate
        % Not all individuals can be evaluated
        % For example, if there is only 5 fitness evaluations in the 
        % current environment, the functions only return the fitness 
        % values of the first 5 individual.
        off_fits = pro.GetFits(off);
        
        % selection
        num = length(off_fits);
        cmp = find(off_fits > fits(num));
        pop(cmp, :) = off(cmp, :);
        fits(cmp) = off_fits(cmp);
        
        % environment changes
        % If all the fitness evaluations in the current envrionment are
        % consumed, the environment will change after executing the 
        % CheckChange function.
        if pro.CheckChange(pop, fits)
            fits = pro.GetFits(pop);
        end
    end
    
    [peak, allpeak] = pro.GetPeak();
    speak = sum(peak, 2);  % the number of peaks found by the algorithm in all environments
    sallpeak = sum(allpeak); % the number of actual peaks in all environments
end