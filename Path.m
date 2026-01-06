classdef Path
    properties
        coefficients
        path
    end
    methods
        function obj = Path(coefficients,path)
            obj.coefficients = coefficients;
            obj.path = path;
            obj = obj.optimise(); % optimise the given path
        end
        function obj = optimise(obj)
            l = max(max(obj.path));
            ref = zeros(l); % preallocated memory for adjacency matrix
            % fill the adjacency matrix with info from the connected nodes
            for p = obj.path.'
                ref(p(1),p(2)) = 1;
                ref(p(2),p(1)) = 1;
            end
            % preallocates memory for the following
            seen = zeros(1,l); % a list of seen elements to identify closed cycles
            seeds = zeros(1,l); % a list to identify the starts of open paths
            for i = 1:l
                count = sum(ismember(obj.path,i),'all');
                if count == 0 % to avoid errors with missing elements
                    seen(i) = 1;
                elseif mod(count,2) == 1 % elements that appear once must be a seed
                    seeds(i) = 1;
                end
            end
            
            % finds all open paths
            L = sum(seeds)/2;
            newpaths = zeros(L,2); % preallocate memory for shortened path
            j = 1;
            for i = 1:l
                t0 = 0;
                t1 = i;
                if seeds(t1) == 1 % start of an open path
                    seeds(t1) = 0; % remove it from seeds to avoid double counting
                    seen(t1) = 1;
                    while seeds(t1) == 0
                        t2 = find(ref(t1,:)==1); % find elements adjacent to t1 (current position)
                        if size(t2) == 1 & t0 ~=0 % to avoid errors with two element paths
                        elseif ref(t1,t2(1)) == 1 && t2(1) ~= t0 % ensures you don't go backwards
                            t0 = t1;
                            t1 = t2(1);
                            seen(t1) = 1;
                        elseif ref(t1,t2(2)) == 1
                            t0 = t1;
                            t1 = t2(2);
                            seen(t1) = 1;
                        end
                    end
                    seeds(t1) = 0; % remove it from seeds to avoid double counting
                    newpaths(j,:) = [i,t1];
                    j = j + 1;
                end
            end

            obj.path = newpaths; % new reduced path
            
            % finds all closed cycles
            for i = 1:l
                if seen(i) == 0 % element of a closed cycle
                    t0 = 0;
                    t1 = i;
                    while seen(t1) == 0 % avoids double counting elements in the cycle
                        seen(t1) = 1;
                        t2 = find(ref(t1,:) == 1); % find elements adjacent to t1 (current position)
                        if size(t2) == 1 % to avoid errors with two element cycles
                            seen(t2) = 1;
                        elseif ref(t1,t2(1)) == 1 && t2(1) ~= t0 % ensures you don't go backwards
                            t0 = t1;
                            t1 = t2(1);
                        else
                            t0 = t1;
                            t1 = t2(2);
                        end
                    end
                    obj.coefficients(2) = obj.coefficients(2) + 1; % number of unknots
                end
            end
        end
    end
end