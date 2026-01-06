classdef Knot
    properties
        diagram
        crossings
        components
        linknumbers
        writhe
        kbracket
        jones
        paths
    end
    methods
        function obj = Knot(diagram,varargin)
            obj.diagram = diagram;
            obj.crossings = size(diagram,1);
            if size(diagram,1)>0 % not the unknot
                if nargin == 1
                    obj = obj.w(); % calculate the writhe
                    obj = obj.kauffman(); % calulcate the kauffman bracket and jones polynomial
                    if abs(subs(obj.jones,1))~=1 % check for multiple components
                        obj = obj.lk();
                    end
                else                             % info from connected sum
                    obj.writhe = varargin{1}(1);
                    obj.kbracket = varargin{1}(2);
                    obj.jones = varargin{1}(3);
                    if size(varargin{1},2) == 4 % hopf sum
                        obj.linknumbers = varargin{1}{4};
                    end
                end
            else % the unknot
                obj.writhe = 0;
                obj.kbracket = 1;
                obj.jones = 1;
            end
        end
        function obj = w(obj)
            obj.writhe = 0;
            for i = 1:obj.crossings
                if abs(obj.diagram(i,2)-obj.diagram(i,4)) == 1
                    mod = 1;
                else
                    mod = -1;
                end
                if obj.diagram(i,2)<obj.diagram(i,4) % negative crossing
                    obj.writhe = obj.writhe - 1*mod;
                else                                 % positive crossing
                    obj.writhe = obj.writhe + 1*mod;
                end
            end
        end
        function obj = lk(obj)
            m = max(max(obj.diagram));
            under = obj.diagram(:,1:2:3);
            over = obj.diagram(:,2:2:4);
            cross = [under; over];
            ref = zeros(m); % preallocate memory
            for c = cross.' % create adjacency matrix
                ref(c(1),c(2)) = 1;
                ref(c(2),c(1)) = 1;
            end
            obj.components = cell(0); % create empty cell array, no preallocation because we don't know how many components there are
            i = 0;
            while i ~= m
                i = i + 1;
                f = find(ref(i,:)==1); 
                if size(f,2) == 0                 % edge case
                    obj.components{end+1} = i;
                elseif size(f,2) == 1             % edge case
                    obj.components{end+1} = i:f;
                    i = f;
                else
                    obj.components{end+1} = i:f(2); % this is the end of the cycle, we know that the arcs increment by 1 each time so we can use this as a shortcut
                    i = f(2);
                end
            end
            l = size(obj.components,2);
            links = cell(0);
            for i = 1:l % loop over each crossing to find which crossings are part of links
                incross = [ismember(under(:,1),obj.components{i}),ismember(over(:,1),obj.components{i})];
                t = mod(sum(ismember(incross,1),2),2).';
                f = find(t==1);
                links{end+1} = f; % store the rows for each link crossing (with the respective component number)
            end
            obj.linknumbers = zeros(l);
            for i = 1:l-1
                for j = i+1:l
                    clinks = intersect(links{i},links{j}); % link crossings containing both components in question
                    linknumber = 0;
                    for link = clinks
                        if abs(obj.diagram(link,2)-obj.diagram(link,4)) == 1
                            if obj.diagram(link,2)<obj.diagram(link,4) % negative crossing
                                linknumber = linknumber-1;
                            else                                       % positive crossing
                                linknumber = linknumber+1;
                            end
                        else                                           % edge cases
                            for component = obj.components
                                if ismember(obj.diagram(link,2),component{1})
                                    L = [component{1}(1),component{1}(end)];
                                end
                            end
                            if isequal([obj.diagram(link,2),obj.diagram(link,4)],L)
                                linknumber = linknumber + 1;
                            else
                                linknumber = linknumber - 1;
                            end
                        end
                    end
                    obj.linknumbers(i,j) = linknumber/2;
                    obj.linknumbers(j,i) = linknumber/2; % fill in the linking number matrix
                end
            end
        end
        function obj = kauffman(obj)
            L = size(obj.diagram,1);
            i = 1;
            obj.paths = cell(L,2); % preallocate memory for resolved crossings
            for crossing = obj.diagram.'
                obj.paths{i,1} = Path([1,0],[crossing(1),crossing(2);crossing(3),crossing(4)]); % A resolution 
                obj.paths{i,2} = Path([-1,0],[crossing(1),crossing(4);crossing(2),crossing(3)]); % A^-1 resolution
                i = i+1;
            end
            
            resolved = {obj.paths{1,1},obj.paths{1,2}}; % initialise combination of resolutions
            
            % finds every combination of resolved crossing
            for i = 2:L
                k = 1;
                newresolved = cell(1,2^i); % preallocate memory for new resolved cell array
                for path = resolved
                    for j = 1:2
                        c = path{1}.coefficients + obj.paths{i,j}.coefficients; % new coefficients
                        p = [path{1}.path;obj.paths{i,j}.path];                 % new path
                        newresolved{k} = Path(c,p);
                        k = k + 1;
                    end
                end
                resolved = newresolved;
            end
            
            syms A
            obj.kbracket = 0; % initialise kauffman bracket
            B = -(A^2+A^-2);  % <OUD> = B*<D>
            for i = 1:2^L
                obj.kbracket = obj.kbracket + A^(resolved{i}.coefficients(1)) * B^(resolved{i}.coefficients(2)-1);
            end
            obj.kbracket = simplify(obj.kbracket);
            jpoly = (-A^3)^-obj.writhe*obj.kbracket;
            syms t
            obj.jones = simplify(subs(jpoly,A,t^(-1/4)));
        end
    end

    methods (Static)
        function Sum = ConnectedSum(K1,K2)
            if size(K1.diagram,1)==0 % left identity
                Sum = K2;
            elseif size(K2.diagram,1)==0 % right identity
                Sum = K1;
            else
                augK2 = K2.diagram + K1.crossings*2 - 1; % delete arcs at the end of K1 and the start of K2 and glue them together
                test2 = [augK2(:,1)-augK2(:,3),augK2(:,2)-augK2(:,4)];
                [r,c] = find(abs(test2)~=1);
                augK2(r,c+1+sign(test2(r,c))) = max(max(augK2)) + 1; % fix the end of the knot so it doesnt loop back to itself
    
                augK1 = K1.diagram;
                test1 = [augK1(:,1)-augK1(:,3),augK1(:,2)-augK1(:,4)];
                [r,c] = find(abs(test1)~=1);
                augK1(r,c+1-sign(test1(r,c))) = max(max(augK2)); % fix the end of the knot to remove the arc we deleted earlier
                diagram = [augK1;augK2];
                writhe = K1.writhe + K2.writhe;
                kauffman = simplify(K1.kbracket*K2.kbracket);
                jones = simplify(K1.jones*K2.jones);
                Sum = Knot(diagram,{writhe,kauffman,jones});
            end
        end

        function HopfSum = Hopf(K1,K2,orientation)
            % same general process as the connected sum, with some tweaks
            syms A t
            if orientation == '+' % sets up positive hopf link
                hopf = [2,4,1,3;4,2,3,1];
                hopfwrithe = 2;
                hopfjones = -t^(1/2)*(t^2 + 1);
                hopflinknumber = [0,1;1,0];
            elseif orientation == '-' % sets up negative hopf link
                hopf = [2,3,1,4;4,1,3,2];
                hopfwrithe = -2;
                hopfjones = -(t^2 + 1)/t^(5/2);
                hopflinknumber = [0,-1;-1,0];
            else
                error('Please enter a valid orientation, + or -.') % error if an invalid orientation is chosen
            end
            hopfbracket = -(A^8 + 1)/A^4;

            id1 = [true,false,false,false;false,true,false,true]; % section of hopf that belongs to K1
            id2 = flip(id1,1); % that belongs to K2

            
            if size(K1.diagram,1)>0 % non-identity
                m1 = max(max(K1.diagram));
                augK1 = K1.diagram;
                test1 = [augK1(:,1)-augK1(:,3),augK1(:,2)-augK1(:,4)];
                [r,c] = find(abs(test1)~=1);
                augK1(r,c+1+sign(test1(r,c))) = m1 + 1;
            else % identity
                augK1 = K1.diagram;
                m1 = 0;
            end
        
            hopf(id1) = hopf(id1) + m1;
            
            if size(K2.diagram,1)>0 % non-identity
                m2 = max(max(K2.diagram));
                augK2 = K2.diagram + m1 + 2;
                test2 = [augK2(:,1)-augK2(:,3),augK2(:,2)-augK2(:,4)];
                [r,c] = find(abs(test2)~=1);
                augK2(r,c+1+sign(test2(r,c))) = m1 + m2 + 3;
            else % identity
                augK2 = K2.diagram;
                m2 = 0;
            end
            
            hopf(id2) = hopf(id2) + m1 + m2;
            hopf(2,3) = hopf(2,3) + m1; % small adjustment to hopf link crossings
        
            diagram = [augK1;hopf;augK2];
            writhe = K1.writhe + hopfwrithe + K2.writhe;
            kauffman = K1.kbracket*hopfbracket*K2.kbracket;
            jones = K1.jones*hopfjones*K2.jones;
            linknumber = hopflinknumber;
            HopfSum = Knot(diagram,{writhe,kauffman,jones,linknumber});
        end
    end
end